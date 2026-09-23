#!/usr/bin/env python3
"""Compile current QE routines with gfortran; no DFPT data or QE build needed.

Requires NumPy and gfortran (or --fc). Numerical production routines are read
from this checkout on every run. Only disk I/O and external state are stubbed.
"""
import argparse
import ctypes as ct
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def routine(path, name):
    text = path.read_text()
    match = re.search(rf'^\s*subroutine {name}\b.*?^\s*end subroutine {name}\b[^\n]*',
                      text, re.I | re.M | re.S)
    if match is None:
        raise AssertionError(f'Cannot locate production routine {name} in {path}')
    return match.group() + '\n'


def run(args, cwd):
    subprocess.run([str(x) for x in args], cwd=cwd, check=True)


def farray(x, dtype=None):
    return np.array(x, dtype=dtype, order='F')


def invoke(fun, *args):
    fun(*[x.ctypes.data_as(ct.c_void_p) if isinstance(x, np.ndarray)
          else ct.byref(ct.c_int(x)) for x in args])


def matrix(phi):
    n = phi.shape[2]
    return phi.transpose(0, 2, 1, 3).reshape(3*n, 3*n, order='F')


def tensor(mat):
    n = mat.shape[0]//3
    return farray(mat.reshape(3, n, 3, n, order='F').transpose(0, 2, 1, 3))


def check_qstar(lib):
    rng = np.random.default_rng(39201)
    n = 3
    x = rng.normal(size=(3*n, 3*n)) + 1j*rng.normal(size=(3*n, 3*n))
    dyn = farray(x + x.conj().T)
    ops = np.zeros((3, 3, 48), dtype='i4', order='F')
    ops[:, :, 0] = np.diag([-1, 1, 1])
    inv = np.ones(48, dtype='i4')
    perm = np.ones((48, n), dtype='i4', order='F')
    perm[0] = [1, 3, 2]
    rtau = np.zeros((3, 48, n), order='F')
    rtau[0, 0] = [0, -1, -1]
    phi = tensor(dyn)
    for anti in (0, 1):
        for qx in (0.137, 0.5):
            q = farray([qx, 0.231, -0.319])
            out = np.zeros_like(phi)
            invoke(lib.rotate_and_add_dyn_, phi, out, n, 1, ops, inv, perm, rtau, q, anti)
            transform = np.zeros_like(dyn)
            for a, target in enumerate(perm[0]-1):
                phase = np.exp(-2j*np.pi*np.dot(q, rtau[:, 0, a]))
                transform[3*target:3*target+3, 3*a:3*a+3] = phase*ops[:, :, 0]
            expected = transform @ (dyn.conj() if anti else dyn) @ transform.conj().T
            np.testing.assert_allclose(matrix(out), expected, atol=2e-14, rtol=0)
            np.testing.assert_allclose(matrix(out), matrix(out).conj().T, atol=2e-14, rtol=0)
            np.testing.assert_allclose(np.linalg.eigvalsh(matrix(out)), np.linalg.eigvalsh(dyn), atol=2e-14, rtol=0)
            np.testing.assert_array_equal(phi, tensor(dyn))
    print('PASS q rotation: unitary/antiunitary, atom permutation, phase, boundary, Hermiticity, eigenvalues')

    # A noncentrosymmetric magnetic point group, with no pure time reversal.
    signs = [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]
    base_ops = [np.diag(s) for s in signs]
    base_tr = [0, 1, 1, 0]
    positions = np.array([[0., 0., 0.], [.5, .5, .5]])
    n = 2
    x = rng.normal(size=(6, 6)) + 1j*rng.normal(size=(6, 6))
    source = farray(x+x.conj().T)

    def star(order, q, source):
        ops = np.zeros((3, 3, 48), dtype='i4', order='F')
        tr = np.zeros(48, dtype='i4')
        inv = np.ones(48, dtype='i4')
        perm = np.zeros((48, n), dtype='i4', order='F')
        rtau = np.zeros((3, 48, n), order='F')
        for k, old in enumerate(order):
            ops[:, :, k] = base_ops[old]; tr[k] = base_tr[old]; inv[k] = k+1
            perm[k] = [1, 2]
            rtau[:, k, :] = (positions @ base_ops[old].T-positions).T
        count = np.zeros(1, dtype='i4'); qs = np.zeros((3, 96), order='F')
        output = np.zeros((3, 3, n, n, 96), dtype='complex128', order='F')
        before = np.zeros_like(source); after = np.zeros_like(source)
        eye = farray(np.eye(3)); work = source.copy(order='F')
        invoke(lib.audit_run_, n, 4, 1, tr, ops, inv, perm, rtau, eye, eye,
               farray(q), work, 0, count, qs, output, before, after)
        return {tuple(np.round(qs[:, i], 10)): matrix(output[:, :, :, :, i])
                for i in range(count[0])}, after

    # Generic q exercises star averaging; axial q also connects q to -q.
    for q in ([.137, .231, -.319], [.137, 0., 0.]):
        original, projected = star([0, 1, 2, 3], q, source)
        reordered, _ = star([0, 3, 1, 2], q, source)
        again, _ = star([0, 1, 2, 3], q, projected)
        assert original.keys() == reordered.keys() == again.keys()
        for point, value in original.items():
            np.testing.assert_allclose(value, reordered[point], atol=2e-14, rtol=0)
            np.testing.assert_allclose(value, again[point], atol=2e-14, rtol=0)
            np.testing.assert_allclose(value, value.conj().T, atol=2e-14, rtol=0)
            minus = tuple(-v for v in point)
            np.testing.assert_allclose(value.conj(), original[minus], atol=2e-14, rtol=0)
    print('PASS q-star: operation reordering, q/-q constraint, reality and idempotence')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fc', default='gfortran')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='qe-magnetic-unit-') as tmp:
        work = Path(tmp)
        common = [args.fc, '-O0', '-g', '-fcheck=all', '-ffree-line-length-none']
        run(common + [ROOT/'Modules/kind.f90', ROOT/'PHonon/PH/ph_nc_response_conv.f90',
                      HERE/'nc_response_probe.f90', '-o', work/'nc.x'], work)
        run([work/'nc.x'], work)
        compact = ROOT/'PHonon/PH/rotate_pattern_add.f90'
        (work/'compact.f90').write_text(routine(compact, 'compact_dyn') + routine(compact, 'scompact_dyn'))
        sources = [ROOT/'Modules/kind.f90', ROOT/'Modules/constants.f90', HERE/'qstar_probe.f90',
                   work/'compact.f90', ROOT/'PHonon/PH/trntnsc.f90',
                   ROOT/'PHonon/PH/rotate_and_add_dyn.f90', ROOT/'PHonon/PH/q2qstar_ph.f90',
                   ROOT/'LR_Modules/star_q.f90', ROOT/'PW/src/eqvect.f90']
        run(common + ['-shared', '-fPIC'] + sources + ['-o', work/'probe.so'], work)
        check_qstar(ct.CDLL(str(work/'probe.so')))
        # Test the actual export block against the spin-space identity shift.
        source = (ROOT/'PHonon/PH/solve_linter.f90').read_text()
        start = source.index('           if(lmetq0) then\n              ! A Fermi-energy shift')
        end = source.index('           endif', start) + len('           endif')
        (work/'export_shift.inc').write_text(source[start:end] + '\n')
        run(common + ['-I', work, HERE/'gamma_stubs.f90', ROOT/'LR_Modules/apply_dpot.f90',
                      HERE/'gamma_probe.f90', '-o', work/'gamma.x'], work)
        run([work/'gamma.x'], work)
    print('ALL QE MAGNETIC RESPONSE UNIT CHECKS PASSED')


if __name__ == '__main__':
    main()
