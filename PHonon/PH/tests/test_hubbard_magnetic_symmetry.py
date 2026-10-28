#!/usr/bin/env python3
"""Exercise production PH Hubbard symmetrizers with independent spinor oracles.

Run: python3 PHonon/PH/tests/test_hubbard_magnetic_symmetry.py
The compiler artifacts live in a temporary directory and are removed afterwards.
Only small linear-algebra problems are used; no material calculation is run.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
J = np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=complex)


@dataclass
class Operation:
    orbital: np.ndarray
    spin: np.ndarray  # PW comp_dspinldau: already J U* for antiunitary ops.
    anti: bool
    atoms: np.ndarray
    translation: np.ndarray
    modes: np.ndarray
    cartesian: np.ndarray


def identity(npe: int) -> Operation:
    return Operation(np.eye(3), np.eye(2, dtype=complex), False,
                     np.arange(2), np.zeros((3, 2)),
                     np.eye(npe, dtype=complex), np.eye(3, dtype=int))


def full(blocks: np.ndarray) -> np.ndarray:
    """PW order: orbital is fastest; spin blocks are uu,ud,du,dd."""
    return np.block([[blocks[:, :, 0], blocks[:, :, 1]],
                     [blocks[:, :, 2], blocks[:, :, 3]]])


def blocks(matrix: np.ndarray) -> np.ndarray:
    return np.stack([matrix[:3, :3], matrix[:3, 3:],
                     matrix[3:, :3], matrix[3:, 3:]], axis=2)


def adjoint(data: np.ndarray) -> np.ndarray:
    result = np.empty_like(data)
    for ia in range(2):
        for ip in range(data.shape[-1]):
            result[:, :, :, ia, ip] = blocks(full(data[:, :, :, ia, ip]).conj().T)
    return result


def spinor_inverse(op: Operation) -> np.ndarray:
    # Derive from c'=A c* (anti) or c'=A c (unitary), rather than copying
    # hubbard_rotate_nc's index expression. Anti inversion is c=A^T c'^*.
    spin = op.spin.T if op.anti else op.spin.conj().T
    return np.kron(spin, op.orbital.T)


def operation_oracle(data: np.ndarray, op: Operation, q: np.ndarray,
                     mode_matrix: np.ndarray | None = None) -> np.ndarray:
    coeffs = op.modes if mode_matrix is None else mode_matrix
    inv = spinor_inverse(op)
    output = np.zeros_like(data)
    for ia, mapped in enumerate(op.atoms):
        phase = np.exp(2j * np.pi * q @ op.translation[:, ia])
        for ip in range(data.shape[-1]):
            response = np.zeros((6, 6), complex)
            for jp in range(data.shape[-1]):
                coefficient = phase * coeffs[jp, ip]
                matrix = full(data[:, :, :, mapped, jp])
                if op.anti:
                    matrix = matrix.conj()
                    coefficient = coefficient.conjugate()
                # N_ab=c_a* c_b; apply the inverse spinor basis to both slots.
                response += coefficient * inv.conj() @ matrix @ inv.T
            output[:, :, :, ia, ip] = blocks(response)
    return output


def serialize(array: np.ndarray) -> str:
    values = np.asarray(array).ravel(order="F")
    if np.iscomplexobj(array):
        return " ".join(f"({z.real:.17e},{z.imag:.17e})" for z in values)
    return " ".join(str(z) for z in values)


class MagneticHubbardSymmetry(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        compiler = shutil.which(os.environ.get("FC", "gfortran"))
        if compiler is None:
            raise unittest.SkipTest("gfortran is required for the production probe")
        cls.directory = tempfile.TemporaryDirectory(prefix="qe_hubbard_magnetic_")
        cls.executable = Path(cls.directory.name) / "probe.x"
        subprocess.run([compiler, "-O0", "-g", "-fcheck=all",
                        "-ffree-line-length-none", "-o", str(cls.executable),
                        str(HERE / "hubbard_magnetic_stubs.f90"),
                        str(ROOT / "LR_Modules/hubbard_nc_response.f90"),
                        str(ROOT / "PHonon/PH/sym_dns.f90"),
                        str(ROOT / "PHonon/PH/syme_dns.f90"),
                        str(HERE / "hubbard_magnetic_probe.f90")],
                       cwd=cls.directory.name, check=True, capture_output=True, text=True)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.directory.cleanup()

    def setUp(self) -> None:
        self.rng = np.random.default_rng(20261027)

    def sample(self, npe: int) -> np.ndarray:
        shape = (3, 3, 4, 2, npe)
        return self.rng.normal(size=shape) + 1j * self.rng.normal(size=shape)

    def run_probe(self, data: np.ndarray, ops: list[Operation], q=None,
                  minus_op: int = 0, tmq=None, electric: bool = False,
                  nsymq: int | None = None) -> np.ndarray:
        npe = data.shape[-1]
        q = np.zeros(3) if q is None else np.asarray(q)
        tmq = np.eye(npe, dtype=complex) if tmq is None else tmq
        count = len(ops) if nsymq is None else nsymq
        lines = [f"{2 if electric else 1} {npe} {count} {len(ops)} {minus_op}", serialize(q)]
        for op in ops:
            lines += [f"{int(op.anti)} " + serialize(op.atoms + 1),
                      serialize(op.orbital), serialize(op.spin),
                      serialize(op.translation), serialize(op.cartesian),
                      serialize(op.modes)]
        lines += [serialize(tmq), serialize(data)]
        run = subprocess.run([str(self.executable)], input="\n".join(lines) + "\n",
                             text=True, capture_output=True, check=True)
        numbers = np.fromstring(run.stdout, sep=" ").reshape(-1, 2)
        return (numbers[:, 0] + 1j * numbers[:, 1]).reshape(data.shape, order="F")

    def close(self, actual, expected):
        np.testing.assert_allclose(actual, expected, atol=8e-14, rtol=8e-14)

    def c2_anti(self, npe: int, complex_modes: bool = True) -> Operation:
        # C2x*Theta has A=J U* with U=-i sigma_x; its anti action squares to +1.
        spatial_spin = -1j * np.array([[0.0, 1.0], [1.0, 0.0]])
        op = identity(npe)
        op.orbital = np.diag([1., -1., -1.])
        op.cartesian = op.orbital.astype(int)
        op.spin = J @ spatial_spin.conj()
        op.anti = True
        op.atoms = np.array([1, 0])
        if complex_modes:
            a = .43
            op.modes = np.exp(.61j) * np.array([[np.cos(a), np.sin(a)],
                                                [np.sin(a), -np.cos(a)]])
            op.translation[0, :] = .375
        else:
            op.modes = np.eye(npe, dtype=complex)
        return op

    def test_antiunitary_complex_modes_phase_and_atom_exchange(self):
        data = self.sample(2)
        op = self.c2_anti(2)
        q = np.array([.37, .13, .19])
        transformed = operation_oracle(data, op, q)
        expected = (data + transformed) / 2
        self.close(self.run_probe(data, [identity(2), op], q=q), expected)
        # This operation closes even with complex mode gauge and nonzero phase.
        self.close(operation_oracle(transformed, op, q), data)
        self.close(self.run_probe(expected, [identity(2), op], q=q), expected)
        # A transpose at fixed q would erase this genuinely non-Hermitian case.
        self.assertGreater(np.max(np.abs(expected - adjoint(expected))), .1)

    def test_general_inverse_spin_rotation_from_wavefunctions(self):
        op = self.c2_anti(2)
        axis = np.array([.2, .7, -.3]); axis /= np.linalg.norm(axis)
        sigma = np.array([[axis[2], axis[0]-1j*axis[1]],
                          [axis[0]+1j*axis[1], -axis[2]]])
        u = np.cos(.38)*np.eye(2) - 1j*np.sin(.38)*sigma
        op.spin = J @ u.conj()
        angle = .41
        op.orbital = np.array([[np.cos(angle), -np.sin(angle), 0.],
                               [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]])
        data = self.sample(2)
        q = np.array([.37, .13, .19])
        expected = (data + operation_oracle(data, op, q)) / 2
        self.close(self.run_probe(data, [identity(2), op], q=q), expected)
        # Independent rank-one wavefunction check of the anti inverse.
        c = self.rng.normal(size=6) + 1j*self.rng.normal(size=6)
        forward = np.kron(op.spin, op.orbital) @ c.conj()
        recovered = spinor_inverse(op) @ forward.conj()
        self.close(recovered, c)
        # Antiunitary pullback must conjugate phase and mode coefficients too.
        wrong = op.modes.conj().copy()
        wrong_op = Operation(op.orbital, op.spin, op.anti, op.atoms,
                             -op.translation, wrong, op.cartesian)
        self.assertGreater(np.max(np.abs(operation_oracle(data, wrong_op, q)
                                        - operation_oracle(data, op, q))), .1)

    def test_unitary_inverse_direction_regression(self):
        data = self.sample(2)
        op = self.c2_anti(2)
        op.anti = False
        op.spin = np.diag(np.exp(1j*np.array([.31, -.31])))
        a = .37
        op.orbital = np.array([[np.cos(a), -np.sin(a), 0.],
                               [np.sin(a), np.cos(a), 0.], [0., 0., 1.]])
        q = np.array([.27, -.11, .13])
        expected = (data + operation_oracle(data, op, q)) / 2
        self.close(self.run_probe(data, [identity(2), op], q=q), expected)

    def test_minus_q_and_antiunitary_are_independent_operations(self):
        data = self.sample(2)
        op = self.c2_anti(2)
        q = np.array([.37, .13, .19])
        # The generalized combination is not a normal magnetic PH setup
        # (pure time reversal is disabled there), but distinguishes both maps.
        transformed = operation_oracle(data, op, q, op.modes)
        expected = (data + adjoint(transformed)) / 2
        actual = self.run_probe(data, [identity(2), op], q=q, minus_op=2,
                                tmq=op.modes, nsymq=1)
        self.close(actual, expected)
        # Neither dropping the anti map nor replacing adjoint by conjugation works.
        self.assertGreater(np.max(np.abs(expected - (data + transformed.conj())/2)), .1)

    def test_minus_q_identity_hermitian_projection(self):
        data = self.sample(2)
        expected = (data + adjoint(data))/2
        self.close(self.run_probe(data, [identity(2)], minus_op=1), expected)
        self.close(self.run_probe(expected, [identity(2)], minus_op=1), expected)

    def test_electric_field_antiunitary_hermitian_projection(self):
        data = self.sample(3)
        op = self.c2_anti(3, complex_modes=False)
        # Electric-vector coefficients are s(ip,jp), hence transpose into
        # operation_oracle's response-index convention (jp,ip).
        op.modes = op.cartesian.T.astype(complex)
        hermitian = (data + adjoint(data))/2
        expected = (hermitian + operation_oracle(hermitian, op, np.zeros(3)))/2
        result = self.run_probe(data, [identity(3), op], electric=True)
        self.close(result, expected)
        self.close(result, adjoint(result))
        self.close(self.run_probe(result, [identity(3), op], electric=True), result)


if __name__ == "__main__":
    unittest.main(verbosity=2)
