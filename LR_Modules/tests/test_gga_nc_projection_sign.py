#!/usr/bin/env python3
"""Audit NC GGA projection signs and rotations against the discrete SCF map.

The probe compiles the actual PW compute_rho and selected LR dgradcorr.
Its FFT and constant quadratic GGA coefficients reuse the complementary
test_gga_nc_variation probe. No material calculation is performed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np

from test_gga_nc_variation import derivative, sample, serial

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

PROGRAM = """program probe
  use kinds
  use fft_types
  use gc_lr
  use noncollin_module
  implicit none
  type(fft_type_descriptor) :: dfft
  integer :: n,k,is
  real(dp) :: q(3),a,c
  real(dp),allocatable :: rho(:,:),spins(:,:),grho(:,:,:),zero(:,:,:),ss(:,:,:),g(:,:)
  complex(dp),allocatable :: drho(:,:),dv(:,:)
  read(*,*) n,q,a,c,ux,lsign
  dfft%nnr=n; dfft%ngm=n
  allocate(rho(n,4),spins(n,2),grho(3,n,2),zero(n,2,2),ss(n,2,2),g(3,n))
  allocate(drho(n,4),dv(n,4),gmag(3,n,4),vsgga(n),segni(n))
  read(*,*) rho
  read(*,*) drho
  read(*,*) vsgga
  call compute_rho(rho,spins,segni,n)
  grho=0;zero=0;g=0;dv=0;gmag=0;ss=c
  ss(:,1,1)=a;ss(:,2,2)=a
  call dgradcorr(dfft,rho,grho,zero,zero,zero,ss,q,drho,4,2,g,dv)
  do k=1,n
    write(*,'(3ES25.16)') segni(k),spins(k,:)
  enddo
  do is=1,4
    do k=1,n
      write(*,'(2ES25.16)') real(dv(k,is)),aimag(dv(k,is))
    enddo
  enddo
end program
"""


def potential(rho, axis, signed, cell_count=1):
    mag = np.linalg.norm(rho[:, 1:], axis=1)
    sign = np.where(rho[:, 1:] @ axis >= 0, 1., -1.) if signed else np.ones(len(rho))
    spins = np.column_stack(((rho[:, 0] + sign * mag) / 2,
                             (rho[:, 0] - sign * mag) / 2))
    lap = derivative(derivative(spins, cell_count), cell_count)
    vsp = -.7 * lap - .23 * lap[:, ::-1]
    value = np.column_stack((vsp.sum(axis=1) / 2,
        (sign * (vsp[:, 0] - vsp[:, 1]) / (2 * mag))[:, None] * rho[:, 1:]))
    return value, vsp, sign, spins


def finite_variation(rho, change, axis, signed, q, h):
    if q == 0:
        value = np.zeros_like(change)
        for direction, factor in ((change.real, 1), (change.imag, 1j)):
            value += factor * (potential(rho + h * direction, axis, signed)[0]
                - potential(rho - h * direction, axis, signed)[0]) / (2 * h)
        return value
    assert q == .25
    n = len(rho)
    x = 2 * np.pi * np.arange(4 * n) / n
    base = np.tile(rho, (4, 1))
    direction = (np.tile(change, (4, 1)) * np.exp(1j * q * x)[:, None]).real
    value = (potential(base + h * direction, axis, signed, 4)[0]
        - potential(base - h * direction, axis, signed, 4)[0]) / (2 * h)
    return 2 * np.mean((value * np.exp(-1j * q * x)[:, None]).reshape(4, n, 4), axis=0)


def build(source, directory):
    original = (HERE / 'gga_nc_variation_stubs.f90').read_text()
    stubs = original.split('program probe\n')[0]
    anchor = 'module noncollin_module\n'
    assert stubs.count(anchor) == 1
    stubs = stubs.replace(anchor, anchor + '  use kinds\n  logical :: lsign=.false.\n'
                          '  real(dp) :: ux(3)=0\n')
    temporary_source = directory / 'projection_probe.f90'
    temporary_source.write_text(stubs + PROGRAM)
    exe = directory / 'projection_probe.x'
    subprocess.run([shutil.which('gfortran'), '-cpp', '-O0', '-g', '-fcheck=all',
        '-ffree-line-length-none', str(temporary_source), str(ROOT / 'PW/src/compute_rho.f90'),
        str(source), '-o', str(exe)], cwd=directory, check=True, capture_output=True, text=True)
    return exe


def run(exe, rho, change, axis, signed, q):
    n = len(rho)
    _, vsp, sign, spins = potential(rho, axis, signed)
    first = f'{n} {q} 0 0 .7 .23 ' + serial(axis) + (' T' if signed else ' F')
    data = '\n'.join((first, serial(rho), serial(change), serial(vsp[:, 0] - vsp[:, 1]))) + '\n'
    result = subprocess.run([str(exe)], input=data, text=True, capture_output=True, check=True)
    values = np.fromstring(result.stdout, sep=' ')
    projection = values[:3 * n].reshape(n, 3)
    response = values[3 * n:].reshape(-1, 2)
    response = (response[:, 0] + 1j * response[:, 1]).reshape(n, 4, order='F')
    error = float(np.max(np.abs(projection - np.column_stack((sign, spins)))))
    return response, error


def rotation():
    axis = np.array([.3, .5, .8])
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    cross = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    theta = .913
    return np.eye(3) + np.sin(theta) * cross + (1 - np.cos(theta)) * (cross @ cross)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT / 'LR_Modules/dgradcorr.f90')
    parser.add_argument('--report', type=Path)
    parser.add_argument('--require-pass', action='store_true')
    args = parser.parse_args()
    records = []
    covariances = []
    rotate = rotation()
    with tempfile.TemporaryDirectory(prefix='qe-gga-projection-sign-') as temporary:
        exe = build(args.source, Path(temporary))
        for n, q in ((17, 0.), (17, .25), (65, .25)):
            rho, change = sample(n, textured=True)
            rho[:, 1] += .07  # Keep every sampled sign branch away from zero.
            for label, axis, signed in (
                ('unsigned', np.array([0., 0., 1.]), False),
                ('positive', np.array([0., 0., 1.]), True),
                ('negative', np.array([0., 0., -1.]), True),
                ('spatially-changing', np.array([1., 0., 0.]), True),
            ):
                response, projection_error = run(exe, rho, change, axis, signed, q)
                first = finite_variation(rho, change, axis, signed, q, 2e-5)
                second = finite_variation(rho, change, axis, signed, q, 1e-5)
                fd = (4 * second - first) / 3
                relative = float(np.linalg.norm(response - fd) / np.linalg.norm(fd))
                records.append(dict(n=n, q=q, sign_case=label,
                    min_abs_m_dot_axis=float(np.min(np.abs(rho[:, 1:] @ axis))),
                    production_compute_rho_max_error=projection_error,
                    relative_error=relative,
                    finite_difference_halving_change=float(np.linalg.norm(second - first) / np.linalg.norm(fd)),
                    pass_=bool(relative < 1e-7 and projection_error < 1e-14)))
                rotated_rho = rho.copy()
                rotated_change = change.copy()
                rotated_rho[:, 1:] = rho[:, 1:] @ rotate.T
                rotated_change[:, 1:] = change[:, 1:] @ rotate.T
                rotated_response, rotated_projection_error = run(exe, rotated_rho, rotated_change,
                    rotate @ axis, signed, q)
                expected = response.copy()
                expected[:, 1:] = response[:, 1:] @ rotate.T
                covariance = float(np.linalg.norm(rotated_response - expected) / np.linalg.norm(expected))
                covariances.append(dict(n=n, q=q, sign_case=label, relative_error=covariance,
                    production_compute_rho_max_error=rotated_projection_error,
                    pass_=bool(covariance < 1e-12 and rotated_projection_error < 1e-14)))
    report = dict(source=str(args.source.resolve()),
        source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
        compute_rho_sha256=hashlib.sha256((ROOT / 'PW/src/compute_rho.f90').read_bytes()).hexdigest(),
        scope='Quadratic discrete GGA, actual PW projection, fixed sign branches, complex q amplitudes and global spin covariance',
        records=records, spin_covariance=covariances,
        all_pass=all(row['pass_'] for row in records + covariances))
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    if args.require_pass and not report['all_pass']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
