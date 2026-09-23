# Magnetic DFPT release review — 2026-09-23

This review compares the local QE 7.3.1 sources with remote branch
`qe-7.3.1-dfpt-u-noncolin` at `69c6dbae977b506d1b7c5e0004f77774d3484996`.
The release branch is `qe-7.3.1-dfpt-u-noncolin-magsym`.

## Checks executed for this publication

| Check | Result |
|---|---|
| `ph.x` and `q2qstar.x` build/link with Intel 2021.3, Intel MPI 2021.3, HDF5 1.12 | Passed |
| All callers of `q2qstar_ph` / `rotate_and_add_dyn` | Reviewed: consistent explicit `t_rev` and `star_q1` order |
| Current-source q-star numerical probes | Passed: unitary/antiunitary, atom permutation, nonzero translation phase, zone boundary, Hermiticity, eigenvalues, operation reordering, q/-q constraint and idempotence |
| Current-source NC residual probe | Passed: NaN padding excluded, partition and sample-count invariance, magnetic-axis rotation, irrep duplication, finite/inner-solve stop gates |
| Current-source Gamma potential-shift probe | 84 configurations passed; largest error below `1e-12` |
| `LR_Modules/tests` Hubbard response | Passed; existing DFPT+U restrictions retained |
| Paired Perturbo release regression | 159 executable runs met assertions, including 22 intentional invalid-input/cache rejections; 1/4 MPI coverage |

The small source probes used gfortran 15 with bounds/runtime checks and NumPy.
They compile production routines with minimal state/I/O stubs, not replacement
implementations of the numerical transformations. Logs and hashes are in
[`PHonon/PH/tests/release-20260923`](../PH/tests/release-20260923).
The paired Perturbo repository contains its detailed run records under
`porting/release-20260923` and its release review under `porting`.

Review also found that Perturbo executable targets did not depend on the linked
QE archives. Both makefiles now include those prerequisites. Make dry runs
confirm that changing `liblrmod.a` forces both executables to relink. This matters
for the Gamma direct-state fix, which lives in the QE response library.

## Material evidence and limits

The earlier Fe comparison used NC/PBE+SOC without U, identical density and
Wannier basis, a 12³ DFPT electronic mesh, a 4³ q mesh, `tr2_ph=1e-16`, and
`nmix_ph=20`. Magnetic and unitary-only runs had respectively 16/8 operations
(8/0 antiunitary). The full 64-point maximum frequency difference was
`7.84288501e-4 cm^-1`, below the 0.1 cm^-1 comparison target. The raw reality
residual was zero; the scaled Hermiticity residual was at most `2.87e-16`.
The paired EPR had identical electronic bands and a complex Cartesian EPC
relative difference of `5.16470744e-6`.

These are **historical reported results**: the large validation working directory
was deleted at the owner's request before this publication review. A fresh full
material DFPT run was not performed here. The standalone tests and publication
regressions above were rerun and their compact evidence is retained.

The Fe comparison demonstrates numerical consistency within the tested setup,
not physical convergence: Gamma before ASR remained near -75 to -76 cm^-1,
different ASR choices affected the interpolated spectrum, and electronic-mesh
and Wannier convergence still require production checks. Real-material WFPT,
arbitrary magnetic directions, noncentrosymmetric magnetic materials, GPU,
USPP/PAW and meta-GGA are not claimed as newly validated end-to-end paths.

See [migration and build instructions](../../README_MAGNETIC_DFPT.md) before
reusing old response data or caches. Passing the checks above is not a claim
that every possible configuration of QE is error-free.
