# Magnetic noncollinear SOC+U release, 2026-10-28

The later [WFPT spin/SOC/U extension](https://github.com/M1KaiLe/perturbo-github/blob/perturbo-3.0.1-ALL/docs/WFPT_SPIN_U_20261028.md) uses the same PH
response export. This update also writes `prefix.dnsbare_pattern_qN` for
scalar/collinear +U after `sym_dns_wrapper`; qe2pert reads it together with
`prefix.dnsscf_qN`. Rebuild PH and regenerate these responses for that path.

Use this QE branch with the matching
[Perturbo 3.0.1-ALL branch](https://github.com/M1KaiLe/perturbo-github/tree/perturbo-3.0.1-ALL).

## Supported calculation

Magnetic noncollinear SOC+U can retain antiunitary magnetic operations:
`no_t_rev` may be omitted or set to `.false.`. Setting it to `.true.` still
selects the unitary subgroup for comparisons. The supported Hubbard scope is
fully relativistic norm-conserving pseudopotentials, Dudarev U only,
`HUBBARD (atomic)`, one Hubbard manifold per species and one band group.
Metals use smearing. Start PH responses with `recover=.false.`.

PH exports the dynamical matrices, ordinary `dvscf`, and Hubbard occupation
responses. Ordinary qe2pert reconstructs the Hubbard contribution to EPC from
these responses and the matching electronic/Wannier data. PH's built-in
`electron_phonon` option remains unsupported with U; leave it unset.
USPP/PAW, U+V, independent Hund terms, ortho-atomic projectors and
meta-GGA phonons are outside this release's SOC+U scope.

## Changes

- PH Hubbard symmetry applies antiunitary conjugation to the complete
  orbital-spin pullback, including the mode matrix and lattice phase.
- The nonlocal dynamical matrix reloads the direct spinor when a magnetic
  pool has only one physical k point. An auxiliary time-reversed wavefunction
  left in the shared buffer previously caused a large pool-count dependence.
- The shared GGA response takes the discrete FFT derivative of the projected
  spin-density response, matching the operation used by the SCF potential.
- Save-file initialization restores the PW GGA spin axis before rebuilding
  the potential; absent XML magnetization and angle fields are initialized.
- The paired converter distinguishes antiunitary symmetry from the separate
  q/−q adjoint relation and checks their metadata independently.

The occupation order is `uu,ud,du,dd`, with orbital-spin bra indices first.
QE's antiunitary spin matrix already contains the time-reversal spin factor.
The real-space Hermiticity relation `N(-q)=N(q)^dagger` is independent of
electronic time reversal. In qe2pert the total q reversal is the exclusive-or
of the antiunitary flag and the additional q-reality flag; both matrix
operations are applied when both flags are set.

## Multi-atom coverage

Occupation responses have independent Hubbard-site and displaced-site axes:
`dns(m1,m2,spin,site,Cartesian_direction,displaced_atom)`.
The species index selects U and orbital dimension, while the atom-specific
offset selects projector columns. Displacing a non-Hubbard atom can induce
occupation changes and Hubbard-potential feedback at other sites.

The production-routine tests cover three atoms with two different U values
and orbital dimensions plus a non-Hubbard atom, all nine displacements and
the full 9×9 Hubbard force-constant matrix. Reordering both atoms and packed
projectors through all six permutations preserves EPC. Four-site file tests
use d/non-Hubbard/p/d sites, twelve complex displacement modes, and nonzero
cross-site responses. Separate magnetic tests cover atom permutations,
nonzero lattice phases, complex modes and all antiunitary/q-reality branches.

Twenty multi-atom, symmetry and derivative unit tests passed on October 27.
The force-derivative errors in the finite-dimensional production-routine
probes were at most 9.02e-13 Ry/bohr; full bare+SCF force-constant errors were
at most 5.12e-13 Ry/bohr². These probes supply independent complex wavefunction
responses. Real-material self-consistency is covered by the Fe calculations
below.

## Material validation and remaining limits

The Fe controls used fully relativistic NC pseudopotentials, atomic projectors
and U=3 eV. Electronic grids and real-space FFT resolution were matched when
comparing primitive cells and supercells.

| Check | Result |
| --- | --- |
| Two-atom Fe, Gamma PH full magnetic group versus unitary subgroup | Dynamical-matrix relative difference 7.7627e-5; maximum frequency difference 0.023956 cm⁻¹ |
| Primitive Fe ordinary EPC, 4³ k and 4³ q grids | Complex EPC relative difference 1.6611e-5; antiunitary subset 3.1032e-5; IFC difference 2.3677e-5; all 67 predefined checks passed |
| PW92 SOC+U, primitive finite-q PH versus four-atom real displacements | Full complex occupation-response difference 0.30905%; force/PH difference 0.01772% with the final production PH; energy/force curvature difference 0.07679% |
| Final production PW92 SOC+U PH full/unitary comparison | Dynamical-matrix relative difference 1.09225e-5; Cartesian Hubbard response 6.13012e-6 |

The two-atom Gamma example has inversion-forbidden first-order Hubbard
occupation changes. Nonzero occupation feedback was checked with the four-atom
cosine displacement `(h,0,-h,0)`, corresponding to primitive `q=(0,0,1/4)`.
PW92 with and without U passed independent energy, force and phonon controls.
It was a functional-consistency control with the same pseudopotential held
fixed, not a new converged material parametrization.

The initial signed-PBE finite-displacement comparison failed. The shared GGA
and save-file defects were repaired, but the existing local spin-label rule
remains non-smooth near label changes. An isolated fixed-label PBE+U control
passed occupation and force checks at both small displacement steps; its
energy/force curvature differences of 9.30% and 3.56% failed acceptance.
Fixed-label, alternative-solver and alternative-functional candidates are
excluded from production. Passing magnetic-group comparisons does not remove
this specific finite-difference limitation.

The full material EPC chain used primitive Fe. Multi-species compounds,
chemically inequivalent Hubbard sites and antiferromagnetic primitive cells
have production-routine coverage but await material-level end-to-end tests.

## Build and reproduce small checks

Build QE with the usual site configuration, then rebuild the matching
qe2pert/Perturbo against those QE archives. Keep the converter in the QE
directory so the production-routine tests can locate both source trees.
Use a serial make for qe2pert and one MPI rank per pool at runtime.

With gfortran and NumPy, run from the QE root:

```bash
python3 -B PHonon/PH/tests/test_hubbard_magnetic_symmetry.py
python3 -B PHonon/PH/tests/test_dynmat_single_k.py
python3 -B PHonon/PH/tests/test_hubbard_metal_weights.py
python3 -B LR_Modules/tests/test_gga_nc_variation.py --require-pass
python3 -B LR_Modules/tests/test_gga_nc_projection_sign.py --require-pass
python3 -B PHonon/PH/tests/dfptu_response_diagnosis/test_qexsd_copy_magnetic.py --report /tmp/qexsd-magnetic-report.json
make -C LR_Modules/tests check
```

The XML test also compiles a retained small pre-fix routine and requires that
negative control to fail. The paired Perturbo release lists its EPC, q-star,
format and multi-atom tests; the format tests additionally require h5py.

This source release contains production fixes, portable small tests and
documentation. Material scratch data, executable files, local compiler
settings, scheduler records and experimental candidate builds remain in the
local archives and are not included in the Git commits.
