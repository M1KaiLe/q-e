# QE 7.3.1 magnetic DFPT branch

Branch: `qe-7.3.1-dfpt-u-noncolin-magsym`.
Use with [Perturbo 3.0.1-ALL](https://github.com/M1KaiLe/perturbo-github/tree/perturbo-3.0.1-ALL).

This update retains antiunitary magnetic operations through the static DFPT
dynamical-matrix export and the matching qe2pert/Perturbo workflow. Its new
magnetic-group scope is **norm-conserving, no Hubbard U, LDA/GGA, noncollinear
magnetism with SOC**. `no_t_rev=.true.` remains supported. Existing DFPT+U
restrictions are retained; USPP/PAW and meta-GGA are not newly validated.

## Numerical changes

- q-star generation and both dynamical-matrix rotation paths share the same
  operation order and `t_rev`. Antiunitary rotations conjugate the source
  matrix, while translation phases use the target q without a second
  conjugation. The static Hessian relation `D(-q)=D(q)*` is distinct from
  electronic time-reversal symmetry and from the MBC response.
- `q2qstar.x` reconstructs magnetic operations from XML local moments and
  preserves polar tensors. `--unitary-magnetic` retains the unitary-only mode.
- At Gamma, the Fermi-energy shift is subtracted from the scalar Pauli
  component, not from the three magnetic-field components. Shared fine/coarse
  grid storage receives the shift only once.
- Magnetic `localdos` explicitly reads the direct wavefunction even when a
  pool holds only one physical k point; it must not reuse the auxiliary
  time-reversed spinor left by a previous response solve.
- NC/no-U magnetic response checks the scalar and three-component magnetic
  residual separately on physical FFT points. Each residual is an RMS divided
  by `max(1, RMS(input), RMS(output))`; the maximum squared residual over
  channels and perturbations must be below `tr2_ph`. Both Sternheimer branches
  must converge. This removes dependence on padding and irrep multiplicity.
  The convergence criterion differs from older builds; reusing the same
  numerical threshold does not imply the same accuracy.
- Standard output prints compact iteration summaries. `verbosity='high'`
  enables per-perturbation residual diagnostics. No new required input exists.

## Build and small tests

Build QE with the site's usual compiler/MPI/library configuration (`make pw ph`),
then build Perturbo against these libraries with the same toolchain. A change of
compiler or module interface requires a clean rebuild. qe2pert must be built
serially. Both Perturbo makefiles track linked QE archives to trigger relinking.

The numerical unit checks compile the **current** source with gfortran, NumPy,
and temporary I/O stubs. They do not require saved DFPT data:

```bash
python PHonon/PH/tests/test_magnetic_response.py --fc gfortran
make -C LR_Modules/tests check
```

Checks include antiunitary/unitary rotations, atom permutation, translation
phases, a zone-boundary case, Hermiticity, eigenvalues, operation reordering,
q/-q constraints, residual invariance to FFT partition/padding, and 84 Gamma
potential-shift configurations. See [the release record](PHonon/Doc/RELEASE_REVIEW_20260923.md)
for the tests actually executed and the limits of material validation.

## Data migration

The paired qe2pert checks the complete raw dynamical-matrix mesh for
`D(-q)=D(q)*` before IFC ASR or discarding imaginary components; scaled errors
above `1e-8` stop conversion. The EPR stores spatial rotations separately from
antiunitary flags; Perturbo uses their combined momentum action for reduction.

If only q-star export was wrong, regenerate `.dyn` from **complete** saved
`PARTIAL_DYN` and displacement patterns, using a new output directory:

```bash
python PHonon/PH/tools/regenerate_magnetic_dyn.py \
  --dyn-dir /path/to/original/dyn \
  --phsave /path/to/original/PREFIX.phsave \
  --prefix PREFIX --output /path/to/new/dyn
```

Add `--unitary-magnetic` when matching a `no_t_rev=.true.` source calculation.
If the original response suffered from the Gamma potential/direct-state or
convergence defects, recompute that response first: re-export cannot fix it.
Regenerate affected IFC, EPR, setup grids, screening tables and scattering
caches. Restarting an old interrupted response is not part of this validation.

The previous large `PHonon/PH/validation` working directory was removed at the
owner's request. It is not required by these small tests, and historical material
summaries are explicitly distinguished from newly executed release checks.
