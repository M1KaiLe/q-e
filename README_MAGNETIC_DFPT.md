# QE 7.3.1 magnetic DFPT branch

Branch: `qe-7.3.1-dfpt-u-noncolin-magsym`.
Use with [Perturbo 3.0.1-ALL](https://github.com/M1KaiLe/perturbo-github/tree/perturbo-3.0.1-ALL).

The September release retained antiunitary magnetic operations through the
static DFPT dynamical-matrix export and matching qe2pert/Perturbo workflow for
**norm-conserving, no Hubbard U, LDA/GGA, noncollinear magnetism with SOC**.

The October 27 update also implements antiunitary operations in the experimental
noncollinear SOC+U response and ordinary qe2pert EPC path. Magnetic calculations
no longer require `no_t_rev=.true.`; that unitary-only setting remains available.
The +U scope is still norm-conserving, fully relativistic pseudopotentials,
Dudarev U, `atomic` projectors, one Hubbard channel per species, and one band
group. USPP/PAW, meta-GGA, general Hund couplings, recovery, and WFPT+U are not
enabled by this change. See the
[mathematical audit](PHonon/Doc/DFPTU_RELEASE_20261028.md) for conventions,
tests, and the distinction between algebraic checks and material validation.
The [October validation record](PHonon/Doc/DFPTU_RELEASE_20261028.md)
reports the actual Fe group comparisons, EPC results, independent finite
displacements, and remaining acceptance limits.
The mathematical probes, Fe magnetic-group comparisons, and EPC comparisons
pass. The initial PBE finite-displacement campaign had a 2.53% occupation-response
error and force constants with opposite signs; those results are preserved in
the historical record. The subsequent
[response diagnosis](PHonon/Doc/DFPTU_RELEASE_20261028.md) repaired
shared GGA/read-in defects and established independent energy, force and phonon
agreement for PW92-LSDA SOC with and without U. For U=3 eV, the occupation-response
error is 0.30893%, the force/PH difference is 0.06012%, and energy/force agreement
is 0.07679%. The final production PH full/unitary comparison passes with complex
dynamical-matrix and Cartesian Hubbard-response relative differences of
1.09225e-5 and 6.13012e-6. The isolated PBE fixed-label control passes both
steps' occupation and force checks, while its energy/force differences of
9.30% and 3.56% fail acceptance. The default signed-PBE label rule remains
non-smooth; fixed-label and solver candidates remain outside production.

For +U EPC, export the ordinary PH `dvscf` and Hubbard responses, then use the
paired qe2pert/Perturbo path. PH's built-in `electron_phonon` mode still rejects
+U. Metals use smearing, and response calculations start with `recover=.false.`.
The [SOC+U release notes](PHonon/Doc/DFPTU_RELEASE_20261028.md)
record the supported settings, multi-atom coverage and material validation.
Generated calculation data and experimental builds are retained separately
from this source release.

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
- Hubbard spinor response symmetrization conjugates the complete antiunitary
  pullback, including the displacement matrix and translation phase. The spin
  representation already includes the time-reversal spin factor; the separate
  q/-q occupation relation uses the combined orbital-spin adjoint.
- The nonlocal bare dynamical matrix reloads the direct spinor even when a pool
  contains only one physical k point. Hubbard initialization can leave an
  auxiliary time-reversed spinor in the shared wavefunction cache; reusing it
  made the result depend severely on the pool count.
- The shared GGA kernel applies the discrete FFT derivative to the projected
  spin-density response, matching the SCF gradient operation.
- Save-file readers restore the PW GGA spin axis before reconstructing the
  potential, and initialize absent XML magnetization/angle fields explicitly.

## Build and small tests

Build QE with the site's usual compiler/MPI/library configuration (`make pw ph`),
then build Perturbo against these libraries with the same toolchain. A change of
compiler or module interface requires a clean rebuild. qe2pert must be built
serially. Both Perturbo makefiles track linked QE archives to trigger relinking.

The numerical unit checks compile the **current** source with gfortran, NumPy,
and temporary I/O stubs. They do not require saved DFPT data:

```bash
python PHonon/PH/tests/test_magnetic_response.py --fc gfortran
python PHonon/PH/tests/test_hubbard_magnetic_symmetry.py
python PHonon/PH/tests/test_dynmat_single_k.py
python PHonon/PH/tests/test_hubbard_metal_weights.py
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
