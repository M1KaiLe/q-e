#!/usr/bin/env python3
"""Compare the production two-branch metal weights with a spectral derivative.

The oracle diagonalizes a complete four-cell spin/orbital Hamiltonian at +/-h.
The response instead calls QE's actual orthogonalize, solves its shifted finite
matrix exactly, and calls the actual Hubbard Kramers partner. This checks the
occupation weights, branch transpose and signs jointly, without material SCF or
an assumed time-reversal/inversion symmetry of the magnetic Hamiltonian.
One control also compares the self-consistent Dudarev response with independent
nonlinear mean-field finite differences at fixed electron number. These probes
do not test the plane-wave CG solver, moving projectors or QE's SCF mixer.
"""
from pathlib import Path
import math
import os
import shutil
import subprocess
import tempfile
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
J = np.kron(np.array([[0., 1.], [-1., 0.]]), np.eye(2))


def occupation(energy, sigma, scheme):
    x = -np.asarray(energy) / sigma
    if scheme == -99:
        return 1 / (1 + np.exp(-x))
    if scheme == -1:
        x = x - 1 / np.sqrt(2)
        return .5 * np.vectorize(math.erf)(x) + np.exp(-x*x) / np.sqrt(2*np.pi) + .5
    return .5 * np.vectorize(math.erfc)(-x)


def matrix_lines(a):
    return [f"{z.real:.17e} {z.imag:.17e}" for z in a.ravel(order="F")]


class HubbardMetalWeights(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which(os.environ.get("FC", "gfortran"))
        if compiler is None:
            raise unittest.SkipTest("Fortran compiler unavailable")
        cls.temporary = tempfile.TemporaryDirectory(prefix="qe-metal-weights-")
        cls.addClassCleanup(cls.temporary.cleanup)
        build = Path(cls.temporary.name)
        # Split the small driver out so its USE sees the production module.
        stubs, driver = (HERE / "hubbard_metal_weights_stubs.f90").read_text().split("PROGRAM probe\n", 1)
        (build / "stubs.f90").write_text(stubs)
        (build / "probe.f90").write_text("PROGRAM probe\n" + driver)
        cls.binary = build / "probe.x"
        subprocess.run([compiler, "-cpp", "-O0", "-fcheck=all", "-ffree-line-length-none",
                        str(build / "stubs.f90"),
                        str(ROOT / "LR_Modules/hubbard_nc_response.f90"),
                        str(ROOT / "LR_Modules/orthogonalize.f90"),
                        str(ROOT / "Modules/wgauss.f90"),
                        str(ROOT / "Modules/w0gauss.f90"),
                        str(build / "probe.f90"), "-o", str(cls.binary)],
                       cwd=build, check=True, text=True, capture_output=True)

    def run_probe(self, action, entries):
        result = subprocess.run([str(self.binary), action], input="\n".join(entries)+"\n",
                                text=True, capture_output=True, check=True)
        values = np.fromstring(result.stdout, sep=" ").reshape(-1, 2)
        return (values[:, 0]+1j*values[:, 1]).reshape((4, 4), order="F")

    def model(self, mixed_occupied=False, near_degenerate=False):
        rng = np.random.default_rng(271027)
        random = lambda shape: rng.normal(size=shape)+1j*rng.normal(size=shape)
        vectors = np.array([np.linalg.qr(random((4, 4)))[0] for _ in range(4)])
        energies = np.array([[-.18, -.055, .018, .073],
                             [-.15, -.031, .029, .082],
                             [-.13, -.041, .011, .091],
                             [-.17, -.028, .035, .067]])
        occupied = np.full(4, 4)
        if mixed_occupied:
            energies[[0, 2], 3] = 1.5
            occupied[[0, 2]] = 3
        if near_degenerate:
            energies[1, 1] = energies[0, 1] + 3e-6
        perturbation = .02 * random((4, 4, 4))
        return energies, vectors, perturbation, occupied

    def response(self, model, scheme=-1, sigma=.04, alpha=.8):
        energies, vectors, perturbation, occupied = model

        def one_side(e0, e1, u0, u1, potential, n0, n1):
            raw = potential @ u0
            entries = [f"{scheme} {sigma} {alpha} {n0} {n1}",
                       " ".join(map(str, e0)), " ".join(map(str, e1))]
            rhs = self.run_probe("weight", entries+matrix_lines(u1)+matrix_lines(raw))
            denominator = e1[:, None] - e0[None, :]
            denominator[:n1, :] += alpha
            dpsi = u1 @ ((u1.conj().T @ rhs) / denominator)
            return u0.conj() @ dpsi.T

        direct = np.zeros((4, 4), complex)
        auxiliary = np.zeros_like(direct)
        for k in range(4):
            kq, mk, mkq = (k+1) % 4, (-k) % 4, (-k-1) % 4
            direct += one_side(energies[k], energies[kq], vectors[k], vectors[kq],
                               perturbation[k], occupied[k], occupied[kq]) / 4
            auxiliary += one_side(energies[mk], energies[mkq],
                                  J @ vectors[mk].conj(), J @ vectors[mkq].conj(),
                                  J @ perturbation[mkq].T @ J.T,
                                  occupied[mk], occupied[mkq]) / 4
        partner = self.run_probe("partner", matrix_lines(auxiliary))
        return direct+partner, direct, auxiliary

    def full_hamiltonians(self, model):
        energies, vectors, perturbation, _ = model
        h0 = np.zeros((16, 16), complex)
        dh = np.zeros_like(h0)
        for k in range(4):
            start = slice(4*k, 4*k+4)
            target = slice(4*((k+1) % 4), 4*((k+1) % 4)+4)
            h0[start, start] = (vectors[k]*energies[k]) @ vectors[k].conj().T
            dh[target, start] = perturbation[k]
            dh[start, target] = perturbation[k].conj().T
        return h0, dh

    def finite_difference(self, model, scheme=-1, sigma=.04, step=2e-6):
        h0, dh = self.full_hamiltonians(model)

        def density(h):
            e, u = np.linalg.eigh(h)
            return (u*occupation(e, sigma, scheme)) @ u.conj().T

        # The perturbation connects distinct crystal momenta, so its first-order
        # total electron-number change vanishes: no first-order chemical shift.
        derivative = (density(h0+step*dh)-density(h0-step*dh))/(2*step)
        result = np.zeros((4, 4), complex)
        for k in range(4):
            start = slice(4*k, 4*k+4)
            target = slice(4*((k+1) % 4), 4*((k+1) % 4)+4)
            result += derivative[target, start].T / 4
        return result

    def self_consistent_response(self, model, strength=.06, scheme=-1, sigma=.04):
        independent, _, _ = self.response(model, scheme, sigma)
        feedback = np.empty((16, 16), complex)
        for element in range(16):
            dns = np.zeros((4, 4), complex)
            dns.flat[element] = 1
            dv = self.run_probe("dudarev", [str(strength)]+matrix_lines(dns))
            response, _, _ = self.response(
                (model[0], model[1], np.repeat(dv[None, :, :], 4, axis=0), model[3]),
                scheme, sigma)
            feedback[:, element] = response.ravel()
        result = np.linalg.solve(np.eye(16)-feedback, independent.ravel()).reshape(4, 4)
        self.assertLess(np.linalg.norm(result.ravel()-independent.ravel()-feedback@result.ravel()), 1e-13)
        self.assertGreater(np.linalg.norm(result-independent)/np.linalg.norm(result), .05)
        return result

    def self_consistent_finite_difference(self, model, strength=.06, scheme=-1,
                                         sigma=.04, step=2e-4):
        h0, dh = self.full_hamiltonians(model)
        fourier = np.exp(2j*np.pi*np.outer(np.arange(4), np.arange(4))/4)/2
        transform = np.kron(fourier, np.eye(4))
        h0, dh = [transform @ h @ transform.conj().T for h in (h0, dh)]
        e0, u0 = np.linalg.eigh(h0)
        initial = (u0*occupation(e0, sigma, scheme)) @ u0.conj().T
        target_electrons = np.trace(initial).real
        blocks = [slice(4*l, 4*l+4) for l in range(4)]

        def interaction(local_density):
            # Independent derivative of U/2 Tr[rho-rho^2] in ordinary density
            # convention. It deliberately does not call the production kernel.
            result = np.zeros((16, 16), complex)
            for l, block in enumerate(blocks):
                result[block, block] = strength*(.5*np.eye(4)-local_density[l])
            return result

        reference = np.array([initial[block, block] for block in blocks])
        external = h0-interaction(reference)

        def scf(amplitude):
            local = reference.copy()
            for iteration in range(400):
                e, u = np.linalg.eigh(external+amplitude*dh+interaction(local))
                lower, upper = e.min()-2, e.max()+2
                for _ in range(70):
                    chemical = .5*(lower+upper)
                    if np.sum(occupation(e-chemical, sigma, scheme)) < target_electrons:
                        lower = chemical
                    else:
                        upper = chemical
                density = (u*occupation(e-chemical, sigma, scheme)) @ u.conj().T
                output = np.array([density[block, block] for block in blocks])
                error = np.max(np.abs(output-local))
                if error < 2e-13:
                    self.assertLess(abs(np.trace(density).real-target_electrons), 2e-12)
                    return output, iteration+1, float(error)
                local = .6*local+.4*output
            self.fail(f"finite-matrix SCF did not converge: {error}")

        plus, nplus, errplus = scf(step)
        minus, nminus, errminus = scf(-step)
        # Transform full complex site occupations N=rho^T to q=2*pi/4.
        derivative = (plus-minus).transpose(0, 2, 1)/(2*step)
        response = np.einsum("l,lab->ab", np.exp(-2j*np.pi*np.arange(4)/4), derivative)/4
        self.last_scf = {"step": step, "iterations_plus_minus": [nplus, nminus],
                         "maximum_density_residual_plus_minus": [errplus, errminus],
                         "target_electron_number": float(target_electrons)}
        return response

    def assert_response(self, actual, expected, tolerance=3e-8):
        self.assertLess(np.linalg.norm(actual-expected) / np.linalg.norm(expected), tolerance)

    def test_magnetic_complex_response_matches_independent_fd(self):
        model = self.model()
        for scheme in (-1, 0, -99):
            with self.subTest(smearing=scheme):
                response, direct, aux = self.response(model, scheme)
                expected = self.finite_difference(model, scheme)
                self.assert_response(response, expected)
                # Both errors remain visible in this complex magnetic model.
                self.assertGreater(np.linalg.norm(direct-expected)/np.linalg.norm(expected), .1)
                self.assertGreater(np.linalg.norm(direct+aux.T-expected)/np.linalg.norm(expected), .1)

    def test_occupied_subspace_shift_cancels(self):
        model = self.model()
        first, _, _ = self.response(model, alpha=.3)
        second, _, _ = self.response(model, alpha=1.7)
        self.assert_response(first, second, 2e-13)

    def test_direct_and_auxiliary_occupied_counts_can_differ(self):
        model = self.model(mixed_occupied=True)
        actual, _, _ = self.response(model)
        self.assert_response(actual, self.finite_difference(model))

    def test_near_degenerate_metal_limit_is_controlled(self):
        model = self.model(near_degenerate=True)
        actual, _, _ = self.response(model)
        self.assert_response(actual, self.finite_difference(model), 2e-7)

    def test_self_consistent_dudarev_feedback_matches_full_mean_field_fd(self):
        model = self.model()
        actual = self.self_consistent_response(model)
        finite = []
        for step in (2e-4, 1e-4):
            with self.subTest(step=step):
                expected = self.self_consistent_finite_difference(model, step=step)
                finite.append(expected)
                self.assert_response(actual, expected, 2e-7)
        self.last_scf = {"comparison": "Richardson extrapolation of the two converged central differences"}
        self.assert_response(actual, (4*finite[1]-finite[0])/3, 3e-8)


if __name__ == "__main__":
    unittest.main(verbosity=2)
