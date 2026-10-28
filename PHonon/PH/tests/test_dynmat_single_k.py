#!/usr/bin/env python3
"""Compile production dynmat_us and check the magnetic single-k cache boundary.

For a single plane wave the expectation of a translated separable nonlocal
projector is independent of its position: its phase cancels between bra and
ket. The exact force constant is zero, including a complex spin potential.
Leaving the auxiliary time-reversed spinor in the shared cache must not change
that result. The legacy condition is compiled separately as a negative control.
"""
from pathlib import Path
import os
import re
import shutil
import subprocess
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "dynmat_us.f90"


class SingleKNonlocalCurvature(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which(os.environ.get("FC", "gfortran"))
        if compiler is None:
            raise unittest.SkipTest("Fortran compiler unavailable")
        cls.temporary = tempfile.TemporaryDirectory(prefix="qe-single-k-")
        cls.addClassCleanup(cls.temporary.cleanup)
        build = Path(cls.temporary.name)
        text = SOURCE.read_text()
        legacy, count = re.subn(
            r"IF\s*\(nksq\s*>\s*1\s*\.OR\.\s*\(noncolin\s*\.AND\.\s*domag\)\)",
            "IF (nksq > 1)", text, count=1, flags=re.IGNORECASE)
        if count != 1:
            raise AssertionError("expected the magnetic direct-spinor reload guard")
        legacy_source = build / "legacy_dynmat_us.f90"
        legacy_source.write_text(legacy)
        cls.executables = []
        for source, name in [(SOURCE, "current"), (legacy_source, "legacy")]:
            binary = build / name
            command = [compiler, "-O0", "-fcheck=all", "-ffree-line-length-none",
                       str(HERE / "dynmat_single_k_stubs.f90"), str(source),
                       "-o", str(binary)]
            subprocess.run(command, cwd=build, check=True, capture_output=True, text=True)
            cls.executables.append(binary)

    def run_probe(self, executable, scenario):
        result = subprocess.run([str(executable), str(scenario)], check=True,
                                capture_output=True, text=True)
        reads, curvature = result.stdout.split()
        return int(reads), float(curvature)

    def test_single_k_reloads_direct_spinor(self):
        reads, curvature = self.run_probe(self.executables[0], 1)
        self.assertEqual(reads, 1)
        self.assertLess(curvature, 2e-16)

    def test_legacy_cache_bug_is_detected(self):
        reads, curvature = self.run_probe(self.executables[1], 1)
        self.assertEqual(reads, 0)
        self.assertGreater(curvature, 1e-2)

    def test_multiple_k_and_nonmagnetic_cache_remain_valid(self):
        for executable in self.executables:
            for scenario, expected_reads in [(2, 2), (3, 0)]:
                with self.subTest(executable=executable.name, scenario=scenario):
                    reads, curvature = self.run_probe(executable, scenario)
                    self.assertEqual(reads, expected_reads)
                    self.assertLess(curvature, 2e-16)


if __name__ == "__main__":
    unittest.main(verbosity=2)
