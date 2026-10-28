#!/usr/bin/env python3
"""Compile the actual XML-species copy routine and exercise optional defaults.

The tested subroutine and its two derived types are extracted verbatim from
production files. Only the kinds/module scaffolding is stubbed. The same
driver must fail against the retained pre-fix routine and pass after the fix.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def extract(text, pattern):
    match = re.search(pattern, text, re.I | re.M | re.S)
    if not match:
        raise RuntimeError("Could not extract " + pattern)
    return match.group() + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--compiler", default="gfortran")
    args = parser.parse_args()
    source = ROOT / "Modules/qexsd_copy.f90"
    types_source = ROOT / "Modules/qes_types_module.f90"
    baseline = HERE / "evidence/qexsd-copy-magnetic-baseline.f90"
    driver = HERE / "qexsd_copy_magnetic_driver.f90"
    routine = extract(source.read_text(),
                      r"^  SUBROUTINE qexsd_copy_atomic_species\b.*?^  END SUBROUTINE qexsd_copy_atomic_species[^\n]*$")
    types = "".join(extract(types_source.read_text(),
                           r"^  TYPE :: " + name + r"\b.*?^  END TYPE " + name + r"[^\n]*$")
                    for name in ("species_type", "atomic_species_type"))
    stub = ("MODULE kinds\n  IMPLICIT NONE\n  INTEGER, PARAMETER :: dp=SELECTED_REAL_KIND(14,200)\n"
            "END MODULE kinds\nMODULE qes_types_module\n  USE kinds, ONLY: dp\n  IMPLICIT NONE\n" +
            types + "END MODULE qes_types_module\n")
    wrapper = "MODULE actual_copy\n  USE kinds, ONLY: dp\n  IMPLICIT NONE\nCONTAINS\n"
    runs = {}
    with tempfile.TemporaryDirectory(prefix="qe-qexsd-copy-") as temporary:
        work = Path(temporary)
        (work / "stubs.f90").write_text(stub)
        for name, text in (("baseline", baseline.read_text()), ("production", routine)):
            tested = wrapper + text + "END MODULE actual_copy\n"
            path = work / (name + ".f90")
            path.write_text(tested)
            exe = work / (name + ".x")
            command = [args.compiler, "-O0", "-g", "-fcheck=all", "-Wall", "-ffree-line-length-none",
                       "-J", str(work), str(work / "stubs.f90"), str(path), str(driver), "-o", str(exe)]
            built = subprocess.run(command, cwd=work, text=True, capture_output=True, check=True)
            process = subprocess.run([str(exe)], cwd=work, text=True, capture_output=True, check=False)
            cases = {line.split(maxsplit=1)[1]: line.startswith("PASS ")
                     for line in process.stdout.splitlines() if line.startswith(("PASS ", "FAIL "))}
            runs[name] = {"returncode": process.returncode, "cases": cases,
                          "stdout": process.stdout, "stderr": process.stderr,
                          "compile_stderr": built.stderr, "compile_command": command,
                          "extracted_routine_sha256": sha(text.encode()),
                          "compiled_source_sha256": sha(tested.encode()),
                          "binary_sha256": sha(exe.read_bytes())}
    success = (runs["production"]["returncode"] == 0 and len(runs["production"]["cases"]) == 8 and
               all(runs["production"]["cases"].values()) and runs["baseline"]["returncode"] != 0 and
               runs["baseline"]["cases"].get("second_read_missing_fields") is False)
    result = {"pass": success, "scope": __doc__, "runs": runs,
              "source_sha256": {str(path): sha(path.read_bytes())
                                for path in (source, types_source, baseline, driver, Path(__file__))},
              "compiler_version": subprocess.check_output([args.compiler, "--version"], text=True).splitlines()[0],
              "types_sha256": sha(types.encode()), "stub_sha256": sha(stub.encode())}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2) + "\n")
    args.report.with_suffix(".production.f90").write_text(routine)
    args.report.with_suffix(".types.f90").write_text(stub)
    print(json.dumps({"pass": success, "baseline": runs["baseline"]["cases"],
                      "production": runs["production"]["cases"]}, indent=2))
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
