#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src/gs/dc/dg_wpw_lda_hartree.f90"

assert SRC.exists(), "missing dg_wpw_lda_hartree module"
text = SRC.read_text().lower()
structures = (ROOT / "src/common/structures.f90").read_text().lower()
dcdft = (ROOT / "src/gs/dc/dcdft.f90").read_text().lower()
salmon_xc = (ROOT / "src/xc/salmon_xc.f90").read_text().lower()

for token in (
    "module dg_wpw_lda_hartree",
    "subroutine validate_core_ownership",
    "subroutine integrate_core_lda_terms",
    "subroutine hartree_energy_global",
    "subroutine update_wpw_lda_hartree",
    "subroutine update_wpw_owned_lda_hartree",
    "use hartree_sub, only: hartree",
    "use salmon_xc, only: exchange_correlation",
    "call hartree(",
    "call exchange_correlation(",
):
    assert token in text, f"missing contract token: {token}"

assert "srg_tot" in structures, "DC total-system state does not retain the XC send/recv grid"
assert "dc%srg_tot" in dcdft and "finalize_sendrecv_grid_storage(dc%srg_tot)" in dcdft, \
    "DC total-system XC communication state lacks initialization/finalization ownership"
assert "any(shape(eexc_tmp) /= mg%num)" in salmon_xc and "deallocate(eexc_tmp)" in salmon_xc, \
    "XC scratch must be reallocated when fragment and total-grid shapes differ"

# Two fragments see the full grid, but each point has exactly one core owner.
core = np.array([
    [True, True, True, False, False, False],
    [False, False, False, True, True, True],
])
assert np.all(core.sum(axis=0) == 1)

rho = np.array([0.4, 0.7, 1.1, 0.8, 0.5, 0.3])
vxc = -0.6 * np.cbrt(rho)
exc = -0.45 * rho ** (4.0 / 3.0)
weights = np.array([0.2, 0.3, 0.25, 0.35, 0.15, 0.4])

# Buffer/halo values are deliberately enormous; the mask must eliminate them.
rho_frag = np.tile(rho, (2, 1))
vxc_frag = np.tile(vxc, (2, 1))
exc_frag = np.tile(exc, (2, 1))
rho_frag[~core] = 1.0e90
vxc_frag[~core] = -1.0e90
exc_frag[~core] = 1.0e90

exc_fragment = np.sum(np.where(core, exc_frag * weights, 0.0))
nvxc_fragment = np.sum(np.where(core, rho_frag * vxc_frag * weights, 0.0))
assert np.isclose(exc_fragment, np.sum(exc * weights))
assert np.isclose(nvxc_fragment, np.sum(rho * vxc * weights))

vh = np.array([0.2, 0.5, 0.7, 0.6, 0.3, 0.1])
eh = 0.5 * np.sum(rho * vh * weights)
assert np.isclose(eh, 0.5 * np.dot(weights, rho * vh))

# Negative ownership controls: one hole and one double owner must fail.
hole = core.copy(); hole[:, 2] = False
double = core.copy(); double[:, 2] = True
assert np.any(hole.sum(axis=0) != 1)
assert np.any(double.sum(axis=0) != 1)

print("PASS DG Wannier+PW fragment-core LDA and global Hartree contract")

fc = shutil.which("gfortran")
assert fc is not None, "gfortran is required"
fixture = ROOT / "tests/dg/fixtures/dg_wpw_lda_hartree_driver.F90"
with tempfile.TemporaryDirectory(prefix="dg-wpw-lda-") as td:
    td = Path(td)
    exe = Path(td) / "check"
    fixture_text = fixture.read_text()
    split_at = fixture_text.index("program dg_wpw_lda_hartree_driver")
    stubs = td / "stubs.F90"; driver = td / "driver.F90"
    stubs.write_text(fixture_text[:split_at])
    driver.write_text(fixture_text[split_at:])
    result = subprocess.run(
        [fc, str(stubs), str(SRC), str(driver), "-o", str(exe)],
        cwd=td, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    assert result.returncode == 0, result.stdout
    result = subprocess.run([str(exe)], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert result.returncode == 0, result.stdout
    assert "PASS DG WPW LDA Hartree Fortran driver" in result.stdout
