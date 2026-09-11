"""Fresh projector bound metadata must not read nps before calc_nps."""
from pathlib import Path
import re

root = Path(__file__).resolve().parents[2]
source = (root / "src/atom/pp/prep_pp.f90").read_text().lower()
initial = source.split("if (.not. allocated(ppg%jxyz_max)) then", 1)[1].split("end if", 1)[0]
assert not re.search(r"ppg%jxyz_min\s*=\s*ppg%nps", initial), (
    "fresh projector bounds read uninitialized nps and corrupt DC cache ownership fingerprints"
)
assert re.search(r"ppg%jxyz_min\s*=\s*0\b", initial)
assert re.search(r"ppg%jxyz_max\s*=\s*0\b", initial)
main = (root / "src/gs/main_dft.f90").read_text().lower()
assert "call hash_alloc_integer_rank2(hash,grid%jxyz_min)" in main
assert "call hash_alloc_integer_rank2(hash,grid%jxyz_max)" in main
print("deterministic initial projector-bound metadata: PASS")
