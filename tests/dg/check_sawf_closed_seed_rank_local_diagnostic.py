from pathlib import Path

root = Path(__file__).resolve().parents[2]
source = (root / "src/gs/dc/lcfo_flux.f90").read_text().lower()

start = source.index("call build_sawf_closed_fragment_seed_basis(", source.index("subroutine prepare_sawf_closed_seed_eigensystem"))
stop = source.index("if(failure/=0) call lcfo_sawf_fatal('sawf closed-seed basis construction failed')", start)
block = source[start:stop]

diagnostic = "[fatal] sawf closed-seed basis rank="
assert diagnostic in block, "missing rank-local closed-seed construction diagnostic"
assert block.index(diagnostic) < block.index("call comm_get_max(failure,dc%icomm_tot)"), (
    "closed-seed failure must be diagnosed on the failing rank before collective reduction"
)
assert "trim(message)" in block[block.index(diagnostic):], "rank-local diagnostic must include the detailed reason"

print("PASS closed-seed failures are diagnosed locally before reduction")
