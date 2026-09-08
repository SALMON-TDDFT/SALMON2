from pathlib import Path


source = (Path(__file__).resolve().parents[2] / "src/gs/dc/lcfo_flux.f90").read_text()
body = source.split("subroutine write_wannier90_global_trace_file", 1)[1].split(
    "end subroutine write_wannier90_global_trace_file", 1
)[0]

assert "sawf_explicit_basis_active" in body
assert "sawf_explicit_buffer" in body
assert "closed SAWF buffer basis count mismatch" in body
assert body.index("sawf_explicit_basis_active") < body.index("binfile_bfb")
assert "matmul" not in body.split("sawf_explicit_basis_active", 1)[1].split("else", 1)[0]
assert "call lcfo_sawf_rank_local_fatal" in body
assert " stop " not in " ".join(body.lower().split())
assert "iostat=" in body.lower()
assert "iomsg=" in body.lower()
normalized = "".join(body.split())
assert "any(nxyz_buffer_seed(1:3)<size(stencil%coef_nab,1))" in normalized
assert body.index("finite-difference stencil radius") < body.index("do face=1,6")
assert "close(iunit,iostat=io_status,iomsg=io_message)" in normalized
assert "cannot close output" in body

print("PASS global BPW traces use the symmetry-closed in-memory buffer basis")
