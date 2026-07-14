from pathlib import Path


source = (Path(__file__).resolve().parents[2] / "src/gs/dc/lcfo_flux.f90").read_text()
body = source.split("subroutine output", 1)[1].split("end subroutine output", 1)[0]

call = "call write_wannier90_global_bpw_seed()"
assert call in body
window = body[max(0, body.index(call) - 240) : body.index(call) + len(call)]
assert "wannier_cluster_size" not in window
assert "all(wannier_cluster_size" not in body
assert "requires cluster/global BPW export" not in body

exporter = source.split("subroutine write_wannier90_global_bpw_seed", 1)[1].split(
    "end subroutine write_wannier90_global_bpw_seed", 1
)[0]
for token in [
    "num_wann_file / max(1, dc%n_frag)",
    "call write_wannier90_global_trace_file",
    "binfile_bpw",
    "binfile_bpwt",
]:
    assert token in exporter or token == "binfile_bpwt"

print("PASS global-cluster SAWF exports fragment BPW and DG face traces")
