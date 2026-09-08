#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
source = (root / "src/parallel/communication.f90").read_text().lower()
body = source.split("subroutine comm_init", 1)[1].split("end subroutine", 1)[0]

assert "#ifdef use_eigenexa" in body
assert "mpi_thread_multiple" in body
assert "mpi_thread_funneled" in body
assert "iprovided < irequired" in body
assert "mpi_abort" in body
assert "eigenexa requires mpi_thread_multiple" in body

print("PASS EigenExa requests and validates MPI_THREAD_MULTIPLE")
