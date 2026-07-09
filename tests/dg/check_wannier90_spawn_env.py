#!/usr/bin/env python3
from pathlib import Path


root = Path(__file__).resolve().parents[2]
src = root / "src" / "gs" / "dc" / "lcfo_flux.f90"
text = src.read_text()

required = [
    "mpi_clean_env_prefix()",
    "-u OMPI_COMM_WORLD_SIZE",
    "-u OMPI_COMM_WORLD_RANK",
    "-u PMIX_NAMESPACE",
    "-u PMIX_RANK",
]

missing = [item for item in required if item not in text]
if missing:
    raise SystemExit("missing W90 spawn MPI-env cleanup markers: " + ", ".join(missing))

if "command_line = trim(change_dir_command)//trim(mpi_clean_env_prefix())//' '//trim(wannier_command)" not in text:
    raise SystemExit("Wannier90 spawn command does not use mpi_clean_env_prefix()")

print("Wannier90 spawn clears inherited MPI/PMIx environment")
