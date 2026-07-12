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

spawn = text.split("subroutine run_wannier90_seed_files", 1)[1].split(
    "end subroutine run_wannier90_seed_files", 1
)[0]
if "trim(mpi_clean_env_prefix())" not in spawn or "trim(resolved_command)" not in spawn:
    raise SystemExit("Wannier90 spawn command does not use mpi_clean_env_prefix() and the resolved command")
if "get_environment_variable" in spawn:
    raise SystemExit("Wannier90 spawn re-resolves the command instead of using the resolved value")
resolver = text.split("logical function dc_lcfo_wannier_import_only_requested", 1)[1].split(
    "end function dc_lcfo_wannier_import_only_requested", 1
)[0]
environment = resolver.index("get_environment_variable('SALMON_WANNIER90_COMMAND'")
selection = resolver.index("select_wannier90_command")
cache = resolver.index("cache_resolved_wannier90_command")
if not environment < selection < cache:
    raise SystemExit("Wannier90 command resolver does not preserve namelist -> environment -> compiled precedence")
if text.count("get_environment_variable('SALMON_WANNIER90_COMMAND'") != 1:
    raise SystemExit("Wannier90 command is resolved more than once")
if "stop " in spawn.lower() or "call comm_bcast" not in spawn or "call lcfo_sawf_fatal" not in spawn:
    raise SystemExit("Wannier90 external-command failure is not collective-safe")

print("Wannier90 spawn preserves command precedence and clears inherited MPI/PMIx environment")
