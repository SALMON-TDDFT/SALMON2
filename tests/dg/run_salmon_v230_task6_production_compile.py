#!/usr/bin/env python3
"""Compile staged Task 5/6 integration in an isolated production overlay."""

from argparse import ArgumentParser
from hashlib import sha256
from pathlib import Path
import os
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
TASK7_STAGES = {
    # Official stage 2 has the pre-DG calc_vlocal_fragment_dcdft interface and
    # cannot compile existing dcdft_soi callers.  Stage 3 is clean and matches
    # those already integrated DG interfaces; this substitution is overlay-only.
    "src/gs/dc/dcdft.f90": "3",
    "src/gs/main_dft.f90": "2",
    # Stage 3 calls the already integrated DG solve/potential/energy interfaces;
    # official stage 2 cannot type-check against those modules.
    "src/gs/scf_iteration_dft.f90": "3",
    "src/io/inputoutput.f90": "2",
    "src/rt/main_tddft.f90": "2",
}
W90_ARCHIVE = Path(
    "/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/"
    "build-mpi-eigenexa/wannier90/src/v3.1.0.tar.gz"
)


def source_digest(source: Path) -> tuple[int, str]:
    digest = sha256()
    count = 0
    for item in sorted(source.rglob("*")):
        relative = item.relative_to(source).as_posix().encode()
        if item.is_symlink():
            payload = os.readlink(item).encode()
            kind = b"L"
        elif item.is_file():
            payload = item.read_bytes()
            kind = b"F"
        else:
            continue
        digest.update(kind + b"\0" + relative + b"\0" + payload + b"\0")
        count += 1
    return count, digest.hexdigest()


parser = ArgumentParser()
parser.add_argument("--jobs", type=int, default=4)
parser.add_argument("--keep", action="store_true")
args = parser.parse_args()

holder = tempfile.TemporaryDirectory(prefix="salmon-task6-production-")
workspace = Path(holder.name)
source = workspace / "source"
build = workspace / "build"
shutil.copytree(ROOT, source, ignore=shutil.ignore_patterns(".git"), symlinks=True)

print(f"repository={ROOT}")
print(f"HEAD={subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()}")
print(f"MERGE_HEAD={subprocess.check_output(['git','rev-parse','MERGE_HEAD'],cwd=ROOT,text=True).strip()}")
print("overlay_recipe=copy current worktree excluding .git; replace only Task7 paths with documented clean index stages")
for relative, selected_stage in TASK7_STAGES.items():
    stage = subprocess.check_output(["git", "ls-files", "-u", "--", relative], cwd=ROOT, text=True)
    blob = next(line.split()[1] for line in stage.splitlines() if line.split()[2] == selected_stage)
    (source / relative).write_bytes(subprocess.check_output(["git", "cat-file", "blob", blob], cwd=ROOT))
    print(f"task7_substitution={relative} stage={selected_stage} blob={blob}")

# The DG stage-3 dcdft carries the required fragment-potential API but predates
# v2.3 init_dft's srg/unfold outputs. Compose only that mechanical upstream API
# at the Task7 overlay boundary; the real conflicted file remains untouched.
dcdft_overlay = source / "src/gs/dc/dcdft.f90"
dcdft_text = dcdft_overlay.read_text()
dcdft_text = dcdft_text.replace(
    "      type(s_stencil) :: stencil_dummy\n      type(s_ofile) :: ofile_dummy\n",
    "      type(s_stencil) :: stencil_dummy\n      type(s_sendrecv_grid) :: srg_dummy\n"
    "      type(s_ofile) :: ofile_dummy\n      type(s_unfold) :: unfold_dummy\n",
    1,
)
dcdft_text = dcdft_text.replace(
    "      & stencil_dummy,dc%fg_tot,dc%poisson_tot,dc%srg_tot,dc%srg_scalar_tot,ofile_dummy)\n"
    "      deallocate(dc%system_tot%rocc)\n",
    "      & stencil_dummy,dc%fg_tot,dc%poisson_tot,srg_dummy,dc%srg_scalar_tot,ofile_dummy,unfold_dummy)\n"
    "      deallocate(dc%system_tot%rocc)\n      call dealloc_cache(srg_dummy)\n",
    1,
)
dcdft_overlay.write_text(dcdft_text)
print(f"task7_overlay_api_composition=src/gs/dc/dcdft.f90 sha256={sha256(dcdft_text.encode()).hexdigest()}")

file_count, tree_digest = source_digest(source)
print(f"overlay_source_files={file_count}")
print(f"overlay_source_sha256={tree_digest}")
print(f"overlay_workspace={workspace}")

configure = [
    "cmake", "-S", str(source), "-B", str(build),
    "-DCMAKE_BUILD_TYPE=Debug",
    "-DUSE_MPI=ON", "-DUSE_SCALAPACK=ON", "-DUSE_EIGENEXA=OFF",
    "-DUSE_WANNIER90=ON", "-DWANNIER90_COMMS=serial",
]
print("configure:", " ".join(configure), flush=True)
subprocess.run(configure, check=True)
if not W90_ARCHIVE.is_file():
    raise SystemExit(f"cached Wannier90 archive is missing: {W90_ARCHIVE}")
w90_destination = build / "wannier90/src/v3.1.0.tar.gz"
w90_destination.parent.mkdir(parents=True, exist_ok=True)
shutil.copy2(W90_ARCHIVE, w90_destination)
print(f"wannier90_archive={W90_ARCHIVE}")
build_command = ["cmake", "--build", str(build), "--target", "salmon", "-j", str(args.jobs)]
print("build:", " ".join(build_command), flush=True)
try:
    subprocess.run(build_command, check=True)
finally:
    if args.keep:
        print(f"kept_overlay={workspace}")
        holder.cleanup = lambda: None
