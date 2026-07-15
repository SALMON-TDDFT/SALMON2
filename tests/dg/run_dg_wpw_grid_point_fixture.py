#!/usr/bin/env python3
from pathlib import Path
import shlex
import subprocess

ROOT=Path(__file__).resolve().parents[2]
BUILD=ROOT/"build-mpi-eigenexa"
OBJ=Path("/tmp/test_dg_wpw_grid_point_adapter.o")
EXE=Path("/tmp/test_dg_wpw_grid_point_adapter")
subprocess.run(["/opt/homebrew/bin/mpifort","-I",str(BUILD),"-c",
                str(ROOT/"tests/dg/test_dg_wpw_grid_point_adapter.f90"),"-o",str(OBJ)],check=True)
link_dir=BUILD/"src"
args=shlex.split((link_dir/"CMakeFiles/salmon.dir/link.txt").read_text())
args.insert(1,str(OBJ)); args.remove("CMakeFiles/salmon.dir/main.f90.o")
args[args.index("-o")+1]=str(EXE)
subprocess.run(args,cwd=link_dir,check=True)
subprocess.run([str(EXE)],check=True)
