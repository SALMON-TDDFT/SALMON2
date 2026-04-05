#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
script="$repo_root/tests/tools/run_primitive_nonorth_stress_fd.sh"
cmake_lists="$repo_root/tests/source/CMakeLists.txt"

grep -Fq 'check_primitive_nonorth_stress_fd_tool.sh' "$cmake_lists"
test -x "$script"
grep -Fq 'primitive nonorthogonal stress finite-difference helper' "$script"
grep -Fq -- '--mode MODE' "$script"
grep -Fq 'scale' "$script"
grep -Fq 'shear' "$script"
grep -Fq "stress_output_level = 'high'" "$script"
grep -Fq 'nproc_rgrid(1) = 1' "$script"
grep -Fq 'nproc_rgrid(2) = 1' "$script"
grep -Fq 'nproc_rgrid(3) = 1' "$script"
