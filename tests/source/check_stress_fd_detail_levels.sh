#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
global_f90="$repo_root/src/io/salmon_global.f90"
input_f90="$repo_root/src/io/inputoutput.f90"
write_f90="$repo_root/src/io/write.f90"

grep -Fq "character(6)   :: stress_fd_detail     = 'high'" "$global_f90"
grep -Fq "stress_fd_detail    = 'high'" "$input_f90"
grep -Fq "stress_fd_detail == 'middle' .or. stress_fd_detail == 'high'" "$write_f90"
grep -Fq "if(stress_fd_detail == 'high') then" "$write_f90"
grep -Fq "stress_fd_detail must be 'low', 'middle', or 'high' (aliases: 'A', 'B', 'C')" "$input_f90"
grep -Fq "case('A')" "$input_f90"
grep -Fq "case('B')" "$input_f90"
grep -Fq "case('C')" "$input_f90"
grep -Fq "case('low','middle','high')" "$input_f90"
