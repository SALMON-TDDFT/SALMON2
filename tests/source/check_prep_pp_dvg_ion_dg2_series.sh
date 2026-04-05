#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
prep_pp_f90="$repo_root/src/atom/pp/prep_pp.f90"

grep -Fq "(g2sq*r1)**4/1680d0" "$prep_pp_f90"
! rg -F -q "(g2sq*r1)**4/2520d0" "$prep_pp_f90"
