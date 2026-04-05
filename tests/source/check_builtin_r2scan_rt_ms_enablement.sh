#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)

xc_file="$repo_root/src/xc/salmon_xc.f90"
rt_init="$repo_root/src/rt/initialization_rt.f90"
rt_step="$repo_root/src/rt/time_evolution_step.f90"
gs_init="$repo_root/src/gs/initialization_dft.f90"
gs_scf="$repo_root/src/gs/scf_iteration.f90"
cmake_lists="$repo_root/tests/source/CMakeLists.txt"

if grep -Fq 'stop "r2SCAN supports only GS calculations"' "$xc_file"; then
  echo "unexpected GS-only r2scan guard still present" >&2
  exit 1
fi

grep -Fq 'system%xc_payload' "$xc_file"

rt_step_xc_calls=$(grep -F 'call exchange_correlation(' "$rt_step" | wc -l | tr -d ' ')
[ "$rt_step_xc_calls" -ge 2 ]

grep -Fq 'call exchange_correlation(' "$rt_init"
grep -Fq 'call exchange_correlation(' "$gs_init"
grep -Fq 'call exchange_correlation(' "$gs_scf"

if grep -Fq 'type(s_xc_operator_payload) :: xc_payload' "$rt_init"; then
  echo "unexpected temporary xc_payload in $rt_init" >&2
  exit 1
fi
if grep -Fq 'type(s_xc_operator_payload) :: xc_payload' "$rt_step"; then
  echo "unexpected temporary xc_payload in $rt_step" >&2
  exit 1
fi
if grep -Fq 'type(s_xc_operator_payload) :: xc_payload' "$gs_init"; then
  echo "unexpected temporary xc_payload in $gs_init" >&2
  exit 1
fi
if grep -Fq 'type(s_xc_operator_payload) :: xc_payload' "$gs_scf"; then
  echo "unexpected temporary xc_payload in $gs_scf" >&2
  exit 1
fi
if grep -Fq 'copy_xc_operator_payload' "$rt_init"; then
  echo "unexpected payload copy helper use in $rt_init" >&2
  exit 1
fi
if grep -Fq 'copy_xc_operator_payload' "$rt_step"; then
  echo "unexpected payload copy helper use in $rt_step" >&2
  exit 1
fi
if grep -Fq 'copy_xc_operator_payload' "$gs_init"; then
  echo "unexpected payload copy helper use in $gs_init" >&2
  exit 1
fi
if grep -Fq 'copy_xc_operator_payload' "$gs_scf"; then
  echo "unexpected payload copy helper use in $gs_scf" >&2
  exit 1
fi
if grep -Fq 'xc_payload=' "$rt_init"; then
  echo "unexpected explicit xc_payload argument in $rt_init" >&2
  exit 1
fi
if grep -Fq 'xc_payload=' "$rt_step"; then
  echo "unexpected explicit xc_payload argument in $rt_step" >&2
  exit 1
fi
if grep -Fq 'xc_payload=' "$gs_init"; then
  echo "unexpected explicit xc_payload argument in $gs_init" >&2
  exit 1
fi
if grep -Fq 'xc_payload=' "$gs_scf"; then
  echo "unexpected explicit xc_payload argument in $gs_scf" >&2
  exit 1
fi

grep -Fq 'check_builtin_r2scan_rt_ms_enablement.sh' "$cmake_lists"
