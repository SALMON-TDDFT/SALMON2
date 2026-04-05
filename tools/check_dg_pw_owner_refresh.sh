#!/bin/sh
set -eu

src="${1:-/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_plane_wave.f90}"

if [ ! -f "$src" ]; then
  echo "missing source: $src" >&2
  exit 1
fi

if ! perl -0ne 'exit 0 if /subroutine compact_plane_wave_basis.*allocated\(dg_frag%coef_pw_owner\).*dg_frag%owned_coef_pw_start = 0.*dg_frag%owned_coef_pw_end = -1.*allocate\(dg_frag%coef_pw_owner\(dg_frag%n_plane_waves\)\).*dg_frag%coef_pw_owner\(i\) = min\(dg_frag%isize - 1, \(\(i - 1\) \* dg_frag%isize\) \/ dg_frag%n_plane_waves\).*end subroutine compact_plane_wave_basis/s; exit 1' "$src"; then
  echo "compact_plane_wave_basis does not rebuild PW owner metadata" >&2
  exit 1
fi

echo "ok: compact_plane_wave_basis rebuilds PW owner metadata"
