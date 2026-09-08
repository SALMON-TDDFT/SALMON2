#!/bin/sh
set -eu

src="${1:-/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/rt/dg/rt_dg_density_reconstruct.f90}"

if [ ! -f "$src" ]; then
  echo "missing source: $src" >&2
  exit 1
fi

if ! rg -q "if \(dg_frag%is_frag_root\) then" "$src"; then
  echo "missing fragment-root guard in density reconstruction" >&2
  exit 1
fi

if ! perl -0ne 'exit 0 if /if \(dg_frag%is_frag_root\) then.*local_fragments_overlap_rank_box/s; exit 1' "$src"; then
  echo "send-side overlap setup is not fragment-root restricted" >&2
  exit 1
fi

if ! perl -0ne 'exit 0 if /if \(dg_frag%is_frag_root\) then.*do ifrag = dg_frag%ifrag_start, dg_frag%ifrag_end/s; exit 1' "$src"; then
  echo "fragment contribution loop is not fragment-root restricted" >&2
  exit 1
fi

echo "ok: density reconstruction is fragment-root restricted"
