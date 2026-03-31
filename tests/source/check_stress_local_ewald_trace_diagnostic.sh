#!/bin/sh
set -eu

repo_root=$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)
write_f90="$repo_root/src/io/write.f90"

if ! rg -F -q 'Local/Ewald trace diagnostics [Hartree]' "$write_f90"; then
  echo "missing local/ewald trace diagnostics heading" >&2
  exit 1
fi

if ! rg -F -q 'virial_loc_lr_residual = stress_tensor_trace(system%stress_loc_lr_grad + system%stress_loc_lr_diag) &' "$write_f90"; then
  echo "missing local LR virial residual definition" >&2
  exit 1
fi

if ! rg -F -q '* system%det_a + system%stress_loc_lr_energy' "$write_f90"; then
  echo "missing flipped local LR virial sign" >&2
  exit 1
fi

for needle in \
  'Tr(loc)*V' \
  'Tr(loc_grad)*V' \
  'Tr(loc_diag)*V' \
  'Tr(loc_sr_grad)*V' \
  'Tr(loc_lr_grad)*V' \
  'Tr(loc_lr)*V + E_lr' \
  'Tr(ewa)*V' \
  'Tr(ewa_G_grad)*V' \
  'Tr(ewa_G_diag+self)*V' \
  'Tr(ewa_R)*V' \
  'E_ion_ion - (E_ewa_G+E_ewa_R)' \
  'Tr(loc+ewa)*V'
do
  if ! rg -F -q "$needle" "$write_f90"; then
    echo "missing local/ewald trace diagnostic line: $needle" >&2
    exit 1
  fi
done
