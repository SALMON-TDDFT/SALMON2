#!/bin/zsh
set -eu

repo_root=${0:A:h:h}
for target in \
  "$repo_root/src/rt/dg/rt_dg_fragment.f90" \
  "$repo_root/src/rt/dg/rt_dg_fragment_soi.f90"; do
  if rg -n "allocate\\(dg_frag%coef\\(dg_frag%n_mat_max,\\s*dg_frag%nstate_tot" "$target" >/dev/null; then
    print -u2 "DG coefficient storage still allocates full n_mat_max x nstate_tot in ${target:t}."
    exit 1
  fi
done

if ! rg -n "local_coef_count" "$repo_root/src/rt/dg/rt_dg_fragment_types.f90" >/dev/null; then
  print -u2 "DG fragment type does not expose local coefficient storage metadata."
  exit 1
fi

if rg -n "allocate\\(dg_frag%(H_mat|H_mat_c|H_mat_kinetic|S_mat|S_mat_c|S_mat_prop|S_mat_prop_c|momentum_mat|momentum_mat_c)\\(.*dg_frag%n_mat_max" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian.f90" \
  "$repo_root/src/rt/dg/rt_dg_fragment_basis_update.f90" \
  "$repo_root/src/rt/dg/rt_dg_plane_wave.f90" \
  "$repo_root/src/rt/dg/rt_dg_fragment_soi.f90" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90" \
  "$repo_root/src/rt/dg/rt_dg_fragment_basis_update_soi.f90" >/dev/null; then
  print -u2 "DG-SOI operator setup still performs full dense n_mat_max x n_mat_max allocation."
  exit 1
fi

if rg -n "allocate\\((H0c|M)\\(n_tot, n_tot\\)" "$repo_root/src/rt/dg/rt_dg_integrator_derivative.f90" >/dev/null; then
  print -u2 "DG derivative still has reachable full dense n_tot x n_tot work arrays."
  exit 1
fi

if ! rg -n "S_mat_blocks_c|S_mat_prop_blocks_c" "$repo_root/src/rt/dg/rt_dg_fragment_types.f90" >/dev/null; then
  print -u2 "SOI overlap has no complex block storage; imaginary Hermitian overlap would be lost."
  exit 1
fi

if ! rg -n "S_mat_blocks_c\\(iblk\\)%val.*cmplx\\(S_blocks_re" "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90" >/dev/null; then
  print -u2 "SOI overlap does not preserve the complex Hermitian overlap in block storage."
  exit 1
fi

if rg -n "call diagonalize_mixed_basis\\(" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian.f90" \
  "$repo_root/src/rt/dg/rt_dg_fragment_basis_update.f90" >/dev/null; then
  if ! rg -n "Mixed basis dense EVP skipped" "$repo_root/src/rt/dg/rt_dg_plane_wave.f90" >/dev/null; then
    print -u2 "Plane-wave mixed basis can still enter the dense EVP path."
    exit 1
  fi
  skip_line=$(rg -n "Mixed basis dense EVP skipped" "$repo_root/src/rt/dg/rt_dg_plane_wave.f90" | head -1 | cut -d: -f1)
  dense_line=$(rg -n "allocate\\(H_work\\(n_total, n_total\\), S_work\\(n_total, n_total\\)" "$repo_root/src/rt/dg/rt_dg_plane_wave.f90" | head -1 | cut -d: -f1)
  if [[ -n "$dense_line" && "$dense_line" -le "$skip_line" ]]; then
    print -u2 "Plane-wave mixed basis dense work arrays are allocated before the block/raw path guard."
    exit 1
  fi
fi

for target in \
  "$repo_root/src/rt/dg/rt_dg_integrator_rk.f90" \
  "$repo_root/src/rt/dg/rt_dg_integrator_aetrs.f90" \
  "$repo_root/src/rt/dg/rt_dg_integrator_derivative.f90"; do
  if rg -n "n\\s*=\\s*dg_frag%n_mat_max|n_frag\\s*=\\s*dg_frag%n_mat_max" "$target" >/dev/null; then
    print -u2 "DG coefficient propagation still uses global n_mat_max as the local coefficient extent in ${target:t}."
    exit 1
  fi
done

dg_return_line=$(rg -n "DG-RT initialization returns before conventional eigen-energy" "$repo_root/src/rt/initialization_rt.f90" | head -1 | cut -d: -f1)
eigen_line=$(rg -n "call calc_eigen_energy\\(energy,spsi_in" "$repo_root/src/rt/initialization_rt.f90" | head -1 | cut -d: -f1)
if [[ -z "$dg_return_line" || -z "$eigen_line" || "$dg_return_line" -ge "$eigen_line" ]]; then
  print -u2 "DG-RT initialization can still enter the conventional calc_eigen_energy path after DC read."
  exit 1
fi

print "DG local coefficient storage check passed."
