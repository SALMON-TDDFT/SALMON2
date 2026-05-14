#!/bin/zsh
set -eu

repo_root=${0:A:h:h}

for target in \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian.f90" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90"; do
  if rg -n "local_self_blocks_only\\s*=\\s*dg_frag%parallel_mode_orbital" "$target" >/dev/null; then
    print -u2 "DG Galerkin volume blocks are still coupled implicitly to orbital mode in ${target:t}."
    exit 1
  fi

  if ! rg -n "subroutine init_matrix_blocks\\(dg_frag, blocks, block_map, n_blocks, diagonal_only\\)" "$target" >/dev/null; then
    print -u2 "DG matrix block setup does not expose an explicit local-vs-boundary policy in ${target:t}."
    exit 1
  fi

  if ! rg -n "subroutine init_momentum_blocks\\(dg_frag, diagonal_only\\)" "$target" >/dev/null; then
    print -u2 "DG momentum block setup does not expose an explicit local-vs-boundary policy in ${target:t}."
    exit 1
  fi
done

if ! rg -n "diagonal_blocks_only \\.and\\. ifrag_col /= dg_frag%ifrag_group" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian.f90" >/dev/null; then
  print -u2 "Orbital H blocks must be row-local while keeping neighbor columns for boundary/flow coupling."
  exit 1
fi

if ! rg -n "subroutine add_dg_surface_flux_blocks" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian.f90" >/dev/null; then
  print -u2 "DG Hamiltonian is missing the explicit Eq. (4) boundary/flow surface term."
  exit 1
fi

if ! rg -n "subroutine add_dg_surface_momentum_blocks" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian.f90" >/dev/null; then
  print -u2 "DG current/momentum is missing the explicit boundary normal-flow term."
  exit 1
fi

if ! rg -n "basis_functions_buffer_soi\\.bin" "$repo_root/src/gs/dc/lcfo_soi.f90" >/dev/null; then
  print -u2 "DC SOI must export a complex buffer basis for DG surface flux stencils."
  exit 1
fi

if ! rg -n "basis_functions_buffer_soi\\.bin" "$repo_root/src/rt/dg/rt_dg_fragment_soi.f90" >/dev/null; then
  print -u2 "RT SOI must read the complex buffer basis, not only core basis_functions_soi.bin."
  exit 1
fi

if ! rg -n "phi_frag_spinor_c" "$repo_root/src/rt/dg/rt_dg_fragment_types.f90" >/dev/null || \
   ! rg -n "phi_frag_spinor_c" "$repo_root/src/rt/dg/rt_dg_fragment_soi.f90" >/dev/null; then
  print -u2 "RT SOI must keep the two-component spinor basis, not collapse it to phi_frag_c component 1."
  exit 1
fi

if rg -n "using spin-1 basis only|spin-1 basis only" "$repo_root/src/rt/dg/rt_dg_fragment_soi.f90" >/dev/null; then
  print -u2 "RT SOI still contains a spinor-component discard path."
  exit 1
fi

if ! rg -n "subroutine add_dg_surface_flux_blocks_soi" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90" >/dev/null; then
  print -u2 "SOI Hamiltonian is missing the explicit Eq. (4) boundary/flow surface term."
  exit 1
fi

if ! rg -n "subroutine add_dg_surface_momentum_blocks_soi" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90" >/dev/null; then
  print -u2 "SOI current/momentum is missing the explicit boundary normal-flow term."
  exit 1
fi

if ! rg -n "momentum_blocks_c" "$repo_root/src/rt/dg/rt_dg_fragment_types.f90" >/dev/null || \
   ! rg -n "momentum_blocks_c" "$repo_root/src/rt/dg/rt_dg_fragment_ops.f90" >/dev/null; then
  print -u2 "SOI current/momentum must preserve complex spinor momentum blocks through apply paths."
  exit 1
fi

if rg -n "call reduce_complex_momentum_blocks" "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90" >/dev/null; then
  print -u2 "SOI momentum must not use the global scratch/reduce path in orbital mode."
  exit 1
fi

if ! rg -n "n_same == 2 \\.and\\. n_adjacent == 1" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian.f90" >/dev/null || \
   ! rg -n "n_same == 2 \\.and\\. n_adjacent == 1" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90" >/dev/null; then
  print -u2 "DG momentum boundary blocks must be face-neighbor only, not edge/corner neighbors."
  exit 1
fi

if ! rg -n "if \\(\\.not\\. dg_frag%parallel_mode_orbital\\) return" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian.f90" >/dev/null; then
  print -u2 "DG surface flux must not silently use the orbital-only face evaluator in legacy real-space mode."
  exit 1
fi

if ! rg -n "call init_matrix_blocks\\(dg_frag, dg_frag%H_mat_blocks, dg_frag%H_block_map, dg_frag%n_H_blocks,.*" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian.f90" >/dev/null || \
   ! rg -n "diagonal_only=\\.false\\.\\)" "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian.f90" >/dev/null; then
  print -u2 "DG Hamiltonian must allocate sparse neighbor blocks only for the explicit boundary/flow term."
  exit 1
fi

if rg -n "call init_matrix_blocks\\(dg_frag, dg_frag%S_mat_blocks, dg_frag%S_block_map, dg_frag%n_S_blocks,.*diagonal_only=\\.false\\.\\)" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian.f90" >/dev/null; then
  print -u2 "DG overlap S must remain block diagonal."
  exit 1
fi

if rg -n "call init_momentum_blocks\\(dg_frag\\)" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian.f90" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90" >/dev/null; then
  print -u2 "DG momentum setup must explicitly request local volume blocks."
  exit 1
fi

if ! rg -n "call init_momentum_blocks\\(dg_frag, diagonal_only=\\(\\.not\\. dg_frag%parallel_mode_orbital\\)\\)" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian.f90" >/dev/null || \
   ! rg -n "call init_momentum_blocks\\(dg_frag, diagonal_only=\\(\\.not\\. dg_frag%parallel_mode_orbital\\)\\)" \
  "$repo_root/src/rt/dg/rt_dg_fragment_hamiltonian_soi.f90" >/dev/null; then
  print -u2 "DG momentum must allocate face-neighbor blocks only in orbital boundary-flow mode."
  exit 1
fi

print "DG Galerkin block policy check passed."
