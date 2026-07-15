#!/usr/bin/env python3
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
operator=ROOT/'src/common/dg_wpw_matrix_free_operator.f90'
adapter=ROOT/'src/rt/dg/rt_dg_wpw_matrix_free_adapter.f90'
layout=ROOT/'src/rt/dg/rt_dg_wpw_column_layout.f90'
builder=ROOT/'src/rt/dg/rt_dg_wpw_sparse_builder.f90'
main=(ROOT/'src/gs/main_dft.f90').read_text().lower()
dense=(ROOT/'src/gs/dc/dg_wannier_pw_scf.f90').read_text().lower()

assert operator.exists(), 'missing matrix-free WPW operator contract'
assert layout.exists(), 'missing distributed PW-row layout'
assert adapter.exists(), 'missing RT block-operator adapter'
assert builder.exists(), 'production route is incomplete without direct windowed sparse builder'
op=operator.read_text().lower(); sub=adapter.read_text().lower(); dist=layout.read_text().lower()
for token in ('abstract interface','apply_h_batch','apply_s_batch','class(*), intent(inout) :: context',
              'xw_owned', 'xp_owned'):
    assert token in op, f'missing matrix-free operator token: {token}'
for token in ('apply_h_wpw_distributed', 'apply_s_wpw_distributed',
              'apply_matrix_blocks_batch_compact', 'apply_wp_owned_columns',
              'apply_pp_owned_rows', 'fetch_rows_from_owners',
              'reduce_w_partial_to_owners', 'rank-local', 'collective'):
    assert token in sub, f'missing block-adapter token: {token}'
assert 'apply_overlap_operator_batch' not in sub, \
    'production adapter must not use the global-replicated overlap apply'
for forbidden in ('prepare_local_fragment_pw_blocks',
                  'allocate(dg_frag%fp_local_pw_ids(n_pw))'):
    assert forbidden not in sub + dist, \
        f'production route uses transitional all-PW layout: {forbidden}'
for token in ('s_dg_wpw_column_layout', 'windowed_kg', 'owned_column_ids', 'wpw_column_owner'):
    assert token in dist, f'missing production (K,G) layout token: {token}'
for forbidden in ('allocate(h_global','allocate(s_global','complex(8)::h_global','complex(8)::s_global'):
    assert forbidden not in sub, f'production solver contains global dense allocation: {forbidden}'
assert 'run_dg_wpw_fixed_scf' not in main, 'dense oracle must never be called from main_dft'
assert 'small-system mathematical oracle' in dense
print('PASS matrix-free DG WPW SCF route contract')
