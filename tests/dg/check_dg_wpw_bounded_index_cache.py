#!/usr/bin/env python3
from pathlib import Path

source = (Path(__file__).resolve().parents[2] / 'src/common/dg_wpw_bounded_operator.f90').read_text().lower()
apply_start = source.index('subroutine apply_blocks')
apply_end = source.index('end subroutine', apply_start)
apply_body = source[apply_start:apply_end]
init_start = source.index('subroutine initialize_dg_wpw_bounded_operator')
init_end = source.index('end subroutine', init_start)
init_body = source[init_start:init_end]

assert 'find_id(' not in apply_body, 'bounded H/S apply repeats stable-ID linear searches'
assert 'find_id(' not in init_body, 'bounded operator initialization must not linearly scan IDs per sparse entry'
assert 'find_id_sorted(' in init_body, 'bounded operator initialization must use logarithmic endpoint lookup'
assert 'strictly_increasing(owned_w)' in init_body and 'strictly_increasing(required_p)' in init_body, \
    'binary endpoint lookup requires collective validation of all stable-ID namespaces'
for name in ('ww_ri', 'ww_ci', 'wp_wi', 'wp_pi', 'pp_ri', 'pp_ci'):
    assert name in source, f'missing bounded operator endpoint cache: {name}'
    assert f'candidate%{name}' in source, f'endpoint cache is not initialized transactionally: {name}'
for name in ('ww_h_dense', 'ww_s_dense', 'wp_h_dense', 'wp_s_dense', 'pp_h_dense', 'pp_s_dense'):
    assert name in source, f'missing rank-local bounded dense cache: {name}'
    assert f'candidate%{name}' in init_body, f'dense cache is not initialized transactionally: {name}'
for name in ('ww_dense', 'wp_dense', 'pp_dense'):
    assert f'matmul({name}' in apply_body.replace(' ', ''), f'bounded apply does not batch {name} over vectors'
for name in ('ww', 'wp', 'pp'):
    assert f'do i=1,size({name})' not in apply_body.replace(' ', ''), \
        f'bounded apply retains scalar sparse-entry loop for {name}'

print('PASS bounded WPW H/S apply uses epoch-local endpoint index caches')
