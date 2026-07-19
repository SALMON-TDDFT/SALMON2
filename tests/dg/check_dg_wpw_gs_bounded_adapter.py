#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'src/common/dg_wpw_bounded_operator.f90'
EXCHANGE = ROOT / 'src/common/dg_wpw_owner_exchange.f90'

assert SOURCE.exists(), 'missing GS/RT-neutral bounded WPW operator adapter'
assert EXCHANGE.exists(), 'owner-targeted WPW exchange has not been moved to the common layer'

source = SOURCE.read_text().lower()
exchange = EXCHANGE.read_text().lower()
for token in (
    's_dg_wpw_bounded_operator',
    'initialize_dg_wpw_bounded_operator',
    'apply_h_dg_wpw_bounded',
    'apply_s_dg_wpw_bounded',
    'global_gram_dg_wpw_bounded',
    'owned_w_ids',
    'owned_p_ids',
    'required_w_ids',
    'required_p_ids',
    'operator_epoch',
    'layout_fingerprint',
    'metric_convention',
    'operator_convention',
    'release_dg_wpw_bounded_operator',
):
    assert token in source, f'missing bounded operator contract: {token}'

for forbidden in (
    'use rt_',
    's_dg_fragment_rt',
    'allocate(h_global',
    'allocate(s_global',
    'global_owner',
    'x_global',
):
    assert forbidden not in source + exchange, f'non-bounded/common-layer dependency detected: {forbidden}'

assert 'fetch_rows_from_owners' in exchange
assert 'reduce_w_partial_to_owners' in exchange
assert 'release_dg_wpw_owner_schedule' in exchange and 'mpi_comm_free' in exchange
print('PASS GS-neutral bounded WPW adapter source contract')
