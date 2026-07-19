#!/usr/bin/env python3
from pathlib import Path
import re
ROOT=Path(__file__).resolve().parents[2]
text=(ROOT/'src/gs/dc/lcfo_flux.f90').read_text().lower()
assert re.search(r"initialize_dg_wpw_volume_accumulator\s*\([^\n]*\n?[^\)]*size\(owned_w\)\s*,\s*nps\s*,\s*nps",text), \
  'piecewise DG volume quadrature must stage owned W rows and every support P row'
assert 'wpw_nonowned_candidates_pending' in text, 'missing explicit non-owned candidate staging state'
assert 'if(volume_info==0)wpw_nonowned_candidates_pending=.true.' in text and \
       'if(candidate_info==0)wpw_nonowned_candidates_pending=.false.' in text, \
  'non-owned candidate staging is not held until successful routing'
publish=text[text.index('subroutine publish_wpw_core_volume_candidates'):
             text.index('end subroutine publish_wpw_core_volume_candidates')]
assert 'nw=size(wpw_owned_w_ids)' in publish.replace(' ', ''), \
  'WP volume publication must enumerate fragment-owned W rows only'
assert 'wp_w(iwp)=wpw_owned_w_ids(iw)' in publish.replace(' ', ''), \
  'WP volume candidates must use the stable ID of each owned W row'
assert 'wp_p(iwp)=wpw_context%support_column_ids(ip)' in publish, \
  'WP volume candidates do not preserve the stable ID of each support P row'
assert 'pp_r(ipp)=wpw_context%support_column_ids(ip)' in publish, \
  'PP volume candidates do not preserve the stable ID of each support P row'
print('PASS non-owned P candidates remain staged until R5 routing')
