#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
text = (ROOT / 'src/gs/dc/lcfo_flux.f90').read_text().lower()

def body(name):
    start = text.index(f'subroutine {name}')
    end = text.index(f'end subroutine {name}', start)
    return text[start:end]

volume_halo = body('prepare_wpw_volume_halo')
quadrature = body('assemble_wpw_core_volume_quadrature')
face_trace = body('prepare_wpw_canonical_face_trace_provider')

assert 'call mpi_allreduce(local_failure,global_failure' in volume_halo, \
    'volume halo enters root exchange without a fragment-root collective preflight'
assert 'if(pack_info/=0)return' not in volume_halo, \
    'volume halo can return locally before the fragment-root exchange'
assert quadrature.count('call mpi_allreduce(local_failure,global_failure') >= 2, \
    'volume quadrature lacks collective preflight before fragment reduction'
assert 'if(point_info/=0)return' not in quadrature, \
    'volume quadrature can return locally before fragment reduction'
assert 'if(point_info/=0)then;coverage(point)=0;cycle;endif' in face_trace and \
       'if(point_info/=0.or.any(column_ids/=p_ids))then;coverage(point)=0;cycle;endif' in face_trace, \
    'face trace evaluation can return locally before fragment reduction'
print('PASS WPW preparation failures follow a collective fail-closed schedule')
