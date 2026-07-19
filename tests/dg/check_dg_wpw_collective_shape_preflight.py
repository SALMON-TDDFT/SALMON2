#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

def body(path, name):
    text = (ROOT / path).read_text().lower()
    start = text.index(f'subroutine {name}')
    end = text.index(f'end subroutine {name}', start)
    return text[start:end]

side = body('src/rt/dg/rt_dg_wpw_face_side_halo.f90', 'reduce_dg_wpw_face_side_parts')
trace = body('src/rt/dg/rt_dg_wpw_trace_halo_provider.f90', 'reduce_dg_wpw_face_trace_parts')
for name, routine in [('face side', side), ('face trace', trace)]:
    assert 'reference_shape' in routine, f'{name} reduction has no fixed-size dimension handshake'
    shape_bcast = routine.index('call mpi_bcast(reference_shape')
    preflight = routine.index('call mpi_allreduce(local_bad,global_bad', shape_bcast)
    payload_bcast = routine.index('call mpi_bcast(reference_grid')
    assert shape_bcast < preflight < payload_bcast, \
        f'{name} reduction communicates rank-dependent payload counts before shape preflight'
print('PASS face reductions validate fixed-size dimensions before payload collectives')
