#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
text = (ROOT / 'src/rt/dg/rt_dg_wpw_production_builder.f90').read_text().lower()
start = text.index('subroutine scan_and_route_dg_wpw_canonical_faces')
end = text.index('end subroutine scan_and_route_dg_wpw_canonical_faces', start)
body = text[start:end]
preflight = body.index('call mpi_allreduce(local_bad,global_bad')
assert preflight < body.index('allocate(staged'), \
    'batch scanner allocates from unvalidated rank-local extents before collective preflight'
assert preflight < body.index('call assemble_wpw_canonical_face_grid'), \
    'batch scanner slices unvalidated face arrays before collective preflight'
print('PASS canonical-face batch scanner validates collectively before slicing')
