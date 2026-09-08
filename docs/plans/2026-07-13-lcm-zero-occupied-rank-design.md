# LCM zero-occupied-rank fix design

## Problem

The Local Chern Marker work arrays retain a dummy orbital column through
`max(1,nocc_local)`.  The occupied-orbital copy routines infer their loop bound
from that allocated extent, so a rank with no occupied orbitals copies one
dummy orbital and may read a wavefunction orbital it does not own.

## Design

Keep the existing dummy-column allocation because downstream BLAS and
communication code relies on nonzero extents.  Pass the actual
`nocc_local` value into each occupied-orbital copy routine and use it as the
only copy-loop bound.  Apply the same contract to the scalar/collinear and SOI
spinor implementations.

## Verification

Add a source-level regression check that establishes both call sites pass the
actual local count and both copy routines use that explicit count instead of
an allocated dummy extent.  Verify the check fails before the production
change and passes afterward, then rebuild the MPI/EigenExa configuration.

