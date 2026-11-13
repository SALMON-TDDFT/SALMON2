# Public SALMON VR desktop bridge validation

2026-11-13, macOS arm64. No SALMON equations or output update order changed.

Public artifact: SALMON-TDDFT/SALMON-VR main, salmon-vr-mac.tar.
SHA256 (matches published Git LFS pointer):
`516da7e747398cfbf00cd9406ea03469479f6324d90ebce87640bdbef1cf7058`.
Official Mac manual: https://salmon-tddft.jp/download/macos_manual.pdf
Publication: https://doi.org/10.1007/s10055-026-01319-0

The public archive's executable lacked an execute permission. User execute
permission was added to the temporary extracted executable for this test.
No quarantine, signature or OS security settings were changed.

Exported four density-difference CUBEs from the existing TB-mBJ output validation,
native steps 0, 8, 16, 24 and dt=0.05 atomic time units. Generated matched all-zero
Ac_x CUBEs explicitly marked as display placeholders. Configuration and application
were placed side-by-side as required by the public manual. The UI displayed the
placeholder note and dt=0.0012094422 fs/native step. After a forward frame action,
it displayed `Current Time Step: 000001 (Interpolated)` with positive/negative
density isosurfaces and atoms. No Exception/Error/NaN/Infinity was found in
Player.log during this check. This verifies this small dataset and app build,
not arbitrary data, large-memory operation, Windows, recording or live streaming.

Unit checks cover explicit zero selection, values/grid preservation, physical
step metadata, mismatched grids, duplicate steps, refusal to overwrite, and
optional app copying without launch. The full Python suite passed 77 tests against the repaired, installed Mac wheel.
