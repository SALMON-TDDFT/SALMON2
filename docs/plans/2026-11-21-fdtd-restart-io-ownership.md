# Common Maxwell restart/helper file-unit ownership

### Legacy EM restart/helper file-unit ownership

The common_maxwell text/binary restart readers and writers now open their own
local newunit, including per-rank data and metadata files and the remapping read
path. txtfile_copy_em independently owns its input and output units instead of
borrowing700/701. Existing helper signatures retain ifn for source compatibility;
it no longer selects an actual file unit in these helpers. Paths, formats,
metadata, array order, copy trimming and communication/reconstruction operations
are unchanged. No additional startup or per-step collective is introduced.

The existing I/O fixture has an explicit --fdtd-restart-io mode: synthetic EM
arrays only, no SCF/GS data. GNU MPI1/2 tests preserve host units600/601/700/701
while exercising text values, real/complex1–5D binary arrays, single/all modes,
same-layout and cyclically exchanged target-domain reads, inactive payloads and
text copying. Control/occupied outputs and prior-runtime outputs match byte for
byte. This tests the helpers and does not certify every complete legacy batch
restart configuration or arbitrary rank counts. Normal LD+FDTD/direct-resume
field, point, plane, shape and checkpoint data also match prior controls.

Linear-response/spectral/ASE/ART and TTM output units are still separate work.
This change isolates file-unit ownership; missing/corrupt restart data, IO device
failure and general allocation failure retain their existing handling. No general
IO-error recovery, Fujitsu execution or multi-node certification is claimed.
