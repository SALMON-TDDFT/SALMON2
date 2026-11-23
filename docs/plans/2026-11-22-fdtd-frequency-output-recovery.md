# Batch frequency-resolved plane-map I/O

Bounded R3 continuation. Native E/H currently rejects these diagnostics with15;
keep that scope unchanged rather than bypassing it for this test.

- Batch frequency-map header/row/close IOSTAT errors warn and continue subsequent
  windows/planes/energies and cleanup. Keep Fourier sums,MPI ordering,unit
  conversions and normal formatting. No new collective or global selector.
- Missing/partial maps are not fabricated or retried after finalization. Other
  output families,required restarts and staging allocation remain separate.
- Reuse --fdtd-spectrum-io --plane-map-faults:two energies,three planes,two windows,
  real directory/read-only faults,caller-unit retention and map-unit closure.
  Focused mode avoids redundant time snapshots;no new test module/CI physics job.

Record:`salmonpy-fdtd-frequency-output-recovery-20261122/summary.json` in canonical
calc2026. Source/wheel GNU Mac MPI1/2:all12 normal maps match prior bytes;8 faults
retain other11 maps and time series. Caller units600/601/602 and MPI program
completion pass. Both suites679 run,one optional skip,no failures. Source/object
runtime and host-probe audits,tracked-only build and sdist→wheel pass. No GS or
optional CI dispatch. Close-only faults were not injected;permissions measured
on a non-root Mac. Native spectral stepping/restart is explicitly not claimed.
