# Fast-PR HSE symmetry-reduced k / spatial SCF restart

Extend the explicit hse_kspace bridge to compute2. Retain existing full-q CPU
exchange and representative weights, one/two k ranks and no orbital split.
Validate local Cartesian boxes against global owned cell/spacing/point count;
do not bypass the full-grid identity check. Recursive exchange refresh uses
the existing full-grid symmetry expansion and local representative sources.

Migration remains mixed-repartition3: accepted numerical rows/k slices and SCF
mixture move directly; only rank-specific geometric phase caches are retired
as in pure-k migration. No density/exchange/SCF replay. Check both changed-k
and changed-space layouts, pre-PBE transition, exact cut interiors, existing
1e-9/1e-10 continuation, native late failure and protected different-input owners,
previous full-k/Gamma and reduced-k results. Source/object audits, unit suite,
source-to-wheel/public sample, JP/EN/AI/list. No routine CI. Packed/tiled SCF
formats remain separate: current built-in SCF generates dense ACE.

Accepted: source2↔4,2→1,4→8,changed-axis4→4,pre-PBE2→4; physical cut
interiors exact, bootstrap0, continuation within unchanged1e-9/1e-10. Protected
different-input owners and failure/continue/retry exact. Previous reduced-k,
full-k, Gamma mixed and ordinary calculation results exact. Source/installed637
tests each636pass/1skip; wheel4→2 and public reduced-k sample4→2 pass. Source/
object audits and source/wheel/docs/header/native byte parity verified.
Campaign-specific source/wheel GS export explicitly authorized. Evidence:
calc2026/salmonpy-hse-symspace-record-20261120. No CI dispatched.
