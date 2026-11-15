# Native field sampling follow-up (F4b)

Code audit after explicit packet transport; this document is not implementation.

Current evidence from coupled.py:

- NativeEHField exposes sampled E at n*dt and integrated A at (n+1/2)*dt
  through its electronic-clock convention. sp_eh_step returns E only. No Python
  point-sampled H accessor exists. Native saved H is centered at E time; it must
  not be mistaken for the raw Yee H endpoint.
- MultiscaleField exposes A at n*dt and its stored E at (n-1)*dt.
- ResponseSession packet mode currently accepts j at integer local n*dt only.
  The legacy driver shifts electronic timing by the field adapter's offset.
- CoupledStepper currently calls step(A); reject packet sessions at construction
  before any field advance until an explicit packet scheduler exists.

Next implementation boundaries:

1. Explicit sample records for available A/E, preserving per-component absolute
   timestamps and physical sample positions. Never advertise absent H as zeros.
2. An explicit user-selected temporal/spatial sampling policy. Require history
   and a stencil covering the requested time/point; no extrapolation or silent
   half-step retagging. Native legacy A update order must remain unchanged.
3. Resolve local versus field clock mapping jointly with return-current timing;
   simply packaging E/A into a FieldPacket is insufficient to close the loop.
4. Add native H sampling only through a documented native accessor, with Yee
   staggering/time metadata and native/reference tests. Preserve runtime source
   inventories and provenance if native code changes.
5. Connect an explicit packet coupling driver; keep current-density semantics.
   Isolated dipole differentiation, macro volume, and self-field choice remain
   explicit separate operations, not native-source assumptions.

Acceptance: old A-only driver parity, manufactured time/space interpolation,
missing-stencil/port/clock rejection before material updates, one-/two-field-rank
sample equality, and bounded local MPI tests. No new physical convergence sweep.
