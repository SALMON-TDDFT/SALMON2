# One-way molecular dipole radiation implementation plan

Goal: distribute physical molecular charge dipoles into explicit FDTD macro cells
and feed finite-difference polarization current to existing NativeEHField.
Approved direction: one-way first; no local-field/self-field feedback or automatic
vacuum-box normalization. No new native numerical code.

1. Add dipole_source.py: immutable DipoleMap with molecule/point IDs, nonnegative
   column-normalized weights and positive explicit volumes, all atomic units.
   polarization(p)=W p/V; current=(Pnew-Pold)/dt. Physical charge dipoles only.
   Test conservation with unequal volumes, static offsets, superposition, strict
   shapes/finite values, immutable ownership and dt/weight rejection first.
2. Add drive_native_dipoles: collective field-group only, fresh NativeEHField,
   validate dt/IDs/volumes/capacity and finite trajectory before stepping. Dipoles
   at n*dt yield current at (n+1/2)*dt. First interval must be constant because
   native starts with zero current. Caller explicitly duplicates p0 to delay a
   recorded molecular waveform by dt. No silent interpolation/padding.
3. Run native one-way regression with a zero-amplitude plane wave and vacuum grid.
   Compare first nonzero injection against native -4*pi*dt*S*S^T j update;
   check a remote observation later becomes nonzero. Compare field-rank1/2 if
   available within four-core budget. Use explicit existing audited library.
   Feed recorded C2H2 total charge dipoles as a separate end-to-end smoke test.
4. Public API inventory, explicit wheel file list, paired guide and runnable
   lightweight example; AI JSON synchronization. Full source and installed wheel
   tests, meaningful native checks. No GitHub Actions, commit or push.
