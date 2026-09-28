# PBE0

[日本語の実装・使用方法と完全な入力例](../pbe0-implementation-ja.md)

Use `xc='pbe0'` in `&functional`: 25% unscreened Fock exchange, 75% PBE exchange, and full PBE correlation. No rVV10 is added. PBEh(40) and HSE06 retain their existing definitions.

The PBEh Coulomb convention and controls are reused, including `pbeh_coulomb_radius`, `exx_mlwf_*`, and source-support ACE. A finite Coulomb radius is a numerical approximation, not part of the PBE0 definition. The radius is in the selected input length unit; zero selects the existing half-shortest-BvK-side limit.

```
&functional
 xc='pbe0'
 pbeh_coulomb_radius=4d0
/
```

DC-GS may use the PBE warm-up (`exx_pre_scf_threshold`) and unlocalized fragment exchange (`yn_exx_dc_mlwf='n'`). RT requires a compatible PBE0 seed: functional name, exchange fraction, Coulomb radius, and other saved functional parameters are checked. A PBEh(40) seed is not a PBE0 seed.

The 32 H2 comparison uses the same mesh, fragment buffers, .999 source support, cutoff 4 bohr, dt=.05 au and 7000 steps as PBEh(40). Each functional has its own converged DC GS and its own impulse/zero-field runs.
