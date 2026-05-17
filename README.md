# SALMON2-SBE — SALMON fork with SBE dephasing

This repository is a fork of the original SALMON project, an open-source
software package for ab-initio quantum-mechanical calculations of light-matter
interactions.

For the upstream SALMON project, documentation, and publications, see:

http://salmon-tddft.jp/

## Fork features

This fork adds a phenomenological decoherence/dephasing time to the
Semiconductor Bloch Equations (SBE) solver. The implementation introduces the
`&sbe` namelist parameter `t2_sbe_fs`, which damps off-diagonal density-matrix
elements during SBE time propagation.

### `t2_sbe_fs` SBE dephasing parameter

`&sbe` now accepts:

```fortran
&sbe
  ! ... existing SBE parameters ...
  t2_sbe_fs = 50.0d0
/
```

Parameter details:

- **Name:** `t2_sbe_fs`
- **Namelist:** `&sbe`
- **Units:** femtoseconds (`fs`), regardless of the global `unit_system`
- **Default:** `1.0d10`, which effectively disables dephasing
- **Active range:** positive values below `1.0d9` enable the dephasing term
- **Physical effect:** off-diagonal density-matrix elements decay with the
  characteristic time `T2`

The SBE Bloch solver converts `t2_sbe_fs` from femtoseconds to atomic units and
adds the T2 term to the same fourth-order Taylor propagation operator used for
the Bloch equation. For the intended SBE regime, for example `dt = 0.001 fs` and
`T2 = 10 fs`, `dt / T2 = 1.0e-4`, so the Taylor approximation to the coherence
damping is effectively the exponential decay

```text
rho_ij -> rho_ij * exp(-dt / T2),  i != j
```

This leaves diagonal occupations unchanged and applies the relaxation only to
coherences. Keep `dt` much smaller than `T2`; this is already satisfied by the
typical SBE setup above.

### Input propagation

The new parameter is:

- declared globally as `t2_sbe_fs`;
- read from the `&sbe` namelist;
- initialized to a no-dephasing default;
- broadcast to all MPI ranks before SBE propagation.

## Minimal SBE input example

```fortran
&calculation
  theory = 'sbe'
/

&sbe
  num_sbe = 1
  ! Other required SBE system parameters must be supplied as usual.
  ! Set a finite T2 to enable dephasing:
  t2_sbe_fs = 50.0d0
/
```

Set `t2_sbe_fs` back to the default value, or any value greater than or equal to
`1.0d9`, to recover the original no-dephasing behavior.

## License

SALMON is available under Apache License version 2.0.

    Copyright 2017-2025 SALMON developers

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
