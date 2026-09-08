# Optical Vortex In `single_scale_maxwell_tddft`

## Scope

This note summarizes the current optical-vortex implementation added to SALMON for `single_scale_maxwell_tddft`.

## Positioning Of The Present Model

This implementation should be regarded as a simplified optical-vortex model for single-scale boundary injection, not as a strict Maxwell-consistent Laguerre-Gaussian vector potential.

- The OAM phase factor `exp(i l phi)` is included.
- A transverse vector potential with specified polarization is included.
- The radial envelope is not the standard Laguerre-Gaussian envelope. It is a `sin^2` cutoff window.
- `A_z` is not included.
- `w(z)`, `R(z)`, and Gouy phase are not included.

So the present implementation is best interpreted as:

- an approximate transverse vortex vector potential,
- with finite-radius apodization,
- intended for practical single-scale injection.

Standard OAM/LG beam background:

- Allen et al., Phys. Rev. A 45, 8185 (1992)
  - [PubMed](https://pubmed.ncbi.nlm.nih.gov/9906912/)

For paraxial vector-potential beam theory and longitudinal corrections:

- Lax, Louisell, McKnight, Phys. Rev. A 11, 1365 (1975)
  - [APS](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.11.1365)
- Davis, Phys. Rev. A 19, 1177 (1979)
  - [APS](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.19.1177)

Relevant source files:

- Input variables: [src/io/salmon_global.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/io/salmon_global.f90)
- Input parsing/checks: [src/io/inputoutput.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/io/inputoutput.f90)
- Optical-vortex generator: [src/maxwell/optical_vortex_field.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/maxwell/optical_vortex_field.f90)
- Boundary injection: [src/maxwell/fdtd_coulomb_gauge.f90](/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/src/maxwell/fdtd_coulomb_gauge.f90)

## Required Input Parameters

### `&singlescale`

- `method_singlescale = '3d'`
  - Optical vortex is supported only for `3d`.
- `yn_optical_vortex = 'y'`
  - Enables optical-vortex field generation.
- `optical_vortex_charge`
  - Topological charge `l`.
  - `0` means no OAM phase.
- `optical_vortex_polarization`
  - Polarization state.
  - Allowed values: `linear_x`, `linear_y`, `left_circular`, `right_circular`
- `optical_vortex_radius`
  - Vortex radius `R`.
  - The field is set to zero for `r >= R`.
- `optical_vortex_center_x`
  - Focus center `x`.
  - If left as a large negative sentinel, box center is used.
- `optical_vortex_center_y`
  - Focus center `y`.
  - If left as a large negative sentinel, box center is used.

### `&emfield`

The temporal envelope and field strength reuse the existing `&emfield` variables.

- `ae_shape1`
  - Currently supported: `Acos2`, `Acos3`, `Acos4`, `Acos6`, `Acos8`
- `omega1`
  - Carrier angular frequency
- `tw1`
  - Pulse width
- `phi_cep1`
  - CEP
- `t1_start`
  - Pulse start offset
- `I_wcm2_1`
  - Intensity
- `E_amplitude1`
  - Direct field amplitude used when `I_wcm2_1 < 0`

## Physical Meaning

- Propagation direction: `z`
- Transverse plane: `x-y`
- Generated components: `A_x(x,y,t)`, `A_y(x,y,t)`
- `A_z = 0`
- Retarded-time shift: `t -> t - z/c`
- Radial envelope: `sin^2` inside radius `R`
- Azimuthal phase: `exp(i l phi)`

## Actual Formula Used In The Code

The implementation itself is in [src/maxwell/optical_vortex_field.f90]

### 1. Field amplitude

When `I_wcm2_1 >= 0`:

```text
f0 = 5.338d-9 * sqrt(I_wcm2_1)
```

When `I_wcm2_1 < 0`:

```text
f0 = E_amplitude1
```

### 2. Beam center

```text
xc = optical_vortex_center_x
yc = optical_vortex_center_y
```

If left unspecified:

```text
xc = 0.5 * Nx * hx
yc = 0.5 * Ny * hy
```

### 3. Radius and azimuth

```text
r = sqrt((x - xc)^2 + (y - yc)^2)
phi = atan2(y - yc, x - xc)
```

### 4. Temporal envelope

```text
tt = t - tw1/2 - t1_start
```

For `ae_shape1 = AcosN`:

```text
env_t = cos(pi * tt / tw1)^N
```

If

```text
|tt| >= tw1/2
```

then `A = 0`.

### 5. Spatial envelope

```text
window(r) = sin[(pi/2) * (1 - r/R)]^2
```

If

```text
r >= R
```

then `A = 0`.

### 6. Radial factor

```text
radial_factor =
  1                  (l = 0)
  (r / R)^|l|        (l /= 0)
```

### 7. Polarization vectors

`linear_x`

```text
zpol = (1, 0, 0)
```

`linear_y`

```text
zpol = (0, 1, 0)
```

`left_circular`

```text
zpol = (1/sqrt(2), i/sqrt(2), 0)
```

`right_circular`

```text
zpol = (1/sqrt(2), -i/sqrt(2), 0)
```

### 8. Phase

```text
phase = omega1 * tt + 2*pi*phi_cep1 + l * phi
```

### 9. Final vector potential

The code uses

```text
A(x,y,t)
= -(f0 / omega1)
   * env_t
   * window(r)
   * radial_factor
   * Im[zpol * exp(i * phase)]
```

In component form:

```text
A_ext(1:3)
= -(f0 / omega1)
   * cos(pi * tt / tw1)^N
   * sin[(pi/2)(1-r/R)]^2
   * (r/R)^|l|
   * Im[zpol * exp(i phase)]
```

For `l = 0`, `(r/R)^|l| = 1`.

## Simple Interpretation

- `l = 0`, `left_circular/right_circular`
  - Circularly polarized beam with a finite-radius `sin^2` aperture
- `l != 0`
  - Approximate optical vortex with OAM phase `exp(i l phi)`
- `A_z`
  - Not included
- Gauge side
  - The field is injected directly as vector potential in the single-scale boundary update

## Current Limitations

- Only for `single_scale_maxwell_tddft`
- Only for `method_singlescale='3d'`
- `ae_shape1` limited to `Acos2/3/4/6/8`
- `A_z = 0`
- This is not a strict exact LG Maxwell solution
- Implemented field is an approximate combination of
  - `AcosN` temporal envelope
  - `sin^2` radial window
  - `exp(i l phi)` azimuthal phase
  - chosen polarization vector

## Difference From A Standard LG Beam

Compared with a standard Laguerre-Gaussian beam description:

- the radial envelope is replaced by `sin^2` cutoff,
- Gouy phase is omitted,
- curvature radius `R(z)` is omitted,
- beam-waist evolution `w(z)` is omitted,
- longitudinal vector-potential correction is omitted.

Accordingly, this implementation should be described as a simplified single-scale vortex-field model, not as a full LG vector-potential solution.

## Example

```fortran
&emfield
  ae_shape1 = 'Acos2'
  omega1 = 0.05696d0
  tw1 = 330.0d0
  phi_cep1 = 0.0d0
  t1_start = 0.0d0
  I_wcm2_1 = 1.0d10
/

&singlescale
  method_singlescale = '3d'
  yn_optical_vortex = 'y'
  optical_vortex_charge = 1
  optical_vortex_polarization = 'left_circular'
  optical_vortex_radius = 20.0d0
  optical_vortex_center_x = -1.0d30
  optical_vortex_center_y = -1.0d30
/
```
