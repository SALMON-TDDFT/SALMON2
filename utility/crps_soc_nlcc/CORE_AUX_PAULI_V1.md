# V2 coherent auxiliary-Pauli NLCC model

This is an explicitly specified frozen-core **auxiliary Pauli ensemble** built
from both blocks of the generator's occupied atomic Dirac orbitals. It is not a
Foldy–Wouthuysen reduction and its tau is not the Dirac kinetic-energy density.
Neither a density-only substitute nor von Weizsaecker tau is used.

`write_core_aux=.true.` is the default and writes all V2 core information.
The consuming functional selects its required quantities: LSDA uses rho and
noncollinear MGGA uses rho/tau/current together. Explicit `.false.` is reserved
for historical V1 reproduction. No existing V1 file is rewritten. The new canonical result is
`/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/ni-soc-nlcc-20261024/v2-aux-pauli-validated/Ni/Ni_soc_nlcc_rps.dat`.

## Same-shell construction

Core orbitals are occupied atomic indices with `AtoP==0`. Each supplies two
angular shells, with amplitudes s=F/r and s=G/r, the original occupation q,
and labels:

- Large block: l=lA, lambda=−(kappa+1).
- Small block: l=lA−sign(kappa), lambda=kappa−1.

Each angular magnetic sublevel has occupation q/(2j+1). Atomic F/G normalization
is integral(F²+G²)dr=1. At the partial-core radius each angular amplitude is
replaced inside by `(r/rc)^l [a+b(r/rc)^2+c(r/rc)^4]`, matching its value and
first two derivatives. Atomic derivatives follow the radial Dirac equations;
a cubic Hermite interpolant supplies values/derivatives between native points.
This gives finite origin limits, including small-component angular shells.

All regular seeds are scaled together by w=sqrt(rhoPC/n_seed), preserving every
legacy density sample. Derivatives include the derivative of w. The density
derivative is the same symmetric five-point rule used by the existing writer
(with its even-origin and endpoint rules); primitive derivative convergence is
therefore verified numerically rather than presumed exact between samples.
Normalized seed ratios avoid squaring underflow in tiny densities.

For each final radial amplitude s and derivative s', full angular sums give:

```
rho = sum q*s²/(4*pi)
tau = sum q*[s'² + l*(l+1)*(s/r)²]/(4*pi)       # no factor 1/2
C   = sum q*lambda*s²/(8*pi*r)
J_i^a = C*epsilon_iab*rhat_b
```

Reference-gauge magnetization, spin-vector tau, and charge current vanish;
spin current does not generally vanish. At r=0, s/r uses s' for l=1 and zero
for higher l; l=0 has zero centrifugal term and lambda=0. C(0)=0.

Ni supplies ten blocks from five occupied core shells, atomic core occupancy12.
The smoothed partial-core density integrates to **2.047762334288 electrons**,
not12: the historical NLCC smoothing does not conserve the atomic core charge.
The metadata occupancy is an atomic configuration label, not a normalization
constraint on each auxiliary block after common scaling.

## Explicit far-boundary continuation

The native atom's finite radial support ends before all positive lobes of the
legacy density interpolant. The first missing positive-density seed occurs at
19.11 Å; the last nonzero native seed anchor is18.69 Å. All affected density
samples are retained. Relative angular amplitudes continue from that anchor as
bound-Dirac asymptotic exp(−decay*r)/r, with derivatives from this same continued
function, before the common density scale is applied.

This path requires r>10 Å and integrated remaining charge<1e−60 electrons,
otherwise generation stops. The selected Ni affected charge is
4.45052501672791e−69 electrons. This also satisfies the tighter1e−68-electron
charge bound. Integrated tau over the affected region is6.86445e−68 a0^-2;
integrated absolute spin-current coefficient is1.34329e−74 a0^-1. This is a finite-atomic-domain/interpolation
artifact, not described as IEEE density underflow. Original failed attempts
remain in sibling `v2-aux-pauli-20261025*` directories. No density cutoff or
support truncation is introduced. The existing atomic solver still reports
its ordinary IEEE underflow flag; all final output entries are finite.

## Exact format

The normal table uses marker `KY_SOC_NLCC_V2`; all numeric header/channel/radial
rows retain the V1 format and units. After the radial rows, before
`END_KY_SOC_NLCC_V2`, the block is:

```
CORE_AUX_PAULI_V1
PAULI_CLOSED_SHELL_V1 ATOMIC_UNITS TAU_NO_HALF
nrows nblocks atomic_core_occupancy
orbital_id l lambda occupation                  # one row per block
r_bohr rho_au tau_au C_au (s_au,sprime_au)_blocks # one row per radial sample
```

The block uses a0^-3 density, a0^-5 tau, a0^-4 C, a0^-3/2 amplitude,
a0^-5/2 amplitude derivative. Conversion uses a0=0.529177210903 Å.
The old SALMON constant differs slightly, so cross-unit comparisons must
account for that documented difference. All rows through20 Å are supplied.
These are reference-gauge primitives/fields; the noncollinear MGGA caller must
apply its consistent electromagnetic/spin-frame transformation to frozen core
and valence together.

## Verification

`test_core_aux.py` independently evaluates spin spherical harmonics using
Clebsch–Gordan coefficients and numerical Cartesian derivatives for l0/1/2,
both allowed j branches. It checks density, tau, spin-current tensor sign and
normalization, zero magnetization and zero charge current. It also validates
the actual Ni V2 file against its V1 control: identical numeric potentials,
waves, density and density derivatives; finite fields; same-primitive closure;
tau>=|grad(rho)|²/(4*rho), up to roundoff; complete endpoint support.

Primitive fourth-order derivative relative max residuals at stride4/2/1 are
0.01798 /0.002285 /0.0003125, respectively. The largest residual is near the
partial-core matching radius, with clear grid convergence. The selected Ni
maxima are tau=1.81027041655 a0^-5 and |C|=0.00110835356506 a0^-4. The minimum
tau−tauVW is −6.91e−77, roundoff in the far tail.

New-generator V1 controls reproduce the selected Ni, W and Se files byte for
byte. Source and all three linked application objects passed the PZHEEV
prohibition checker before linking. Source/executable/table hashes, compiler,
source snapshots and independent test output are in the canonical final
`manifest.json`. No material GS data was collected; cloud synchronization is
not confirmed. The previous radial ghost/transferability limitations still
apply because V2 preserves the ionic potentials and projectors unchanged.
