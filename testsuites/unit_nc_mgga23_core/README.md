# Auxiliary-Pauli frozen-core NLCC tests

The production bridge accepts nonzero cores only with complete V2 provenance,
joint reference density, kinetic density and matrix current. Historical
`tau_nlcc` is never treated as orbital-derived core tau. The model is a regular
auxiliary-Pauli frozen core; it is not an exact relativistic Dirac CKED model.

In the ionic-rest reference gauge, `Ncore4=(rho,0,0,0)` and
`Tcore4=(tau,0,0,0)`, with tau defined **without** a factor one half.
The reference current retains all spin components, with zero charge component.
For the prescribed velocity-gauge phase `q=-Ac`, the bridge adds

```
Jcore(:,i) = Jref(:,i) + q(i)*Ncore
Tcore = Tref + sum_i [2*q(i)*Jref(:,i) + q(i)**2*Ncore]
```

Core and valence fields are combined before density gradients, Laplacians and
functional evaluation. The gKS action differentiates valence orbitals only.
The occupied-state band expectation therefore contracts with **valence**
density, tau and current, while the energy sees the augmented fields.
The explicit integrated core response follows by chain rule:

```
dE/dAc(i) = integral [-2*dtau.(Jref(:,i)-Ac(i)*Ncore) - dj(:,i).Ncore]
```

`test_core.f90` compares the phase prescription with explicit phase-transformed
core spinor gradients and Laplacians; checks all three field derivatives against
energy differences; co-transforms both core and valence under U(1) and local
SU(2); tests every scalar helper input with NaN/Infinity; and checks zero core.

`test_core_bridge.f90` uses the actual full SALMON bridge and Taylor propagation
on deterministic synthetic spinors. It checks scalar and spin-current reference
cores, orbital energy differences, valence-only band contraction, Hermiticity,
all three explicit `Ac` derivatives, zero-core regression and invalid metadata,
negative tau, missing tau and invalid layouts. Global spin conservation is tested
only for a scalar reference core without spin current or ionic SOC; a fixed
spin-current background need not conserve valence spin alone.

Pure tests:

```sh
python3 testsuites/unit_nc_mgga23_core/run.py \
  --checker /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/tools/check_forbidden_eigensolvers.py
```

Add `--build-dir BUILD --libs LIBRARIES` for the serial full bridge tests.
Add `--mpi --fc mpifort` with an MPI build for two-rank orbital decomposition.
Every build/link retains source and application-object prohibition audits.
Pure GNU Fortran tests use bounds checks and invalid/zero/overflow traps.
The bridge uses actual compiled SALMON objects, so rebuild the selected source
before running. These are numerical component tests, not physical GS collection,
a Ni demagnetization benchmark, or proof of pseudopotential transferability.
