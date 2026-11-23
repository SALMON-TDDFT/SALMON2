# Non-diagonal background energy accounting

The field observer now calls the independently tested Fortran accounting kernel using the actual background interpolation operator. The current observation buffer has18 components; the last is the separate constitutive exchange.

Let K be the symmetric off-diagonal background interpolation operator used by the field update, e0/e1 the old/new electric edge vectors, P the selected edge region, and V the cell volume. The additional electric energy assigned to the region is

    U_cross(e) = V/(8 pi) e^T P K e.

The outward constitutive exchange in one accepted step is

    Q_cross = V/(8 pi) [e0^T P K e1 - e1^T P K e0].

Then ΔU_cross + Q_cross equals V/(4 pi) times the region sum of E_mid dot K ΔE. For the entire domain Q_cross vanishes by symmetry. Across a region boundary it is generally nonzero and must be reported separately from physical E cross H flux. Local assigned cross energy need not be positive on its own.

`sp_fdtd_cross_energy` consumes E_old, E_new, K E_old and K E_new with a common Fortran (3,nx,ny,nz) edge layout. It returns old/new cross energies and outward constitutive exchange. It only sums the owned intersection with the observation box, performs scalar finite checks and publishes no output on failure. It performs no MPI, I/O, interpolation or field update.

The independent manufactured test uses a symmetric nonlocal matrix, checks nonzero regional exchange, cancellation on the full domain, additive spatial partitions, the midpoint-work identity and NaN/Inf rejection. This establishes the algebra, not the correctness of the solver interpolation or physical convergence.

Integration obtains K E_old/new from eh_tensor_off at pole0. Scratch is allocated before halo communication and failures are agreed within the field group. Admission adds one field-group max, allocation adds one field-group max, two operator applications add six halo exchanges per active material each, followed by the existing single scalar-vector reduction. Diagonal-only observation retains its previous communication count. No application-wide collective is added. Subpixel and Ohmic observations remain separate scopes.
