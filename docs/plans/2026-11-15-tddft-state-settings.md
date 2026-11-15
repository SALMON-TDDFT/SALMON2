# TDDFT settings ownership: fixed timestep

Continue the approved independent-context work without GitHub Actions. First extract
fixed timestep from the common global input into each native RT state. Keep the
existing input path as the default and add an optional construction-time value to
the internal context, with finite/positive checks before allocation or output.
Preserve all fast HSE arithmetic and update order. No global timestep swapping.

Test two states with dt=0.001 and 0.0015 against separate processes whose input dt
matches each state. Deliberately change the global dt during interleaved steps;
check currents, energies, output times and full synthetic state traces, reverse
order and close/reopen. Keep same-dt coverage and the unchanged fast-HSE full-driver
pulse comparison. Audit sources/objects, run local tests and sdist-to-wheel checks,
and synchronize Japanese/English/AI documentation.

This stage does not allow different grids, XC, k meshes or adaptive timesteps.
Investigation found module-owned XC work arrays in salmon_xc and modx/mody/modz
optimization tables. These need explicit ownership before heterogeneous grids are
safe. No public Python ABI or general native failure atomicity is claimed.
