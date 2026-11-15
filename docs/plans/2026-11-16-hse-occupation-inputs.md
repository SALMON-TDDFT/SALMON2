# HSE occupation input isolation

Route all exx_native nelec/temperature/nstate decisions through the initialized system occupation snapshot, with a legacy fallback for manually constructed/uninitialized systems. Keep arithmetic, branch priority and execution order. Other HSE switches, spin flags, shared Wannier/spatial state and outer TDDFT guards are not generalized here.

Red/green native owner test: A nelec4/temp-1/nstate2; B nelec2/temp-2/nstate3 (RT target capacity2), separate geometry and screening; poison globals99/.5/0 before binding. Compare refresh/action/energy across interleave, reference states, solo processes and close/reopen in ACE/full-action modes, OMP1/2. Verify original GS fixed/finite-temperature parity and full TDDFT regression at supported common inputs. Retain TDDFT configuration guard because non-HSE output consumers still use legacy inputs. Update guides, source inventory and wheel evidence; no remote CI.
