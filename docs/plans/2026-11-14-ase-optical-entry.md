# ASE optical entry implementation plan

**Goal:** Infer and check neutral electron counts from pseudopotentials and connect ASE Atoms to independent A-driven TDDFT states.

**Approved design:** Keep the optional ASE adapter thin. `Pseudopotential.valence_electrons` reads KY's fourth header field (matching Read_PS_KY) or NC-UPF's PP_HEADER z_valence. Reject absent, nonfinite, noninteger or impossible charges rather than round or guess. `dft_from_atoms` infers electrons when omitted and rejects explicit counts inconsistent with a neutral cell. Grid, states and kmesh remain required, since no universal converged defaults exist. Calculator uses the same validation. Core DFT stays explicit for compatibility.

`tddft_from_atoms(atoms, directory, *, executable, pseudopotentials, states, grid, kmesh, ...)` runs the existing DFT and GS-to-TDDFT loader, returns its context-managed initial state, and keeps the GS result in directory. No new propagator, implicit runtime acquisition, changed update order or E-driven crystal API. Failure retains GS diagnostics. Native one-rank, periodic, unpolarized, fixed-ion restrictions remain.

**Alternatives:** Inferring from atomic numbers would count core electrons incorrectly. Automatic grid/state defaults would obscure convergence choices. A second propagation wrapper would duplicate the existing crystal API. All are intentionally avoided.

**Implementation:**
1. Add failing unit tests in test_dft.py and test_ase.py for header parsing, count mismatch, preserved structure, preflight rejection and GS-to-initial composition. Run them with optional ASE installed in an external venv.
2. Implement the valence property and the two ASE entry changes; allow Calculator electron inference. Re-run focused tests, including Calculator cache/failure behavior.
3. Add an executable Si sample and integration test: use the existing starter inputs/runtime/reference, check charge current tolerances and replay independent states from one initial state. One worker at a time, one thread, 4 GiB/4-core user budget; resource measurement status separate from calculation.
4. Update Japanese ASE guidance, add English guidance, synchronize AI JSON and sample catalog. Run complete Python suite and source-distribution-to-wheel packaging check because a guide is added to wheel installation.
5. Record measured outcomes and limits; inspect final diff. No new force/XC/spin/isolated-system or warm-start claims.
