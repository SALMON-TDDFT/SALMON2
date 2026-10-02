# SR-only HSE without MLWF localization

For static-ion Gamma HSE with uniformly doubly occupied orbitals, positive SR tolerance, full WF support (norm fraction 0 or 1, radius 0), source ACE, no pair screening and no WF snapshots, use canonical occupied orbitals directly. Spatial decomposition and screened kernel remain enabled. Other paths retain localization. A rejected canonical source ACE stops rather than falling back to full exchange.

Validation: 8 MPI x 2 OMP Si8x1x1, SR1e-2, dt0.16, 16steps: all ranks exit0, finite data, ACE rejects0, norm error1.09e-6. Max current-table difference1.36e-17, energy-table difference9.95e-13 Ha vs prior fixed binary. Action invariance tested for SR1e-2 and1e-3 with active neighborhood FFT on8ranks. No significant speed gain demonstrated; rank peak RSS193.17->185.36MiB in single runs.

check_canonical_route.py checks a completed SALMON output for bypass, absence of localization output, active SR and no ACE rejection. This change does not include experimental multi-step ACE reuse or change its default.
