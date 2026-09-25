# Si 4³ strong pump–probe: equilibrium-model issue under analysis

**Latest assessment (2026-09-25):** gamma=.004 has a substantial negative equilibrium response near2.06eV that survives cubic and exponential window changes. Physical interpretation of its pump-induced peak changes is suspended pending model diagnosis. See [gamma diagnostic](gamma_diagnostic/README.md). All propagation is stopped; older progress entries below are historical.

Fixed alpha=.2, beta=0; gamma is the restoring coefficient, not damping. Acos2 pump: omega=.2 au (5.44 eV), width60 au, intensity1e13 W/cm², z. Probe80 au; ±1e-4 vector-potential steps. Full shifted4³ mesh, same real-space grid and method-specific converged ground states as the completed impulse comparison. Runtime and raw files: calculations/si_hse_native/k4_pump_probe.

At dt=.08, gamma=.001 pump-only and ±probe fail around17.2fs with large XC fields and catastrophic electron-number loss. The unpumped trajectory completes. No physical spectra are inferred from these failed pumped runs. Gamma=.004 pump completes23.22fs; its probe/reference calculations are still running. Half-dt pump diagnostics for both coefficients and half-amplitude probes at gamma=.004 are queued after the initial8 cases. These discriminate numerical error and nonlinear probe contamination; successful completion is not by itself long-time/model validation.

## HSE changes under validation

Full-k linearly polarized Acos2 with optional impulse is enabled for HSE Taylor4/ACE and fullEXX. The nonlocal pseudopotential current phases now use endpoint A consistently with kinetic current; propagation uses midpoint A. Exchange is rebuilt from the propagated occupied density matrix. A common spatially uniform A cancels in the momentum-transfer difference; no separate fixed ground-state exchange is substituted. The self-consistent exchange contribution to the total commutator current cancels when summed over the same occupied density matrix; ACE is checked against fullEXX to assess its approximation.

Short artificial pulse (1.28au, dt=.08/.04/.02): ACE/full relative current difference1.94e-6; timestep refinement ratio3.947. These tests passed but are not validation of long-time spectra at production dt=.16. Probe and zero-field cases and energy-work diagnostics remain required.

Pulse RT restart is explicitly rejected pending pulse-parameter metadata support. Pulse checkpoints now use version3/4 so changing the input to impulse cannot load them as legacy version1/2. The metadata patch has been rebuilt; native checkpoint rejection tests remain pending. Existing impulse restarts retain their previous formats.

No HSE long pump–probe production job has started. Do not report the two-method comparison complete.

### Additional checks completed

Rebuilt metadata versioning passes both negative native checks: pulse restart rejected, and pulse checkpoint reinterpreted as impulse rejected. Pulse+probe and zero-field16step runs complete. Zero-field max|J|8.49e-12 au, max energy drift2.80e-12 Ha. Energy-work mismatch for the short pulse decreases5.15e-7 →4.36e-8 →4.43e-9 Ha with dt=.08,.04,.02 (endpoint-average J integrated against delta A). Actual width60au pump pilots at dt=.16/.08 are queued after TDCDFT refinements, MPI8OMP1, through80au. Their completed traces must be compared before long HSE production.

### TDCDFT refinement complete (2026-09-24)

Gamma=.004 pump and ±probe plus equilibrium reference all finish23.22fs. Half-dt pump also finishes. Relative L2 current difference(.08 vs .04) is0.0829%. Probe amplitudes±1e-4 vs±5e-5 yield0.0338% relative post-probe current difference and0.0122% relative Im epsilon difference in2–6eV. These are finite-time sensitivity checks, not k/grid convergence. Gamma=.001 still breaks near17.5fs at half dt; halving dt alone does not fix it.

The preliminary optical-range plot retains negative spectral features and observation-window dependence. Below1.5eV the raw transforms have large signed artifacts/features; CSV retains the full range. A simple exciton bleaching/shift interpretation is premature. HSE comparison remains outstanding.

HSE production runner waits for both80au pump pilots to finish and requires current relative dt difference<0.1% plus logged norm error<1e-5. On pass it sequentially runs6000steps dt=.16 (960au) for pump, ±probe, equilibrium at MPI8OMP1, with48GiB/6hr limits per case. Monitor hse_production_status.json; do not duplicate. Passing short pilots does not establish long-time convergence, which must be assessed on the resulting trajectories.

### First full HSE comparison (2026-09-25)

All four HSE trajectories completed6000steps. Each took about2h10m; maximum summed RSS1.42GiB. comparison_preliminary.png/pdf and comparison_metrics.json compare identical post-probe windows. At880au the response at the *unperturbed maximum's energy* falls from86.92 to52.18 at4.30eV in HSE, and259.51 to27.32 at2.76eV in fixed-alpha TDCDFT(gamma=.004). These are fixed-energy reductions, not matched excited-state peak heights or fitted peak shifts. The signed2.5–4.5eV integrals drop45.92→15.96(HSE) and93.21→26.40(TDCDFT). Additional signed structures and observation-window dependence preclude a simple single-peak interpretation. The two reference maxima should not be assumed to represent the same exciton.

HSE even/odd probe-current ratio is4.62%; central differencing removes the even component, but odd nonlinear contamination still requires an amplitude check. Native ±5e-5 HSE runs are therefore launched sequentially under run_hse_halfprobe.py, monitored in hse_halfprobe_status.json. Comparison remains preliminary until these complete. TDCDFT amplitude refinement was already checked.

## Analysis after user-requested stop (2026-09-25)

User stopped additional HSE half-amplitude calculations. The halfprobe controller and MPI group were terminated and their status marked stopped_by_user. No partial halfprobe data enter this analysis. The heartbeat was paused. Use the completed ±1e-4 central derivative as the primary HSE result; amplitude refinement is not established for HSE and is not required to describe the present qualitative comparison.

At the common880au post-probe window, HSE unpumped maxima near3.68,3.88,4.30eV change to a substantially weaker, restructured spectrum. The local maximum in4.1–4.5eV is4.30eV/86.92 before and4.26eV/59.30 after pump (about32% lower local maximum). The response at fixed4.30eV is40% lower. The apparent−.04eV displacement is only the location of a finite-window maximum, not a precision exciton binding-energy shift. Signed area in3.6–4.1eV falls20.33→4.87 and4.1–4.5eV falls23.74→11.27.

TDCDFT(gamma=.004) is qualitatively less uniform: the dominant2.76eV structure falls259.51→27.32 at fixed energy, while the3.30eV/101.97 local maximum is replaced by a3.37eV/148.33 structure and a3.64eV/98.09 structure develops. Thus saying all peaks bleach would be incorrect. The signed3.1–3.6eV area still decreases28.27→18.32, while3.6–4.1eV rises2.74→10.68; peak height alone obscures the redistribution.

Both show strong pump-induced spectral restructuring, but they do not agree feature-by-feature. Comparison is at equal external pump, not equal excited carrier density (not measured here). These are different baseline functionals and TDCDFT gamma=.004 changes its equilibrium spectrum relative to the prior gamma=.001 comparison. Large signed lobes and observation-window dependence remain. Do not infer a unique linewidth/lifetime from these non-isolated peaks, or label negative lobes optical gain without further analysis. The strongest defensible result is reduced response around original major features together with new spectral structures, with a much stronger redistribution in this TDCDFT model. No new dynamics are launched.

## Consolidated TDCDFT assessment (2026-09-25)

The [validation note](../../tdcdft-validation.md#2026-09-25-update-fixed-alpha-strong-excitation-and-equilibrium-response-failure)
now consolidates the gamma scan, negative unpumped response, beta scan and revised
physical interpretation. The beta scan was subsequently authorized and completed:
beta=.001/.004 reduces the original negative peak but does not remove it; beta=.016
creates a larger negative structure near2.67 eV and a growing late-time Axc envelope.
See [beta results](beta_scan/README.md). The equilibrium issue remains unresolved;
normal completion is not physical validation. No new strong-pump calculations were
performed for beta and the heartbeat remains paused.
