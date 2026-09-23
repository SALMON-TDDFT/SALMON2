# ELF pump-probe calculation plan

Goal: determine whether the exciton-like optical peak changes after the same Si pump under the new fully dynamic ELF-alpha model.

Use the authorized existing TDCDFT branch and experimental Si setup (4³ k, no convergence scan). Reuse the prior pump-probe delay and observation window: pump duration60 au, omega.2 au, intensity1e13 W/cm², probe at80 au, nt12000, dt.08. Keep alpha0=.2 and ELF stride10 active throughout. No new SALMON behavior or unrequested damping.

1. Prepare reproducible inputs and isolated calculation directories for pump only, pump±probe, half-size pump±probe, and unpumped probe. Start with pump only, positive probe1e-4, unpumped probe. Propagate from GS to verify pre-probe equality. Store sparse checkpoints and all current/alpha output; check completion, electron number and finite arrays before extracting spectra.
2. If trajectories are valid, use [J(pump+eta)-J(pump-eta)]/(2eta), compare against eta/2 central difference and one-sided pump subtraction. Extract finite-window differential epsilon with the existing tested Fourier sign/window convention. Compare against a fresh unpumped ELF response using the same probe time/window. Retain alpha response to the probe.
3. Compare peak position, height, signed2–4 eV integral, apparent FWHM and windows600/720/880 au. Energy sampling.01 eV is not spectral resolution. Do not infer exciton binding or lifetime from coarse finite-window peaks.
4. If native full feedback fails at long times or lacks a linear-probe limit, document the failure and do not label a truncated/invalid transform a converged peak. Examine the alpha/current/field evolution; any frozen-alpha auxiliary calculation must be explicitly labeled partial response, not silently substituted.
5. Verify central-difference analysis with synthetic even/odd probe terms, run existing Fourier tests, review the scientific comparison, save inputs, spectra, metrics and plots, and commit locally. No push.

## User steering: include beta/gamma stabilization

The initial beta=gamma=0 pump-only and pump+probe runs fail norm at6.50 fs. The user points out that this requires the Proca stabilizers, not rejection of the ELF alpha law. Continue the original goal by connecting existing ELF measurements to the existing normalized Proca oscillator, without altering the LRC polarization mode.

For tdcdft='proca', tdcdft_screening='elf', use a''+beta a'+gamma a=alpha(t)J, explicitly a variable-alpha extension of Williams-Ullrich Eq26. This differs from differentiating a'=-alpha P (which would produce an extra -alpha'P term); no such term is added to the Proca mode. Existing checkpoint saved_mode distinguishes these closures. Keep alpha0=.2, stride10, no upper clipping. Reject direct a2/a0 in ELF mode, retaining explicit normalized coefficients.

Test acceptance, oscillator recurrence with changing alpha and nonzero beta/gamma, normalization, odd restart, freeze, changed gamma and invalid alpha. Preserve legacy tests. Start beta=0 as in the paper; screen gamma1e-4,4e-4,1e-3 on short strong-pump extensions, then verify stable candidates over the full window and compare unpumped spectra. Do not treat literature2D thresholds as Si constants. Complete pump-only, central +/- probes and half-amplitude checks with the chosen common stabilizers. Check gamma sensitivity if it materially changes the optical peak. Paper low-frequency auxiliary resonance must remain distinguishable from the optical peak.


## Completed outcome

ELF-Proca support and serial/MPI recurrence/restart tests are complete. Short gamma
screen and full gamma=.001 central/half-probe comparison complete; gamma=.0004
unpumped reference also complete. Finite-window optical probe amplitude agreement
is0.2122%, but the strongly signed, window-dependent pumped response does not
support a quantitative exciton bleaching/shift claim. Results, raw compressed
traces and reproduction scripts are in ../results/si-elf-pump-probe/proca/README.md.
No k convergence scan, no silently frozen alpha, no fit of beta/gamma to bleaching.
