# Polarization-consistent instantaneous screening

User approved testing E_xc=alpha(t)P, beta=gamma=0, same 4³ k laser setup.
Keep 'instant' acceleration closure for controlled comparison; add 'polarization' mode.
Use the existing instantaneous estimator and P=-integral j, with no temporal average.
Integrate a'=-alpha P with a second-order explicit Taylor step:
  a_next = a_now - dt*alpha_now*P + dt²*alpha_now*j/2
           - dt*(alpha_now-alpha_previous)*P/2.
Backward alpha derivative is first-order, adequate in the dt² term for global second-order
accuracy for smooth alpha. An abrupt alpha change needs timestep sensitivity testing.
For constant alpha and trapezoidal P this reproduces the existing undamped LRC update with
consistent initialization. No extra restart state: previous alpha is already saved.
Centered output E agrees with alpha P to time-discretization error during coefficient changes;
after alpha becomes constant, the offset vanishes to roundoff, unlike the old closure.

1. Failing unit tests: varying smooth alpha with known integral; constant alpha equivalence;
   sudden alpha change followed by no current must not leave a persistent field offset.
2. Implement pure advance helper; route only new mode through it. Require beta=gamma=0.
   Preserve old mode and fixed-alpha outputs. Restart mode checking rejects closure changes.
3. Serial/MPI builds, analytic and Si integration tests including restart and invalid settings.
4. Same 12000-step strong laser run at 4³ k, matching existing inputs except mode. Then a
   half-dt run over the same duration if needed to check observed behavior. Compare a, E, j,
   norm, alpha and E-alphaP. Keep original data intact.
5. Review, documentation, scientific results and commit; no push.
