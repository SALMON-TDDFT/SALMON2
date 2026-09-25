# Beta scan: completed

Three unpumped trajectories completed normally (487–497 s each), alpha=.2, gamma=.004; beta=0 reused. Same dt=.08, full4³, 960 au duration, probe at80 au. No strong-pump runs were added.

Increasing beta reduces the negative minimum in the original1.7–2.4 eV band, but does not remove negative response. The beta=.016 trace develops a stronger negative lobe near2.67 eV and a growing late-time Axc envelope. Its positive maximum increases, rather than simply broadening and weakening. Normal process termination does not establish physical stability.

beta=0: minimum over 1.5–4 eV = -173.384 at 2.06 eV
beta=0.001: minimum over 1.5–4 eV = -139.207 at 2.06 eV
beta=0.004: minimum over 1.5–4 eV = -82.281 at 2.09 eV
beta=0.016: minimum over 1.5–4 eV = -211.160 at 2.67 eV

The narrow-band metrics.json must not be used alone to select beta: they miss the new negative structure outside2.4 eV. comparison.png displays this failure. Peak heights use the same880 au cubic window; finite-window spectra of the growing response are not stationary absorption spectra. These runs do not support beta alone as a remedy for the equilibrium-model issue. Further timestep/model feedback analysis would be needed to distinguish the origin of growth; no additional runs are launched here.
