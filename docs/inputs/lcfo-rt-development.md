# DC LCFO initial states and real-space RT

The fixed-LCFO-subspace RT path has been removed. DC-LCFO remains available for ground-state diagonalization and for constructing initial real-space orbitals. RT propagates the reconstructed orbitals on the real-space grid without retaining or projecting onto the LCFO basis.

The former `yn_hse_lcfo_rt`, direct-WF, LCFO ACE/U interval, seed and diagnostic input options are no longer accepted. Remove them from old inputs. See [HSE inputs](hse.md) for the current real-space exchange controls. Historical reports and plans describe earlier experimental versions and are not current execution instructions.
