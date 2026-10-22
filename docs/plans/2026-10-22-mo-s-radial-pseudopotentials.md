# Mo / S radial-potential exploration

User authorized Mo and S creation, one-channel Rps variations, radial-shape comparison.
Subsequent explicit user approval: fix unused/unoccupied Mo 4f atomic solve without changing neutral configuration; verify W/Se reproducibility.

Use current CRPS SOC/NLCC generator; Mo12 (5s/4p/4d), S6 (3s/3p/3d), LDA/PZ inherited functional.
Fix core radius0.6/0.4 Angstrom respectively. Sweep each s/p/d pair2.0,2.2,2.4,2.6,2.8,3.0,3.2 bohr, other channels2.6.
Preserve failures; require unique c2 root, finite complete output, correct channel/valence metadata.
Compare full radial curves, origin detail and curvature; form one combined shape candidate per species.
All outputs/plots/inputs/provenance in canonical OneDrive/calc2026/mo-s-rps-scan-20261022.
No scattering/ghost/solid acceptance, no existing binaries/jobs or standard pseudopotentials replaced.

Generator regression: original Mo fails classical turning point before5s; patched Mo generates complete finite12-valence file.
Bound unused orbitals retained to preserve historical continuum reference energy; unbound nonreference states excluded,
wavefunctions zeroed and next orbital energy seed restored. W/Se byte-identical after repair.
