# UEG density-matrix distance implementation plan

Goal: evaluate the approved unclipped alpha(t)=0.2*d(t)/d(0) offline every10 steps on existing Si trajectories.
Architecture: Fourier-transform the periodic occupied orbitals on the existing real-space grid; construct a homogeneous free-electron reference on the same shifted supercell momentum mesh with the same electron number. Fractionally fill the last degenerate Fermi shell. Match mean flow by convex interpolation of adjacent integer momentum translations. Compute Hilbert-Schmidt distance via occupied Gram matrices and plane-wave diagonal weights, without storing a full density matrix. Preserve reference purity and finite-mesh limitations explicitly.
Tech stack: Python, NumPy, existing checkpoint reader, Matplotlib. Existing TDCDFT branch and trajectories; no production changes or k scan.

1. Test FFT normalization/momentum indexing, Fermi shell count/isotropy, shift number/current conservation and positivity, HS identity against explicit density matrices, unitary invariance and unclipped alpha.
2. Implement core in docs/results/si-ueg-distance/ueg.py; run tests first failing then passing.
3. Reference matching: primary target canonical momentum is measured SALMON electron current / density minus effective A/c (classical plus XC). Also evaluate orbital-FFT momentum matching and resting-reference controls, as current includes finite-difference/nonlocal-pseudopotential terms. Report errors and reference mixing/purity.
4. Analyze all3 x161 snapshots. Record per-k fixed-rank lower bound: propagated occupied subspace has16 states at each k, whereas a Fermi sea need not. Do not claim this diagnostic guarantees metallic alpha=0 along a reachable trajectory.
5. Compare with stored MLWF/invariant spreads; inspect plots, obtain independent numerical review, document findings, test and commit locally.
