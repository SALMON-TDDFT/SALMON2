# Occupation recovery — 2026-11-16

The temporary worktree and previous artifacts were unavailable after reboot. The occupation snapshot, native probe, original-solver comparisons, package entries and guides were reconstructed from conversation records atop 0a5d6e04. This is a reconstruction, not recovery of byte-identical uncommitted files.

Fresh outputs live under `/Users/otobetoshihito/.codex/outputs/salmonpy-state-inputs-20261116`: occupation-build.json/reference-build.json, occupation-probe/result.json, occupation-full/result.json, recovery-lifecycle/result.json, recovery-dist, recovery-source-tests.log and recovery-installed-tests.log. Apple Silicon/GNU/OpenBLAS/static Libxc, serial CPU, two build workers. Original comparison passed for six PZ/HSE fixed/zero/finite-temperature cases, with exact output/restart hashes and residuals. Two occupation states with different electron counts/temperatures ignore poisoned globals and match unchanged original ne2mu. Four retained lifecycle cases pass. Python source and installed suites each run435 tests. No remote CI or physical convergence claim.

HSE, DC and parallel global consumers remain at this recovery checkpoint; subsequent HSE work is separate. Historical /tmp links in earlier validation records no longer identify available artifacts. The persistent checkout replaces the missing temporary checkout; no upstream merge/push.
