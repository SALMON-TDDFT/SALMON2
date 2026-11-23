Distribution notes recorded 2026-11-24; current instructions are in SI_INTRO.md / SI_INTRO.en.md.

# Si starter: ground state to current

[日本語](SI_INTRO.md)

Compute the fixed Si SCF input and compare 21 samples of A/c-driven charge current against the reference. The current electron-number-bearing worker protocol was revalidated with a newly built calculation program on 2026-11-23. Inputs, potential, expected results and tolerances were unchanged.

## Required distributions

Obtain an updated SALMONpy wheel for your OS, CPU and Python, plus an outer runtime ZIP listed in the current `runtime_artifacts` of `provenance.json`. NumPy is the only Python dependency of this example.

The new locally validated distribution was built on Apple Silicon macOS. It has not been uploaded to a public Release. Obtain the matching wheel and ZIP from the distributor; do not invent a download URL. Use the wheel-only [LD/SBE test](README.en.md) if the matching new ZIP is unavailable.

The old `runtime-2026.11.14` Release lacks the current electron-number reply. Its hashes and origins remain under `historical_runtime_artifacts` and related fields, but the current entry rejects it. New Intel Mac and Linux distributions still need validation. Windows DFT/TDDFT distribution is not established.

