# Public fast-PR symmetry and protected mixed connection

Goal: Expose the verified scalar/full-q/internal-k1 factor-off Taylor profile
through explicit NativeFastHSEBackend(...,allow_symmetry=True), preserving old
constructor behavior and C ABI, and verify native E/H/LD/SBE/TDDFT/HSE mixtures.

Architecture: retain sp_fast_hse_init default-denied; add sp_fast_hse_init_symmetry
with identical arguments and sp_fast_hse_symmetry_protocol()==1. Both call one
initializer with an explicit permission flag. Status34 means incompatible A,
no mutation, ready state retained. Python requires unanimous finite status34 and
this new protocol before raising retryable InputError; unknown/old cases fail.
The opt-in changes backend identity and group-size admission, not input matrices.

1. Fail-first constructor/lazy/profile/identity and field-rejection recovery tests;
   native query missing in old library. Implement C/Python entry and bounded
   profile admission before native lifetime acquisition.
2. Direct public MPI1-group tests: nonzero A, forbidden direction/retry, save/load
   and staged/protected restore, independent groups with different settings.
3. Ordinary mixed sample adds --hse-symmetry and explicit --hse-reference-bundle;
   never convert arbitrary GS payloads. Reference is supplied full-k equivalent,
   serial fast-PR worker. Include flag/reference input in experiment identity.
   Validate two field-local LD + SBE + PBE TDDFT2 + symmetry HSE1, field1/2.
4. Protected whole-step late/field-admission failure and protected restore stage
   failure preserve all original field/LD/electronic bytes and ordinary streams;
   retry/direct resumed A/E/j matches continuous without elapsed replay.
5. Keep conventional HSE, learned/PT-CN symmetry/internal-k2 rejected. JP/EN/AI,
   sample help/tests/checklist; source/installed tests and sdist-wheel. Local-only,
   warning-only resource advice, persistent evidence and skip-ci commit/push.


Measured scope extension/limitation: phase0 exposed pre-existing Libxc5.2.3
near-zero-gradient sensitivity (pure C diagnostic reproduces it); do not weaken
tolerances or patch the library. Phase.23 is explicitly selected only for synthetic
connectivity tests. Preserve failed results and put the low-gradient treatment at
the top of the remaining numerical checklist; no robust physical-accuracy claim.
