# Collective retained-GS array admission verification

The existing unit contract first reproduced local TypeError for list/dict
quantities and native access for differing supported rank selections. Validation
now agrees supported strings before native array calls. Invalid requests reject
collectively with InputError; ordinary numerical copies are unchanged.

628 unit methods pass, one skipped, using the full dependency environment and
SHA-verified existing wheel LD/SBE kernels. No new library was built.
Synthetic MPI4 (two groups of two) covers PZ orbital partition and Libxc-PBE
spatial partition. One rank supplies a list while its peer supplies density;
then the ranks request density versus occupations. Both fail collectively.
After each failure, resources/current and original GS orbital bytes are exact;
a corrected read succeeds. Subsequent currents match the independent file route
and direct8+8 resume, with zero replay and no SCF calls. The orbital current SHA
also matches the pre-change f6005fb9 validation. No physical convergence claim.

The first MPI attempt without local socket access and the first source suite
without installed kernels were environment failures; preserved diagnostic logs
are separate from successful calculation checks. Repeating with local MPI access
and the already validated kernels resolved them. Source/docs/AI checks pass.

Records (logs/hashes/results only, no collected GS arrays):
/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/calc2026/salmonpy-gs-array-admission-20261119

Local record writes verified; cloud sync not established. No CI dispatched.
Full SCF continuation remains unimplemented; see the SCF checkpoint follow-up plan.
