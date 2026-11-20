# Group5: TDDFT file scopes and common native call wrapper

Reuse the existing native scope signature, with protocol2 declaring additional
TDDFT path coverage. Retained GS accepts1/2; TDDFT selects scope only for2 and
keeps scoped-chdir for older runtimes. Route RT core and fast-PR adapter file
operations through resolve_path without changing arithmetic or checkpoint bytes.

Consolidate common directory entry/collective admission/finally cleanup in the
existing _native_group module, then use it for GS, standard legacy/shared/registry
TDDFT and fast-PR legacy/shared/registry wrappers. Preserve main-thread checks;
remove the redundant HSE registry _call override. Do not add a new module.
Native FDTD wrapper remains a separate next slice.

Verify wrapper tests before/after, ordinary/full/reduced HSE current histories,
interleaved models and direct resume against audited previous controls. Add an
explicit no-chdir guard only to new-runtime groups. Validate old protocol1
fallback and collective failure cleanup. Build/audit source+objects; source/wheel
and JP/EN/AI/list parity. No routine CI. Collected GS fixture copies require
campaign-specific authorization before the real bundle-based tests.

## Accepted local checks

Source and installed637 tests each636pass/1skip. Old wrapper no-chdir checks fail;
new variants, GS→RT capability selection and protocol1 TDDFT fallback pass.
MPI2/group, three groups: PBE/TB-mBJ heterogeneous and shared standard vs audited
old controls have exact currents/direct resume; full/reduced learned/PT-CN HSE
registry current/energy/history and shared-HSE controls exact. All new groups
forbid chdir, native direct replay0. Existing-GS ordinary/HSE bridge and old-HSE
fallback pass; no SCF, interleaved import/source-copy removal/independent child
continuation/direct resume exact, measured current/energy route difference0.
Wheel matrices and bridge controls pass. Stream numerical writers and HSE
checkpoint numerical payload declarations unchanged. Source/object audits and
prepared-source/wheel/JP/EN/AI parity pass. Campaign GS copies/checkpoints were
explicitly authorized; source data were reused, not newly computed. Evidence:
calc2026/salmonpy-td-directory-record-20261120. No CI; cloud sync unconfirmed.
