# Shared HSE in normal mixed sessions

Opt in with NativeFastHSEBackend(...,shared_model=True) and layout
points_per_group={material_id:N}. Defaults preserve one point/group and existing
identity/output behavior. Optional backend.open_group(initial,dt,context) returns
a session-local factory exposing new_state(point_id),close(); None uses existing
create_state. Resource lifetime is ResponseSession-owned, not cached on a reused
backend. User factories remain optional and are not restricted to built-ins.

Layer planning chunks points of each material in encounter order into groups of
up toN; ranks_per_point retains its existing value as the group's rank width.
Preserve all cells/point/material IDs. Never auto change N or silently drop points.
Known singleton native backends reject unsupported grouping before launch.

Checkpoint helper keeps old filenames/schema for one binary state per rank, but
uses rank+point filenames and native-direct/4 for multi-binary groups. All paths
(admission/save/fresh restore/protected stage) use one naming function. Preserve
old readers for versions1–3; do not imply cross-topology resume. Update identity
admission for version4, including explicit legacy/no-identity mode.

Verify tests first, source/installed and packaging; ordinary mixed sample includes
LD/SBE/PBE and two HSE points in one group, compared with separate groups and
independent serial responses. Continuous vs protected direct split/resume has no
replay. Test staging failure in a later point preserves all owners and permits
retry. Update JP/EN/AI/examples and latest checklist. No Actions.
