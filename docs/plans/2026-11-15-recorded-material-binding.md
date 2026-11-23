# Recorded endpoint-model binding

Approved M5 adapter continuation. Scope: serial endpoint-A models (configured SBE
and explicitly recorded user models), one MPI rank per material group. Existing
bind_material numerical checking remains. No implicit LD E conversion or bypass
of managed TDDFT provider rules. No native changes or Actions.

1. Tests first: bind_recorded_material copies/validates identity; identity record
   drift before prepare rejected; response_identity builds point-run/1 in integer
   transport order; incompatible shared IDs rejected. Known ConfiguredMaterial
   TDDFT wrappers must not bypass bind_material's managed-worker restriction.
2. Extend MaterialBinding with optional immutable material_record, default None.
   New public bind_recorded_material(initial,key) uses existing endpoint adapter.
   Public response_identity(bindings,point_map,run_id) constructs validated record,
   no model/state creation or MPI calls. Session cross-checks provided metadata
   against recorded bindings before callbacks, keeping legacy bindings allowed.
3. Add custom model + configured SBE local/MPI regression, exact direct-response
   parity; paired docs/AI JSON and sample. Full tests, sdist->wheel, installed sample.
   TDDFT/LD dedicated provider factories and restart identity remain separate work.
