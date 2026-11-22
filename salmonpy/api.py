"""Explicit public names; importing this inventory does not import components.

Listed APIs are preview APIs covered by API.md/API.en.md, not a stable C ABI.
Returned state/result objects expose the documented methods and properties;
private names and incidental imports are not public APIs.
"""
from types import MappingProxyType

__all__ = ['PUBLIC_API', 'native_contract', 'query_host_abi']
PUBLIC_API = MappingProxyType({'salmonpy.subpixel': ('SubpixelSmoothing','effective_material'), 'salmonpy.optical_materials': ('OpticalMaterial','list_materials','load_material','wavelength_to_frequency'),
 'salmonpy.analysis': ('SurfaceFlux', 'surface_flux', 'PortSpectrum', 'absorbed_energy_density', 'reflection_transmission', 'Table', 'read_table', 'Spectrum', 'spectrum', 'OpticalSpectrum',
 'OpticalSensitivity', 'KCurrentContributions', 'optical_sensitivity', 'k_current_contributions', 'LongWavelengthCorrection', 'long_wavelength_correction', 'optical_conductivity', 'dielectric_function', 'dielectric_tensor', 'waveplate_diagnostic', 'optical_linearity', 'optical_sum_rule'),
 'salmonpy.api': ('PUBLIC_API', 'native_contract', 'query_host_abi'),
 'salmonpy.doctor': ('diagnose', 'probe_native_protocols'),
 'salmonpy.convergence': ('plan_grid_scan','run_grid_scan','plan_dt_scan','run_dt_scan','run_current_grid_scan'),
 'salmonpy.gs_input': ('GSInputMetadata','describe_gs_input'),
 'salmonpy.native_steps': ('native_step_targets',),
 'salmonpy.stepping': ('StepTarget','StepRecoveryError','guarded_step'),
 'salmonpy.memory': ('MemoryEstimate','MemoryReport','MemoryBudgetWarning','process_memory','memory_advice'),
 'salmonpy.checkpointing': ('CheckpointCapability', 'check_checkpoint_compatibility', 'CheckpointTarget', 'execute_checkpoint', 'CheckpointState', 'restore_checkpoint', 'AtomicRestoreError'),
 'salmonpy.ase': ('BOHR_ANGSTROM', 'explicit_path',
                  'atoms_from_initial',
                  'validate_structure',
                  'to_ase_band_structure',
                  'ASEBandRun',
                  'run_bands',
                  'dft_from_atoms',
                  'tddft_from_atoms'),
 'salmonpy.ase_calculator': ('SalmonCalculator',),
 'salmonpy.bands': ('HARTREE_EV', 'BandStructure',
                    'BandPath',
                    'ExplicitKPoints',
                    'BandRun',
                    'BandInitialState'),
 'salmonpy.point_runs': ('PointRun', 'PointState', 'open_point_arrays'),
 'salmonpy.initial_states': ('tddft_from_catalog',),
 'salmonpy.builtin_materials': ('builtin_material_registry', 'material_from_catalog', 'ConfiguredMaterial', 'ConfiguredState', 'ConfiguredGroundState'),
 'salmonpy.c_material': ('CMaterial', 'CABIError'),
 'salmonpy.config': ('MaterialRegistry', 'MaterialCatalog', 'resolve_materials', 'load_materials'),
 'salmonpy.fast_hse_mpi': ('NativeFastHSEBackend', 'NativeFastHSEModel', 'NativeFastHSERegistry'),
 'salmonpy.fast_hse': ('FastHSEMaterialProvider', 'FastHSEResponseBackend'),
 'salmonpy.layers': ('NativeLayerAssembly', 'assemble_native_layers'),
 'salmonpy.coupled': ('NativeEHField', 'MultiscaleField', 'CoupledStepper', 'PacketCoupledStepper'),
 'salmonpy.crystal': ('AResponseState', 'AInitialState', 'CrystalResponse', 'respond'),
 'salmonpy.native_gs': ('NativeGSBackend','NativeGSSession','NativeGSState','NativeGSArray','NativeGSInfo','NativeSCFSnapshot'),
 'salmonpy.dft': ('Pseudopotential', 'SCFSettings', 'SymmetryOperations', 'DFT', 'GroundStateResult'),
 'salmonpy.distributed': ('SolverGroup',
                          'ParallelPlan',
                          'ParallelContext',
                          'NativeParallelContext'),
 'salmonpy.embedded_tddft': ('NativeTDDFTInitialState',
                             'NativeTDDFTMaterialProvider',
                             'NativeTDDFTBackend',
                             'NativeTDDFTModel',
                             'NativeTDDFTRegistry'),
 'salmonpy.errors': ('SalmonpyError',
                     'InputError',
                     'InputTypeError',
                     'NumericalError',
                     'RuntimeFailure',
                     'SnapshotError',
                     'SCFReviewRequired'),
 'salmonpy.fdtd': ('PlaneMonitor', 'PlaneData', 'read_native_planes', 'join_plane_segments', 'combine_plane_samples', 'EnergyBalance', 'energy_balance', 'Boundaries',
                   'Dielectric',
                   'Impulse',
                   'PlaneWave',
                   'FDTD',
                   'FDTDResult',
                   'FDTDMPI',
                   'SalmonExecutable'),
 'salmonpy.isolated_tddft': ('IsolatedRuntime', 'load_isolated_runtime', 'DipoleObservables', 'IsolatedTDDFTInitialState', 'IsolatedTDDFTState'),
 'salmonpy.dipole_source': ('DipoleMap', 'drive_native_dipoles', 'MolecularDrive', 'select_molecular_field'),
 'salmonpy.field_sampling': ('SampledField', 'sample_packet'),
 'salmonpy.field_input': ('FieldValue', 'FieldPacket', 'FieldContract', 'FieldAdapter'),
 'salmonpy.fields': ('FieldSamples', 'AcosPulse', 'sample_pulses'),
 'salmonpy.gs_native': ('GammaRestartConverter', 'SerialRestartConverter'),
 'salmonpy.gs_archive': ('GroundStateArchive', 'save_ground_state'),
 'salmonpy.geometry': ('vector', 'Grid', 'Shape', 'Geometry'),
 'salmonpy.ld': ('LDModelConfig', 'LDTensorModelConfig', 'LDStateConfig', 'Response', 'LorentzDrude', 'TensorLorentzDrude', 'LDState', 'FitReport', 'fit_lorentz_drude', 'fit_lorentz_drude_tensor'),
 'salmonpy.material_response': ('SolverContext',
                                'MaterialBinding',
                                'packetize_A_binding',
                                'ResponseSession',
                                'LDMaterialProvider',
                                'LDResponseBackend',
                                'TDDFTPreparedState',
                                'TDDFTMaterialProvider',
                                'TDDFTResponseBackend'),
 'salmonpy.materials': ('MaterialCheck', 'bind_material', 'check_material', 'bind_recorded_material', 'response_identity', 'bind_configured_ld_endpoints', 'bind_configured_tddft'),
 'salmonpy.maxwell': ('SPEED_OF_LIGHT', 'FieldStateConfig',
                      'Transverse',
                      'Longitudinal',
                      'ThinFilm',
                      'Multiscale1D',
                      'CoupledResponse',
                      'couple'),
 'salmonpy.preparation': ('prepare_inputs', 'prepare_fdtd', 'prepare_dft', 'prepare_tddft', 'verify_inputs', 'preflight_inputs', 'check_job'),
 'salmonpy.outputs': ('CubeGrid',
                      'read_cube',
                      'DecomposedCurrent',
                      'read_decomposed_current',
                      'OutputArchive',
                      'save_arrays'),
 'salmonpy.runtime': ('RuntimeUnavailable', 'RuntimeValidation', 'prepare_runtime_job',
                      'host_platform',
                      'inspect_runtime',
                      'load_runtime',
                      'load_validated_build',
                      'install_runtime',
                      'build_runtime'),
 'salmonpy.salmon_vr': ('VRExport', 'export_cubes'),
 'salmonpy.sbe': ('SBEInitialState', 'SBEState'),
 'salmonpy.snapshots': ('StateSnapshot',),
 'salmonpy.structure': ('StructureInput', 'GridSelection', 'grid_from_spacing'),
 'salmonpy.tddft': ('OutputSaveError', 'TDDFTInitialState', 'TDDFTState'),
 'salmonpy.timing': ('SCHEMA', 'InputPort',
                     'OutputPort',
                     'TimeSpec',
                     'crystal_default',
                     'matter_timing',
                     'charge_sign',
                     'check_loop')})



# One declaration supplies documentation and ctypes bindings. Tokens describe
# C value/pointer types, not Python-owned buffers or Fortran derived objects.
_HOST_OPERATIONS = (
    ('sp_attach', ('int64',), 'int', 'solver group', 'host MPI initialized; valid intracommunicator; no attachment',
     'duplicates communicator; host retains original', '0 success; -1 handle width/range; otherwise MPI status'),
    ('sp_detach', (), 'int', 'solver group', 'attached; all owned states closed',
     'releases duplicated communicators; never finalizes host MPI', '0 success; otherwise MPI status'),
    ('sp_info', ('int*', 'int*'), 'void', 'local', 'attached',
     'writes rank and size to caller buffers', 'no returned status; valid output pointers required'),
    ('sp_barrier', (), 'void', 'solver group', 'attached',
     'synchronizes attached solver group', 'no returned status; MPI failure recovery not guaranteed'),
    ('sp_sum', ('double',), 'double', 'solver group', 'attached',
     'returns sum on every solver rank; does not modify solver state', 'no returned status; MPI failure recovery not guaranteed'),
    ('sp_library', ('char*', 'int', 'int*'), 'void', 'local', 'valid buffers; capacity >= 1; query before attach allowed',
     'writes truncated NUL-terminated MPI library text and initialized/not-finalized flag',
     'no returned status; invalid pointer/capacity is outside contract'),
)


def _host_operations():
    return {name: {'arguments': list(args), 'returns': result,
                   'collective': collective, 'requires': requires,
                   'effect': effect, 'failure': failure,
                   'thread': 'main thread, sequential',
                   'buffer_lifetime': 'caller owns buffers; not retained after return'}
            for name, args, result, collective, requires, effect, failure in _HOST_OPERATIONS}


def _bind_host_operations(library):
    """Bind the complete common operation set without executing native calls.

    Resolve all symbols first so missing functions cannot leave partially bound
    declarations. A ctypes declaration is not proof of a binary's implementation;
    the authenticated build, host descriptor and host probes remain required.
    """
    import ctypes as c
    from .errors import RuntimeFailure
    functions = []
    for name, args, result, *_ in _HOST_OPERATIONS:
        function = getattr(library, name, None)
        if function is None:
            raise RuntimeFailure('runtime lacks required host operation ' + name)
        functions.append((function, args, result))
    types = {'int64': c.c_int64, 'int': c.c_int, 'double': c.c_double,
             'int*': c.POINTER(c.c_int), 'char*': c.POINTER(c.c_char), 'void': None}
    for function, args, result in functions:
        function.argtypes = [types[arg] for arg in args]
        function.restype = types[result]



def _state_operations(component, route):
    """Current required state entries only; optional diagnostics stay separate."""
    if component not in ('tddft', 'fast_hse'):
        return {}
    prefix = 'sp_tddft_' if component == 'tddft' else 'sp_fast_hse_'
    if route == 'single':
        specs = ({'init': ('string', 'double', 'int', 'vector3', 'decomposition5', 'vector3'),
                  'step': ('vector3', 'vector3', 'vector3'), 'close': ()}
                 if component == 'tddft' else
                 {'init': ('string', 'double', 'int', 'vector3'),
                  'step': ('vector3', 'vector3', 'vector3', 'double*'),
                  'checkpoint': ('string', 'int'),
                  'info': ('int*', 'vector3', 'vector3'), 'close': ()})
    else:
        specs = {'handle_protocol': (), 'model_protocol': (),
                 'handle_info': ('int64', 'int*', 'vector3', 'vector3'),
                 'handle_checkpoint': ('int64', 'string', 'int'),
                 'handle_transaction': ('int64', 'int'),
                 'handle_close': ('int64',), 'handle_label': ('int64', 'buffer', 'int')}
        if component == 'tddft':
            specs.update(handle_init=('string', 'double', 'int', 'vector3', 'decomposition5', 'int64*', 'vector3'),
                         model_new=('int64*', 'vector3'),
                         handle_step=('int64', 'vector3', 'vector3', 'vector3'))
        else:
            specs.update(handle_init=('string', 'double', 'int', 'int', 'int64*', 'vector3'),
                         model_new=('double', 'int', 'int64*', 'vector3'),
                         handle_step=('int64', 'vector3', 'vector3', 'vector3', 'double*'),
                         handle_history=('int64', 'history8'))
        if route == 'registry':
            specs.update(registry_protocol=(), model_close=('int64',))
            if component == 'tddft':
                specs.update(model_register=('string', 'double', 'int', 'vector3', 'decomposition5', 'int64*', 'int64*', 'vector3'),
                             model_spawn=('int64', 'int64*', 'vector3'))
            else:
                specs.update(model_register=('string', 'double', 'int', 'int', 'int64*', 'int64*', 'vector3'),
                             model_spawn=('int64', 'double', 'int', 'int64*', 'vector3'))
    report = {}
    for name, args in specs.items():
        query = name.endswith('_protocol')
        entry = {'arguments': list(args), 'returns': 'int',
                 'collective': 'local' if query else 'solver group',
                 'requires': 'none' if query else 'live attachment; matching call order; correct model I/O directory',
                 'buffer_lifetime': 'caller-owned; not retained after return',
                 'thread': 'main thread, sequential',
                 'failure': 'protocol version' if query else '0 success; nonzero operation-specific status; outputs valid only on success'}
        if name.endswith('checkpoint'):
            entry['modes'] = {0:'save',1:'fresh load',2:'stage',3:'commit',4:'cancel'}
            entry['path'] = 'NUL-terminated rank-local path; caller authenticates input/runtime/topology compatibility'
            entry['atomicity'] = 'stage/commit/cancel participate in caller-coordinated restore; raw load alone is not whole-session atomic'
        if name.endswith('transaction'):
            entry['modes'] = ({0:'capture',1:'rollback and close streams',2:'release',3:'reopen streams'}
                              if component == 'tddft' else
                              {0:'capture',1:'rollback and close streams',2:'reopen/flush/release; retain on failure',3:'reopen streams'})
            entry['atomicity'] = 'caller restores only owned stream bytes between rollback and reopen'
        if name.endswith('step'):
            entry['inputs'] = 'old and new Cartesian A/c endpoints, atomic units; old equals accepted endpoint'
            entry['outputs'] = ('SALMON matter J=-j' if component=='tddft' and route=='single' else 'charge-current density j')
            if component=='fast_hse':entry['energy']='total energy, Hartree'
            entry['atomicity'] = 'only documented numerical/transaction rejection paths preserve accepted state; not arbitrary MPI/process failures'
        if name in ('model_register','model_spawn','model_new','handle_init'):
            entry['ownership'] = 'returned positive state token is owning and generation-checked; close once; never copy as an independent state'
        if name == 'model_register':entry['model_token']='negative opaque generation token; distinct from state token'
        if name == 'model_close':entry['effect']='retire allocator; existing children retain model until last state closes'
        if name == 'handle_close':entry['effect']='consumes state token; stale token rejected even after slot reuse; does not detach MPI'
        if name == 'handle_info':entry['effect']='read accepted step, A/c and charge j without propagation'
        if name == 'handle_label':entry['capacity']='bytes including NUL; insufficient capacity returns13; no truncated successful label'
        names = {
            'init': ['path','dt','capacity']+(['initial_A','decomposition'] if component=='tddft' else [])+['j'],
            'step': ['old_A','new_A','j']+(['energy'] if component=='fast_hse' else []),
            'info': ['step','accepted_A','j'], 'checkpoint': ['path','mode'], 'close': [],
            'handle_init': ['path','dt','capacity']+(['initial_A','decomposition'] if component=='tddft' else ['allow_symmetry'])+['state','j'],
            'model_new': (['dt','capacity'] if component=='fast_hse' else [])+['state','j'],
            'model_register': ['path','dt','capacity']+(['initial_A','decomposition'] if component=='tddft' else ['allow_symmetry'])+['model','state','j'],
            'model_spawn': ['model']+(['dt','capacity'] if component=='fast_hse' else [])+['state','j'],
            'model_close': ['model'],
            'handle_step': ['state','old_A','new_A','j']+(['energy'] if component=='fast_hse' else []),
            'handle_info': ['state','step','accepted_A','j'], 'handle_history': ['state','history'],
            'handle_checkpoint': ['state','path','mode'], 'handle_transaction': ['state','mode'],
            'handle_close': ['state'], 'handle_label': ['state','label','capacity'],
        }.get(name, [])
        if component=='tddft' and route=='single':names=['J_matter' if n=='j' else n for n in names]
        if len(names)!=len(args):raise AssertionError('state parameter declaration differs: '+name)
        entry['parameters']=[{'name':parameter,'type':kind,
                              'direction':'out' if parameter in ('j','J_matter','energy','step','accepted_A','history','label') or kind=='int64*' else 'in'}
                             for parameter,kind in zip(names,args)]
        if 'decomposition' in names:entry['decomposition_order']=['k','orbital','grid_x','grid_y','grid_z']
        report[prefix+name]=entry
    return report


def _bind_state_operations(library, component, *, route):
    """Bind required current signatures atomically with respect to symbol lookup."""
    import ctypes as c
    import numpy as np
    from .errors import InputError
    specs=_state_operations(component,route)
    if not specs:raise InputError('no built-in state operation contract for '+str(component))
    functions=[]
    for name,entry in specs.items():
        function=getattr(library,name,None)
        if function is None:raise InputError('runtime lacks required state operation '+name)
        functions.append((function,entry['arguments']))
    types={'string':c.c_char_p,'double':c.c_double,'int':c.c_int,'int64':c.c_int64,
           'int*':c.POINTER(c.c_int),'int64*':c.POINTER(c.c_int64),
           'double*':c.POINTER(c.c_double),'buffer':c.c_void_p,
           'vector3':np.ctypeslib.ndpointer(dtype=np.float64,shape=(3,),flags='C_CONTIGUOUS'),
           'decomposition5':np.ctypeslib.ndpointer(dtype=np.int32,shape=(5,),flags='C_CONTIGUOUS'),
           'history8':np.ctypeslib.ndpointer(dtype=np.int32,shape=(8,),flags='C_CONTIGUOUS')}
    for function,args in functions:
        function.argtypes=[types[arg] for arg in args];function.restype=c.c_int


def native_contract(component, *, route):
    """Describe a built-in execution contract without loading a native library.

    Explicit names: tddft/fast_hse with single/shared/registry; fdtd with single.
    This is declared documentation, not a runtime capability check or a stable
    C ABI promise. Unknown/custom components return unverified, never a refusal
    to execute them. Their backend owns its separate contract.
    """
    report = {'schema': 'salmonpy.native-contract/1', 'component': component,
              'route': route, 'known': False, 'runtime_verified': False,
              'scope': 'declared interface only; no library loading or execution',
              'limitations': ['Runtime capabilities and checkpoint compatibility require separate queries.']}
    definitions = {
        'tddft': ('native_context.h', 'standard_tddft_handles.h', 'standard_tddft_models.h'),
        'fast_hse': ('native_context.h', 'fast_hse_handles.h',),
        'fdtd': ('native_context.h',)}
    if (component not in definitions or route not in ('single', 'shared', 'registry')
            or (component == 'fdtd' and route != 'single')):
        report['status'] = 'unverified'
        report['limitations'].append('Use the custom backend contract; this description imposes no execution restriction.')
        return report
    report.update(known=True, status='experimental', abi_stability='experimental',
                  host_abi_version=1,host_abi_query='sp_host_abi',
                  host_operations=_host_operations(),
                  state_operations=_state_operations(component,route),
                  host_thread='main thread, sequential calls',
                  host_owns_mpi_initialization_and_finalization=True,
                  attachment='one native attachment per process; close owned states before detach',
                  collective_order='collective operations must match on the solver communicator; query locality is entry-point specific',
                  multiple_live_states=route != 'single',
                  model_settings_owner=('single field lifetime' if component == 'fdtd' else
                                        'explicit single-state settings' if route == 'single' else
                                        'explicit model'),
                  headers=['salmonpy/native_adapters/'+name for name in definitions[component]],
                  failure_guarantee='Only declared numerical/transaction failures; not arbitrary STOP, allocation or MPI failure recovery')
    if component != 'fdtd':
        report['boundary'] = 'Cartesian A/c and physical charge-current density j, atomic units'
        report['native_current'] = ('SALMON matter J=-j; Python adapter returns charge j'
                                    if component == 'tddft' and route == 'single' else
                                    'physical charge-current density j')
    else:
        report['boundary'] = 'Native E/H and mapped A/current coupling; follow COUPLED time/normalization contract'
    return report


def query_host_abi(library):
    """Validate an explicitly loaded, trusted built-in library's current host ABI.

    Does not load libraries, initialize MPI, attach, compile or execute a solver.
    Unknown user backends keep their own contracts instead of using this query.
    Per-operation protocols and MPI-library identity are separate requirements.
    """
    import ctypes as c
    from .errors import RuntimeFailure
    probe=getattr(library,'sp_host_abi',None)
    if probe is None:raise RuntimeFailure('runtime lacks current host ABI; explicitly rebuild or use a matching runtime')
    probe.argtypes=[c.c_int,c.POINTER(c.c_int64)];probe.restype=c.c_int
    words=(c.c_int64*10)()
    status=probe(len(words),words)
    if status:raise RuntimeFailure('host ABI query failed with status '+str(status))
    values=tuple(words)
    expected=(1,32,64,64,32,64,1,1,1,3)
    if values!=expected:raise RuntimeFailure('unsupported host ABI descriptor: '+str(values))
    if (c.sizeof(c.c_int)*8,c.sizeof(c.c_double)*8,c.sizeof(c.c_int64)*8)!=values[1:4]:
        raise RuntimeFailure('Python and runtime host ABI type widths differ')
    return {'schema':'salmonpy.host-abi/1','version':values[0],
            'primitive_bits':{'c_int':values[1],'c_double':values[2],'int64':values[3],
                              'fortran_integer':values[4],'real8':values[5]},
            'communicator':'MPI_Fint widened to int64; never raw MPI_Comm',
            'MPI_ownership':'host','thread_policy':'main thread, sequential',
            'directory_protocol_required':values[9],'initializes_MPI':False,
            'scope':'host descriptor only; operation/profile/runtime checks remain separate'}
