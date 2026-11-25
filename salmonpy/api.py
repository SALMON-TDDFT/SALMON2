"""Explicit public names; importing this inventory does not import components.

Listed APIs are preview APIs covered by API.md/API.en.md, not a stable C ABI.
Returned state/result objects expose the documented methods and properties;
private names and incidental imports are not public APIs.
"""
from types import MappingProxyType

__all__ = ['PUBLIC_API', 'native_contract', 'query_host_abi', 'unit_contract']
PUBLIC_API = MappingProxyType({'salmonpy.quickstart': ('run',), 'salmonpy.subpixel': ('SubpixelSmoothing','effective_material'), 'salmonpy.optical_materials': ('OpticalMaterial','list_materials','load_material','wavelength_to_frequency'),
 'salmonpy.analysis': ('SurfaceFlux', 'surface_flux', 'PortSpectrum', 'absorbed_energy_density', 'reflection_transmission', 'Table', 'read_table', 'Spectrum', 'spectrum', 'OpticalSpectrum',
 'OpticalSensitivity', 'KCurrentContributions', 'optical_sensitivity', 'k_current_contributions', 'LongWavelengthCorrection', 'long_wavelength_correction', 'optical_conductivity', 'dielectric_function', 'dielectric_tensor', 'waveplate_diagnostic', 'optical_linearity', 'optical_sum_rule'),
 'salmonpy.api': ('PUBLIC_API', 'native_contract', 'query_host_abi', 'unit_contract'),
 'salmonpy.manual': ('load','search','task_details'),
 'salmonpy.doctor': ('diagnose', 'probe_native_protocols'),
 'salmonpy.convergence': ('plan_grid_scan','run_grid_scan','plan_dt_scan','run_dt_scan','run_current_grid_scan'),
 'salmonpy.gs_input': ('GSInputMetadata','describe_gs_input'),
 'salmonpy.native_steps': ('native_step_targets',),
 'salmonpy.stepping': ('StepTarget','StepRecoveryError','guarded_step'),
 'salmonpy.memory': ('MemoryEstimate','MemoryReport','MemoryBudgetWarning','process_memory','memory_advice'),
 'salmonpy.checkpointing': ('CheckpointCapability', 'check_checkpoint_compatibility', 'CheckpointTarget', 'execute_checkpoint', 'CheckpointState', 'restore_checkpoint', 'AtomicRestoreError', 'CheckpointSchedule'),
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
 'salmonpy.pseudopotentials': ('PseudopotentialReference', 'silicon_reference', 'list_references', 'catalog_reference', 'download'),
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
 'salmonpy.isolated_tddft': ('IsolatedRuntime', 'load_isolated_runtime', 'DipoleObservables', 'ChargeDensitySamples', 'OrbitalMagneticResponse', 'IsolatedTDDFTInitialState', 'IsolatedTDDFTState', 'IsolatedSavepoint'),
 'salmonpy.dipole_source': ('DipoleMap', 'drive_native_dipoles', 'MolecularDrive', 'select_molecular_field', 'MolecularFrame', 'MolecularFeedback'),
 'salmonpy.near_field': ('ElectricPotentialSeries','NearFieldSeries','NearFieldFrame','couple_path_samples','GridFieldSeries'),
 'salmonpy.field_sampling': ('SampledField', 'sample_packet'),
 'salmonpy.field_input': ('FieldValue', 'FieldPacket', 'FieldContract', 'FieldAdapter'),
 'salmonpy.fields': ('FieldSamples', 'AcosPulse', 'sample_pulses'),
 'salmonpy.gs_native': ('GammaRestartConverter', 'SerialRestartConverter'),
 'salmonpy.gs_initial_guess': ('GSInitialGuess',),
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
 'salmonpy.outputs': ('SpinDynamics', 'read_spin_dynamics', 'CubeGrid',
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



def _field_operations():
    specs = {
        'init': [('path','string','in'),('points','int','in'),('indices','indices','in')],
        'step': [('points','int','in'),('J','matrix','in'),('E','matrix','out')],
        'close': [],
        'current_mask': [('points','int','in'),('mask','mask','out')],
        'step_owned': [('points','int','in'),('local_points','int','in'),('J','array','in'),('A','array','in'),('E','array','out')],
        'step_performance': [('points','int','in'),('local_points','int','in'),('owned','int','in'),('J','array','in'),('old_A','array','in'),('A','array','out'),('E','array','out')],
        'step_checked': [('points','int','in'),('J','array','in'),('A','array','in'),('E','array','out')],
        'step_edges': [('points','int','in'),('edges','int','in'),('indices','indices','in'),('J','matrix','in'),('A','matrix','in'),('E','matrix','out')],
        'edges_read': [('edges','int','in'),('indices','indices','in'),('E','matrix','out')],
        'fields_read': [('edges','int','in'),('indices','indices','in'),('fields','array','out')],
        'edge_geometry': [('counts','mask','out'),('origins','matrix','out'),('spacing','array','out')],
        'sample_h': [('points','int','in'),('H','matrix','out')],
        'box_energy_step': [('box','indices','in'),('values','buffer','out')],
        'plane_read': [('axis','int','in'),('offset','int','in'),('values_count','int','in'),('info','buffer','out'),('values','buffer','out')],
        'step_indices': [('values','mask','out')],
        'checkpoint_protocol': [], 'step_transaction_protocol': [],
        'checkpoint': [('path','string','in'),('mode','int','in')],
        'step_transaction': [('mode','int','in')],
    }
    result = {}
    for name,parameters in specs.items():
        query=name.endswith('_protocol')
        local=query or name in ('current_mask','plane_read','edges_read','fields_read','edge_geometry')
        entry={'parameters':[{'name':n,'type':t,'direction':d} for n,t,d in parameters],
               'arguments':[t for _,t,_ in parameters], 'returns':'int',
               'collective':'local' if local else 'solver group',
               'required_at_initialization':name in ('init','step','close'),
               'requires':'none' if query else 'main thread, sequential; live field attachment; matching collective order',
               'ownership':'caller buffers, not retained; one field lifetime per process',
               'failure':'protocol version' if query else '0 success; nonzero operation-specific status; consume outputs only on success'}
        if name in ('init','step','step_checked','step_owned','sample_h'):
            entry['layout']='C contiguous (points,3), equivalent to Fortran (3,points); int32 indices, float64 fields'
        if name=='init':entry['indices']='one-based global grid triples in requested point order; Python converts zero-based indices'
        if name in ('step','step_checked','step_owned'):
            entry['current']='physical charge-current density j in atomic units (C parameter named J); previous half-step; no volume multiplier or sign conversion'
            entry['atomicity']='use declared protected-step transaction/coordinator; raw call is not universal rollback'
        if name=='step_edges':
            entry['routing']='rank-local lists: only locally owned cells consumed; different list sizes allowed; duplicates add'
            entry['layout']='indices int32 C(edges,3), zero-based global offsets; J float64 C(edges,3) on E components, no averaging/volume/sign conversion'
            entry['atomicity']='uses the existing checked E/LD/H update and recovery gates'
        if name=='edges_read':entry['layout']='same zero-based edge offsets; only owned E components returned, zero elsewhere; local, read-only'
        if name=='fields_read':entry['layout']='float64 C(2,N,3), raw E and half-step-centered H on distinct Yee component grids; only owned values, zero elsewhere; local read-only'
        if name=='edge_geometry':entry['layout']='global counts int32[3], E-component origins float64 C(3,3) [component,axis], spacing float64[3]; bohr'
        if name=='step_owned':entry['routing']='J rows follow ascending current-mask rows; local count differs by field rank; A/E retain all mapped rows'
        if name=='step_performance':
            entry['requires']='performance mode; caller guarantees matching collective order, current ownership and valid data; no rollback/retry guarantee'
            entry['layout']='float64 C(points,3); compact owned J rows when owned=1; old_A input, next A/E outputs; no finite scans or validation collectives'
        if name=='step_indices':entry['layout']='int32[2]: accepted step and staged candidate step; -1 means absent; read-only'
        if name=='current_mask':entry['layout']='int32[points]; rank-local constitutive ownership, overlaps permitted'
        if name=='box_energy_step':
            entry['requires']='accepted step with live unconsumed protected-step backup; identical zero-based inclusive box on the field group; read-only, no I/O'
            entry['layout']='box int32 C(2,3); output float64[19]; modified/physical energy, LD/source work, six-face Poynting flux, separate constitutive exchange and Ohmic work'
        if name=='sample_h':entry['clock']='raw macro-centered H at (n+1/2)dt; accepted step >=1; no time advance'
        if name=='plane_read':
            entry['layout']='info int32[21]; values float64[values_count]; first call count0 queries metadata; axis1..3, local-grid offset'
            entry['requires']='committed usable field; no pending restore; sampled step; no native MPI communication'
        if name=='checkpoint':
            entry['modes']={0:'save',1:'fresh load',2:'stage',3:'commit',4:'cancel'}
            entry['ownership']='distinct path per rank; caller validates input/runtime/topology and publishes multi-file save atomically'
            entry['atomicity']='coordinated stage/vote/commit/cancel required for whole-session protected restore'
        if name=='step_transaction':
            entry['modes']={0:'capture',1:'rollback',2:'release',3:'no-op'}
            entry['scope']='supported field/LD numerical state only; unsupported output modes return15; no arbitrary MPI/process/allocator recovery'
        result['sp_eh_'+name]=entry
    return result


def _bind_field_operations(library, names=('sp_eh_init','sp_eh_step','sp_eh_close')):
    """Bind requested FDTD entries, leaving optional operations opt-in."""
    import ctypes as c
    import numpy as np
    from .errors import RuntimeFailure
    specs=_field_operations();functions=[]
    for name in names:
        if name not in specs:raise RuntimeFailure('unknown built-in field operation '+str(name))
        fn=getattr(library,name,None)
        if fn is None:raise RuntimeFailure('runtime lacks required field operation '+name)
        functions.append((fn,specs[name]['arguments']))
    types={'string':c.c_char_p,'int':c.c_int,'buffer':c.c_void_p,
           'indices':np.ctypeslib.ndpointer(dtype=np.int32,ndim=2,flags='C_CONTIGUOUS'),
           'mask':np.ctypeslib.ndpointer(dtype=np.int32,ndim=1,flags='C_CONTIGUOUS'),
           'matrix':np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS'),
           'array':np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')}
    for fn,args in functions:
        fn.argtypes=[types[a] for a in args];fn.restype=c.c_int



_GS_CORE=('protocol','create','info','array_info','read','close','resources')
_GS_SCF=('scf_protocol','start','continue','scf_info','scf_snapshot','scf_restore','scf_snapshot_close')


def _gs_operations():
    specs={
        'protocol':[], 'create':[('path','string'),('state','int64*')],
        'array_guess_protocol':[], 'create_from_arrays':[('path','string'),('steps','int'),('reader','buffer'),('cookie','buffer'),('state','int64*')],
        'guess_protocol':[], 'create_from_gs':[('path','string'),('source','int64'),('steps','int'),('state','int64*')],
        'import':[('path','string'),('state','int64*')],
        'info':[('state','int64'),('stats','buffer'),('values','buffer')],
        'array_info':[('state','int64'),('quantity','int'),('ndim','int*'),('shape','buffer'),('lower','buffer')],
        'read':[('state','int64'),('quantity','int'),('real','buffer'),('imag','buffer'),('length','int64')],
        'close':[('state','int64')], 'resources':[('counts','buffer')],
        'scf_protocol':[], 'start':[('path','string'),('steps','int'),('state','int64*')],
        'continue':[('state','int64'),('steps','int')],
        'scf_spin_field':[('state','int64'),('field','buffer')],
        'scf_spin_info':[('state','int64'),('values','buffer')],
        'scf_adaptive_mixing':[('state','int64'),('enabled','int'),('min_rate','double'),('growth','double'),('patience','int'),('cooldown','int')],
        'scf_mixing_info':[('state','int64'),('values','buffer'),('counts','buffer')],
        'scf_set_mixing':[('state','int64'),('rate','double')],
        'scf_extend_limit':[('state','int64'),('new_limit','int')], 'scf_info':[('state','int64'),('values','buffer')],
        'scf_snapshot':[('state','int64'),('snapshot','int64*')],
        'scf_restore':[('state','int64'),('snapshot','int64')],
        'scf_snapshot_close':[('state','int64'),('snapshot','int64')],
        'scf_bootstrap':[('path','string'),('state','int64*')],
        'scf_disk':[('state','int64'),('path','string'),('mode','int')],
        'scf_candidate':[('state','int64'),('snapshot','int64'),('mode','int')],
        'scf_format_signature':[('kinds','buffer'),('bytes','buffer')],
    }
    for axis in ('k','o','r','m'):
        specs['scf_'+axis+'_repartition']=[('state','int64'),('directory','string'),('old_size','int'),('bounds','buffer')]
    result={}
    for name,parameters in specs.items():
        local=name in ('protocol','guess_protocol','array_guess_protocol','scf_protocol','resources','scf_format_signature')
        entry={'arguments':[t for _,t in parameters], 'returns':'int',
               'parameters':[{'name':n,'type':t,'direction':'out' if n in ('stats','values','ndim','shape','lower','real','imag','counts','kinds','bytes') or t=='int64*' else 'in'} for n,t in parameters],
               'collective':'local' if local else 'solver group',
               'requires':'main thread, sequential; live attachment and matching token/order for state calls',
               'ownership':'caller owns buffers; not retained; state/snapshot tokens own resources until explicit close',
               'failure':'query version' if name.endswith('protocol') else '0 success; nonzero operation-specific status; outputs valid only on success'}
        if name in ('create','start','import','scf_bootstrap'):
            entry['input']='immutable input/assets; call inside selected model I/O directory'
            entry['effect']={'create':'initialize and run SCF; retain last iterate even if unconverged',
                             'start':'initialize once, run requested main-iteration budget, retain mixing/history',
                             'import':'read existing GS; no SCF or GS writes; convergence/forces unavailable',
                             'scf_bootstrap':'prepare static resources for restart, no initial orbitals/density/mixing or SCF iterations'}[name]
        if name=='create_from_gs':entry['effect']='copy wavefunctions/occupations by value; fresh density/mixing/HSE; steps0 full SCF,positive budget; same session/layout/model except positions; no source mutation or GS I/O'
        if name=='create_from_arrays':entry['effect']='synchronous rank-local array reader into a fresh SCF owner; callback bounds/shape and complex128/float64 layout are in standard_gs.h; steps0 full SCF,positive budget; numerical finite/sampling validation in Fortran; caller authenticates declared model and owns reader through return; no pointers retained or intermediate GS files'
        if name=='info':entry['layout']='stats int32[4]: iterations/converged/atoms/bands; values float64[3]: Hartree energy/residual/tolerance; import uses unknown convergence and profile-specific absent values'
        if name=='array_info':entry['layout']='ndim int32; shape/lower int64[7]; quantity1 density,2 occupation,3 eigenvalue,4 force,5 orbitals,6 SCF history'
        if name=='read':entry['layout']='float64 real/imag arrays of length; Fortran order; use array_info first; local arrays can include halos'
        if name=='resources':entry['layout']='int32[3]: states/live child communicators/slot high-water; rank-local'
        if name=='scf_extend_limit':entry['scope']='increase only the owned SCF total-iteration cap; no SCF, wave/density/mixing reset or input-model mutation; serialized gs_input retains the new runtime cap'
        if name=='scf_info':entry['layout']='int32[6]: paused/iterations/total cap/initializations/executed main iterations/live snapshots'
        if name in ('scf_snapshot','scf_restore','scf_snapshot_close'):
            entry['scope']='same live GS owner; numerical/mixing/HSE history only; no portable files or text-output rewind; restore does not rewind lifetime execution/init counters'
        if name=='scf_candidate':
            entry['modes']={2:'stage snapshot',3:'publish',4:'discard'}
            entry['scope']='one pending numerical payload; snapshot0 for publish/discard; host votes before publication; no arbitrary allocation/process/MPI rollback guarantee'
        if name=='scf_disk':
            entry['modes']={0:'save new rank-local file',2:'stage disk read'}
            entry['scope']='host authenticates current runtime/input/topology; publish with candidate(state,0,3); no multi-file atomic publication'
        if name.endswith('repartition'):
            entry['bounds']='old ranks in order; k/orbital2, spatial6, mixed10 int32 values per old rank; explicit layout-specific capability query required'
        if name=='scf_format_signature':entry['layout']='int32 kinds[8], int8 bytes[128]; query only, no GS initialization'
        result['sp_gs_'+name]=entry
    return result


def _bind_gs_operations(library,names):
    """Resolve requested GS entries before binding; optional APIs remain opt-in."""
    import ctypes as c
    from .errors import InputError
    specs=_gs_operations();functions=[]
    for name in names:
        if name not in specs:raise InputError('unknown built-in GS operation '+str(name))
        fn=getattr(library,name,None)
        if fn is None:raise InputError('runtime lacks required GS operation '+name)
        functions.append((fn,specs[name]['arguments']))
    types={'double':c.c_double,'string':c.c_char_p,'int':c.c_int,'int64':c.c_int64,'buffer':c.c_void_p,
           'int*':c.POINTER(c.c_int),'int64*':c.POINTER(c.c_int64)}
    for fn,args in functions:
        fn.argtypes=[types[a] for a in args];fn.restype=c.c_int


def native_contract(component, *, route):
    """Describe a built-in execution contract without loading a native library.

    Explicit names: tddft/fast_hse with single/shared/registry; fdtd/gs with single.
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
        'fdtd': ('native_context.h','native_fdtd.h'),
        'gs': ('native_context.h','standard_gs.h')}
    if (component not in definitions or route not in ('single', 'shared', 'registry')
            or (component in ('fdtd','gs') and route != 'single')):
        report['status'] = 'unverified'
        report['limitations'].append('Use the custom backend contract; this description imposes no execution restriction.')
        return report
    report.update(known=True, status='experimental', abi_stability='experimental',
                  host_abi_version=1,host_abi_query='sp_host_abi',
                  host_operations=_host_operations(),
                  state_operations=_state_operations(component,route),
                  field_operations=_field_operations() if component=='fdtd' else {},
                  gs_operations=_gs_operations() if component=='gs' else {},
                  host_thread='main thread, sequential calls',
                  host_owns_mpi_initialization_and_finalization=True,
                  attachment='one native attachment per process; close owned states before detach',
                  collective_order='collective operations must match on the solver communicator; query locality is entry-point specific',
                  multiple_live_states=component=='gs' or route != 'single',
                  model_settings_owner=('single field lifetime' if component == 'fdtd' else
                                        'explicit GS state' if component=='gs' else
                                        'explicit single-state settings' if route == 'single' else
                                        'explicit model'),
                  headers=['salmonpy/native_adapters/'+name for name in definitions[component]],
                  failure_guarantee='Only declared numerical/transaction failures; not arbitrary STOP, allocation or MPI failure recovery')
    if component=='gs':
        report['boundary']='atomic units; Hartree energy, Cartesian forces; array_info defines local Fortran layout'
    elif component != 'fdtd':
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


def unit_contract():
    """Describe existing conversion conventions; never changes numeric inputs.

    The legacy solver input constants and current Python presentation constants
    deliberately remain distinct to preserve established input/results.
    """
    from .bands import HARTREE_EV
    from .structure import _BOHR_ANGSTROM
    return {'schema':'salmonpy.units/1',
            'boundary':'atomic units; no implicit conversion',
            'solver_input':{'hartree_eV':27.21138505,'origin':'salmon_global au_energy_ev; retained input convention'},
            'geometry':{'bohr_angstrom':_BOHR_ANGSTROM,'origin':'SALMON geometry convention; shared by ASE adapter'},
            'output':{'hartree_eV':HARTREE_EV,'origin':'phys_constants/CODATA energy presentation'},
            'forces':'Hartree/bohr at boundary; ASE eV/angstrom uses output hartree_eV / geometry bohr_angstrom',
            'scope':'declared current conventions, not inspection of a custom runtime or an automatic unit conversion'}
