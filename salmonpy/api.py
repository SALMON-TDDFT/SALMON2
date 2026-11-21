"""Explicit public names; importing this inventory does not import components.

Listed APIs are preview APIs covered by API.md/API.en.md, not a stable C ABI.
Returned state/result objects expose the documented methods and properties;
private names and incidental imports are not public APIs.
"""
from types import MappingProxyType

__all__ = ['PUBLIC_API', 'native_contract']
PUBLIC_API = MappingProxyType({'salmonpy.subpixel': ('SubpixelSmoothing','effective_material'), 'salmonpy.optical_materials': ('OpticalMaterial','list_materials','load_material','wavelength_to_frequency'),
 'salmonpy.analysis': ('PortSpectrum', 'absorbed_energy_density', 'reflection_transmission', 'Table', 'read_table', 'Spectrum', 'spectrum', 'OpticalSpectrum',
 'OpticalSensitivity', 'KCurrentContributions', 'optical_sensitivity', 'k_current_contributions', 'LongWavelengthCorrection', 'long_wavelength_correction', 'optical_conductivity', 'dielectric_function', 'dielectric_tensor', 'waveplate_diagnostic', 'optical_linearity', 'optical_sum_rule'),
 'salmonpy.api': ('PUBLIC_API', 'native_contract'),
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
 'salmonpy.fdtd': ('EnergyBalance', 'energy_balance', 'Boundaries',
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


def native_contract(component, *, route):
    """Describe a built-in execution contract without loading a native library.

    Explicit names: tddft/fast_hse with legacy/shared/registry; fdtd with legacy.
    This is declared documentation, not a runtime capability check or a stable
    C ABI promise. Unknown/custom components return unverified, never a refusal
    to execute them. Their backend owns its separate contract.
    """
    report = {'schema': 'salmonpy.native-contract/1', 'component': component,
              'route': route, 'known': False, 'runtime_verified': False,
              'scope': 'declared interface only; no library loading or execution',
              'limitations': ['Runtime capabilities and checkpoint compatibility require separate queries.']}
    definitions = {
        'tddft': ('standard_tddft_handles.h', 'standard_tddft_models.h'),
        'fast_hse': ('fast_hse_handles.h',),
        'fdtd': ()}
    if (component not in definitions or route not in ('legacy', 'shared', 'registry')
            or (component == 'fdtd' and route != 'legacy')):
        report['status'] = 'unverified'
        report['limitations'].append('Use the custom backend contract; this description imposes no execution restriction.')
        return report
    report.update(known=True, status='experimental', abi_stability='experimental',
                  host_thread='main thread, sequential calls',
                  host_owns_mpi_initialization_and_finalization=True,
                  attachment='one native attachment per process; close owned states before detach',
                  collective_order='collective operations must match on the solver communicator; query locality is entry-point specific',
                  multiple_live_states=route != 'legacy',
                  model_settings_owner=('single field lifetime' if component == 'fdtd' else
                                        'process compatibility settings' if route == 'legacy' else
                                        'explicit model'),
                  headers=['salmonpy/native_adapters/'+name for name in definitions[component]],
                  failure_guarantee='Only declared numerical/transaction failures; not arbitrary STOP, allocation or MPI failure recovery')
    if component != 'fdtd':
        report['boundary'] = 'Cartesian A/c and physical charge-current density j, atomic units'
        report['native_current'] = ('SALMON matter J=-j; Python adapter returns charge j'
                                    if component == 'tddft' and route == 'legacy' else
                                    'physical charge-current density j')
    else:
        report['boundary'] = 'Native E/H and mapped A/current coupling; follow COUPLED time/normalization contract'
    return report
