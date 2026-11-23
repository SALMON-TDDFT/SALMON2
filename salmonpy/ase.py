"""Optional ASE bridge for structures, bands and A-driven TDDFT initial states.

ASE and f90nml are imported on use; native execution stays in the physics modules.
"""
from .errors import InputError, InputTypeError, RuntimeFailure
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from .bands import BandInitialState, BandStructure, ExplicitKPoints, HARTREE_EV, _input_groups

# SALMON inputoutput.f90 length conversion; do not silently use ASE's CODATA default.
BOHR_ANGSTROM = 0.52917721067


def _ase_types():
    try:
        from ase import Atoms
        from ase.dft.kpoints import BandPath
        from ase.spectrum.band_structure import BandStructure as ASEBandStructure
    except ImportError as exc:
        raise ImportError('ASE bridge requires the optional ase dependency; install ase and f90nml') from exc
    return Atoms, BandPath, ASEBandStructure


def _copy_path(path):
    _,PathType,_=_ase_types()
    if not isinstance(path,PathType): raise InputTypeError('expected ase.dft.kpoints.BandPath')
    if path.cell.rank!=3 or not np.isfinite(path.cell.array).all(): raise InputError('full periodic path cell required')
    return PathType(cell=path.cell.copy(),kpts=np.array(path.kpts,copy=True),
                    special_points={k:np.array(v,copy=True) for k,v in path.special_points.items()},path=path.path)


def explicit_path(path):
    """Copy every ASE fractional k point, preserving order and discontinuities."""
    return ExplicitKPoints(_copy_path(path).kpts)


def atoms_from_initial(initial):
    """Read original GS input geometry into a new ASE Atoms (Å).

    Does not infer binary checkpoint provenance. Supports one inline or relative
    file coordinate source, Cartesian or reduced, and SALMON au/A_eV_fs units.
    """
    Atoms,_,_=_ase_types()
    try:
        import f90nml
    except ImportError as exc:
        raise ImportError('GS geometry loading requires optional f90nml') from exc
    if not initial._snapshot.exists(): raise RuntimeFailure('initial state is closed')
    groups=_input_groups(initial._input)
    nml=f90nml.reads('\n'.join('&'+name+body+'/' for name,body in groups.items()
                             if name not in ('atomic_coor','atomic_red_coor')))
    system=nml.get('system',{})
    if system.get('yn_periodic')!='y' or system.get('spin','unpolarized')!='unpolarized':
        raise InputError('periodic unpolarized initial state required')
    units=nml.get('units',{}).get('unit_system','au')
    if units not in ('au','a.u.','A_eV_fs'): raise InputError('unsupported unit_system')
    scale=BOHR_ANGSTROM if units in ('au','a.u.') else 1.
    vectors=np.array([system.get('al_vec'+str(i),[0.,0.,0.]) for i in range(1,4)],dtype=float)
    if vectors.shape!=(3,3): raise InputError('three complete lattice vectors required')
    if np.any(vectors-np.diag(np.diag(vectors))):
        cell=vectors
    else:
        al=np.array(system.get('al',[0.,0.,0.]),dtype=float)
        if al.shape!=(3,): raise InputError('complete al vector required')
        if np.prod(al)==0:
            if np.prod(np.diag(vectors))!=0: al=np.diag(vectors)
            else:
                grid=nml.get('rgrid',{})
                al=np.array(grid.get('dl',[0.,0.,0.]))*np.array(grid.get('num_rgrid',[0,0,0]))
        cell=np.diag(al)
    cell=cell*scale
    if not np.isfinite(cell).all() or abs(np.linalg.det(cell))<1e-12: raise InputError('invalid cell')
    sources=[]
    for group,field,reduced in [('atomic_coor','file_atom_coor',False),('atomic_red_coor','file_atom_red_coor',True)]:
        if group in groups: sources.append((groups[group],reduced))
        filename=system.get(field,'none')
        if filename!='none':
            p=Path(filename)
            if p.is_absolute() or '..' in p.parts or '\\' in filename: raise InputError('coordinate path must be bundle-relative')
            sources.append(((initial._snapshot/p).read_text(),reduced))
    if len(sources)!=1: raise InputError('exactly one explicit GS coordinate source required')
    import shlex
    positions=[];species=[]
    pseudo=nml.get('pseudo',{});numbers=pseudo.get('izatom',[])
    if not isinstance(numbers,list): numbers=[numbers]
    start=pseudo.start_index.get('izatom',[1])[0]
    for line in sources[0][0].splitlines():
        lexer=shlex.shlex(line,posix=True);lexer.whitespace_split=True;lexer.commenters='!'
        fields=list(lexer)
        if not fields: continue
        if len(fields)!=5: raise InputError('expected atom label, three coordinates and species index')
        try:
            positions.append([float(x.replace('d','e').replace('D','E')) for x in fields[1:4]])
            idx=int(fields[4])-start
            if idx<0 or idx>=len(numbers) or numbers[idx] is None: raise InputError('undefined species')
            z=int(numbers[idx])
            if z!=numbers[idx] or not 1<=z<=118: raise InputError('invalid atomic number')
            species.append(z)
        except (ValueError,TypeError,IndexError) as exc:
            raise InputError('invalid atomic coordinate/species record') from exc
    if len(species)!=system.get('natom') or not species: raise InputError('natom does not match coordinate records')
    positions=np.array(positions)
    if not np.isfinite(positions).all(): raise InputError('nonfinite atomic coordinates')
    kwargs={'scaled_positions':positions} if sources[0][1] else {'positions':positions*scale}
    return Atoms(numbers=species,cell=cell,pbc=True,**kwargs)


def validate_structure(initial,atoms,*,tolerance_angstrom=1e-8):
    """Require same ordered species, cell and positions (modulo lattice translations)."""
    Atoms,_,_=_ase_types()
    if not isinstance(atoms,Atoms): raise InputTypeError('expected ASE Atoms')
    if not np.isfinite(tolerance_angstrom) or not 0<tolerance_angstrom<=1e-5:
        raise InputError('tolerance_angstrom must be in (0, 1e-5]')
    reference=atoms_from_initial(initial)
    if not np.all(atoms.pbc): raise InputError('all three periodic boundaries required')
    if not np.array_equal(atoms.numbers,reference.numbers): raise InputError('ordered species differ from GS')
    if not np.allclose(atoms.cell.array,reference.cell.array,rtol=0,atol=tolerance_angstrom):
        raise InputError('cell differs from GS; prepare a matching ground state')
    delta=(atoms.positions-reference.positions) @ np.linalg.inv(reference.cell.array)
    delta-=np.rint(delta)
    if not np.isfinite(delta).all() or np.max(np.abs(delta @ reference.cell.array))>tolerance_angstrom:
        raise InputError('atomic positions differ from GS; prepare a matching ground state')


def to_ase_band_structure(bands,path,*,reference_ev):
    """Export copies in ASE (spin,k,band) eV convention; keep energy zero explicit."""
    _,_,ASEBandStructure=_ase_types()
    if not isinstance(bands,BandStructure): raise InputTypeError('expected SALMONpy BandStructure')
    if not np.isscalar(reference_ev) or not np.isfinite(reference_ev): raise InputError('finite reference_ev required')
    path=_copy_path(path)
    # Express actual Cartesian k in the supplied ASE basis. Propagate the same
    # 1e-8 Å componentwise cell tolerance used by preflight, rather than imposing
    # a tighter reciprocal tolerance after an expensive native calculation.
    actual_fractional=bands.k_cartesian @ path.cell.array.T/(2*np.pi*BOHR_ANGSTROM)
    error_bound=1e-8*np.sum(np.abs(bands.k_cartesian),axis=1)[:,None]/(2*np.pi*BOHR_ANGSTROM)+1e-12
    if actual_fractional.shape!=path.kpts.shape or not np.all(np.abs(actual_fractional-path.kpts)<=error_bound):
        raise InputError('ASE cell/k points differ from computed Cartesian k points')
    if bands.energies.shape[2] not in (1,2): raise InputError('ASE supports one or two spin channels')
    return ASEBandStructure(path,bands.energies.transpose(2,0,1).copy()*HARTREE_EV,reference=float(reference_ev))


@dataclass(frozen=True)
class ASEBandRun:
    native: object
    band_structure: object

    def plot(self,**kwargs):
        return self.band_structure.plot(**kwargs)


def run_bands(initial,atoms,path,directory,*,reference_ev,**run_options):
    """Validate ASE inputs, run independent native bands, save an ASE JSON result.

    run_options are BandInitialState.run options (bands, tolerance, etc.). No
    Calculator is attached, and neither input Atoms nor path is changed.
    """
    if not isinstance(initial,BandInitialState): raise InputTypeError('expected BandInitialState')
    if not np.isscalar(reference_ev) or not np.isfinite(reference_ev): raise InputError('finite reference_ev required')
    atoms=atoms.copy();path=_copy_path(path)
    validate_structure(initial,atoms)
    if not np.allclose(path.cell.array,atoms_from_initial(initial).cell.array,rtol=0,atol=1e-8):
        raise InputError('path cell differs from GS; no implicit cell transformation')
    native=initial.run(explicit_path(path),directory,**run_options)
    result=to_ase_band_structure(native.bands,path,reference_ev=reference_ev)
    result.write(str(native.directory/'ase-bands.json'))
    from ase.io.jsonio import encode
    (native.directory/'ase-input.json').write_text(encode({'atoms':atoms,'path':path,'reference_ev':float(reference_ev),
         'bohr_angstrom':BOHR_ANGSTROM,'hartree_ev':HARTREE_EV}))
    return ASEBandRun(native,result)


def dft_from_atoms(atoms,*,pseudopotentials,states,kmesh,grid=None,
                   grid_spacing=None,grid_spacing_unit=None,electrons=None,xc='pz',cval=None):
    """Copy ASE geometry into an independent DFT model; no calculation is launched.

    ASE cell/positions use Å. Electrons default to the neutral pseudopotential
    valence sum; an explicit count must match. Choose grid counts or grid_spacing with an explicit unit. States and kmesh
    remain explicit; grid rounding is not a convergence check. Charged/magnetic states are rejected.
    """
    from .dft import DFT, Pseudopotential
    from .sbe import _integer
    Atoms,_,_=_ase_types()
    if not isinstance(atoms,Atoms): raise InputTypeError('expected ASE Atoms')
    from .structure import _validate_occupancy
    _validate_occupancy(atoms)
    atoms=atoms.copy()
    if not np.all(atoms.pbc): raise InputError('three-dimensional periodic boundary conditions required')
    if grid_spacing is not None:
        if grid is not None:
            raise InputError('specify grid or grid_spacing, not both')
        from .structure import grid_from_spacing
        grid = grid_from_spacing(atoms.cell.array, spacing=grid_spacing,
                                 unit=grid_spacing_unit).counts
    elif grid is None or grid_spacing_unit is not None:
        raise InputError('specify grid counts or grid_spacing with grid_spacing_unit')

    if np.any(atoms.get_initial_magnetic_moments()!=0) or np.any(atoms.get_initial_charges()!=0):
        raise InputError('ASE initial charge/spin settings are unsupported; neutral closed-shell model required')
    pseudopotentials = tuple(pseudopotentials)
    if any(not isinstance(p, Pseudopotential) for p in pseudopotentials):
        raise InputTypeError('expected Pseudopotential objects')
    by_number = {p.atomic_number: p for p in pseudopotentials}
    if len(by_number) != len(pseudopotentials) or set(by_number) != set(atoms.numbers):
        raise InputError('exactly one pseudopotential per atomic species required')
    valences = {z: p.valence_electrons for z, p in by_number.items()}
    neutral = sum(valences[z] for z in atoms.numbers)
    if electrons is not None and _integer(electrons, 'electrons', 2, 2**31-1) != neutral:
        raise InputError(f'electrons must equal neutral pseudopotential valence sum ({neutral}); charged ASE systems are unsupported')
    electrons = neutral
    return DFT(atoms.numbers,atoms.cell.array/BOHR_ANGSTROM,atoms.get_scaled_positions(wrap=False),
               pseudopotentials,electrons=electrons,states=states,grid=grid,kmesh=kmesh,xc=xc,cval=cval)


def tddft_from_atoms(atoms, directory, *, executable, pseudopotentials, states,
                     kmesh, grid=None, grid_spacing=None, grid_spacing_unit=None,
                     electrons=None, xc='pz', scf=None,
                     max_steps=1000, threads=1, timeout=600.):
    """Run GS in a new directory and return a context-managed TDDFT initial state.

    Uses the existing DFT/GS loader without changing propagation. The audited
    executable must support both DFT and the TDDFT worker. GS files remain in
    directory, also on failure. Each initial.new_state creates an independent
    A/c-driven calculation. No runtime acquisition or compilation is implicit.
    """
    from .dft import SCFSettings
    from .sbe import _integer
    from .ld import _positive
    max_steps = _integer(max_steps, 'max_steps', 1, 2**31-2)
    threads = _integer(threads, 'threads', 1, 2**31-1)
    timeout = _positive(timeout, 'timeout')
    scf = SCFSettings() if scf is None else scf
    model = dft_from_atoms(atoms, pseudopotentials=pseudopotentials, states=states,
                           grid=grid, grid_spacing=grid_spacing, grid_spacing_unit=grid_spacing_unit,
                           kmesh=kmesh, electrons=electrons, xc=xc)
    gs = model.run(executable, directory, scf=scf, threads=threads, timeout=timeout)
    return gs.tddft_initial_state(executable=executable, max_steps=max_steps,
                                  threads=threads, timeout=timeout)
