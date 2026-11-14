"""Strict DFT-generated serial Gamma GS/1→GS/2 adapter; no solver execution."""
import hashlib,json,os,re,shutil,subprocess,tempfile
from pathlib import Path
import numpy as np
from .errors import InputError,InputTypeError,RuntimeFailure
from .gs_native import _digest,_environment


def _eigenvalues(text,nb,nk=1):
    lines=[line.strip() for line in text.splitlines() if line.strip()]
    if (len(lines)!=3+nk*(nb+1) or lines[:3]!=['#esp: single-particle energies (eigen energies)',
        '#occ: occupation numbers, io: orbital index','# 1:io, 2:esp[a.u.], 3:occ']):
        raise InputError('unpolarized a.u. eigenvalue table required')
    eigen=np.empty((nk,nb));occ=np.empty((nk,nb))
    try:
        for k in range(nk):
            offset=3+k*(nb+1)
            if re.fullmatch(r'k=\s*'+str(k+1)+r',\s*spin=\s*1',lines[offset]) is None:raise InputError('invalid eigenvalue k/spin order')
            for b,line in enumerate(lines[offset+1:offset+1+nb]):
                parts=line.split()
                if len(parts)!=3 or int(parts[0])!=b+1:raise InputError('invalid eigenvalue band index')
                eigen[k,b]=float(parts[1].replace('D','e').replace('d','e'))
                occ[k,b]=float(parts[2].replace('D','e').replace('d','e'))
    except (ValueError,OverflowError) as exc:raise InputError('invalid eigenvalue numbers') from exc
    if not np.isfinite(eigen).all() or not np.isfinite(occ).all():raise InputError('nonfinite eigenvalue table')
    return eigen,occ.ravel()


def _kpoints(text,nk,cell):
    lines=[s.strip() for s in text.splitlines() if s.strip()]
    if (len(lines)!=nk+7 or lines[0]!='# k-point distribution' or 'Reduced coordinate' not in lines[2]
            or '[none]' not in lines[4] or not lines[-2].startswith('# coefficients (2*pi/a [a.u.])')):
        raise InputError('orthogonal native k-point distribution required')
    points=np.empty((nk,3));weights=np.empty(nk)
    try:
        for k,line in enumerate(lines[5:5+nk]):
            parts=line.split()
            if len(parts)!=5 or int(parts[0])!=k+1:raise InputError('invalid k-point index/order')
            points[k]=[float(v) for v in parts[1:4]];weights[k]=float(parts[4])
        reciprocal=np.array([float(v) for v in lines[-1].split()])
    except (ValueError,OverflowError) as exc:raise InputError('invalid k-point distribution values') from exc
    if reciprocal.shape!=(3,) or not np.isfinite(reciprocal).all() or not np.allclose(reciprocal,2*np.pi/np.diag(cell),rtol=1e-14,atol=0):
        raise InputError('k-point reciprocal basis differs from cell')
    if not np.isfinite(points).all() or not np.isfinite(weights).all() or np.any(weights<0) or not np.isclose(weights.sum(),1.,rtol=0,atol=1e-12):
        raise InputError('invalid native k-point coordinates/weights')
    return points,weights


def _model_and_input(source,record,multi_k=False):
    from .dft import DFT,Pseudopotential,SCFSettings
    from .bands import _input_groups
    if not multi_k and record.get('kmesh')!=[1,1,1]:raise InputError('native export supports Gamma mesh only')
    potentials=[];refs=[]
    for i,p in enumerate(record['pseudopotentials']):
        candidates=[source.directory/'native'/f'pseudo{i+1}{suffix}' for suffix in ('_rps.dat','.upf')]
        candidates=[p for p in candidates if p.is_file()]
        if len(candidates)!=1:raise InputError('one DFT-generated potential file per species required')
        path=candidates[0];potential=Pseudopotential.from_file(path,atomic_number=p['atomic_number'],local_channel=p['local_channel'])
        if potential.sha256!=p['sha256']:raise InputError('potential bytes differ from GS metadata')
        potentials.append(potential);refs.append(dict(atomic_number=potential.atomic_number,filename=path.name,
            sha256=potential.sha256,source='verified native GS/1'))
    model=DFT(record['numbers'],record['cell_bohr'],record['fractional_positions'],potentials,
              record['electrons'],record['states'],record['grid'],record['kmesh'],record['xc'])
    if multi_k and (not np.array_equal(model.cell,np.diag(np.diag(model.cell))) or np.any(np.diag(model.cell)<=0)):
        raise InputError('multi-k native export requires aligned orthogonal cell')
    settings=SCFSettings(record['max_iterations'],record['tolerance'])
    expected=_input_groups(model.input_text(settings));actual=_input_groups((source.directory/'native/gs.inp').read_text())
    if set(actual)!=set(expected):raise InputError('only DFT-generated GS input groups supported')
    for name in expected:
        if name!='control' and actual[name].strip()!=expected[name].strip():raise InputError('GS metadata/input mismatch: '+name)
    control=expected['control'].strip()
    warm=control.replace("yn_restart='n'", "yn_restart='y', read_gs_restart_data='wfn', yn_reset_step_restart='y', directory_read_data='./initial_guess/'",1)
    if actual['control'].strip() not in (control,warm):raise InputError('unsupported native GS control settings')
    if _digest(source.directory/'native/gs.inp')!=record['input_sha256']:raise InputError('GS input hash differs from metadata')
    return model,refs


def export_gamma(converter,result,directory,*,orbital_block=16,timeout=300.):
    from .dft import GroundStateResult
    from .gs_archive import save_ground_state,_integer
    from .ld import _positive
    if not isinstance(result,GroundStateResult):raise InputTypeError('GroundStateResult required')
    orbital_block=_integer(orbital_block,'orbital_block');timeout=_positive(timeout,'timeout')
    if converter.exporter is None or converter.exporter_sha256 is None:raise RuntimeFailure('converter has no probed native exporter; explicitly build a matching new converter')
    if _digest(converter.exporter)!=converter.exporter_sha256:raise InputError('native exporter hash changed')
    build=json.loads(converter.build_record.read_text())
    if _digest(build['target_executable'])!=converter.target_sha256:raise InputError('paired runtime hash changed')
    source=GroundStateResult.open(result.directory);mp=source.directory/'ground-state.json'
    raw_manifest=mp.read_bytes();manifest_hash=hashlib.sha256(raw_manifest).hexdigest();record=json.loads(raw_manifest)
    if record.get('iterations')!=source.iterations or record.get('energy_hartree')!=source.energy_hartree:
        raise InputError('GS metadata changed during export validation')
    if record.get('executable_sha256')!=converter.target_sha256:raise InputError('GS and native exporter require the same runtime')
    try:model,refs=_model_and_input(source,record,converter._multi_k)
    except (KeyError,TypeError) as exc:raise InputError('incomplete DFT-generated GS metadata') from exc
    eigenpath=source.directory/'native/gs_eigen.data'
    if 'native/gs_eigen.data' not in record['files'] or _digest(eigenpath)!=record['files']['native/gs_eigen.data']:raise InputError('missing or changed native eigenvalue output')
    nk=int(np.prod(model.kmesh))
    points=np.zeros((1,3));weights=np.ones(1)
    if converter._multi_k and nk>1:
        kp=source.directory/'native/gs_k.data'
        if 'native/gs_k.data' not in record['files'] or _digest(kp)!=record['files']['native/gs_k.data']:raise InputError('missing or changed native k-point output')
        points,weights=_kpoints(kp.read_text(),nk,model.cell)
        mx,my,mz=model.kmesh
        expected=np.array([[(ix+.5)/mx-.5,(iy+.5)/my-.5,(iz+.5)/mz-.5]
                           for iz in range(mz) for iy in range(my) for ix in range(mx)])
        if not np.allclose(points,expected,rtol=0,atol=1e-14) or not np.allclose(weights,1./nk,rtol=1e-14,atol=0):
            raise InputError('k-point output differs from DFT-generated source mesh')
    eigen,eigen_occ=_eigenvalues(eigenpath.read_text(),model.states,nk)
    root=Path(directory).expanduser()
    if root.resolve().is_relative_to(source.directory):raise InputError('portable archive must be outside native GS result')
    if root.exists() or root.is_symlink():raise FileExistsError(root)
    environment=dict(os.environ);profile=_environment(environment)
    mappings=[]
    with tempfile.TemporaryDirectory(prefix='salmonpy-native-gs-export-') as tmp:
        work=Path(tmp);native=work/'native_input';native.mkdir()
        try:
            for name in ('info.bin','occupation.bin','wfn.bin'):
                relative='native/data_for_restart/'+name
                shutil.copy2(source.directory/relative,native/name)
                if _digest(native/name)!=record['files'][relative]:raise InputError('native GS changed while copying')
            nx,ny,nz=model.grid;nb=model.states;grid=(nz,ny,nx);ng=nx*ny*nz
            (work/'request.txt').write_text(f'{nx} {ny} {nz} {nb}\n{model.electrons}\n{source.iterations}\n'+(f'{nk}\n' if converter._multi_k else ''))
            if nk>1:(work/'portable_weights.bin').write_bytes(np.asarray(weights,dtype='<f8').tobytes())
            with (work/'exporter.log').open('w') as log:
                completed=subprocess.run([str(converter.exporter)],cwd=work,env=environment,stdout=log,stderr=subprocess.STDOUT,timeout=timeout,check=False)
            if completed.returncode:raise RuntimeFailure('native GS export failed: '+(work/'exporter.log').read_text()[-4000:])
            for name,count in [('portable_real.bin',ng*nb*nk),('portable_imag.bin',ng*nb*nk),('portable_occupations.bin',nb*nk)]:
                if (work/name).stat().st_size!=8*count:raise RuntimeFailure('invalid exported stream size')
            kind=(work/'native_kind.txt').read_text().strip()
            if kind not in ('real','complex'):raise RuntimeFailure('unknown exported native orbital kind')
            occ=np.fromfile(work/'portable_occupations.bin',dtype='<f8')
            if not np.array_equal(occ,eigen_occ):raise InputError('native occupation and eigenvalue table disagree')
            occ=occ.reshape(nk,nb)
            real=np.memmap(work/'portable_real.bin',dtype='<f8',mode='r',shape=(nk,nb,*grid));mappings.append(real)
            imag=np.memmap(work/'portable_imag.bin',dtype='<f8',mode='r',shape=(nk,nb,*grid));mappings.append(imag)
            waves=np.memmap(work/'complex-map.bin',dtype=np.complex128,mode='w+',shape=(nk,nb,*grid));mappings.append(waves)
            density=np.zeros(grid)
            try:
                with np.errstate(over='raise',invalid='raise'):
                    for k in range(nk):
                        for b in range(nb):
                            waves[k,b].real=real[k,b];waves[k,b].imag=imag[k,b]
                            density+=(weights[k]*occ[k,b])*(real[k,b]*real[k,b]+imag[k,b]*imag[k,b])
            except FloatingPointError as exc:raise InputError('nonfinite reconstructed density') from exc
            waves.flush()
            provenance=dict(schema='salmonpy.gs-native-export/1',source_manifest_sha256=manifest_hash,
                target_sha256=converter.target_sha256,exporter_sha256=converter.exporter_sha256,
                density_source='reconstructed_from_native_orbitals',eigenvalue_source='native/gs_eigen.data',
                native_real_orbital=kind=='real',local_channels=[dict(atomic_number=p.atomic_number,local_channel=p.local_channel) for p in model.pseudopotentials],
                record_environment=dict(profile))
            if nk>1:provenance['kpoint_source']='native/gs_k.data'
            return save_ground_state(root,wavefunctions=waves,density=density,eigenvalues=eigen,
                occupations=occ,kpoints=points,weights=weights,cell_bohr=model.cell,
                fractional_positions=model.fractional_positions,atomic_numbers=list(model.numbers),electrons=model.electrons,
                energy_hartree=source.energy_hartree,xc=model.xc,pseudopotentials=refs,orbital_block=orbital_block,provenance=provenance)
        finally:
            # Explicitly release mapping handles before TemporaryDirectory cleanup.
            for mapping in reversed(mappings):mapping._mmap.close()
