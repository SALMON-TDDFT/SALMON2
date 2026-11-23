"""Synthetic orthogonal multi-k native→GS/2→TDDFT parity. No SCF/actual GS collection.

python tools/test_gs_multik_runtime.py CONVERTER_BUILD NEW_OUTPUT
Uses existing SHA-validated runtime/probe, builds nothing. Fixtures deliberately
include synthetic complete-result metadata only to exercise the API contract.
"""
import argparse,hashlib,json,os,shutil,subprocess
from pathlib import Path
import numpy as np
from salmonpy.gs_native import SerialRestartConverter
from salmonpy.dft import DFT,Pseudopotential,SCFSettings,GroundStateResult
from salmonpy.fdtd import SalmonExecutable


def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('build',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    record=json.loads(a.build.read_text());converter=SerialRestartConverter.from_build(a.build)
    engine=SalmonExecutable(record['target_executable'],record['target_sha256']);probe=Path(record['probe'])
    if digest(probe)!=record['probe_sha256']:raise ValueError('probe hash changed')
    root=a.output.resolve();root.mkdir(parents=True,exist_ok=False);env=dict(os.environ)
    package=Path(__import__('salmonpy').__file__).resolve().parent
    pp=Pseudopotential.from_file(package/'examples/si_intro_data/Si_rps.dat',atomic_number=14,local_channel=2)
    summaries=[]
    for nk,mesh in [(2,(2,1,1)),(8,(2,2,2)),(3,(3,1,1)),(6,(3,2,1))]:
        kind='complex'
        case=root/('k'+str(nk));native=case/'source/native';restart=native/'data_for_restart';restart.mkdir(parents=True)
        mode='write-gs-norm'
        with (case/'probe.log').open('w') as log:
            subprocess.run([str(probe),mode,str(restart)+os.sep,str(nk)],stdout=log,stderr=subprocess.STDOUT,env=env,check=True,timeout=60)
        model=DFT([14],np.eye(3)*8,[[.3,.4,.45]],[pp],4,2,(8,8,8),mesh,'pz');settings=SCFSettings()
        (native/'gs.inp').write_text(model.input_text(settings));(native/'pseudo1_rps.dat').write_bytes(pp.data)
        (native/'gs_status.dat').write_text('SALMONPY_GS_V1\n1 1\n1e-10 1e-9 -1.0\nrho_dne\n')
        eigen=['#esp: single-particle energies (eigen energies)','#occ: occupation numbers, io: orbital index','# 1:io, 2:esp[a.u.], 3:occ']
        for k in range(nk):eigen += [f'k= {k+1}, spin= 1',f'1 {-.5+.01*k:.16e} 2',f'2 {.2+.01*k:.16e} 2']
        (native/'gs_eigen.data').write_text('\n'.join(eigen)+'\n')
        points=np.array([[(ix+.5)/mesh[0]-.5,(iy+.5)/mesh[1]-.5,(iz+.5)/mesh[2]-.5]
                         for iz in range(mesh[2]) for iy in range(mesh[1]) for ix in range(mesh[0])])
        # Simulate native vec_k/reciprocal conversion and limited E23.15 text precision.
        reciprocal=2*np.pi/np.diag(model.cell);reduced=(points*reciprocal)/reciprocal
        kt=['# k-point distribution','# ik: k-point index','# kx,ky,kz: Reduced coordinate of k-points',
            '# wk: Weight of k-point','# 1:ik[none] 2:kx[none] 3:ky[none] 4:kz[none] 5:wk[none]']
        kt += [str(k+1)+' '+' '.join(format(float(v),'.15g') for v in (*reduced[k],1/nk)) for k in range(nk)]
        kt += ['# coefficients (2*pi/a [a.u.]) in kx, ky, kz',' '.join(format(v,'.15g') for v in reciprocal)]
        (native/'gs_k.data').write_text('\n'.join(kt)+'\n')
        gs=dict(model._restart_signature(),schema='salmonpy.gs/1',status='complete',fractional_positions=model.fractional_positions.tolist(),
                executable_sha256=engine.sha256,energy_hartree=-1.,iterations=1,tolerance=settings.tolerance,max_iterations=settings.max_iterations,
                input_sha256=digest(native/'gs.inp'),files={x.relative_to(case/'source').as_posix():digest(x) for x in native.rglob('*') if x.is_file()})
        (case/'source/ground-state.json').write_text(json.dumps(gs,indent=2)+'\n')
        source=GroundStateResult.open(case/'source');archive=source.export_archive(case/'archive',converter=converter,orbital_block=1)
        converted=converter.write_bundle(archive,case/'reimported',executable=engine,pseudopotentials=[pp])
        check='read-gs-norm'
        with (case/'probe.log').open('a') as log:
            subprocess.run([str(probe),check,str(converted/'data_for_restart')+os.sep,str(nk)],stdout=log,stderr=subprocess.STDOUT,env=env,check=True,timeout=60)
        values=np.array([[0,0,.0001*n] for n in range(1,7)])
        with source.tddft_initial_state(executable=engine,max_steps=12,timeout=60) as original,archive.tddft_initial_state(
             converter=converter,executable=engine,pseudopotentials=[pp],max_steps=12,timeout=60) as imported:
            with original.new_state(dt=.001) as one,imported.new_state(dt=.001) as two:
                first=[];second=[]
                for index,A in enumerate(values):
                    first.append(one.step(A));second.append(two.step(A))
                    if index==1:token=two.snapshot()
                left=np.array(first);right=np.array(second);np.testing.assert_allclose(left,right,rtol=1e-12,atol=1e-14)
                imported.close();two.restore(token)
                replay=np.array([two.step(A) for A in values[2:]])
                assert replay.tobytes()==right[2:].tobytes()
                for label,state in [('original',one),('imported',two)]:shutil.copy2(Path(state._work.name)/'run/worker.log',case/(label+'-worker.log'))
        density_electrons=float(archive.array('density').sum()*512/512)
        if not np.isclose(density_electrons,4.,atol=1e-12,rtol=0):raise ValueError('reconstructed density integral differs')
        summaries.append(dict(nk=nk,kmesh=mesh,native_kind=kind,current_rtol=1e-12,current_atol=1e-14,native_vs_exported_allclose=True,unchanged_native_reader_wave_parity=True,bitwise_native_vs_exported_currents=left.tobytes()==right.tobytes(),
            bitwise_replay=True,max_current_difference=float(np.max(abs(left-right))),density_integral=density_electrons,
            occupations=archive.array('occupations').tolist(),currents=right.tolist(),provenance=dict(archive.metadata['provenance'])))
    # Independent explicit-native file_kw baseline for a nonuniform-weight import.
    from salmonpy.gs_archive import GroundStateArchive,save_ground_state
    from salmonpy.bands import _input_groups
    base=GroundStateArchive.open(root/'k2/archive');weighted_root=root/'weighted';weighted_root.mkdir()
    shutil.copytree(root/'k2/source',weighted_root/'source')
    native=weighted_root/'source/native';input_file=native/'gs.inp';groups=_input_groups(input_file.read_text())
    groups['kgrid']="\n file_kw='reference-kpoints.dat'\n"
    input_file.write_text('\n'.join('&'+k+v+'/' for k,v in groups.items())+'\n')
    weights=np.array([.25,.75]);points=base.array('kpoints')
    lines=['2']+[str(k+1)+' '+' '.join(format(float(v),'.17e') for v in (*points[k],weights[k])) for k in range(2)]
    (native/'reference-kpoints.dat').write_text('\n'.join(lines)+'\n')
    gs=json.loads((weighted_root/'source/ground-state.json').read_text());gs['input_sha256']=digest(input_file)
    gs['files']={x.relative_to(weighted_root/'source').as_posix():digest(x) for x in native.rglob('*') if x.is_file()}
    (weighted_root/'source/ground-state.json').write_text(json.dumps(gs)+'\n')
    source=GroundStateResult.open(weighted_root/'source');waves=np.stack([base.wavefunctions(k) for k in range(2)])
    occ=base.array('occupations');density=np.einsum('k,kb,kbzyx->zyx',weights,occ,abs(waves)**2)
    weighted_archive=save_ground_state(weighted_root/'archive',wavefunctions=waves,density=density,
        eigenvalues=base.array('eigenvalues'),occupations=occ,kpoints=points,weights=weights,cell_bohr=base.array('cell_bohr'),
        fractional_positions=base.array('fractional_positions'),atomic_numbers=base.array('atomic_numbers'),electrons=4.,
        energy_hartree=-1.,xc='pz',pseudopotentials=[dict(atomic_number=14,filename='Si_rps.dat',sha256=pp.sha256,source='synthetic')])
    with source.tddft_initial_state(executable=engine,max_steps=12,timeout=60) as original,weighted_archive.tddft_initial_state(
         converter=converter,executable=engine,pseudopotentials=[pp],max_steps=12,timeout=60) as imported:
        with original.new_state(dt=.001) as one,imported.new_state(dt=.001) as two:
            left=np.array([one.step(A) for A in values]);right=np.array([two.step(A) for A in values])
            np.testing.assert_allclose(left,right,rtol=1e-12,atol=1e-14)
    summaries.append(dict(nk=2,weights=weights.tolist(),scope='independent explicit-native file_kw vs portable nonuniform import',
        native_vs_exported_allclose=True,bitwise_native_vs_exported_currents=left.tobytes()==right.tobytes(),
        max_current_difference=float(np.max(abs(left-right))),currents=right.tolist()))
    # MappingProxy/tuple nested provenance values are serialized by explicit conversion.
    def plain(x):
        from collections.abc import Mapping
        if isinstance(x,Mapping):return {k:plain(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [plain(v) for v in x]
        return x
    result=dict(scope='synthetic multi-k native records/metadata; actual unchanged native runtime, NOT SCF/physical convergence',
                target_sha256=engine.sha256,exporter_sha256=converter.exporter_sha256,cases=summaries)
    (root/'result.json').write_text(json.dumps(plain(result),indent=2)+'\n');print(json.dumps(plain(result),indent=2))

if __name__=='__main__':main()
