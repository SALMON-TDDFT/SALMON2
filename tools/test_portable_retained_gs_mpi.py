"""Synthetic GS/2→retained GS→TDDFT vs independent file route; no SCF/real GS.

mpiexec -n 4 python tools/test_portable_retained_gs_mpi.py NATIVE_BUILD CONVERTER_BUILD NEW_OUTPUT --split orbital
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import warnings
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from salmonpy.gs_archive import GroundStateArchive,save_ground_state
from salmonpy.gs_native import GammaRestartConverter,SerialRestartConverter
from salmonpy.native_gs import NativeGSBackend
from salmonpy.embedded_tddft import NativeTDDFTBackend,NativeTDDFTMaterialProvider
from salmonpy.fdtd import SalmonExecutable
from salmonpy.dft import Pseudopotential
from salmonpy.errors import InputError


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('converter',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--split',choices=('orbital','spatial','k'),default='orbital')
    p.add_argument('--xc',choices=('pz','libxc_pz','libxc_pbe'),default='pz')
    p.add_argument('--gamma',action='store_true')
    p.add_argument('--nonuniform',action='store_true')
    a=p.parse_args()
    if a.gamma and a.nonuniform:p.error('nonuniform test uses two k points')
    if a.gamma and a.split=='k':p.error('this test fixture needs two k points for its two-rank k split')
    from mpi4py import MPI
    world=MPI.COMM_WORLD
    if world.size!=4:raise ValueError('use four ranks (two groups of two)')
    group=world.rank//2;comm=world.Split(group,world.rank)
    dec={'orbital':(1,2,1,1,1),'spatial':(1,1,2,1,1),'k':(2,1,1,1,1)}[a.split]
    package=Path(__file__).resolve().parents[1]/'salmonpy'
    ppfile=package/'examples/si_intro_data/Si_rps.dat'
    if not ppfile.exists():ppfile=Path(__file__).resolve().parents[2]/'examples/si_intro_data/Si_rps.dat'
    pp=Pseudopotential.from_file(ppfile,atomic_number=14,local_channel=2)
    converter=(GammaRestartConverter if a.gamma else SerialRestartConverter).from_build(a.converter)
    receipt=json.loads(a.converter.read_text())
    engine=SalmonExecutable(receipt['target_executable'],receipt['target_sha256'])
    error=None
    if world.rank==0:
        try:
            a.output.mkdir(parents=True,exist_ok=False)
            z,y,x=np.indices((8,8,8))
            u=np.exp(2j*np.pi*(x+2*y+3*z)/8)/np.sqrt(512.)
            waves=np.stack([np.stack([np.ones((8,8,8))/np.sqrt(512.),u]),
                            np.stack([np.ones((8,8,8))/np.sqrt(512.),u])]).astype(complex)
            waves.imag[0,1,0,0,0]=-0.
            if a.gamma:waves=waves[:1]
            nk=len(waves)
            occupations=np.full((nk,2),2.)
            weights=np.array([1.]) if a.gamma else np.array([.25,.75]) if a.nonuniform else np.array([.5,.5])
            density=np.sum(abs(waves)**2*occupations[:,:,None,None,None]*weights[:,None,None,None,None],axis=(0,1))
            save_ground_state(a.output/'archive',wavefunctions=waves,density=density,
                eigenvalues=[[-.5,.2]] if a.gamma else [[-.5,.2],[-.4,.3]],occupations=occupations,
                kpoints=[[0,0,0]] if a.gamma else [[-.2,0,0],[.3,0,0]] if a.nonuniform else [[-.25,0,0],[.25,0,0]],weights=weights,cell_bohr=np.eye(3)*8,fractional_positions=[[.3,.4,.45]],
                atomic_numbers=[14],electrons=4.,energy_hartree=-1.,xc=a.xc,orbital_block=1,
                pseudopotentials=[dict(atomic_number=14,filename='Si_rps.dat',sha256=pp.sha256,source='synthetic orbitals')])
        except Exception as exc:error=str(exc)
    error=world.bcast(error,root=0)
    if error:raise RuntimeError(error)
    root=a.output/str(group)
    if comm.rank==0:root.mkdir()
    comm.Barrier();context=SimpleNamespace(solver_comm=comm,point_ids=(1,),closed=False,_dependents=0)
    archive=GroundStateArchive.open(a.output/'archive')
    def A(i):return 1e-4*np.array([np.sin(.31*i),np.cos(.23*i)-1,np.sin(.17*i)])
    trace=[];captured=[]
    try:
        if group==0:
            with NativeGSBackend(a.build,root,decomposition=dec).open(context) as session:
                with warnings.catch_warnings(record=True) as captured:
                    warnings.simplefilter('always')
                    gs=session.load_archive(archive,root/'converted',converter=converter,executable=engine,
                                            pseudopotentials=[pp],material_id='synthetic')
                assert gs.info.converged is None and gs.calculation is None
                origin=gs.archive_origin
                assert origin['material_id']=='synthetic' and origin['scf_converged'] is None
                rho=gs.array('density');orbitals=gs.array('orbitals')
                starts=[lo-base for lo,base in zip(rho.lower_indices,orbitals.lower_indices)]
                region=tuple(slice(lo,lo+n) for lo,n in zip(starts,rho.values.shape))
                selected=orbitals.values[region+(0,slice(None),slice(None),0)].transpose(4,3,2,1,0)
                ks=orbitals.lower_indices[5]-1;bs=orbitals.lower_indices[4]-1
                xs,ys,zs=(v-1 for v in rho.lower_indices);nx,ny,nz=rho.values.shape
                expected=np.stack([archive.wavefunctions(k,bs,bs+selected.shape[1])[:,zs:zs+nz,ys:ys+ny,xs:xs+nx]
                                   for k in range(ks,ks+selected.shape[0])])
                assert selected.tobytes()==expected.tobytes(), 'native complex bit patterns changed'
                np.testing.assert_array_equal(gs.array('occupations').values[:,:,0].T,archive.array('occupations'))
                before=session.resources();raw_before=orbitals.values.tobytes()
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore',RuntimeWarning)
                    reexport=gs.export_archive(root/'reexport',orbital_block=1)
                assert session.resources()==before
                assert gs.array('orbitals').values.tobytes()==raw_before
                error=None
                if comm.rank==0:
                    try:
                        for k in range(archive.wavefunction_shape[0]):
                            assert reexport.wavefunctions(k).tobytes()==archive.wavefunctions(k).tobytes()
                        for name in ('occupations','eigenvalues','kpoints','weights','cell_bohr','fractional_positions','atomic_numbers'):
                            np.testing.assert_array_equal(reexport.array(name),archive.array(name))
                        assert reexport.metadata['energy_hartree']==archive.metadata['energy_hartree']
                        provenance=reexport.metadata['provenance']
                        assert provenance['schema']=='salmonpy.retained-gs-reexport/1'
                        assert provenance['energy_source']==provenance['eigenvalue_source']==provenance['sampling_source']=='source_archive'
                        assert provenance['converged'] is None
                    except Exception as exc:error=str(exc) or type(exc).__name__
                error=comm.bcast(error,root=0)
                if error:raise AssertionError(error)
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore',RuntimeWarning)
                    second_gs=session.load_archive(root/'reexport',root/'reimported',converter=converter,
                        executable=engine,pseudopotentials=[pp])
                assert second_gs.material_id=='synthetic' and second_gs.info.converged is None
                reimport_model=second_gs.tddft_model(dt=.001,max_steps=24);second_gs.close()
                reimport_state=reimport_model.new_state('reimported-response')
                model=gs.tddft_model(dt=.001,max_steps=24);gs.close()
                assert gs.archive_origin==origin
                first=model.new_state('a');sibling=model.new_state('b')
                for i in range(1,9):
                    trace.append(first.step(A(i)));sibling.step(-A(i))
                    np.testing.assert_array_equal(reimport_state.step(A(i)),trace[-1])
                first.checkpoint(root/f'state-{comm.rank}.bin')
                before=session.resources();current=first.current.copy()
                bad=Pseudopotential(14,2,'_rps.dat',b'deliberately wrong bytes')
                caught=False
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore',RuntimeWarning)
                    try:session.load_archive(archive,root/'retry',converter=converter,executable=engine,
                                             pseudopotentials=[bad],material_id='retry')
                    except InputError:caught=True
                assert all(comm.allgather(caught))
                assert session.resources()==before
                np.testing.assert_array_equal(first.current,current)
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore',RuntimeWarning)
                    other=session.load_archive(archive,root/'retry',converter=converter,executable=engine,
                                               pseudopotentials=[pp],material_id='retry')
                # Re-export while RT owners are active; fail publication once,
                # preserve all owners and retry the same output path.
                from unittest.mock import patch
                from salmonpy import gs_archive as archive_io
                writer=archive_io._write;calls=0
                def broken(path,data):
                    nonlocal calls
                    calls+=1
                    if calls==3:raise OSError('injected reexport write failure')
                    return writer(path,data)
                before=session.resources();current=first.current.copy();caught=False
                with warnings.catch_warnings(),patch.object(archive_io,'_write',side_effect=broken):
                    warnings.simplefilter('ignore',RuntimeWarning)
                    try:other.export_archive(root/'retry-save')
                    except OSError:caught=True
                assert all(comm.allgather(caught))
                assert session.resources()==before
                np.testing.assert_array_equal(first.current,current)
                if comm.rank==0:assert not (root/'retry-save').exists()
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore',RuntimeWarning)
                    saved=other.export_archive(root/'retry-save')
                assert session.resources()==before
                np.testing.assert_array_equal(first.current,current)
                if comm.rank==0:
                    for k in range(archive.wavefunction_shape[0]):
                        assert saved.wavefunctions(k).tobytes()==archive.wavefunctions(k).tobytes()
                other.close()
                for i in range(9,17):
                    trace.append(first.step(A(i)))
                    np.testing.assert_array_equal(reimport_state.step(A(i)),trace[-1])
                restored=model.new_state('restored');restored.checkpoint(root/f'state-{comm.rank}.bin',restore=True)
                for i in range(9,17):np.testing.assert_array_equal(restored.step(A(i)),trace[i-1])
                independent=model.new_state('independent')
                for i in range(1,17):np.testing.assert_array_equal(independent.step(A(i)),trace[i-1])
                for state in (first,sibling,restored,independent):state.close()
                reimport_state.close();reimport_model.close()
                model.close()
        else:
            error=None
            if comm.rank==0:
                try:converter.write_bundle(archive,root/'reference',executable=engine,pseudopotentials=[pp])
                except Exception as exc:error=str(exc)
            error=comm.bcast(error,root=0)
            if error:raise RuntimeError(error)
            initial=NativeTDDFTMaterialProvider(root/'reference','rt.inp',max_steps=24).prepare(context)
            state=NativeTDDFTBackend(a.build,root,decomposition=dec).create_state(initial,1,.001,context)
            try:
                for i in range(1,17):trace.append(state.advance(state.A.copy(),A(i)))
            finally:state.close()
        traces=world.allgather(np.asarray(trace))
        np.testing.assert_array_equal(traces[0],traces[2])
        if world.rank==0:
            result=dict(status='pass',active_RT_reexport_failure_safe=True,active_RT_reexport_retry_exact=True,reexport_exact=True,reimport_current_exact=True,reexport_preserved_owner=True,
                nonuniform=a.nonuniform,gamma=a.gamma,xc=a.xc,current_sha256=hashlib.sha256(np.asarray(trace).tobytes()).hexdigest(),split=a.split,synthetic_native_orbitals_exact=True,
                portable_to_memory_current_exact=True,independent_exact=True,resume_exact=True,
                conversion_failure_preserved_state=True,retry_exact=True,SCF_calls=0,replayed_steps=0,
                warnings=[str(w.message) for w in captured])
            (a.output/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
    except BaseException:
        import traceback;traceback.print_exc();world.Abort(1)


if __name__=='__main__':main()
