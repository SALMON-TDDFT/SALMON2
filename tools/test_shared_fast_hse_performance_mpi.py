"""Interleaved shared/registered HSE modes vs independent previous-runtime controls.

mpiexec -n 6 python tools/test_shared_fast_hse_performance_mpi.py BUILD BASELINE FIXTURE NEW_OUT
Use --symmetry for a reduced first model. --second-fixture/--second-symmetry
register a different model on the same attachment. Existing GS is linked, never
calculated/copied. No electronic disk checkpoint. Use OMP/BLAS=1 and MPI timeout.
"""
import argparse,ctypes as c,hashlib,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import NativeParallelContext
from salmonpy.fast_hse_mpi import _configure,_history_info
from salmonpy.fast_hse import _input_text
from salmonpy._shared_hse import NativeFastHSEModel,_configure_handles
from salmonpy._hse_registry import NativeFastHSERegistry,_ModelLease,_ModelLibrary,_RegisteredModel
from salmonpy._hse_registry import _configure as configure_registry
p=argparse.ArgumentParser(description=__doc__)
for name in ('build','baseline','fixture','output'):p.add_argument(name)
p.add_argument('--symmetry',action='store_true');p.add_argument('--second-fixture');p.add_argument('--second-symmetry',action='store_true')
a=p.parse_args();world=MPI.COMM_WORLD
if world.size<3 or world.size%3:raise ValueError('three equally sized solver groups required')
n=world.size//3;g=world.rank//n;comm=world.Split(g,world.rank)
record=json.loads(Path(a.build if g==0 else a.baseline).read_text())
assert hashlib.sha256(Path(record['library']).read_bytes()).hexdigest()==record['sha256']
lib=c.CDLL(record['library']);_configure(lib);_configure_handles(lib)
if g==0 and a.second_fixture:configure_registry(lib)
native=NativeParallelContext(SimpleNamespace(solver_comm=comm,closed=False,_dependents=0),record['library'],record['sha256'])
root=Path(a.output)/str(g);models=[];states=[]
registry=NativeFastHSERegistry(SimpleNamespace(),native,lib) if g==0 and a.second_fixture else None

def model_for(tag,fixture,symmetry):
    fixture=Path(fixture).resolve();work=root/tag
    groups,propagator=_input_text((fixture/'inputfile').read_text())
    if comm.rank==0:
        work.mkdir(parents=True,exist_ok=False)
        groups['parallel']=groups.get('parallel','')+'\n nproc_k='+str(n)+', nproc_ob=1, nproc_rgrid=1,1,1\n'
        (work/'inputfile').write_text('\n'.join('&'+k+v+'/' for k,v in groups.items())+'\n')
        for name in ('pseudo1_rps.dat','sym.dat'):(work/name).write_bytes((fixture/name).read_bytes())
        (work/'data_for_restart').symlink_to(fixture/'data_for_restart',target_is_directory=True)
    comm.Barrier()
    args=(work,.001,16,propagator,tag,())
    if registry is None:
        model=NativeFastHSEModel(native,lib,*args,symmetry=symmetry,execution_mode='performance' if g==0 else 'checked')
    else:
        lease=_ModelLease(registry)
        model=_RegisteredModel(lease,_ModelLibrary(registry),*args,symmetry=symmetry,execution_mode='performance')
        lease.model=model;model._registry_id=len(models);registry._models.add(model)
    models.append(model);return model

if g==0:
    first=model_for('first',a.fixture,a.symmetry)
    second=model_for('second',a.second_fixture,a.second_symmetry) if a.second_fixture else first
    specs=[(first,1,'checked'),(first,1,None),(second,2,'performance')]
    if a.second_fixture:specs.append((second,2,'checked'))
else:
    second=(g==2 and a.second_fixture is not None)
    model=model_for('control',a.second_fixture if second else a.fixture,a.second_symmetry if second else a.symmetry)
    specs=[(model,g,'checked')]
for model,k,mode in specs:states.append(model.new_state(str(len(states)),execution_mode=mode))
if g==0:assert states[1].execution_mode=='performance' and states[0].execution_mode=='checked'
traces=[[] for _ in states];saved=[]
def drive(i,k):
    symmetric=a.second_symmetry if k==2 and a.second_fixture else a.symmetry
    return (np.array([0.,0.,.001*np.sin(.17*i)]) if symmetric else np.array([.001*np.sin(.2*i),.0007*(np.cos(.13*i)-1),.0009*np.sin(.17*i)]))*(1. if k==1 else -.7)
for i in range(1,17):
    for index,(state,(_,k,_)) in enumerate(zip(states,specs)):
        before=[(s.A.copy(),s.current.copy(),s.steps,s.energy) for s in states]
        current=state.step(drive(i,k));history=state._call(lambda:_history_info(state.lib))
        traces[index].append((current.copy(),state.energy,history.copy()))
        for other,s in enumerate(states):
            if other!=index:
                np.testing.assert_array_equal(s.A,before[other][0]);np.testing.assert_array_equal(s.current,before[other][1])
                assert s.steps==before[other][2] and s.energy==before[other][3]
        if i==8:
            assert state._call(state.lib.sp_fast_hse_step_transaction,0)==0
            saved.append((state.energy,state.execution_mode))
for index,(state,(_,k,_)) in enumerate(zip(states,specs)):
    assert state._call(state.lib.sp_fast_hse_step_transaction,1)==0
    step=c.c_int();A=np.zeros(3);j=np.zeros(3)
    assert state._call(state.lib.sp_fast_hse_info,c.byref(step),A,j)==0
    state.steps=step.value;state.A=A;state.current=j;state.energy=saved[index][0]
    assert state.steps==8 and state.execution_mode==saved[index][1]
    assert state._call(state.lib.sp_fast_hse_step_transaction,3)==0
    for i in range(9,17):
        np.testing.assert_array_equal(state.step(drive(i,k)),traces[index][i-1][0])
        assert state.energy==traces[index][i-1][1]
        np.testing.assert_array_equal(state._call(lambda:_history_info(state.lib)),traces[index][i-1][2])
    assert state._call(state.lib.sp_fast_hse_step_transaction,2)==0
old_token=states[0].handle
for state in states:state.close()
if g==0:
    fresh=models[0].new_state('generation-reuse',execution_mode='performance')
    assert fresh.handle!=old_token
    vector=np.ctypeslib.ndpointer(dtype=np.float64,shape=(3,),flags='C_CONTIGUOUS')
    fn=lib.sp_fast_hse_handle_step_performance;fn.argtypes=[c.c_int64,vector,vector,vector,c.POINTER(c.c_double)];fn.restype=c.c_int
    assert models[0]._call(fn,old_token,np.zeros(3),np.zeros(3),np.zeros(3),c.byref(c.c_double()))==12
    np.testing.assert_array_equal(fresh.step(drive(1,1)),traces[0][0][0]);fresh.close()
for model in models:model.close()
if registry is not None:registry.close()
rows=world.gather((g,comm.rank,[k for _,k,_ in specs],traces),root=0)
if world.rank==0:
    controls={k:row[3][index] for row in rows if row[0] in (1,2) for index,k in enumerate(row[2])}
    for group,rank,keys,streams in rows:
        for index,k in enumerate(keys):
            for step,(current,energy,history) in enumerate(streams[index]):
                expected=controls[k][step]
                np.testing.assert_array_equal(current,expected[0]);assert energy==expected[1]
                np.testing.assert_array_equal(history,expected[2])
    arrays={}
    for k,stream in controls.items():
        arrays[f'current{k}']=np.array([x[0] for x in stream]);arrays[f'energy{k}']=np.array([x[1] for x in stream]);arrays[f'history{k}']=np.array([x[2] for x in stream])
    np.savez(Path(a.output)/'traces.npz',**arrays)
    report={'ranks_per_group':n,'registered':bool(a.second_fixture),'symmetry':a.symmetry,'second_symmetry':a.second_symmetry,'steps':16,
        'interleaved_modes_vs_independent_previous_array_equal':True,'siblings_unchanged':True,
        'explicit_memory_restore_array_equal':True,'mode_preserved':True,'stale_generation_rejected':True,
        'GS':'existing read links only; no new GS/copies/electronic disk saves'}
    (Path(a.output)/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
