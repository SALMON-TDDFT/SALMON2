"""Compare interleaved TDDFT checked/performance owners to separate controls.

mpiexec -n 6 python tools/test_shared_tddft_performance_mpi.py BUILD BASELINE FIXTURE NEW_OUT pz
Add --second-xc tbmbj for heterogeneous registry switching. FIXTURE is the small
normalized synthetic input rt.inp/pseudo1_rps.dat/data_for_restart used by
 test_tddft_performance_mpi.py. No GS calculation/copy or electronic disk save.
Use an MPI timeout and OMP/BLAS=1. Physical convergence is outside this probe.
"""
import argparse,ctypes as c,hashlib,json,re
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from mpi4py import MPI
from salmonpy.distributed import NativeParallelContext
from salmonpy._shared_tddft import NativeTDDFTModel,_configure
from salmonpy._tddft_registry import (_configure_registry,NativeTDDFTRegistry,
    _RegistryLease,_RegisteredLibrary,_RegisteredModel)

p=argparse.ArgumentParser(description=__doc__)
for name in ('build','baseline','fixture','output','xc'):p.add_argument(name)
p.add_argument('--second-xc',choices=('pz','tbmbj','libxc_pz','libxc_pbe'))
a=p.parse_args();world=MPI.COMM_WORLD
if world.size<3 or world.size%3:raise ValueError('three equally sized solver groups required')
if a.xc not in ('pz','tbmbj','libxc_pz','libxc_pbe'):raise ValueError('unsupported probe XC')
n=world.size//3;g=world.rank//n;comm=world.Split(g,world.rank)
record=json.loads(Path(a.build if g==0 else a.baseline).read_text())
assert hashlib.sha256(Path(record['library']).read_bytes()).hexdigest()==record['sha256']
lib=c.CDLL(record['library'])
heterogeneous=a.second_xc is not None
if heterogeneous and g==0:_configure_registry(lib)
else:_configure(lib)
context=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
native=NativeParallelContext(context,record['library'],record['sha256'])
root=Path(a.output)/str(g);fixture=Path(a.fixture).resolve()
models=[];states=[]
def make_work(tag,xc):
    work=root/tag
    if comm.rank==0:
        work.mkdir(parents=True,exist_ok=False)
        text=re.sub(r"xc\s*=\s*'[^']+'","xc='"+xc+"'",(fixture/'rt.inp').read_text())
        (work/'rt.inp').write_text(text)
        (work/'pseudo1_rps.dat').write_bytes((fixture/'pseudo1_rps.dat').read_bytes())
        (work/'data_for_restart').symlink_to(fixture/'data_for_restart',target_is_directory=True)
    comm.Barrier();return work

def standard_model(tag,xc):
    work=make_work(tag,xc)
    model=NativeTDDFTModel(native,lib,work,SimpleNamespace(max_steps=16,manifest=(),input_file='rt.inp'),
                          .001,(n,1,1,1,1),tag,Path('.'))
    models.append(model);return model

registry=None
if heterogeneous and g==0:
    registry=NativeTDDFTRegistry(SimpleNamespace(execution_mode='performance'),native,lib)
    def registered_model(tag,xc):
        work=make_work(tag,xc);mt=c.c_int64();st=c.c_int64();j=np.zeros(3)
        code=registry._call(work,lib.sp_tddft_model_register,b'rt.inp',.001,16,np.zeros(3),
                            np.array([n,1,1,1,1],np.int32),c.byref(mt),c.byref(st),j)
        assert code==0
        proxy=_RegisteredLibrary(registry,mt.value,st.value,j)
        model=_RegisteredModel(_RegistryLease(registry),proxy,work,
            SimpleNamespace(max_steps=16,manifest=(),input_file='rt.inp'),.001,(n,1,1,1,1),tag,Path('.'),
            execution_mode='performance')
        model.registry=registry;model.model_token=mt.value;model._initialized=True;model._attempted=True
        registry._models.add(model);models.append(model);return model
    first=registered_model('first',a.xc);second=registered_model('second',a.second_xc)
    specs=[(first,1,'checked'),(first,1,'performance'),(second,2,'checked'),(second,2,None)]
elif g==0:
    model=standard_model('shared',a.xc)
    specs=[(model,1,'checked'),(model,1,'performance'),(model,2,'performance')]
else:
    model=standard_model('control',a.second_xc if heterogeneous and g==2 else a.xc)
    specs=[(model,g,'checked')]
for model,k,mode in specs:
    state=model.new_state(str(len(states)),execution_mode=mode);states.append(state)
assert not(g==0) or states[1].execution_mode=='performance'
traces=[[] for s in states];saved=[]
def drive(i,k):return np.array([.001*np.sin(.2*i),.0007*(np.cos(.13*i)-1),.0009*np.sin(.17*i)])*(1. if k==1 else -.7)
for i in range(1,17):
    for index,(state,(_,k,_)) in enumerate(zip(states,specs)):
        siblings=[s.current.copy() for s in states]
        traces[index].append(state.step(drive(i,k)))
        for j,s in enumerate(states):
            if j!=index:np.testing.assert_array_equal(s.current,siblings[j])
        if i==8:
            assert state._call(state.lib.sp_tddft_step_transaction,0)==0
            saved.append(state.execution_mode)
for index,(state,(_,k,_)) in enumerate(zip(states,specs)):
    assert state._call(state.lib.sp_tddft_step_transaction,1)==0
    state._refresh()
    assert state.steps==8 and state.execution_mode==saved[index]
    assert state._call(state.lib.sp_tddft_step_transaction,3)==0
    for i in range(9,17):np.testing.assert_array_equal(state.step(drive(i,k)),traces[index][i-1])
    assert state._call(state.lib.sp_tddft_step_transaction,2)==0
for state in states:state.close()
if g==0:
    fn=lib.sp_tddft_handle_step_performance
    vector=np.ctypeslib.ndpointer(dtype=np.float64,shape=(3,),flags='C_CONTIGUOUS')
    fn.argtypes=[c.c_int64,vector,vector,vector];fn.restype=c.c_int
    assert models[0]._call(fn,states[0].handle,np.zeros(3),np.zeros(3),np.zeros(3))==12

for model in models:model.close()
if registry is not None:registry.close()
rows=world.gather((g,comm.rank,[k for _,k,_ in specs],np.array(traces)),root=0)
if world.rank==0:
    controls={k:row[3][index] for row in rows if row[0] in (1,2) for index,k in enumerate(row[2])}
    for group,rank,keys,traces in rows:
        for index,k in enumerate(keys):np.testing.assert_array_equal(traces[index],controls[k])
    report={'ranks_per_group':n,'xc':a.xc,'second_xc':a.second_xc,'steps':16,
        'mixed_modes_and_independent_controls_array_equal':True,'siblings_unchanged':True,
        'explicit_memory_restore_array_equal':True,'mode_preserved':True,'stale_token_rejected_locally':True,'GS':'existing symlink; no new GS/copies/disk electron state'}
    np.savez(Path(a.output)/'traces.npz',first=controls[1],second=controls[2])
    (Path(a.output)/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
