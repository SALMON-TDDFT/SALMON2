"""Manual nine-rank LD/SBE/TDDFT protected-update measurement.

Run separate checked/performance MPI jobs; numerical kernels remain Fortran.
Arguments: BUILD EXISTING_GS NEW_OUTPUT MODE POINTS STEPS STRIDE.
Uses read-only GS links solely in this test, never creates or collects GS.
Peak RSS is per process; its sum is not a simultaneous machine peak.
Communication times include waiting and packing, not pure network latency.
"""
import hashlib, json, os, shutil, sys, time, resource
from pathlib import Path
from unittest.mock import patch
import numpy as np
from mpi4py import MPI
from salmonpy.fdtd import FDTD, PlaneWave
from salmonpy.fields import AcosPulse
from salmonpy.geometry import Grid, Geometry, Shape
from salmonpy.ld import LorentzDrude
from salmonpy.distributed import ParallelContext, ParallelPlan, SolverGroup
from salmonpy.material_response import MaterialBinding, ResponseSession
from salmonpy.embedded_tddft import NativeTDDFTMaterialProvider, NativeTDDFTBackend
from salmonpy.coupled import NativeEHField, CoupledStepper
from salmonpy.native_steps import native_step_targets
from salmonpy.materials import bind_material
from salmonpy.examples.mpi_mixed_multiscale import models
from salmonpy.memory import process_memory
world=MPI.COMM_WORLD
record,fixture,output,mode,count,steps,stride=sys.argv[1:]
fixture=Path(fixture).resolve();out=Path(output)
count,steps,stride=map(int,(count,steps,stride));dt=.003
assert mode in ('checked','performance') and world.size==9
assert count in (16,32,64) and steps>0 and 1<=stride<=steps
record_data=json.loads(Path(record).read_text())
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert digest(Path(record_data['library']))==record_data['sha256']
gs=fixture/'data_for_restart';before={str(p.relative_to(gs)):digest(p) for p in gs.rglob('*') if p.is_file()}
def link_copy(source,dest):
 dest.mkdir()
 for p in sorted(source.rglob('*')):
  rel=p.relative_to(source)
  if rel.parts[0]=='data_for_restart':continue
  if p.is_dir():(dest/rel).mkdir(exist_ok=True)
  elif p.is_file():shutil.copy2(p,dest/rel)
 (dest/'data_for_restart').symlink_to((source/'data_for_restart').resolve(),target_is_directory=True)
 return test_manifest(source)
def test_manifest(root):
 records=[]
 for p in root.rglob('*'):
  if p.is_symlink():
   assert p==root/'data_for_restart' and p.resolve()==gs.resolve()
  elif p.is_file():records.append((str(p.relative_to(root)),digest(p)))
 if (root/'data_for_restart').is_symlink():
  records.extend(('data_for_restart/'+name,sha) for name,sha in before.items())
 return tuple(sorted(records))

class ModeSBE:
    backend='SALMON-SBE-velocity'
    def __init__(self,initial):self.initial=initial
    def new_state(self,**kwargs):
        return self.initial.new_state(**kwargs,execution_mode=mode)

class MeasuredComm:
    """Observe Python collectives without adding a collective to the loop."""
    operations={'allgather','bcast','Bcast','Allreduce','allreduce','Alltoallv',
                'Gatherv','Scatterv','gather','scatter','Barrier','barrier'}
    def __init__(self,comm):self.comm=comm;self.calls={};self.seconds={};self.errors=[]
    def __getattr__(self,name):
        method=getattr(self.comm,name)
        if name not in self.operations:return method
        def measured(*args,**kwargs):
            start=None
            try:start=time.perf_counter()
            except Exception as exc:self.errors.append(str(exc))
            result=method(*args,**kwargs)
            try:
                self.calls[name]=self.calls.get(name,0)+1
                if start is not None:self.seconds[name]=self.seconds.get(name,0.)+time.perf_counter()-start
            except Exception as exc:self.errors.append(str(exc))
            return result
        return measured

def peak_memory():
    if os.environ.get('SALMONPY_TEST_RESOURCE_FAILURE')=='1':
        return None,'injected optional resource measurement failure'
    return process_memory()


def output_inventory(root):
    try:
        files=[p for p in root.rglob('*') if p.is_file() and not p.is_symlink()]
        return dict(regular_files=len(files),bytes=sum(p.stat().st_size for p in files),
                    observation_files=sum(p.suffix=='.data' for p in files),
                    observation_bytes=sum(p.stat().st_size for p in files if p.suffix=='.data'),
                    excludes='symlink targets, report.json and launcher log; logical file sizes, not disk traffic'),None
    except Exception as exc:return None,str(exc)


def usage():
    try:
        if os.environ.get('SALMONPY_TEST_RESOURCE_FAILURE')=='1':
            raise OSError('injected optional resource measurement failure')
        value=resource.getrusage(resource.RUSAGE_SELF)
        return dict(user_seconds=value.ru_utime,system_seconds=value.ru_stime,
                    input_blocks=value.ru_inblock,output_blocks=value.ru_oublock),None
    except Exception as exc:return {},str(exc)

try:
    if world.rank==0:out.mkdir(parents=True,exist_ok=False)
    world.Barrier();started=time.perf_counter()
    ids=tuple(range(1,count+1));electronic=count//4
    groups=(SolverGroup('field',(0,1,2,3),role='field'),
            SolverGroup('tddft-left',(4,5),point_ids=ids[:electronic//2]),
            SolverGroup('tddft-right',(6,7),point_ids=ids[electronic//2:electronic]),
            SolverGroup('sbe',(8,),point_ids=ids[electronic:]))
    locations=[(x,y,z) for x in (20,26,32,38) for y in (16,26,36,46) for z in (16,26,36,46)][:count]
    points=dict(zip(ids,locations))
    grid=Grid((64,64,64),(2.,2.,2.))
    geometry=Geometry(grid,(Shape('box',(8.,16.,16.),(25.,-15.,0.),material=1),
                            Shape('box',(8.,16.,16.),(25.,15.,0.),material=2)))
    ld1=LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.])
    ld2=LorentzDrude(plasma_frequency=20.,strength=[.7,.3],resonance=[0.,25.],damping=[2.,3.])
    calc=FDTD(geometry,(ld1,ld2),dt,steps,source=PlaneWave(AcosPulse(.001,20.,.3),0,-36.),
              observations=((25.,-15.,0.),(25.,15.,0.)),observation_stride=stride)
    with patch('salmonpy.embedded_tddft._copy_bundle',side_effect=link_copy),patch('salmonpy.embedded_tddft._manifest',side_effect=test_manifest),ParallelContext.create(world,ParallelPlan(groups)) as context:
        field=session=None;history=[]
        try:
            if context.group.role=='field':
                field=NativeEHField(context,calc,record,out/'field',points=points,grid=(2,2,1),execution_mode=mode)
            provider=NativeTDDFTMaterialProvider(fixture,'rt.inp',max_steps=steps)
            backend=NativeTDDFTBackend(record,out,decomposition=(2,1,1,1,1),shared_model=True,execution_mode=mode)
            tddft=MaterialBinding('synthetic-pz',provider,backend)
            _,sbe=models();bound=bind_material(ModeSBE(sbe),key='synthetic-sbe',identity='two-band-mode-explicit')
            bindings={p:tddft if p<=electronic else bound for p in ids}
            session=ResponseSession(context,bindings,dt=dt,execution_mode=mode)
            coupled=CoupledStepper(session,field);targets=native_step_targets(coupled)
            setup=time.perf_counter()-started
            solver=MeasuredComm(context.solver_comm);context.solver_comm=solver
            leaders=None
            if context.is_leader:leaders=MeasuredComm(context.leaders);context.leaders=leaders
            pre,pre_error=usage();rss_before,rss_basis=peak_memory()
            world.Barrier();begin=time.perf_counter()
            for i in range(steps):
                frame=coupled.step_protected(targets)
                if world.rank==0:history.append(np.stack([frame[k] for k in ('A_over_c','E_at_points','J_charge')]))
            elapsed=time.perf_counter()-begin
            post,post_error=usage();rss_after,rss_basis_after=peak_memory()
            local=dict(rank=world.rank,group=context.group.name,setup_seconds=setup,loop_seconds=elapsed,
                       peak_RSS_bytes=rss_after,peak_RSS_basis=rss_basis_after,pre_loop_peak_RSS_bytes=rss_before,
                       resource_delta={k:post[k]-pre[k] for k in post.keys()&pre.keys()},
                       measurement_errors=[v for v in (pre_error,post_error) if v],
                       solver_calls=solver.calls,solver_seconds=solver.seconds,
                       leaders_calls=leaders.calls if leaders else {},leaders_seconds=leaders.seconds if leaders else {},
                       communication_measurement_errors=solver.errors+(leaders.errors if leaders else []))
            ranks=world.gather(local,root=0)
            communication=session.communication_report()
        finally:
            if session is not None:session.close()
            if field is not None:field.close()
    assert before=={str(p.relative_to(gs)):digest(p) for p in gs.rglob('*') if p.is_file()}
    if world.rank==0:
        values=np.array(history);assert np.isfinite(values).all() and np.max(np.abs(values[:,2]))>0
        np.savez(out/'traces.npz',values=values)
        io,io_error=output_inventory(out)
        metric_errors=[error for r in ranks for error in r['measurement_errors']+r['communication_measurement_errors']]
        metric_errors.extend(r['peak_RSS_basis'] for r in ranks if r['peak_RSS_bytes'] is None)
        if io_error:metric_errors.append(io_error)
        report=dict(schema='salmonpy.mixed-performance-measurement/1',mode=mode,points=count,
                    tddft_points=electronic,sbe_points=count-electronic,steps=steps,dt=dt,
                    grid=[64]*3,observation_stride=stride,ranks=ranks,communication=communication,io=io,
                    simulated_time_au=steps*dt,
                    component_modes=dict(field=mode,LD_in_field=mode,SBE=mode,TDDFT=mode,outer=mode),
                    loop_seconds_max=max(r['loop_seconds'] for r in ranks),
                    summed_rank_peak_RSS_bytes=(sum(r['peak_RSS_bytes'] for r in ranks)
                        if all(r['peak_RSS_bytes'] is not None for r in ranks) else None),
                    calculation_status='passed',resource_measurement_status='partial' if metric_errors else 'passed',
                    resource_measurement_errors=metric_errors,
                    GS_collected=False,GS_bytes_unchanged=True,protected_every_step=True,
                    library_sha256=record_data['sha256'],
                    limitations=['single Mac, nine local ranks; no multi-node claim',
                                 'synthetic PZ/SBE inputs, not physical convergence',
                                 'Python MPI calls only; Fortran-internal collectives excluded',
                                 'sum of rank peaks is not simultaneous physical memory peak'])
        (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        print('mixed measurement PASS',mode,count,steps,report['loop_seconds_max'],flush=True)
except BaseException:
    import traceback
    traceback.print_exc();world.Abort(1)
