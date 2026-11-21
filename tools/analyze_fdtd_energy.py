"""Manual real-solver observation validation, not required CI.

MODE COMPLETED_NATIVE_FIXTURE NEW_OUTPUT --checker WORKSPACE_CHECKER
Does not run SALMON, SCF or MPI and does not copy GS. Explicitly builds/audits a
small C observation reference. LD scope is the fixed isotropic slab fixture
(wp20,w025,gamma2,epsilon_inf1, transverse periodicity, uniform transverse mode);
J is reconstructed from actual native E with the identical extracted LD kernel,
and compared to a native final LD checkpoint if present. Surface mode checks
a fixed z-invariant dielectric grating against its matched vacuum reference. TDDFT scope is native
2-point 1D Weyl outputs. run.json provides dt/steps and (LD) x/area.
Native profiles/observations must be complete and unmodified. No real-material
optical accuracy is certified; see the energy-budget validation plan.
"""
from pathlib import Path
import ctypes as c,json,sys,re
import numpy as np
from salmonpy.analysis import absorbed_energy_density,reflection_transmission
from salmonpy.fdtd import energy_balance
from salmonpy.ld import LorentzDrude
import argparse,hashlib,os,subprocess
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('mode',choices=('ld','tddft','surface'))
parser.add_argument('directory',type=Path,help='completed native validation fixture directory')
parser.add_argument('output',type=Path,help='NEW report directory; no solver is started')
parser.add_argument('--checker',type=Path,required=True)
parser.add_argument('--cc',default='cc')
parser.add_argument('--vacuum',type=Path,help='surface mode: matched vacuum control')
parser.add_argument('--left-index',type=int,default=-8,help='TDDFT native macro-grid port index')
parser.add_argument('--right-index',type=int,default=9)
args=parser.parse_args()
folder=args.directory.resolve();root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
source=Path(__file__).resolve().parents[1]/'tests/native/test_energy_budget.c'
build=root/'reference';build.mkdir()
subprocess.run([sys.executable,str(args.checker),'--source-root',str(source.parent)],check=True)
obj=build/'reference.o'
subprocess.run([args.cc,'-std=c99','-O2','-fPIC','-c',str(source),'-o',str(obj)],check=True)
subprocess.run([sys.executable,str(args.checker),'--object-tree',str(build)],check=True)
shared=build/('reference.dylib' if sys.platform=='darwin' else 'reference.so')
subprocess.run([args.cc,'-shared',str(obj),'-lm','-o',str(shared)],check=True)
(build/'audit.json').write_text(json.dumps({'sha256':hashlib.sha256(shared.read_bytes()).hexdigest(),
    'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'object_audit':'passed'},indent=2)+'\n')
lib=c.CDLL(str(shared))
lib.energy_center_h.argtypes=[c.c_int,c.c_int,c.c_void_p,c.c_void_p]
lib.energy_split_x.argtypes=[c.c_int]+[c.c_void_p]*5
lib.energy_field_sum.argtypes=[c.c_int,c.c_void_p,c.c_double,c.c_void_p]
lib.energy_weighted_sum.argtypes=[c.c_int,c.c_void_p,c.c_void_p,c.c_void_p]
def split(a):
 a=np.ascontiguousarray(a);o=[np.empty((len(a),3)) for _ in range(4)]
 assert lib.energy_split_x(len(a),a.ctypes.data,*(v.ctypes.data for v in o))==0
 return o

def field_sum(a,volume):
 a=np.ascontiguousarray(a);v=c.c_double()
 assert lib.energy_field_sum(len(a),a.ctypes.data,volume,c.byref(v))==0
 return v.value

def weighted(a,w):
 a=np.ascontiguousarray(a);w=np.ascontiguousarray(w);v=c.c_double()
 assert lib.energy_weighted_sum(len(a),a.ctypes.data,w.ctypes.data,c.byref(v))==0
 return v.value

def check(name,t,field,charge,volumes,left,right,delta,area,notes):
 inc,hi,ref,hr=split(left);tra,ht,back,hb=split(right)
 port=reflection_transmission(t,inc,hi,ref,hr,tra,ht,normal=(1,0,0),units='atomic_gaussian')
 work=absorbed_energy_density(t,field,charge,units='atomic_gaussian',current_convention='charge',field_kind='total')
 report=energy_balance(work,volumes,port,material_ids=tuple(name for _ in volumes),transverse_area=area,field_energy_change=delta,units='atomic_gaussian',rtol=.01)
 data={key:getattr(report,key) for key in ('incident_energy','reflected_energy','transmitted_energy','material_work','field_energy_change','residual','relative_residual','within_tolerance')}
 data.update(notes=notes,points=len(volumes),steps=len(t)-1)
 print(name,json.dumps(data,indent=2),flush=True)
 (root/(name+'-budget.json')).write_text(json.dumps(data,indent=2)+'\n')
 np.savez(root/(name+'-observations.npz'),time=t,macro_E=field,charge_j=charge,point_volumes=volumes,
          material_ids=np.array([name]*len(volumes)),normal=np.array([1,0,0]),area=np.array(area),
          field_energy_change=np.array(delta),other_outgoing_energy=np.array(0.),
          incident_E=inc,incident_H=hi,reflected_E=ref,reflected_H=hr,transmitted_E=tra,transmitted_H=ht)
 np.savez(root/(name+'-spectra.npz'),omega=port.omega,R=port.R,T=port.T,valid=port.valid,work=work)
 return data

if args.mode=='surface':
 if args.vacuum is None:raise ValueError('surface mode requires a matched --vacuum control')
 run=json.loads((folder/'run.json').read_text());ny=len(run['y']);dt=run['dt']
 def read_planes(directory):
  rows=[np.loadtxt(directory/f'obs{i+1}_at_point_rt.data') for i in range(4*ny)]
  if any(not np.array_equal(a[:,0],rows[0][:,0]) for a in rows):raise ValueError('plane clocks differ')
  t=np.r_[0.,rows[0][:,0]]
  raw=np.concatenate([np.zeros((1,4*ny,6)),np.stack([a[:,1:] for a in rows],axis=1)])
  raw=np.ascontiguousarray(raw.reshape(len(t),4,ny,6));out=np.empty((len(t),2,ny,6))
  fn=lib.energy_periodic_planes;fn.argtypes=[c.c_int,c.c_int,c.c_void_p,c.c_void_p]
  if fn(len(t),ny,raw.ctypes.data,out.ctypes.data):raise ValueError('native plane collocation failed')
  return t,out
 t,material=read_planes(folder);tv,vacuum=read_planes(args.vacuum)
 if not np.array_equal(t,tv):raise ValueError('reference time differs')
 from salmonpy.analysis import surface_flux
 dA=np.full(ny,run['spacing'][1]*run['z_extent'])
 inc=np.ascontiguousarray(vacuum[:,0,:]);ref=np.ascontiguousarray(material[:,0,:]-inc)
 tra=np.ascontiguousarray(material[:,1,:])
 ports=reflection_transmission(t,inc[:,:,:3],inc[:,:,3:],ref[:,:,:3],ref[:,:,3:],tra[:,:,:3],tra[:,:,3:],
       normal=(1,0,0),units='atomic_gaussian',area_weights=dA)
 def native_field(name):
  import struct
  directory=folder/'data_for_restart'
  info=struct.unpack('=11i',(directory/(name+'_0_inf.bin')).read_bytes())
  lo,hi=info[:5],info[5:10];shape=tuple(b-a+1 for a,b in zip(lo,hi))
  return np.fromfile(directory/(name+'_0.bin'),dtype=np.float64).reshape(shape,order='F')[:,:,:,0,0],lo
 ez,lo=native_field('ez_x');ez=ez+native_field('ez_y')[0]
 hx=native_field('hx_y')[0]+native_field('hx_z')[0]
 hy=native_field('hy_z')[0]+native_field('hy_x')[0]
 nx,ny,nz=run['grid'];dx,dy,dz=run['spacing']
 xx=(np.arange(nx)+lo[0]-.5)*dx;yy=(np.arange(ny)+lo[1]-.5)*dy
 epsilon=np.ones((nx,ny,nz));epsilon[(abs(xx[:,None,None])<10)&(abs(yy[None,:,None])<16)&np.ones((1,1,nz),dtype=bool)]=run['material_epsilon']
 arrays=[np.ascontiguousarray(a) for a in (ez,hx,hy,epsilon)]
 first=int(np.where(xx==run['planes'][1])[0][0]);last=int(np.where(xx==run['planes'][3])[0][0])
 fn=lib.energy_yee_te_end
 fn.argtypes=[c.c_int]*5+[c.c_double]*4+[c.c_void_p]*5
 stats=np.empty(2)
 if fn(nx,ny,nz,first,last,dx,dy,dt,dx*dy*dz,*(a.ctypes.data for a in arrays),stats.ctypes.data):
  raise ValueError('native field energy reconstruction failed')
 budget=energy_balance([0.],[1.],ports,material_ids=('instantaneous-dielectric',),
       transverse_area=run['area'],field_energy_change=stats[0],units='atomic_gaussian')
 # Field averaging is deliberately the incorrect comparison, to expose missed orders.
 average=reflection_transmission(t,inc[:,:,:3].mean(axis=1),inc[:,:,3:].mean(axis=1),
   ref[:,:,:3].mean(axis=1),ref[:,:,3:].mean(axis=1),tra[:,:,:3].mean(axis=1),tra[:,:,3:].mean(axis=1),
   normal=(1,0,0),units='atomic_gaussian')
 report={key:getattr(budget,key) for key in ('incident_energy','reflected_energy','transmitted_energy',
       'material_work','field_energy_change','residual','relative_residual','within_tolerance')}
 report.update(area=run['area'],samples_per_plane=ny,z_invariance=True,
       mean_field_transmitted_energy=average.transmitted_fluence*run['area'],
       field_storage_scope='TE Yee H recentered in time, endpoint half weights; initial fields zero',
       profile='fixed transversely periodic lossless dielectric bar; not general optical accuracy certification')
 (root/'surface-budget.json').write_text(json.dumps(report,indent=2)+'\n')
 np.savez(root/'surface-observations.npz',time=t,macro_E=np.zeros((len(t),1,3)),charge_j=np.zeros((len(t),1,3)),
   point_volumes=np.ones(1),material_ids=np.array(['instantaneous-dielectric']),normal=np.array([1,0,0]),
   area=np.array(run['area']),area_weights=dA,field_energy_change=np.array(stats[0]),other_outgoing_energy=np.array(0.),
   incident_E=inc[:,:,:3],incident_H=inc[:,:,3:],reflected_E=ref[:,:,:3],reflected_H=ref[:,:,3:],
   transmitted_E=tra[:,:,:3],transmitted_H=tra[:,:,3:])
 np.savez(root/'surface-spectra.npz',omega=ports.omega,R=ports.R,T=ports.T,valid=ports.valid,
       mean_field_R=average.R,mean_field_T=average.T)
 print(json.dumps(report,indent=2));sys.exit(0)

if args.mode=='ld':
 case=folder.name
 from salmonpy.bands import _input_groups,_scalar
 body=_input_groups((folder/'inputfile').read_text())['maxwell']
 count=int(_scalar(body,'media_num','0'))
 vacuum=count==0
 if not vacuum:
  expected={'media_num':'1','media_type(1)':'lorentz-drude',
      'omega_p_ld(1)':'20','pole_num_ld(1)':'1','f_ld(1,1:1)':'1',
      'omega_ld(1,1:1)':'25','gamma_ld(1,1:1)':'2','epsilon_em(1)':'1'}
  for key,value in expected.items():
   actual=_scalar(body,key,'1' if key=='epsilon_em(1)' else None)
   if value=='lorentz-drude':
    if actual!=value:raise ValueError('LD reference requires the fixed scalar slab fixture: '+key)
   elif actual is None or float(actual.replace('d','e'))!=float(value):
    raise ValueError('LD reference requires the recorded wp20/w025/gamma2/epsilon_inf1 fixture: '+key)
 run=json.loads((folder/'run.json').read_text());xs=np.array(run['x']);dt=run['dt'];area=run['area']
 tables=[np.loadtxt(folder/f'obs{i+1}_at_point_rt.data') for i in range(len(xs))]
 t=np.r_[0.,tables[0][:,0]]
 assert all(np.array_equal(a[:,0],tables[0][:,0]) for a in tables), 'port times differ'
 assert np.allclose(np.diff(t),dt,rtol=1e-8,atol=1e-10), 'LD reconstruction requires every native step'
 raw=np.concatenate((np.zeros((1,len(xs),6)),np.stack([a[:,1:] for a in tables],axis=1)))
 centered=np.empty((len(t),len(xs)-1,6))
 assert lib.energy_center_h(len(t),len(xs),raw.ctypes.data,centered.ctypes.data)==0
 cx=xs[1:];left=1;right=len(cx)-2
 volume=(xs[1]-xs[0])*area
 delta=field_sum(centered[-1,left:right+1],volume)-field_sum(centered[0,left:right+1],volume)
 matter=(abs(xs)<10) if not vacuum else np.array([i==len(xs)//2 for i in range(len(xs))])
 E=np.ascontiguousarray(raw[:,:,0:3][:,matter]);j=np.zeros_like(E)
 if not vacuum:
  model=LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[25.],damping=[2.])
  for i in range(E.shape[1]):j[:,i]=model.respond(E[:,i],dt=dt).current
 if (folder/'data_for_restart/rjz_ld_0_inf.bin').is_file() and not vacuum:
  import struct
  info=struct.unpack('=11i',(folder/'data_for_restart/rjz_ld_0_inf.bin').read_bytes())
  lo,hi=info[:5],info[5:10]
  native=np.fromfile(folder/'data_for_restart/rjz_ld_0.bin',dtype=np.float64).reshape(
      tuple(b-a+1 for a,b in zip(lo,hi)),order='F')
  dx=xs[1]-xs[0];xind=np.rint(xs[matter]/dx+.5).astype(int)-lo[0]
  # Transverse fixture is uniform; select the origin-aligned y/z nodes.
  np.testing.assert_allclose(j[-1,:,2],native[xind,-lo[1],-lo[2],0,0],rtol=1e-10,atol=1e-18)
  (root/'native-current-parity.json').write_text(json.dumps({'native_final_LD_current_matches':True})+'\n')
 check(case,t,E,j,np.full(E.shape[1],volume),centered[:,left],centered[:,right],delta,area,
       ('Native SALMON E/H, spatially centered H; vacuum control.' if vacuum else 'Native SALMON E/H; H spatially centered along x. LD current reconstructed from observed E with the identical extracted native LD kernel; no heat-only assumption.'))
if args.mode=='ld':
 sys.exit(0)
run=json.loads((folder/'run.json').read_text());dt=run['dt']
files=sorted((folder/'rt_RT_Ac').glob('rt_Ac_*.data'))
assert len(files)==run['steps']+1, 'complete native profiles required'
assert [int(p.stem.rsplit('_',1)[1]) for p in files]==list(range(run['steps']+1)), 'missing native profile'
frames=np.stack([np.loadtxt(p) for p in files]);indices=frames[0,:,0].astype(int)
assert np.array_equal(frames[:,:,0],np.broadcast_to(indices,frames[:,:,0].shape)), 'grid order changed'
t=np.arange(len(frames))*dt
left=int(np.where(indices==args.left_index)[0][0]);right=int(np.where(indices==args.right_index)[0][0]);m=(indices>=1)&(indices<=2)
from salmonpy.bands import _input_groups,_scalar
spacing=float(_scalar(_input_groups((folder/'input.inp').read_text())['multiscale'],'hx_m').replace('d','e'))
weights=np.full(right-left+1,spacing);weights[[0,-1]]=spacing/2
initial=weighted(frames[0,left:right+1,15],weights)
final=weighted(frames[-1,left:right+1,15],weights)
result=check(folder.name,t,np.ascontiguousarray(frames[:,:,6:9][:,m]),np.ascontiguousarray(frames[:,:,12:15][:,m]),
 np.full(2,spacing),frames[:,left,6:12],frames[:,right,6:12],final-initial,1.,
 'Native full TDDFT+Weyl FDTD, 2 macro points, 8 MPI ranks; all E/H/j and field density read directly from native outputs; measured initial field energy subtracted (zero for delayed incidence).')
