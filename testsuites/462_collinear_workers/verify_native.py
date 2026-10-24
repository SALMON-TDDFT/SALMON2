"""Minimal scalar GS/learned RT ownership regression; not converged production acceptance."""
from pathlib import Path
import argparse,os,subprocess,sys,json,time,re,itertools,math
p=argparse.ArgumentParser();p.add_argument('--binary',required=True);p.add_argument('--output',required=True);p.add_argument('--pseudo',required=True);p.add_argument('--threads',type=int,default=2);p.add_argument('--unique-output',action='store_true');a=p.parse_args()
out=Path(a.output).resolve()
if a.unique_output:out=out/('run-'+str(time.time_ns()))
out.mkdir(parents=True,exist_ok=False);exe=Path(a.binary).resolve();pp=Path(a.pseudo).resolve()
atoms=[(0,0,0),(1,1,1),(2,0,2),(0,2,2),(2,2,0),(3,1,3),(1,3,3),(3,3,1)]
sym=[]
for swap in [False,True]:
 for sx,sy in itertools.product([-1,1],repeat=2):
  rot=[[0,0,0],[0,0,0],[0,0,1]];rot[0][int(swap)]=sx;rot[1][1-int(swap)]=sy
  for tr in atoms:
   image={tuple((sum(rot[i][j]*v[j] for j in range(3))+tr[i])%4 for i in range(3)) for v in atoms}
   if image==set(atoms):sym += [' '.join(map(str,rot[i]+[tr[i]/4])) for i in range(3)]
assert len(sym)==96
base="""&calculation theory='THEORY',yn_dc='n' /
&control sysname='Si',yn_restart='n',CONTROL /
&units unit_system='au' /
&system yn_periodic='y',yn_symmetry='yyn',al=10.26d0,10.26d0,10.26d0,
 nstate=32,nelec=32,nelem=1,natom=8 /
&pseudo izatom(1)=14,file_pseudo(1)='PSEUDO',lmax_ps(1)=2,lloc_ps(1)=2 /
&functional xc='hse06',yn_exx_fft_orbits='y',exx_pre_scf_threshold=0d0,
 exx_factor_mode='FACTOR',exx_factor_exact_interval=2,exx_factor_n_teach=2,
 exx_factor_history_frames=3,exx_factor_warmup_steps=12,
 exx_worker_rprime_groups=GROUPS,hse_block_rows=1,hse_fft_layout='strided',yn_hse_profile='y',hse_sr_tolerance=0d0 /
&rgrid num_rgrid=4,4,4 /
&kgrid num_kgrid=2,2,2 /
&parallel nproc_k=2,nproc_ob=1,nproc_rgrid=1,1,1 /
&scf nscf=4,ncg=4,ncg_init=4,nscf_init_no_diagonal=0,method_init_wf='random',threshold=1d-8 /
&tgrid nt=32,dt=.03d0 /
&emfield trans_longi='tr',ae_shape1='impulse',e_impulse=1d-4,epdir_re1=0d0,0d0,1d0 /
&analysis out_rt_energy_step=1,nenergy=10,de=.01d0 /
&propagation propagator='hse_taylor4',yn_predictor_corrector='y' /
&atomic_coor
ATOMS
/
""".replace('PSEUDO',str(pp)).replace('ATOMS','\n'.join(" 'Si' "+' '.join(str(x*2.565) for x in v)+' 1' for v in atoms))
wrapper=out/'rank.py';wrapper.write_text('''import os,sys,subprocess,json,resource
from pathlib import Path
r=int(os.environ['OMPI_COMM_WORLD_RANK'])
with open('stdout.rank-%d'%r,'x') as o,open('stderr.rank-%d'%r,'x') as e:
 q=subprocess.run([sys.argv[1]],stdin=sys.stdin,stdout=o,stderr=e)
Path('exit.rank-%d'%r).write_text(str(q.returncode))
Path('memory.rank-%d.json'%r).write_text(json.dumps({'maxrss':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}))
sys.exit(q.returncode if q.returncode>=0 else 128-q.returncode)
''')
env={k:v for k,v in os.environ.items() if not k.startswith(('SALMON_','OMPI_'))};env.update(OMP_NUM_THREADS=str(a.threads),OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',OMPI_MCA_btl='self,sm',OMPI_MCA_pml='ob1')
report={}
def run(d,text,np):
 d.mkdir(parents=True);(d/'inputfile').write_text(text);(d/'sym.dat').write_text('\n'.join(sym)+'\n')
 with (d/'inputfile').open() as f,(d/'launcher.log').open('x') as log:
  q=subprocess.run(['mpiexec','--oversubscribe','--bind-to','none','-n',str(np),sys.executable,str(wrapper),str(exe)],cwd=d,env=env,stdin=f,stdout=log,stderr=subprocess.STDOUT,timeout=300)
 logs=[(d/('stdout.rank-%d'%r)).read_text() for r in range(np)]
 assert q.returncode==0,(str(d),q.returncode,logs[0][-500:])
 assert all((d/('exit.rank-%d'%r)).read_text()=='0' for r in range(np)),str(d)
 assert 'end SALMON' in logs[0],str(d)
 for r in range(np):
  e=(d/('stderr.rank-%d'%r)).read_text()
  assert not re.search(r'ERROR STOP|SIGSEGV|NaN|Infinity',e,re.I),(str(d),r,e)
 if np>2:
  assert all('EXX_SERVICE_COMPUTE' in x and 'EXX_SERVICE_STOP' in x for x in logs[2:]),str(d)
  teams=[int(x) for log in logs for x in re.findall(r'WAVE_FFT_OMP_THREADS\s+(\d+)',log)]
  assert teams and all(x==a.threads for x in teams),(str(d),teams)
 return logs
for label,np,groups in [('owners2',2,1),('helpers4',4,2),('helpers8',8,4)]:
 gs=out/label/'gs';rt=out/label/'rt'
 text=base.replace('GROUPS',str(groups)).replace('THEORY','dft').replace('CONTROL',"write_gs_restart_data='wfn'").replace('FACTOR','off')
 logs=run(gs,text,np)
 energies=[float(x) for x in re.findall(r'Total Energy=\s*(\S+)',logs[0])];assert len(energies)==4,(label,energies)
 text=base.replace('GROUPS',str(groups)).replace('THEORY','tddft_response').replace('CONTROL',"directory_read_data='"+str(gs/'data_for_restart')+"/'").replace('FACTOR','learned')
 logs=run(rt,text,np)
 flags=re.findall(r'FACTOR_HISTORY step/predicted:\s+\d+\s+([TF])',logs[0]);assert 'T' in flags and 'F' in flags,(label,flags)
 def table(path):
  rows=[]
  for line in path.read_text().splitlines():
   if not line.strip() or line.lstrip().startswith('#'):continue
   try:row=list(map(float,line.split()))
   except ValueError:continue
   assert all(math.isfinite(x) for x in row),(path,line)
   rows.append(row)
  return rows
 cur=table(rt/'Si_rt.data');en=table(rt/'Si_rt_energy.data');assert len(cur)==32 and len(en)==33,(len(cur),len(en))
 eigen=table(gs/'Si_eigen.data');assert eigen
 svd=[]
 for log in logs[:2]:
  values=[]
  for line in log.splitlines():
   if 'FACTOR_SVD_QUALITY n/residual/unit_u/unit_v/smin/smax:' not in line:continue
   vals=list(map(float,line.split(':',1)[1].split()));n=int(vals[0])
   assert n>0 and all(math.isfinite(x) for x in vals)
   assert vals[1]<=1e-10 and max(vals[2:4])<=1e-10*math.sqrt(n)
   assert vals[4]>0 and vals[5]>=vals[4]
   values.append(vals)
  assert len(values)>=6,(label,'SVD records',len(values))
  svd+=values
 report[label]={'MPI':np,'OMP':a.threads,'GS_energy':energies,'eigen':eigen,'current':cur,'energy':en,'predicted':flags.count('T'),'exact':flags.count('F'),'SVD_records':len(svd),'SVD_max_residual':max(x[1] for x in svd)}
 if label!='owners2':
  ref=report['owners2']
  ed=max(abs(x-y) for x,y in zip(energies,ref['GS_energy']));assert ed<1e-6,(label,'GSenergy',ed)
  for key,tol in [('eigen',1e-8),('current',1e-7),('energy',1e-7)]:
   aa=report[label][key];bb=ref[key];assert len(aa)==len(bb)
   diff=max(abs(x-y) for ra,rb in zip(aa,bb) for x,y in zip(ra,rb));assert diff<tol,(label,key,diff)
   report[label][key+'_maxabs']=diff
 (out/'summary.json').write_text(json.dumps(report,indent=2))
for label,old,new in [('reject-no-symmetry',"yn_symmetry='yyn'","yn_symmetry='nnn'"),
                      ('reject-rt-full',"propagator='hse_taylor4'","propagator='hse_taylor4_full'")]:
 d=out/label;d.mkdir()
 text=base.replace('GROUPS','2').replace('THEORY','tddft_response').replace('FACTOR','learned')
 text=text.replace('CONTROL',"directory_read_data='"+str(out/'owners2/gs/data_for_restart')+"/'").replace(old,new)
 (d/'inputfile').write_text(text);(d/'sym.dat').write_text('\n'.join(sym)+'\n')
 with (d/'inputfile').open() as f,(d/'launcher.log').open('x') as log:
  q=subprocess.run(['mpiexec','--oversubscribe','--bind-to','none','-n','4',sys.executable,str(wrapper),str(exe)],cwd=d,env=env,stdin=f,stdout=log,stderr=subprocess.STDOUT,timeout=30)
 assert q.returncode!=0,(label,'unsupported layout accepted')
 reason='FFT orbits require native CPU dense-k hybrid symmetry with all q points' if label=='reject-no-symmetry' else 'Factor history requires native static HSE Taylor4 predictor-corrector'
 assert reason in ''.join(x.read_text() for x in list(d.glob('stdout.rank-*'))+list(d.glob('stderr.rank-*'))),label
print('SI_GS_RT_WORKERS_PASS',str(out))
