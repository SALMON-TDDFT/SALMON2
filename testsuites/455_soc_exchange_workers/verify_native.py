"""Short native SOC owner/helper equivalence, not physical GS acceptance."""
import argparse,json,os,re,subprocess,sys,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--binary',required=True);p.add_argument('--output',required=True);p.add_argument('--mpiexec',default='mpiexec');args=p.parse_args()
root=Path(__file__).resolve().parents[2];out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False);exe=Path(args.binary).resolve()
base="""&calculation theory='dft',yn_dc='n' /
&control sysname='WSe2',yn_restart='n',write_gs_restart_data='no' /
&units unit_system='au' /
&system yn_periodic='y',spin='noncollinear',yn_spinorbit='y',yn_symmetry='yny',
 al=6.2171989500187843d0,10.768504462196411d0,47.621098340569411d0,
 nelem=2,natom=6,nelec=48,nstate=48,temperature=-1d0 /
&pseudo izatom=74,34,
 file_pseudo(1)='WPATH',
 file_pseudo(2)='SEPATH' /
&functional xc='hse06',exx_factor_mode='off',exx_pre_scf_threshold=0d0,hse_block_rows=1,
 hse_fft_layout='strided',yn_hse_profile='y',exx_worker_rprime_groups=GROUPS /
&rgrid num_rgrid=4,6,16 /
&kgrid num_kgrid=4,2,1,dk_shift=0d0,0d0,0d0,yn_gamma_centered='n' /
&parallel nproc_k=4,nproc_ob=1,nproc_rgrid=1,1,1 /
&atomic_red_coor
 'W' 0d0 0d0 .5d0 1
 'W' .5d0 .5d0 .5d0 1
 'Se' 0d0 .3333333333333333d0 .4337301587301587d0 2
 'Se' 0d0 .3333333333333333d0 .5662698412698413d0 2
 'Se' .5d0 .8333333333333333d0 .4337301587301587d0 2
 'Se' .5d0 .8333333333333333d0 .5662698412698413d0 2
/
&scf nscf=2,ncg=4,ncg_init=4,nscf_init_no_diagonal=0,threshold=1d-8,
 method_mixing='simple_dm',mixrate=.2d0,method_init_wf='gauss',method_init_density='wf' /
""".replace('WPATH',str(root/'testsuites/pseudo/W_rps.dat')).replace('SEPATH',str(root/'testsuites/pseudo/Se_rps.dat'))
sym=[]
for tx,ty in [(0.,0.),(.5,.5)]:
 for sx in [1,-1]:
  for sz in [1,-1]:sym.extend([f'{sx} 0 0 {tx}',f'0 1 0 {ty}','0 0 %d 0'%sz])
# Identical SCF-owner rank count and input in every comparison.
wrapper=out/'rank.py';wrapper.write_text("""from pathlib import Path
import os,subprocess,sys,json,resource
r=int(os.environ['OMPI_COMM_WORLD_RANK'])
with open('stdout.rank-%d'%r,'x') as o,open('stderr.rank-%d'%r,'x') as e:
 q=subprocess.run([sys.argv[1]],stdin=sys.stdin,stdout=o,stderr=e)
Path('exit.rank-%d'%r).write_text(str(q.returncode))
peak=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
Path('memory.rank-%d.json'%r).write_text(json.dumps({'maxrss_bytes':peak if sys.platform=='darwin' else peak*1024}))
sys.exit(q.returncode if q.returncode>=0 else 128-q.returncode)
""")
env={k:v for k,v in os.environ.items() if not k.startswith(('SALMON_','OMPI_'))};env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMPI_MCA_btl='self,sm',OMPI_MCA_pml='ob1')
report={}
for label,np,groups in [('owners4',4,1),('helpers8',8,1),('helpers12',12,2)]:
 d=out/label;d.mkdir();(d/'inputfile').write_text(base.replace('GROUPS',str(groups)));(d/'sym.dat').write_text('\n'.join(sym)+'\n');t=time.monotonic()
 with (d/'inputfile').open() as i,(d/'launcher.log').open('x') as log:
  q=subprocess.run([args.mpiexec,'--oversubscribe','--bind-to','none','-n',str(np),sys.executable,str(wrapper),str(exe)],cwd=d,env=env,stdin=i,stdout=log,stderr=subprocess.STDOUT,timeout=180)
 assert q.returncode==0,(label,q.returncode)
 logs=[(d/('stdout.rank-%d'%r)).read_text() for r in range(np)]
 assert all((d/('exit.rank-%d'%r)).read_text()=='0' for r in range(np))
 assert 'end SALMON' in logs[0]
 if np>4:
  assert all('EXX_SERVICE_COMPUTE' in x and 'EXX_SERVICE_STOP' in x for x in logs[4:])
 energy=[float(x) for x in re.findall(r'Total Energy=\s*(\S+)',logs[0])]
 rows=[]
 for line in (d/'WSe2_eigen.data').read_text().splitlines():
  if line.lstrip().startswith('#'):continue
  fields=line.split()
  if fields:
   try:rows.append([float(x) for x in fields])
   except ValueError:pass
 report[label]={'exit':q.returncode,'MPI':np,'groups':groups,'wall_s':time.monotonic()-t,'energy_eV':energy,'eigen_rows':rows,
  'max_rank_rss_bytes':max(json.loads((d/('memory.rank-%d.json'%r)).read_text())['maxrss_bytes'] for r in range(np))}
 if label!='owners4':
  a=report['owners4'];assert len(energy)==len(a['energy_eV']) and energy
  ed=max(abs(x-y) for x,y in zip(energy,a['energy_eV']));assert ed<1e-6,(label,'energy',ed)
  assert len(rows)==len(a['eigen_rows']) and rows
  diff=max(abs(x-y) for ra,rb in zip(rows,a['eigen_rows']) for x,y in zip(ra,rb));assert diff<1e-8,(label,'eigen',diff)
  report[label].update(energy_maxabs_eV=ed,eigen_maxabs=diff)
 (out/'summary.json').write_text(json.dumps(report,indent=2))
# Unsupported explicit layouts must fail collectively before helpers wait.
for label, change in [('reject-reference', "exx_k_algorithm='reference',"),
                      ('reject-groups', ''), ('reject-unfold', '')]:
 d=out/label;d.mkdir()
 text=base.replace('GROUPS','5' if label=='reject-groups' else '1')
 if label=='reject-unfold':text += "&unfolding dm_unfold_option='super' /\n"
 text=text.replace("&functional xc=", "&functional "+change+"xc=")
 (d/'inputfile').write_text(text);(d/'sym.dat').write_text('\n'.join(sym)+'\n')
 with (d/'inputfile').open() as i,(d/'launcher.log').open('x') as log:
  q=subprocess.run([args.mpiexec,'--oversubscribe','--bind-to','none','-n','8',
                    sys.executable,str(wrapper),str(exe)],cwd=d,env=env,stdin=i,
                   stdout=log,stderr=subprocess.STDOUT,timeout=30)
 logs=''.join(x.read_text() for x in d.glob('*rank-*'))
 assert q.returncode!=0 and 'EXX_SERVICE_UNSUPPORTED' in logs,(label,q.returncode)
print('SOC_NATIVE_HELPERS_PASS')
