from pathlib import Path
import argparse,tempfile,subprocess,os,re,math,json,signal
P=argparse.ArgumentParser();P.add_argument('--binary',required=True);P.add_argument('--mpi');P.add_argument('--auto-input',action='store_true');P.add_argument('--np-flag',default='-n');args=P.parse_args()
root=Path(__file__).resolve().parents[2];energies={};results=[]
# One SCF iteration exercises real init_pp/input_pp and MPI broadcasts; this is
# a reader smoke, not a converged atomic/solid accuracy calculation.
for name,species,soc in [('legacy-Si',['Si'],'n'),('W-Se',['W','Se'],'y'),('W-Si',['W','Si'],'y'),('Si-W',['Si','W'],'y'),('Se',['Se'],'y')]:
 with tempfile.TemporaryDirectory(prefix='ky-reader-') as tmp:
  d=Path(tmp);sp=species;z={'W':74,'Se':34,'Si':14};ne={'W':12,'Se':6,'Si':4};n=sum(ne[x] for x in sp);nst=n*2 if soc=='y' else n
  atoms=['W','Si'] if name in ('W-Si','Si-W') else sp
  coord='\n'.join(" '%s' %s %s %s %s"%(x,*(('0','0','0') if j==0 else ('6','6','6')),sp.index(x)+1) for j,x in enumerate(atoms))
  parallel='nproc_k=2' if args.mpi else 'nproc_k=1'
  s=f'''&calculation theory='dft' /
&control sysname='smoke',yn_restart='n',checkpoint_interval=-1 /
&units unit_system='au' /
&system yn_periodic='y',spin='{ 'noncollinear' if soc=='y' else 'unpolarized'}',yn_spinorbit='{soc}',yn_symmetry='nnn',
 al=12d0,12d0,12d0,nelem={len(sp)},natom={len(atoms)},nelec={n},nstate={nst},temperature_k=300d0 /
&pseudo
'''+''.join(f" izatom({i+1})={z[x]},file_pseudo({i+1})='{root}/testsuites/pseudo/{x}_rps.dat',lloc_ps({i+1})={2 if x=='Si' else 0}\n" for i,x in enumerate(sp))+f'''/
&functional xc='pz' /
&rgrid num_rgrid=8,8,8 /
&kgrid num_kgrid=2,2,2 /
&parallel {parallel},nproc_ob=1,nproc_rgrid=1,1,1 /
&atomic_coor
{coord}
/
&scf method_init_density='wf',method_init_wf='random',nscf=1,ncg=1,ncg_init=1,nscf_init_no_diagonal=10 /
'''
  if args.auto_input:
   s=s.replace(f"spin='{ 'noncollinear' if soc=='y' else 'unpolarized'}',yn_spinorbit='{soc}',",'')
   for i,x in enumerate(sp):
    if x!='Si':s=s.replace(f",lloc_ps({i+1})=0",'')
  assert all(len(l)<=256 for l in s.splitlines());(d/'inputfile').write_text(s)
  env={k:v for k,v in os.environ.items() if not k.startswith('SALMON_')};env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
  cmd=[args.binary]
  if args.mpi:cmd=[args.mpi,args.np_flag,'2']+cmd
  with (d/'inputfile').open() as inp:
   proc=subprocess.Popen(cmd,cwd=d,env=env,stdin=inp,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,start_new_session=True)
   try:text,_=proc.communicate(timeout=90)
   except subprocess.TimeoutExpired:
    os.killpg(proc.pid,signal.SIGTERM)
    try:proc.communicate(timeout=10)
    except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.communicate()
    raise RuntimeError(name+' reader smoke timeout')
  assert proc.returncode==0,(name,proc.returncode,text[-4000:]);assert 'end SALMON' in text,(name,text[-2000:])
  e=[float(x.replace('D','E')) for x in re.findall(r'iter=\s*\d+\s+Total Energy=\s*(\S+)',text)];assert e and all(math.isfinite(x) for x in e),(name,text[-2000:])
  if args.auto_input:
   variables=(d/'variables.log').read_text()
   assert f"spin={'noncollinear' if soc=='y' else 'unpolarized'}" in variables
   assert f"yn_spinorbit={soc}" in variables
  energies[name]=e[-1];results.append(dict(name=name,exit=proc.returncode,energy=e[-1],MPI=2 if args.mpi else 1))
assert abs(energies['W-Si']-energies['Si-W'])<=1e-7,energies
print(json.dumps(results,indent=2));print('KY_READER_FULL_SERIAL_MPI_SMOKE_PASS')
