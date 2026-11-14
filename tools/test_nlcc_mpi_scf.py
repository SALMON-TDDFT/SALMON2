"""Compare new serial/MPI NLCC forces with space/orbital/k replicas.

Only scalar reports are exported from disposable native working directories.
"""
import argparse,hashlib,json,os
from pathlib import Path
import numpy as np
from salmonpy.dft import DFT,Pseudopotential,SCFSettings,_read_status,_read_forces
from measure_nlcc_force_work import run_scalar_only,write_json

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('serial_record',type=Path);p.add_argument('mpi_record',type=Path)
p.add_argument('pseudo',type=Path);p.add_argument('output',type=Path)
a=p.parse_args();engines={}
for mode,path in [('serial',a.serial_record),('mpi',a.mpi_record)]:
 record=json.loads(path.read_text());exe=Path(record['executable'])
 if record.get('status')!='built-and-audited' or hashlib.sha256(exe.read_bytes()).hexdigest()!=record['sha256']:
  raise ValueError('audited executable record mismatch')
 if record.get('build_configuration',{}).get('mpi_enabled',False) is not (mode=='mpi'):
  raise ValueError('runtime MPI capability does not match test role')
 engines[mode]=exe
raw=a.pseudo.read_bytes()
if hashlib.sha256(raw).hexdigest()!='39822757f53f36e3bf3bfb779356152a8d3f21199c7db9dd5a931e5d18c45282':raise ValueError('wrong test UPF')
pseudo=Pseudopotential(14,2,'.upf',raw)
base=np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0],[.25,.25,.25],[.25,.75,.75],[.75,.25,.75],[.75,.75,.25]],float)
base[0,0]+=.12/10.26
model=DFT((14,)*8,np.eye(3)*10.26,base,(pseudo,),32,32,(16,)*3,(2,)*3,'libxc_pbe')
a.output.mkdir(parents=True,exist_ok=False)
write_json(a.output/'provenance.json',{'serial':json.loads(a.serial_record.read_text()),'mpi':json.loads(a.mpi_record.read_text()),'pseudo_sha256':pseudo.sha256,'threads':1})
results=[]
cases=[('serial',1,(1,1,1),1,1),('mpi-r',2,(2,1,1),1,1),('mpi-ob',2,(1,1,1),2,1),('mpi-k',2,(1,1,1),1,2),('mpi-r-ob',4,(2,1,1),2,1)]
for name,ranks,rgrid,ob,k in cases:
 text=model.input_text(SCFSettings(max_iterations=600,tolerance=1e-11)).replace("write_gs_restart_data='all'","write_gs_restart_data='no'")
 text+=f'\n&parallel\n nproc_rgrid={rgrid[0]},{rgrid[1]},{rgrid[2]}, nproc_ob={ob}, nproc_k={k}\n/\n'
 mode='serial' if name=='serial' else 'mpi'
 launcher=() if mode=='serial' else ('/opt/homebrew/bin/mpiexec','--bind-to','none','-n',str(ranks))
 dest=a.output/name
 env={key:value for key,value in os.environ.items() if key!='SALMONPY_TDDFT_WORKER'}
 env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
 code=run_scalar_only(engines[mode],text,raw,dest,env,launcher=launcher)
 if code or 'end SALMON' not in (dest/'stdout.log').read_text():raise RuntimeError('MPI test failed: '+name)
 status=_read_status(dest/'gs_status.dat');force=_read_forces(dest/'gs_forces.dat',8)
 row=dict(name=name,ranks=ranks,rgrid=rgrid,orbital_groups=ob,k_groups=k,**status,forces=force.tolist())
 if results:
  row['energy_difference']=row['energy_hartree']-results[0]['energy_hartree']
  row['max_force_difference']=float(np.max(np.abs(force-np.asarray(results[0]['forces']))))
  if abs(row['energy_difference'])>1e-8 or row['max_force_difference']>1e-8:
   write_json(a.output/'results.json',results+[row]);raise AssertionError('MPI replication mismatch: '+str(row))
 results.append(row);write_json(a.output/'results.json',results)
 print(json.dumps({k:v for k,v in row.items() if k!='forces'}),flush=True)
print('NLCC_MPI_SCF_PASS')
