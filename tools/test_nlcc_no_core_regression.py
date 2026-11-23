"""Old/new scalar-only energy/force regressions for KY and core-free NC-UPF."""
import argparse,hashlib,json,os,re
from pathlib import Path
import numpy as np
from salmonpy.dft import DFT,Pseudopotential,SCFSettings,_read_status,_read_forces
from measure_nlcc_force_work import run_scalar_only,write_json
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('before_record',type=Path);p.add_argument('after_record',type=Path)
p.add_argument('upf',type=Path);p.add_argument('ky',type=Path);p.add_argument('output',type=Path)
a=p.parse_args();engines={}
for name,path in [('before',a.before_record),('after',a.after_record)]:
 record=json.loads(path.read_text());exe=Path(record['executable'])
 if record.get('status')!='built-and-audited' or hashlib.sha256(exe.read_bytes()).hexdigest()!=record['sha256']:raise ValueError('audited executable mismatch')
 engines[name]=exe
raw=a.upf.read_bytes()
if hashlib.sha256(raw).hexdigest()!='39822757f53f36e3bf3bfb779356152a8d3f21199c7db9dd5a931e5d18c45282':raise ValueError('UPF mismatch')
text=re.sub(r'<PP_NLCC\b[^>]*>.*?</PP_NLCC>','',raw.decode(),flags=re.S)
text=re.sub(r'core_correction\s*=\s*"[^"]*"','core_correction="F"',text)
ky=a.ky.read_bytes()
if hashlib.sha256(ky).hexdigest()!='8c51fbf4341de256a868ab446bdb3ed19eef1d1dcc3813cfe5214c07213b8c7c':raise ValueError('KY mismatch')
pseudos=[('core-free-upf','libxc_pbe',Pseudopotential(14,2,'.upf',text.encode())),('ky','pz',Pseudopotential(14,2,'_rps.dat',ky))]
base=np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0],[.25,.25,.25],[.25,.75,.75],[.75,.25,.75],[.75,.75,.25]],float);base[0,0]+=.12/10.26
a.output.mkdir(parents=True,exist_ok=False);results={}
for label,xc,pseudo in pseudos:
 model=DFT((14,)*8,np.eye(3)*10.26,base,(pseudo,),32,32,(16,)*3,(2,)*3,xc)
 inp=model.input_text(SCFSettings(max_iterations=600,tolerance=1e-11)).replace("write_gs_restart_data='all'","write_gs_restart_data='no'")
 runs={}
 for name,exe in engines.items():
  dest=a.output/label/name
  env={k:v for k,v in os.environ.items() if k!='SALMONPY_TDDFT_WORKER'};env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
  code=run_scalar_only(exe,inp,pseudo.data,dest,env,pseudo_name=model._pseudo_name(0))
  if code or 'end SALMON' not in (dest/'stdout.log').read_text():raise RuntimeError(label+' '+name+' failed')
  status=_read_status(dest/'gs_status.dat');f=_read_forces(dest/'gs_forces.dat',8)
  runs[name]={'energy':status['energy_hartree'],'force':f.tolist()}
 e=runs['after']['energy']-runs['before']['energy'];f=float(np.max(np.abs(np.array(runs['after']['force'])-runs['before']['force'])))
 results[label]=dict(xc=xc,pseudo_sha256=pseudo.sha256,energy_difference_hartree=e,max_force_difference_hartree_per_bohr=f,byte_equal_force_reports=(a.output/label/'before/gs_forces.dat').read_bytes()==(a.output/label/'after/gs_forces.dat').read_bytes(),runs=runs)
 write_json(a.output/'results.json',results)
 print(label,e,f,results[label]['byte_equal_force_reports'],flush=True)
 if abs(e)>1e-10 or f>1e-10:raise AssertionError('non-NLCC regression')
print('NLCC_NO_CORE_REGRESSION_PASS')
