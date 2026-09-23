"""Native ELF + normalized Proca integration: EXE GS_DIRECTORY [MPI prefix]."""
from pathlib import Path
import os,sys,tempfile,subprocess
import numpy as np
root=Path(__file__).resolve().parents[2];exe=Path(sys.argv[1]).resolve();gs=Path(sys.argv[2]).resolve();launcher=sys.argv[3:]
work=Path(tempfile.mkdtemp(prefix='salmon-elf-proca-test-'));print(work,flush=True)
base=(root/'docs/results/si-elf-feedback/strong.inp').read_text().replace("tdcdft='lrc'","tdcdft='proca'").replace('tdcdft_damping=0.0','tdcdft_damping=0.01').replace('tdcdft_restoring=0','tdcdft_restoring=0.0004').replace('nt = 1600','nt = 80').replace('tw1 = 60','tw1 = 6').replace('checkpoint_interval=200','checkpoint_interval=37').replace("projection_option='gs'","projection_option='no'")
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
def run(name,inp=base,restart=None,fail=None):
 p=work/name;p.mkdir();(p/'restart').symlink_to(restart or gs);(p/'Si_rps.dat').symlink_to(root/'testsuites/pseudo/Si_rps.dat')
 if restart:inp=inp.replace("sysname = 'Si'","sysname = 'Si'\n yn_restart='y'")
 (p/'inputfile').write_text(inp)
 with (p/'outputfile').open('w') as f:result=subprocess.run(launcher+[str(exe)],input=inp,text=True,cwd=p,stdout=f,stderr=subprocess.STDOUT,env=env)
 text=(p/'outputfile').read_text()
 if fail:
  assert result.returncode!=0 and fail in text,(name,text[-1000:]);return p,None,None
 assert result.returncode==0 and 'end SALMON' in text,(name,text[-1000:])
 return p,np.loadtxt(p/'Si_rt.data'),np.loadtxt(p/'Si_rt_xc.data')
p,r,x=run('full')
np.testing.assert_allclose(x[:,7],.2*x[:,8]/x[:,12],atol=1e-14,rtol=0)
assert np.ptp(x[:,7])>1e-6
# Directly verify the declared Proca recurrence, including the variable alpha drive.
dt=.08;beta=.01;gamma=.0004
prediction=((2-gamma*dt**2)*x[1:-1,1:4]-(1-.5*beta*dt)*x[:-2,1:4]+dt**2*x[1:-1,7,None]*r[1:-1,13:16])/(1+.5*beta*dt)
np.testing.assert_allclose(prediction,x[2:,1:4],atol=1e-14,rtol=1e-10)
_,rr,xx=run('resume',restart=p/'checkpoint_rt_000037')
np.testing.assert_allclose(rr,r[37:],atol=1e-11,rtol=1e-8);np.testing.assert_allclose(xx,x[37:],atol=1e-11,rtol=1e-8)
run('changed_gamma',base.replace('restoring=0.0004','restoring=0.0005'),restart=p/'checkpoint_rt_000037',fail='incompatible checkpoint')
_,_,stopped=run('stopped',base.replace('tdcdft_elf_stride=10','tdcdft_elf_stride=10\n tdcdft_screen_stop=3.2'))
np.testing.assert_allclose(stopped[:40],x[:40],atol=1e-12,rtol=1e-9)
assert np.ptp(stopped[39:,7])==0
run('negative_alpha',base.replace('tdcdft_alpha=0.2','tdcdft_alpha=-0.2'),fail='normalized nonnegative alpha')
print('PASS ELF-Proca normalization, variable-alpha oscillator equation, odd restart, freeze, parameter guards',flush=True)
