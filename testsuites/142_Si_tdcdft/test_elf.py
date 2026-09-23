"""ELF feedback integration: EXE GS_DIRECTORY [MPI launcher prefix]."""
from pathlib import Path
import sys,os,subprocess,tempfile,struct,shutil
import numpy as np
root=Path(__file__).resolve().parents[2];exe=Path(sys.argv[1]).resolve();gs=Path(sys.argv[2]).resolve();launcher=sys.argv[3:]
work=Path(tempfile.mkdtemp(prefix='salmon-elf-test-'));print(work,flush=True)
base=(root/'docs/results/si-time-wannier/strong_dense.inp').read_text().replace('tdcdft_alpha=1.0','tdcdft_alpha=0.2').replace("tdcdft_screening='polarization'","tdcdft_screening='elf'\n tdcdft_elf_stride=10").replace('nt = 1600','nt = 80').replace('tw1 = 60','tw1 = 6').replace('checkpoint_interval=10','checkpoint_interval=37').replace("projection_option='gs'","projection_option='no'")
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
def run(name,inp=base,restart=None,fail=None):
 p=work/name;p.mkdir();(p/'restart').symlink_to(restart or gs);(p/'Si_rps.dat').symlink_to(root/'testsuites/pseudo/Si_rps.dat')
 if restart:inp=inp.replace("sysname = 'Si'","sysname = 'Si'\n yn_restart='y'")
 (p/'inputfile').write_text(inp)
 with (p/'outputfile').open('w') as f:r=subprocess.run(launcher+[str(exe)],input=inp,text=True,cwd=p,stdout=f,stderr=subprocess.STDOUT,env=env)
 text=(p/'outputfile').read_text()
 if fail:
  assert r.returncode!=0 and fail in text,(name,text[-1500:]);return p,None
 assert r.returncode==0 and 'end SALMON' in text,(name,text[-1500:])
 return p,np.atleast_2d(np.loadtxt(p/'Si_rt_xc.data'))
p,x=run('full')
assert x.shape==(80,13) and np.all(np.isfinite(x))
np.testing.assert_allclose(x[:,7],.2*x[:,8]/x[:,12],atol=1e-14,rtol=0)
np.testing.assert_allclose(x[:9,7],.2,atol=1e-14,rtol=0)
for step in range(2,81):
 if step%10:assert x[step-1,8]==x[step-2,8]
assert np.ptp(x[:,7])>1e-6
pr,resume=run('resume',restart=p/'checkpoint_rt_000037')
np.testing.assert_allclose(resume,x[37:],atol=1e-11,rtol=1e-8)
_,even=run('even_resume',restart=p/'checkpoint_rt_000074')
_,chain=run('chain_resume',restart=pr/'checkpoint_rt_000074')
np.testing.assert_allclose(even,x[74:],atol=1e-11,rtol=1e-8)
np.testing.assert_allclose(chain,x[74:],atol=1e-11,rtol=1e-8)
run('bad_stride',base.replace('stride=10','stride=0'),fail='stride must be positive')
run('changed_stride',base.replace('stride=10','stride=5'),restart=p/'checkpoint_rt_000037',fail='incompatible checkpoint')
_,stopped=run('stop',base.replace('tdcdft_elf_stride=10','tdcdft_elf_stride=10\n tdcdft_screen_stop=3.2'))
np.testing.assert_allclose(stopped[:40],x[:40],atol=1e-12,rtol=1e-9)
assert np.ptp(stopped[39:,8])==0 and np.ptp(stopped[39:,9:12],axis=0).max()>1e-8
# A consistent synthetic checkpoint verifies alpha above alpha0 is accepted on restart.
cp=work/'above';shutil.copytree(p/'checkpoint_rt_000037',cp)
raw=(cp/'tdcdft.bin').read_bytes();records=[];at=0
while at<len(raw):
 n=struct.unpack('<i',raw[at:at+4])[0];records.append(bytearray(raw[at+4:at+4+n]));at+=n+8
reference=struct.unpack('<d',records[4][4:12])[0]
records[2][-16:]=struct.pack('<dd',1.1*reference,.22)
(cp/'tdcdft.bin').write_bytes(b''.join(struct.pack('<i',len(r))+r+struct.pack('<i',len(r)) for r in records))
_,above=run('above_resume',base.replace('nt = 80','nt = 38'),restart=cp)
assert abs(above[0,7]-.22)<1e-14
# Reject corrupted normalization before dividing by it.
records[4][4:12]=struct.pack('<d',0.)
(cp/'tdcdft.bin').write_bytes(b''.join(struct.pack('<i',len(r))+r+struct.pack('<i',len(r)) for r in records))
run('bad_reference',restart=cp,fail='incompatible checkpoint')
print('PASS ELF normalization, cadence, non-update restart, alpha>alpha0 restart, stop, input guards',flush=True)

# Odd restart is also required for the pre-existing fixed-alpha path.
fixed=base.replace("tdcdft_screening='elf'","tdcdft_screening='none'")
pf,xf=run('fixed',fixed)
_,rf=run('fixed_odd_resume',fixed,restart=pf/'checkpoint_rt_000037')
np.testing.assert_allclose(rf,xf[37:],atol=1e-11,rtol=1e-8)
print('PASS fixed-alpha odd-step restart',flush=True)

# Enable only with a four-rank MPI launcher. Exercise halo exchange and both
# reductions; nonlinear ELF must be formed after summing orbital/k moments.
if os.environ.get('SALMON_ELF_TEST_RSPACE')=='1':
 split=base.replace('nt = 80','nt = 40')+'\n&parallel nproc_k=2, nproc_ob=1, nproc_rgrid=2,1,1 /\n'
 _,xr=run('rspace',split)
 np.testing.assert_allclose(xr,x[:40],atol=1e-11,rtol=1e-8)
 print('PASS combined real-space/k-space MPI decomposition',flush=True)
