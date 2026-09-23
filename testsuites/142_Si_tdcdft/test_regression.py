"""Integration checks. Usage: python test_regression.py EXE GS_DIRECTORY BASELINE_RT [MPI ...]."""
from pathlib import Path
import os
import subprocess
import sys
import tempfile
import numpy as np

exe,gs,baseline=map(lambda x:Path(x).resolve(),sys.argv[1:4])
launcher=sys.argv[4:]
root=Path(__file__).resolve().parents[2]
base=(root/'testsuites/112_bulk_Si_rt_response_dp/inputfile').read_text()
env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1')
work=Path(tempfile.mkdtemp(prefix='salmon-xc-regression-'))
print('Results:',work,flush=True)

def run(name,setting='',steps=180,dt=0.16,impulse=0.01,control='',restart=None,extra='',pulse=False,fail=None):
    path=work/name;path.mkdir()
    (path/'restart').symlink_to(restart or gs,target_is_directory=True)
    (path/'Si_rps.dat').symlink_to(root/'testsuites/pseudo/Si_rps.dat')
    inp=base.replace("xc = 'PZ'", "xc = 'PZ'\n"+setting)
    inp=inp.replace('nt = 180',f'nt = {steps}').replace('dt = 0.16',f'dt = {dt}')
    inp=inp.replace("ae_shape1 = 'impulse'",f"ae_shape1 = 'impulse'\n e_impulse = {impulse}")
    inp=inp.replace("sysname = 'Si'","sysname = 'Si'\n"+control)+extra
    if pulse:
        inp=inp.replace("theory='tddft_response'", "theory='tddft_pulse'")
        inp=inp.replace("ae_shape1 = 'impulse'", "ae_shape1='Acos2'\n tw1=20\n omega1=0.2\n I_wcm2_1=1e12")
    (path/'inputfile').write_text(inp)
    with (path/'outputfile').open('w') as out:
        result=subprocess.run(launcher+[str(exe)],input=inp,text=True,cwd=path,stdout=out,stderr=subprocess.STDOUT,env=env)
    if fail:
        assert result.returncode != 0 and fail in (path/'outputfile').read_text(), name
        return path,None
    if result.returncode: raise RuntimeError(f'{name} failed: {path}/outputfile')
    return path,np.loadtxt(path/'Si_rt.data')

_,off=run('disabled')
ref=np.loadtxt(baseline)
np.testing.assert_allclose(off,ref,atol=1e-12,rtol=1e-10)
_,zero=run('zero',"tdcdft='lrc'\n tdcdft_alpha=0")
np.testing.assert_array_equal(zero,off)
lrc="tdcdft='lrc'\n tdcdft_alpha=0.2"
proca="tdcdft='proca'\n tdcdft_alpha=0.2\n tdcdft_restoring=0.01"
p,full=run('proca',proca,control='checkpoint_interval=90')
_,lr=run('lrc',lrc)
assert np.max(abs(full[:,15]-off[:,15]))>1e-7
assert np.max(abs(lr[:,15]-full[:,15]))>1e-7
_,resume=run('resume',proca,restart=p/'checkpoint_rt_000090',control="yn_restart='y'")
np.testing.assert_allclose(resume[:,13:16],full[90:,13:16],atol=1e-11,rtol=1e-8)
np.testing.assert_allclose(np.loadtxt(work/'resume/Si_response.data'),np.loadtxt(p/'Si_response.data'),atol=1e-9,rtol=1e-8)
_,small=run('small_probe',proca,impulse=0.001)
_,half=run('half_probe',proca,impulse=0.0005)
rel=np.linalg.norm(half[:,15]*2-small[:,15])/np.linalg.norm(small[:,15])
assert rel<0.005,('weak probe linearity',rel)
_,fine=run('half_dt',proca,steps=360,dt=0.08,impulse=0.001)
rel_dt=np.linalg.norm(fine[1::2,15]-small[:,15])/np.linalg.norm(fine[1::2,15])
assert rel_dt<0.08,('time-step sensitivity',rel_dt)
print('PASS disabled baseline, exact alpha=0, finite coupling, restart and spectrum, linearity',rel,'dt',rel_dt,flush=True)

# Laser alpha -> 0 limit must not switch the classical nonlocal-current phase.
_,pz=run('pulse_zero',"tdcdft='lrc'\n tdcdft_alpha=0",pulse=True)
_,pe=run('pulse_epsilon',"tdcdft='lrc'\n tdcdft_alpha=1e-10",pulse=True)
np.testing.assert_allclose(pe[:,13:16],pz[:,13:16],atol=1e-11,rtol=1e-7)
assert abs(np.loadtxt(p/'Si_rt_xc.data')[0,3])>1e-10
run('negative',"tdcdft='lrc'\n tdcdft_alpha=-1",fail='parameters must be nonnegative')
run('nan',"tdcdft='lrc'\n tdcdft_alpha=NaN",fail='parameters must be finite')
run('unknown',"tdcdft='bad'",fail='tdcdft must be')
run('bad_propagator',lrc,extra="\n&propagation propagator='aetrs' /\n",fail='requires middlepoint')
run('mismatched_restart',proca.replace('0.2','0.3'),restart=p/'checkpoint_rt_000090',
    control="yn_restart='y'",fail='TDCDFT')
print('PASS initial half step, pulse zero-coupling continuity, rejected invalid input/restart',flush=True)
