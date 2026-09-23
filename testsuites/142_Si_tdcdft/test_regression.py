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

def run(name,setting='',steps=180,dt=0.16,impulse=0.01,control='',restart=None,extra='',pulse=False,fail=None,unit_fs=False):
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
    if unit_fs:
        length=.52917721067; time=.02418884326505; energy=27.21138505
        inp=inp.replace("unit_system = 'a.u.'", "unit_system = 'A_eV_fs'")
        inp=inp.replace('al = 10.26d0, 10.26d0, 10.26d0','al = '+', '.join([str(10.26*length)]*3))
        inp=inp.replace('dl = 0.855d0, 0.855d0, 0.855d0','dl = '+', '.join([str(.855*length)]*3))
        inp=inp.replace(f'dt = {dt}',f'dt = {dt*time:.17g}')
        inp=inp.replace(f'e_impulse = {impulse}',f'e_impulse = {impulse*energy/length*time:.17g}')
        inp=inp.replace('de = 0.001',f'de = {.001*energy:.17g}')
        begin,coor=inp.split('&atomic_coor',1)
        rows=[]
        for line in coor.splitlines():
            if line.strip().startswith("'Si'"):
                cols=line.split(); cols[1:4]=[str(float(v)*length) for v in cols[1:4]]
                line=' '.join(cols)
            rows.append(line)
        inp=begin+'&atomic_coor'+'\n'.join(rows)
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

# Paper coefficients normalize to the same equation, including old checkpoint metadata.
a2=-20*np.pi
ksp=f"tdcdft='proca'\n tdcdft_a2={a2:.17g}\n tdcdft_a0={a2*.01:.17g}"
pd,direct=run('direct_proca',ksp,control='checkpoint_interval=90')
np.testing.assert_allclose(direct,full,atol=1e-12,rtol=1e-9)
np.testing.assert_allclose(np.loadtxt(pd/'Si_rt_xc.data'),np.loadtxt(p/'Si_rt_xc.data'),atol=1e-12,rtol=1e-9)
_,direct_resume=run('direct_resume',ksp,restart=p/'checkpoint_rt_000090',control="yn_restart='y'")
np.testing.assert_allclose(direct_resume[:,13:16],full[90:,13:16],atol=1e-11,rtol=1e-8)
_,massless=run('direct_massless',f"tdcdft='proca'\n tdcdft_a2={a2:.17g}\n tdcdft_a0=0")
np.testing.assert_allclose(massless,lr,atol=1e-12,rtol=1e-9)
_,positive=run('positive_a2',"tdcdft='proca'\n tdcdft_a2=10\n tdcdft_a0=.25")
assert np.isfinite(positive).all()
for name,setting,reason in [
 ('mixed',ksp+"\n tdcdft_alpha=.2",'do not mix'),
 ('zero_a2',"tdcdft='proca'\n tdcdft_a0=-.2",'a2 must be nonzero'),
 ('bad_ratio',"tdcdft='proca'\n tdcdft_a2=-10\n tdcdft_a0=.25",'require a0/a2'),
 ('nan_a2',"tdcdft='proca'\n tdcdft_a2=NaN\n tdcdft_a0=-.2",'must be finite'),
 ('other_mode',"tdcdft='lrc'\n tdcdft_a2=-10",'a2/a0 require')]:
 run(name,setting,fail=reason)
run('changed_a2',ksp.replace(f'{a2:.17g}','-60'),restart=pd/'checkpoint_rt_000090',
    control="yn_restart='y'",fail='incompatible checkpoint')
print('PASS direct Proca equivalence, legacy restart, massless limit, both signs, invalid coefficients',flush=True)

# a0 is inverse-time squared in the selected unit system; a2 is unchanged.
paper="tdcdft='proca'\n tdcdft_a2=-85\n tdcdft_a0=-.2"
_,paper_au=run('paper_au',paper)
paper_fs=paper.replace('tdcdft_a0=-.2',f'tdcdft_a0={-.2/.02418884326505**2:.17g}')
_,paper_fs_data=run('paper_fs',paper_fs,unit_fs=True)
np.testing.assert_allclose(paper_fs_data[:,13:16]*.02418884326505*.52917721067**2,
    paper_au[:,13:16],atol=1e-11,rtol=1e-8)
print('PASS atomic-unit / A_eV_fs Proca equivalence',flush=True)

# Instantaneous screening: no temporal averaging, zero-strength reference, active restart.
screen="\n tdcdft_screening='instant'\n tdcdft_screen_omega=0.2\n tdcdft_screen_reference=0"
_,fixed_pulse=run('screen_fixed',proca,pulse=True)
p0,screen_zero=run('screen_zero',proca+screen+'\n tdcdft_screen_strength=0',pulse=True)
np.testing.assert_array_equal(screen_zero,fixed_pulse)
ps,screened=run('screen_active',proca+screen,pulse=True,control='checkpoint_interval=90')
xs=np.loadtxt(ps/'Si_rt_xc.data')
assert xs.shape[1]==12 and np.isfinite(xs).all()
assert np.min(xs[:,7])<0.2 and np.all((xs[:,7]>=0)&(xs[:,7]<=0.2))
assert np.max(abs(screened[:,15]-fixed_pulse[:,15]))>1e-10
_,screen_resume=run('screen_resume',proca+screen,pulse=True,restart=ps/'checkpoint_rt_000090',
                   control="yn_restart='y'")
np.testing.assert_allclose(screen_resume[:,13:16],screened[90:,13:16],atol=1e-11,rtol=1e-8)
np.testing.assert_allclose(np.loadtxt(work/'screen_resume/Si_rt_xc.data'),xs[90:],atol=1e-11,rtol=1e-8)
run('screen_changed',proca+screen+'\n tdcdft_screen_strength=2',pulse=True,
    restart=ps/'checkpoint_rt_000090',control="yn_restart='y'",fail='TDCDFT')
run('screen_removed',proca,pulse=True,restart=ps/'checkpoint_rt_000090',control="yn_restart='y'",fail='TDCDFT')
run('screen_bad_omega',proca+screen.replace('omega=0.2','omega=0'),pulse=True,fail='TDCDFT')
run('screen_impulse',proca+screen,fail='TDCDFT')
print('PASS instantaneous screening, zero strength, bounded active feedback, restart and guards',flush=True)

# Build an actual version-1 file from the unchanged first two Fortran records.
import shutil,struct
legacy=work/'version1_checkpoint'
shutil.copytree(p/'checkpoint_rt_000090',legacy)
raw=(legacy/'tdcdft.bin').read_bytes()
endian='<' if struct.unpack('<i',raw[:4])[0]==56 else '>'
first=struct.unpack(endian+'i',raw[:4])[0]
assert first==56
second_offset=first+8
second=struct.unpack(endian+'i',raw[second_offset:second_offset+4])[0]
raw=raw[:4]+struct.pack(endian+'i',1)+raw[8:second_offset+second+8]
(legacy/'tdcdft.bin').write_bytes(raw)
_,legacy_resume=run('version1_resume',proca,restart=legacy,control="yn_restart='y'")
np.testing.assert_allclose(legacy_resume[:,13:16],full[90:,13:16],atol=1e-11,rtol=1e-8)
run('version1_screening',proca+screen,pulse=True,restart=legacy,control="yn_restart='y'",fail='TDCDFT')
print('PASS version-1 fixed restart and rejection of missing screening state',flush=True)

# Polarization closure: same estimator, a'=-alpha P, no damping/restoring.
polar=screen.replace("'instant'","'polarization'")
_,lrc_pulse=run('polar_fixed_reference',lrc,pulse=True)
_,polar_zero=run('polar_zero',lrc+polar+'\n tdcdft_screen_strength=0',pulse=True)
np.testing.assert_allclose(polar_zero,lrc_pulse,atol=1e-12,rtol=1e-9)
pp,polar_current=run('polar_active',lrc+polar,pulse=True,control='checkpoint_interval=90')
xp=np.loadtxt(pp/'Si_rt_xc.data')
assert np.isfinite(xp).all() and xp.shape[1]==12
post=xp[:,0]>20.5
assert np.ptp(xp[post,7])==0
np.testing.assert_allclose(xp[post,4:7],xp[post,7,None]*xp[post,9:12],atol=1e-12,rtol=1e-9)
assert np.max(abs(polar_current[:,15]-lrc_pulse[:,15]))>1e-10
_,polar_resume=run('polar_resume',lrc+polar,pulse=True,restart=pp/'checkpoint_rt_000090',
                   control="yn_restart='y'")
np.testing.assert_allclose(polar_resume[:,13:16],polar_current[90:,13:16],atol=1e-11,rtol=1e-8)
np.testing.assert_allclose(np.loadtxt(work/'polar_resume/Si_rt_xc.data'),xp[90:],atol=1e-11,rtol=1e-8)
run('polar_changed_closure',lrc+screen,pulse=True,restart=pp/'checkpoint_rt_000090',
    control="yn_restart='y'",fail='TDCDFT')
run('polar_restoring',proca+polar,pulse=True,fail='polarization screening requires zero damping and restoring')
print('PASS polarization closure, constant-alpha limit, no post-pulse field offset, restart and guards',flush=True)
