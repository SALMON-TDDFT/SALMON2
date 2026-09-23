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
    inp=inp.replace("sysname = 'Si'","sysname = 'Si'\n"+control)
    if pulse:
        inp=inp.replace("theory='tddft_response'", "theory='tddft_pulse'")
        inp=inp.replace("ae_shape1 = 'impulse'", "ae_shape1='Acos2'\n tw1=20\n omega1=0.2\n I_wcm2_1=1e12")
    if extra: inp=inp.replace("ae_shape1='Acos2'", "ae_shape1='Acos2'\n"+extra)
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

setting="tdcdft='lrc'\n tdcdft_alpha=1\n tdcdft_screening='polarization'\n tdcdft_screen_omega=.2\n tdcdft_screen_stop=20"
probe="\n ae_shape2='impulse', epdir_re2=0,0,1, T1_T2=14\n"
p,j=run('stopped',setting,pulse=True,extra=probe,control='checkpoint_interval=150')
x=np.loadtxt(p/'Si_rt_xc.data');post=x[:,0]>20
assert np.ptp(x[post,7])==0 and np.ptp(x[post,8])==0
assert np.ptp(x[post,11])>1e-8
_,jr=run('resume',setting,pulse=True,extra=probe,restart=p/'checkpoint_rt_000150',control="yn_restart='y'")
np.testing.assert_allclose(jr,j[150:],atol=1e-11,rtol=1e-8)
run('changed_stop',setting.replace('stop=20','stop=21'),pulse=True,extra=probe,restart=p/'checkpoint_rt_000150',control="yn_restart='y'",fail='TDCDFT')
pd,default=run('default',setting.replace('\n tdcdft_screen_stop=20',''),pulse=True,extra=probe,control='checkpoint_interval=150')
np.testing.assert_allclose(j[x[:,0]<=20],default[x[:,0]<=20],atol=1e-12,rtol=1e-10)
assert np.max(abs(default[:,15]-j[:,15]))>1e-8
print('PASS frozen screening under delayed probe, restart and parameter guard',flush=True)

# Version 2 stored the same first three records, with no stop-time record.
import struct,shutil
legacy=work/'v2_checkpoint';shutil.copytree(pd/'checkpoint_rt_000150',legacy)
raw=(legacy/'tdcdft.bin').read_bytes();endian='<' if struct.unpack('<i',raw[:4])[0]==56 else '>'
offset=0
for _ in range(3):
 size=struct.unpack(endian+'i',raw[offset:offset+4])[0];offset+=size+8
raw=raw[:4]+struct.pack(endian+'i',2)+raw[8:offset]
(legacy/'tdcdft.bin').write_bytes(raw)
_,v2=run('v2_resume',setting.replace('\n tdcdft_screen_stop=20',''),pulse=True,extra=probe,restart=legacy,control="yn_restart='y'")
np.testing.assert_allclose(v2,default[150:],atol=1e-11,rtol=1e-8)
run('v2_with_stop',setting,pulse=True,extra=probe,restart=legacy,control="yn_restart='y'",fail='TDCDFT')
print('PASS v2 compatibility and rejection of changed stop',flush=True)
