"""Compare audited before/after SALMON executables on synthetic SBE exports.
Usage: PYTHONPATH=. python tools/test_sbe_entrypoint.py BEFORE AFTER NEW_DIRECTORY
No ground-state calculation or external scientific dataset is used.
"""
from pathlib import Path
import os
import subprocess
import sys
import numpy as np
from salmonpy.sbe import SBEInitialState
from salmonpy.fields import FieldSamples
from salmonpy.crystal import respond

before,after,dest=map(Path,sys.argv[1:]);dest.mkdir(exist_ok=False)
threads=int(os.environ.get('SALMONPY_TEST_THREADS','1'))
if threads<1: raise ValueError('thread count must be positive')
fixture=Path(__file__).resolve().parents[1]/'salmonpy/tests/sbe_reference.npz'
with np.load(fixture) as data:
    nb,nk=data['eigen'].shape
    for order in range(4):
        for vnl in (False,True):
            for active in (2,3):
                dirs=[]
                for tag,exe in [('before',before),('after',after)]:
                    d=dest/f'{order}-{int(vnl)}-{active}-{tag}';d.mkdir(); dirs.append(d)
                    (d/'test_k.data').write_text('# header\n'*5+''.join(f'{k+1} 0 0 0 {data["weights"][k]}\n' for k in range(nk))+'# coefficients (2*pi/a [a.u.]) in kx, ky, kz\n1 1 1\n')
                    text='# eigen\n# occupation\n# 1:io, 2:esp[a.u.], 3:occ\n'
                    for k in range(nk):
                        text+=f'k= {k+1}, spin= 1\n'
                        for b in range(nb): text+=f'{b+1} {data["eigen"][b,k]} 0.123\n'
                    (d/'test_eigen.data').write_text(text)
                    text='# header\n'*3
                    for part,name in enumerate(('momentum','nonlocal_velocity')):
                        if part: text+='# nonlocal\n'
                        for k in range(nk):
                            for b in range(nb):
                                for c in range(nb):
                                    text+=f'{k+1} {b+1} {c+1} '+' '.join(f'{z.real:.17e} {z.imag:.17e}' for z in data[name][b,c,:,k])+'\n'
                    (d/'test_tm.data').write_text(text)
                    inp=f'''&calculation theory='sbe' /
&control sysname='test' /
&units unit_system='au' /
&system yn_periodic='y', al=4,5,6, nelem=1, natom=1, nelec=2, nstate=3 /
&rgrid num_rgrid=8,8,8 /
&kgrid num_kgrid=2,1,1 /
&tgrid dt=.013, nt=40 /
&emfield ae_shape1='Acos2', E_amplitude1=.02, omega1=2, tw1=.52, epdir_re1=1,0,0 /
&sbe gauge_sbe='velocity_gauge', nstate_sbe(1)={active}, norder_correction={order}, yn_vnl_correction='{'y' if vnl else 'n'}' /
&analysis out_projection_step=10 /
&atomic_red_coor
'H' 0 0 0 1
/
'''
                    (d/'input.inp').write_text(inp)
                    result=subprocess.run([str(exe)],input=inp,text=True,cwd=d,env=dict(os.environ,OMP_NUM_THREADS=str(threads),OPENBLAS_NUM_THREADS='1'),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=30)
                    (d/'run.log').write_text(result.stdout)
                    if result.returncode or not (d/'test_sbe_rt.data').exists():
                        raise RuntimeError(f'{d}: {result.stdout[-2000:]}')
                for name in ('test_sbe_rt.data','test_sbe_rt_energy.data','test_sbe_nex.data'):
                    if threads==1:
                        assert (dirs[0]/name).read_bytes()==(dirs[1]/name).read_bytes(),name
                    else:
                        np.testing.assert_allclose(np.loadtxt(dirs[0]/name),np.loadtxt(dirs[1]/name),rtol=3e-13,atol=3e-15)
                initial=SBEInitialState.from_salmon(dirs[1],'test',bands=3,kpoints=2,electrons=2,lattice=data['lattice'],active_bands=active,correction_order=order,nonlocal_correction=vnl)
                rt=np.loadtxt(dirs[1]/'test_sbe_rt.data'); assert len(rt)==40
                drive=FieldSamples(np.arange(41)*.013,np.vstack([np.zeros(3),rt[:,1:4]]),'A')
                actual=respond(initial,drive)
                np.testing.assert_allclose(actual.current[1:],rt[:,13:16],rtol=3e-13,atol=3e-15)
print(f'16 full-application configurations, threads={threads}: '+
      ('output files byte-identical' if threads==1 else 'outputs agree within tolerance')+
      '; Python currents agree')
