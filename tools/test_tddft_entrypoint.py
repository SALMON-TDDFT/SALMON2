"""Native integration using synthetic normalized orbitals, not a converged GS.
PYTHONPATH=. python tools/test_tddft_entrypoint.py BEFORE AFTER NEW_DIRECTORY
Both executable source/object audits must pass before invoking this test.
"""
from pathlib import Path
import hashlib
import os
import shutil
import struct
import subprocess
import sys
import numpy as np
from salmonpy.tddft import TDDFTInitialState
from salmonpy.fdtd import SalmonExecutable
from salmonpy.fields import FieldSamples
from salmonpy.crystal import respond

before,after,dest=map(Path,sys.argv[1:]);dest.mkdir(exist_ok=False)
root=Path(__file__).resolve().parents[1]; bundle=dest/'bundle';bundle.mkdir()
restart=bundle/'data_for_restart';restart.mkdir()
def record(stream,data):
    stream.write(struct.pack('=i',len(data)));stream.write(data);stream.write(struct.pack('=i',len(data)))
with (restart/'info.bin').open('wb') as stream:
    for value in (1,2,0,1,0): record(stream,struct.pack('=i',value))
with (restart/'occupation.bin').open('wb') as stream: record(stream,np.array([2.,2.],dtype='=f8').tobytes())
x=np.arange(8)[:,None,None]
psi=np.empty((8,8,8,2),dtype='=c16',order='F')
psi[:,:,:,0]=1/np.sqrt(8.**3)
psi[:,:,:,1]=np.exp(2j*np.pi*x/8)/np.sqrt(8.**3)
(restart/'wfn.bin').write_bytes(psi.tobytes(order='F'))
runtime_source=os.environ.get('SALMONPY_RUNTIME_SOURCE')
if not runtime_source: raise RuntimeError('Set SALMONPY_RUNTIME_SOURCE to the prepared external runtime source')
shutil.copy2(Path(runtime_source)/'samples/exercise_06_bulkSi_rt/Si_rps.dat',bundle/'Si_rps.dat')
inp='''&calculation
 theory='tddft_pulse'
/
&control
 sysname='fixture', directory_read_data='./data_for_restart/'
/
&units
 unit_system='au'
/
&system
 yn_periodic='y', al=8,8,8, nelem=1, natom=1, nelec=4, nstate=2
/
&pseudo
 file_pseudo(1)='Si_rps.dat', izatom(1)=14, lloc_ps(1)=2
/
&functional
 xc='PZ'
/
&rgrid
 num_rgrid=8,8,8
/
&kgrid
 num_kgrid=1,1,1
/
&tgrid
 dt=.001, nt=8
/
&emfield
 ae_shape1='Acos2', E_amplitude1=.01, omega1=2, tw1=.008, epdir_re1=1,0,0
/
&analysis
 out_rt_energy_step=1
/
&atomic_red_coor
'Si' .5 .5 .5 1
/
'''
xc=os.environ.get('SALMONPY_TEST_XC','PZ')
pc=os.environ.get('SALMONPY_TEST_PC','n')
inp=inp.replace("xc='PZ'",f"xc='{xc}'")+f"\n&propagation\n yn_predictor_corrector='{pc}'\n/\n"
if 'SALMONPY_TEST_CVAL' in os.environ:
    cval=float(os.environ['SALMONPY_TEST_CVAL'])
    if not np.isfinite(cval): raise ValueError('finite cval required')
    inp=inp.replace(f"xc='{xc}'",f"xc='{xc}', cval={cval}")
(bundle/'rt.inp').write_text(inp)
outputs=[]
for tag,exe in [('before',before),('after',after)]:
    d=dest/tag;shutil.copytree(bundle,d)
    result=subprocess.run([str(exe)],input=inp,text=True,cwd=d,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                          env={**os.environ,'OMP_NUM_THREADS':'1'},timeout=60)
    (d/'run.log').write_text(result.stdout)
    if result.returncode or 'end SALMON' not in result.stdout: raise RuntimeError(result.stdout[-6000:])
    outputs.append(d)
for name in ('fixture_rt.data','fixture_rt_energy.data'):
    assert (outputs[0]/name).read_bytes()==(outputs[1]/name).read_bytes(),name
rt=np.loadtxt(outputs[1]/'fixture_rt.data')
print('ordinary before/after output byte-identical',flush=True)
exe=SalmonExecutable(after,hashlib.sha256(after.read_bytes()).hexdigest())
with TDDFTInitialState.from_salmon(bundle,'rt.inp',executable=exe,max_steps=8,timeout=60) as initial:
    drive=FieldSamples(np.arange(9)*.001,np.vstack([np.zeros(3),rt[:,1:4]]),'A')
    result=respond(initial,drive)
    np.testing.assert_allclose(result.current[1:],rt[:,13:16],rtol=1e-12,atol=1e-14)
    with initial.new_state(dt=.001) as a, initial.new_state(dt=.001) as b:
        assert a._proc.pid!=b._proc.pid
        for i,A in enumerate(drive.values[1:],1):
            np.testing.assert_array_equal(a.step(A),result.current[i])
            self_current=b.current.copy()
            with np.testing.assert_raises(ValueError): b.step([np.nan,0,0])
            np.testing.assert_array_equal(b.current,self_current)
            np.testing.assert_array_equal(b.step(A),result.current[i])
        with np.testing.assert_raises(ValueError): a.step([0,0,0])
    # Terminate one worker between steps; the failed step must restore accepted history.
    with initial.new_state(dt=.001) as state:
        state.step(drive.values[1]);saved=state.current.copy();state._proc.kill();state._proc.wait()
        with np.testing.assert_raises(RuntimeError): state.step(drive.values[2])
        assert state.time==.001
        np.testing.assert_array_equal(state.current,saved)
        np.testing.assert_array_equal(state.step(drive.values[2]),result.current[2])
# Unsupported modes/path escapes must fail before physical initialization.
for bad in (inp.replace("sysname='fixture'","sysname='../escape'"),
            inp.replace("sysname='fixture'","sysname='fixture', write_rt_wfn_k='y'"),
            inp.replace('num_kgrid=1,1,1',"num_kgrid=1,1,1, file_kw='/outside.dat'"),
            inp.replace(f"xc='{xc}'", "xc='libxc_pw'")):
    (bundle/'bad.inp').write_text(bad)
    with TDDFTInitialState.from_salmon(bundle,'bad.inp',executable=exe,max_steps=8,timeout=60) as initial:
        with np.testing.assert_raises(RuntimeError): initial.new_state(dt=.001)
print('native TDDFT A→J matches baseline; batch, interleaved states, rejection and recovery pass')
