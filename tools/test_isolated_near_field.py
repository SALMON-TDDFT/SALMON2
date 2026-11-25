"""Explicit C2H2 dipole-source validation; reuse a saved GS, never run SCF.

EXE REFERENCE UNIFORM_VALIDATION NEW_OUTPUT. All new data stay in NEW_OUTPUT.
Run test_isolated_worker.py first for this executable's uniform-field evidence.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import numpy as np
from isolated_worker_driver import IsolatedWorker
from test_isolated_reference import section
from test_isolated_worker import run_direct_disk, run_savepoint


def run(executable, reference, uniform, output):
    executable=Path(executable).resolve(strict=True)
    reference=Path(reference).resolve(strict=True)
    uniform=Path(uniform).resolve(strict=True)
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False)
    binary_hash=hashlib.sha256(executable.read_bytes()).hexdigest()
    base=json.loads((uniform/'report.json').read_text())
    if base['executable_sha256']!=binary_hash:raise ValueError('uniform validation executable differs')
    original=(reference/'first/input.inp').read_text()
    # Input length/time/energy units remain Angstrom/fs/eV.
    contents=section(original,'emfield',""" ae_shape1='Ecos2'
 ae_shape2='none'
 I_wcm2_1=0d0
 tw1=0.020d0
 tw2=0d0
 omega1=3.1d0
 phi_cep1=0.75d0
 num_dipole_source=2
 rad_dipole_source=0.5d0
 cood_dipole_source(:,1)=3d0,0d0,0d0
 cood_dipole_source(:,2)=-3d0,0d0,0d0
 vec_dipole_source(:,1)=0d0,0d0,0.02d0
 vec_dipole_source(:,2)=0d0,0.01d0,0d0""")
    fixture=root/'reference';(fixture/'gs').mkdir(parents=True)
    for asset in ('C_rps.dat','H_rps.dat'):shutil.copy2(reference/'gs'/asset,fixture/'gs'/asset)
    shutil.copytree(reference/'gs/data_for_restart',fixture/'gs/data_for_restart')
    def prepare(folder):
        folder.mkdir(parents=True,exist_ok=False)
        for asset in ('C_rps.dat','H_rps.dat'):shutil.copy2(fixture/'gs'/asset,folder/asset)
        shutil.copytree(fixture/'gs/data_for_restart',folder/'data_for_restart')
        return folder
    batch=prepare(fixture/'first')
    batch_input=contents.replace("directory_read_data = './restart/'", "directory_read_data = './data_for_restart/'")
    # The source example leaves directory_read_data at its default ./restart/.
    (batch/'restart').symlink_to('data_for_restart',target_is_directory=True)
    (batch/'input.inp').write_text(contents)
    env={k:v for k,v in os.environ.items() if not k.startswith('SALMONPY_')}
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    result=subprocess.run([str(executable)],cwd=batch,input=batch_input,text=True,capture_output=True,env=env,timeout=120)
    (batch/'batch.log').write_text(result.stdout+result.stderr)
    if result.returncode or not (batch/'C2H2_rt.data').is_file():
        raise RuntimeError('batch failed: '+str(batch/'batch.log'))
    dt=.001/.02418884326505
    from salmonpy.field_input import FieldPacket,FieldValue
    def step(worker,n,E=(0.,0.,0.)):
        return worker.step_fields(FieldPacket(E=FieldValue(E,time=n*dt,position=(0,0,0),unit='E:au')))
    with IsolatedWorker(executable,prepare(root/'near'),contents,dt=dt,capacity=20,worker_mode='E3N') as a:
        actual=np.array([step(a,n) for n in range(1,21)])
    expected=np.loadtxt(batch/'C2H2_rt.data')[:,13:19]/.52917721067
    np.testing.assert_allclose(actual[:,:6],expected,rtol=0,atol=1.e-10)
    reference_error=float(np.max(abs(actual[:,:6]-expected)))
    if np.max(abs(actual[:,:6]))<1.e-9:raise AssertionError('near-field test has no response')
    zero=section(contents,'emfield'," ae_shape1='none'\n num_dipole_source=0")
    traces=[]
    for mode in ('E3','E3N'):
        with IsolatedWorker(executable,prepare(root/('zero-'+mode)),zero,dt=dt,capacity=20,worker_mode=mode) as w:
            traces.append(np.array([step(w,n,(.001,-.002,.003)) for n in range(1,21)]))
    np.testing.assert_array_equal(*traces)
    run_savepoint(executable,fixture,root/'savepoint',worker_mode='E3N')
    run_direct_disk(executable,fixture,root/'disk',worker_mode='E3N')
    base['native_near_field']=dict(schema='salmonpy.isolated-near-field-validation/1',steps=20,
        reference_max_abs_error=reference_error,zero_source_max_abs_error=0.,
        interleaved_independence=True,savepoint_same_process=True,direct_disk_process_exit_exact=True)
    (root/'report.json').write_text(json.dumps(base,indent=2)+'\n')
    np.savez(root/'near-field.npz',values=actual,expected=expected,dt=dt)
    print(json.dumps(base['native_near_field'],indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable');p.add_argument('reference');p.add_argument('uniform_validation');p.add_argument('output')
    a=p.parse_args();run(a.executable,a.reference,a.uniform_validation,a.output)
