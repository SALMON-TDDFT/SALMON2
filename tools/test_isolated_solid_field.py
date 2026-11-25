"""One-way saved 3D FDTD E -> Fortran Ylm fit -> C2H2 TDDFT validation.

Reuse an existing GS. No SCF or physical-convergence claim. The uniform case is
compared with the SALMON batch reference. Direct save/restart uses a new process.
"""
import argparse
import json
import shutil
from pathlib import Path
import numpy as np
from isolated_worker_driver import IsolatedWorker
from test_isolated_reference import section
from test_isolated_worker import run_savepoint,run_direct_disk
from salmonpy.field_sampling import SampledField
from salmonpy.near_field import ElectricPotentialSeries


def run(executable,reference,validation,fdtd,output):
    exe=Path(executable).resolve(strict=True);reference=Path(reference);root=Path(output)
    root.mkdir(parents=True,exist_ok=False)
    evidence=json.loads(Path(validation).read_text())
    import hashlib
    if evidence['executable_sha256']!=hashlib.sha256(exe.read_bytes()).hexdigest():raise ValueError('validation executable differs')
    original=(reference/'first/input.inp').read_text()
    contents=section(original,'emfield'," ae_shape1='none'\n num_dipole_source=0")
    def fixture(name,series):
        destination=root/name;(destination/'gs').mkdir(parents=True)
        for asset in ('C_rps.dat','H_rps.dat'):shutil.copy2(reference/'gs'/asset,destination/'gs'/asset)
        shutil.copytree(reference/'gs/data_for_restart',destination/'gs/data_for_restart')
        first=destination/'first';first.mkdir();(first/'input.inp').write_text(contents)
        series.write(first/'salmonpy-ylm.dat')
        return destination
    def prepare(fixture,name):
        folder=root/name;folder.mkdir()
        for asset in ('C_rps.dat','H_rps.dat'):shutil.copy2(fixture/'gs'/asset,folder/asset)
        shutil.copytree(fixture/'gs/data_for_restart',folder/'data_for_restart')
        shutil.copy2(fixture/'first/salmonpy-ylm.dat',folder/'salmonpy-ylm.dat')
        (folder/'rt.inp').write_text(contents)
        return folder
    def history(fixture,name,steps,dt):
        folder=prepare(fixture,name)
        with IsolatedWorker(exe,folder,contents,dt=dt,capacity=steps,worker_mode='E3Y') as worker:
            rows=[worker.dipoles.copy()]
            from salmonpy.field_input import FieldPacket,FieldValue
            for n in range(1,steps+1):
                rows.append(worker.step_fields(FieldPacket(E=FieldValue([0,0,0],time=n*dt,position=(0,0,0),unit='E:au'))))
        return np.array(rows),folder
    # Uniform harmonics reproduce the existing full SALMON time history, not a
    # different hold/endpoint approximation of the field.
    batch=np.loadtxt(reference/'first/C2H2_rt.data');dt=.001/.02418884326505
    E=batch[:,10:13]/(27.21138505/.52917721067)
    E=np.vstack((np.zeros(3),E,E[-1]))
    points=np.array([[x,y,z] for x in (-2.,0.,2.) for y in (-2.,0.,2.) for z in (-2.,0.,2.)])
    series=ElectricPotentialSeries.from_samples(SampledField('E',np.arange(22)*dt,points,np.repeat(E[:,None,:],len(points),axis=1),'E:au'),center=(0,0,0),radius=8.,lmax=1)
    uniform=fixture('uniform-reference',series)
    actual,_=history(uniform,'uniform-worker',20,dt)
    expected=batch[:,13:19]/.52917721067
    np.testing.assert_allclose(actual[1:,:6],expected,rtol=0,atol=1.e-10)
    uniform_error=float(np.max(abs(actual[1:,:6]-expected)))
    with np.load(Path(fdtd)/'samples.npz',allow_pickle=False) as data:
        samples=SampledField('E',data['times'],data['positions'],data['values'],'E:au')
    fitted=ElectricPotentialSeries.from_samples(samples,center=(1,1,1),radius=8.,lmax=2)
    real=fixture('fdtd-reference',fitted)
    trace,bundle=history(real,'fdtd-worker',80,.003)
    if not np.isfinite(trace).all():raise AssertionError('nonfinite TDDFT response')
    dipole_max=float(np.max(abs(trace[:,3:6]-trace[0,3:6])))
    if dipole_max<1.e-12:raise AssertionError('no molecular response to the fitted FDTD field')
    # The same saved field sampled only to l=1 is an explicit negative control.
    first_order=ElectricPotentialSeries.from_samples(samples,center=(1,1,1),radius=8.,lmax=1)
    control=fixture('fdtd-dipole-reference',first_order)
    dipole_trace,_=history(control,'fdtd-dipole-worker',80,.003)
    higher_order_effect=float(np.max(abs(trace-dipole_trace)))
    if higher_order_effect<1.e-12:raise AssertionError('quadrupole terms had no detected effect')
    run_savepoint(exe,real,root/'savepoint',worker_mode='E3Y')
    run_direct_disk(exe,real,root/'direct-disk',worker_mode='E3Y')
    evidence['solid_field']=dict(schema='salmonpy.solid-field-validation/1',uniform_max_abs_error=uniform_error,
        savepoint_same_process=True,direct_disk_process_exit_exact=True,quadrupole_response_nonzero=True)
    (root/'report.json').write_text(json.dumps(evidence,indent=2)+'\n')
    np.savez(root/'response.npz',values=trace,dipole_only=dipole_trace,coefficients=fitted.coefficients,
             relative_residual=fitted.relative_residual,dt=.003)
    from salmonpy.isolated_tddft import load_isolated_runtime,IsolatedTDDFTInitialState
    runtime=load_isolated_runtime(exe,root/'report.json')
    with IsolatedTDDFTInitialState.from_salmon(bundle,'rt.inp',runtime=runtime,max_steps=80,electric_potential=True) as model:
        with model.new_state(dt=.003) as state:
            for n in range(1,41):np.testing.assert_array_equal(state.step([0,0,0]).dipole,trace[n,3:6])
            state.save_direct_checkpoint(root/'public-checkpoint')
            expected=np.array([state.step([0,0,0]).dipole for _ in range(40)])
        with model.new_state(dt=.003) as state:
            state.restore_direct_checkpoint(root/'public-checkpoint')
            assert state.worker_statistics()['updates']==0
            resumed=np.array([state.step([0,0,0]).dipole for _ in range(40)])
            np.testing.assert_array_equal(resumed,expected)
    report=dict(uniform_max_abs_error=uniform_error,fdtd_dipole_max=dipole_max,
                higher_order_effect=higher_order_effect,field_relative_residual_max=float(fitted.relative_residual.max()),
                field_relative_residual_final=float(fitted.relative_residual[-1]),public_process_exit_exact=True,
                restored_extra_evolutions=0,steps=80,saved_step=40,physical_convergence=False,
                input_field=Path(fdtd).as_posix(),
                field_time_origin=json.loads((Path(fdtd)/'report.json').read_text()).get('time_origin',0.),
                molecular_initial_state='field-free GS at the start of the explicitly selected field window',mode='one-way saved field; scalar electric potential, no magnetic multipoles or radiation feedback')
    (root/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('executable','reference','validation','fdtd','output'):p.add_argument(name)
    a=p.parse_args();run(a.executable,a.reference,a.validation,a.fdtd,a.output)
