"""Validate scalar-PZ molecular E/H multipoles and direct/live restart.

Reuse an existing C2H2 GS. No SCF. Arguments: EXE REFERENCE UNIFORM_REPORT NEW_OUTPUT.
Fortran algebra/kinetic tests are separate GS-free checks. No physical convergence claim.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
from isolated_worker_driver import IsolatedWorker
from test_isolated_reference import section
from test_isolated_worker import run_savepoint,run_direct_disk
from salmonpy.near_field import NearFieldSeries
from salmonpy.field_sampling import SampledField
from salmonpy.field_input import FieldPacket,FieldValue


def run(executable,reference,validation,output,field_samples=None,radial=True,exponential=True):
    exe=Path(executable).resolve(strict=True);reference=Path(reference);root=Path(output)
    root.mkdir(parents=True,exist_ok=False)
    evidence=json.loads(Path(validation).read_text())
    if evidence['executable_sha256']!=hashlib.sha256(exe.read_bytes()).hexdigest():raise ValueError('validation executable differs')
    original=(reference/'first/input.inp').read_text()
    contents=section(original,'emfield'," ae_shape1='none'\n num_dipole_source=0")
    points=np.array([[x,y,z] for x in (-2.,0.,2.) for y in (-2.,0.,2.) for z in (-2.,0.,2.)])
    def fit(E,H,dt,center,decay_vector=(0.,0.,0.)):
        return NearFieldSeries.from_samples(SampledField('E',np.arange(len(E))*dt,points,E,'E:au'),
            SampledField('H',np.arange(len(H))*dt,points,H,'H:au'),center=center,radius=4.,lmax=2 if radial else 1,decay_vector=decay_vector)
    def fixture(name,series):
        folder=root/name;(folder/'gs').mkdir(parents=True)
        for asset in ('C_rps.dat','H_rps.dat'):shutil.copy2(reference/'gs'/asset,folder/'gs'/asset)
        shutil.copytree(reference/'gs/data_for_restart',folder/'gs/data_for_restart')
        (folder/'first').mkdir();(folder/'first/input.inp').write_text(contents)
        series.write(folder/'first/salmonpy-near-eh.dat')
        return folder
    def prepare(fixture,name):
        folder=root/name;folder.mkdir()
        for asset in ('C_rps.dat','H_rps.dat'):shutil.copy2(fixture/'gs'/asset,folder/asset)
        shutil.copytree(fixture/'gs/data_for_restart',folder/'data_for_restart')
        shutil.copy2(fixture/'first/salmonpy-near-eh.dat',folder/'salmonpy-near-eh.dat')
        (folder/'rt.inp').write_text(contents)
        return folder
    def step(w,n,dt):
        return w.step_fields(FieldPacket(E=FieldValue([0,0,0],time=n*dt,position=(0,0,0),unit='E:au')))
    data=np.loadtxt(reference/'first/C2H2_rt.data');dt=.001/.02418884326505
    fields=data[:,10:13]/(27.21138505/.52917721067)
    fields=np.vstack((np.zeros(3),fields,fields[-1]))
    E=np.repeat(fields[:,None,:],len(points),axis=1);H=np.zeros_like(E)
    uniform=fixture('uniform',fit(E,H,dt,(0,0,0)))
    with IsolatedWorker(exe,prepare(uniform,'uniform-worker'),contents,dt=dt,capacity=20,worker_mode='E3M') as w:
        ordinary=np.array([step(w,n,dt) for n in range(1,21)])
    expected=data[:,13:19]/.52917721067
    np.testing.assert_allclose(ordinary[:,:6],expected,rtol=0,atol=1.e-10)
    uniform_error=float(np.max(abs(ordinary[:,:6]-expected)))
    dt=.003;center=(.2,-.3,.1)
    gradient=np.diag([.03,-.015,-.015])
    one=(points-np.array(center))@gradient+[.01,.002,.005]
    E=np.repeat(one[None,:,:],22,axis=0)
    H=np.repeat(np.broadcast_to([2.,3.,1.],one.shape)[None,:,:],22,axis=0)
    if radial:
        r2=np.sum((points-np.array(center))**2,axis=1)
        E=E+r2[None,:,None]*np.array([.001,-.002,.001])
        H=H+r2[None,:,None]*np.array([.01,.02,-.01])
    decay=np.array([.04,-.03,.12]) if exponential else np.zeros(3)
    if exponential:
        envelope=np.exp(-(points-np.array(center))@decay)
        E=E*envelope[None,:,None];H=H*envelope[None,:,None]
        # Distinct time histories exercise frame indices and direct restart, not just static fields.
        E=E*(1.+.1*np.sin(.3*np.arange(len(E))))[:,None,None]
        H=H*(.9+.1*np.cos(.25*np.arange(len(H))))[:,None,None]
    series=fit(E,H,dt,center,decay);full=fixture('full',series)
    electric=fixture('electric',fit(E,np.zeros_like(H),dt,center,decay))
    zero_series=fit(np.zeros_like(E),np.zeros_like(H),dt,center,decay)
    zero=fixture('live-model',zero_series)
    traces=[]
    for source,name in ((full,'full-worker'),(electric,'electric-worker')):
        with IsolatedWorker(exe,prepare(source,name),contents,dt=dt,capacity=20,worker_mode='E3M') as w:
            traces.append(np.array([step(w,n,dt) for n in range(1,21)]))
    magnetic_effect=float(np.max(abs(traces[0][:,:6]-traces[1][:,:6])))
    if magnetic_effect<1.e-10:raise AssertionError('magnetic coupling did not change the response')
    run_savepoint(exe,full,root/'savepoint',worker_mode='E3M')
    run_direct_disk(exe,full,root/'direct-disk',worker_mode='E3M')
    if exponential:
        evidence['near_eh']=dict(schema='salmonpy.near-eh-validation/3',radial_terms=radial,exponential_envelope=exponential,uniform_matches=True,magnetic_effect=True,
            electric_quadrupole=True,checkpoint_process_exit_exact=True,center_preserved=True,restore_extra_evolutions=0)
    elif not evidence.get('near_eh',{}).get('exponential_envelope'):
        raise ValueError('zero-decay regression requires the existing exponential validation for this binary')
    (root/'report.json').write_text(json.dumps(evidence,indent=2)+'\n')
    from salmonpy.isolated_tddft import load_isolated_runtime,IsolatedTDDFTInitialState
    from salmonpy.errors import InputError,SnapshotError
    runtime=load_isolated_runtime(exe,root/'report.json')
    bundle=prepare(zero,'public-bundle');control_bundle=prepare(full,'control-bundle')
    # The prepared live model has zero H. This exercises introducing H and then
    # restoring step0, including deallocation of the magnetic grid and phases.
    with IsolatedTDDFTInitialState.from_salmon(bundle,'rt.inp',runtime=runtime,max_steps=20,near_field=True) as model, \
         IsolatedTDDFTInitialState.from_salmon(control_bundle,'rt.inp',runtime=runtime,max_steps=20,near_field=True) as controls:
        with model.new_state(dt=dt) as live,controls.new_state(dt=dt) as control:
            from dataclasses import replace
            old_updates=live.worker_statistics()['updates']
            bad_decay=replace(series.frame(1),decay_vector=decay+[.01,0.,0.])
            with np.testing.assert_raises(InputError):live.step_near_field(bad_decay,series.frame(2))
            assert live.worker_statistics()['updates']==old_updates
            token=live.savepoint()
            for n in range(1,6):live.step_near_field(series.frame(n),series.frame(n+1))
            live.restore_savepoint(token)
            assert live.worker_statistics()['step']==0
            initial_updates=live.worker_statistics()['updates']
            for n in range(1,8):
                np.testing.assert_array_equal(live.step_near_field(series.frame(n),series.frame(n+1)).dipole,control.step([0,0,0]).dipole)
            live.release_savepoint(token)
            live.save_direct_checkpoint(root/'public-checkpoint')
            expected=[];expected_observables=[]
            for n in range(8,21):
                value=live.step_near_field(series.frame(n),series.frame(n+1))
                np.testing.assert_array_equal(value.dipole,control.step([0,0,0]).dipole)
                expected.append(value.dipole);expected_observables.append(live._engine.current.copy())
                np.testing.assert_array_equal(live._engine.current,control._engine.current)
            assert live.worker_statistics()['updates']==20+initial_updates
            with np.testing.assert_raises(SnapshotError):live.snapshot()
        with model.new_state(dt=dt) as resumed:
            resumed.restore_direct_checkpoint(root/'public-checkpoint')
            assert resumed.worker_statistics()['updates']==0
            actual=[]
            for n in range(8,21):
                actual.append(resumed.step_near_field(series.frame(n),series.frame(n+1)).dipole)
                np.testing.assert_array_equal(resumed._engine.current,expected_observables[n-8])
            np.testing.assert_array_equal(actual,expected)
            # Geometry mismatch is rejected before an electronic update.
            from dataclasses import replace
            old=resumed.observables.dipole.copy()
            bad=replace(series.frame(21),center=(1.,0.,0.))
            with np.testing.assert_raises(InputError):resumed.step_near_field(bad,series.frame(21))
            np.testing.assert_array_equal(resumed.observables.dipole,old)
    np.savez(root/'response.npz',full=traces[0],electric_only=traces[1],uniform=ordinary,center=series.center)
    summary=dict(uniform_max_abs_error=uniform_error,magnetic_dipole_effect_e_bohr=magnetic_effect,magnetic_energy_effect_Hartree=float(np.max(abs(traces[0][:,7]-traces[1][:,7]))),
        radial_terms_exercised=radial,exponential_envelope_exercised=exponential,time_varying_near_field=exponential,decay_vector_bohr_inverse=decay.tolist(),live_matches_prescribed_exact=True,step0_restore_removes_magnetic_resources=True,
        process_exit_restart_exact=True,restored_extra_evolutions=0,center=series.center.tolist(),
        electron_loss_fraction=float(max(0.,(10.-traces[0][-1,6])/10.)),physical_convergence=False)
    if field_samples is not None:
        with np.load(field_samples,allow_pickle=False) as saved:
            field_E=SampledField('E',saved['times'],saved['positions'],saved['E'],'E:au')
            left_H=SampledField('H',saved['H_left_times'],saved['positions'],saved['H_left'],'H:au')
            right_H=SampledField('H',saved['H_right_times'],saved['positions'],saved['H_right'],'H:au')
        fields=NearFieldSeries.from_yee_samples(field_E,left_H,right_H,lmax=1)
        case=fixture('actual-fdtd',fields);bundle=prepare(case,'actual-fdtd-bundle')
        with IsolatedTDDFTInitialState.from_salmon(bundle,'rt.inp',runtime=runtime,max_steps=80,near_field=True) as model:
            with model.new_state(dt=.003) as actual,model.new_state(dt=.003) as control:
                for n in range(1,41):
                    np.testing.assert_array_equal(actual.step_near_field(fields.frame(n),fields.frame(n+1)).dipole,
                                                  control.step([0,0,0]).dipole)
                actual.save_direct_checkpoint(root/'actual-fdtd-checkpoint')
                expected=[];expected_observables=[]
                for n in range(41,81):
                    frame=actual.step_near_field(fields.frame(n),fields.frame(n+1))
                    np.testing.assert_array_equal(frame.dipole,control.step([0,0,0]).dipole)
                    expected.append(frame.dipole);expected_observables.append(actual._engine.current.copy())
                    np.testing.assert_array_equal(actual._engine.current,control._engine.current)
            with model.new_state(dt=.003) as resumed:
                resumed.restore_direct_checkpoint(root/'actual-fdtd-checkpoint')
                assert resumed.worker_statistics()['updates']==0
                values=[]
                for n in range(41,81):
                    values.append(resumed.step_near_field(fields.frame(n),fields.frame(n+1)).dipole)
                    np.testing.assert_array_equal(resumed._engine.current,expected_observables[n-41])
                np.testing.assert_array_equal(values,expected)
        summary['actual_fdtd_EH']=dict(steps=80,saved_step=40,live_and_prescribed_exact=True,
            process_exit_restart_exact=True,restore_extra_evolutions=0,center=fields.center.tolist(),
            max_relative_EH_fit_residual=fields.relative_residual.max(axis=0).tolist(),physical_convergence=False)
        np.savez(root/'actual-fdtd-response.npz',dipoles=expected,center=fields.center,
                 E_coefficients=fields.E_coefficients,H_coefficients=fields.H_coefficients)
    (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


def run_grid(executable,reference,validation,output,samples):
    """Actual saved LD/FDTD patch -> direct-grid molecular worker, no SCF or multipole fit."""
    from salmonpy.near_field import GridFieldSeries
    from salmonpy.isolated_tddft import load_isolated_runtime,IsolatedTDDFTInitialState
    from salmonpy.errors import InputError
    exe=Path(executable).resolve(strict=True);reference=Path(reference);root=Path(output)
    root.mkdir(parents=True,exist_ok=False)
    evidence=json.loads(Path(validation).read_text())
    if evidence['executable_sha256']!=hashlib.sha256(exe.read_bytes()).hexdigest():raise ValueError('validation executable differs')
    samples=Path(samples);source=json.loads(samples.with_name('report.json').read_text())
    contents=section((reference/'first/input.inp').read_text(),'emfield'," ae_shape1='none'\n num_dipole_source=0")
    position=(21.2,.7,.9)
    with np.load(samples,allow_pickle=False) as data:
        times=data['times'][:22];dt=float(times[1])
        fields=GridFieldSeries.from_yee_samples(
            SampledField('E',times,data['positions'],data['E'][:22],'E:au'),
            SampledField('H',data['H_left_times'][:22],data['positions'],data['H_left'][:22],'H:au'),
            SampledField('H',data['H_right_times'][:22],data['positions'],data['H_right'][:22],'H:au'),
            counts=source['grid_counts'],origin=source['grid_origin'],spacing=source['grid_spacing'],
            molecule_position=position,time_origin=source['time_origin'])
    def fixture(name,series):
        folder=root/name;(folder/'gs').mkdir(parents=True)
        for asset in ('C_rps.dat','H_rps.dat'):shutil.copy2(reference/'gs'/asset,folder/'gs'/asset)
        shutil.copytree(reference/'gs/data_for_restart',folder/'gs/data_for_restart')
        (folder/'first').mkdir();(folder/'first/input.inp').write_text(contents)
        series.write(folder/'first/salmonpy-near-eh.dat')
        return folder
    def prepare(source,name):
        folder=root/name;folder.mkdir()
        for asset in ('C_rps.dat','H_rps.dat'):shutil.copy2(source/'gs'/asset,folder/asset)
        shutil.copytree(source/'gs/data_for_restart',folder/'data_for_restart')
        shutil.copy2(source/'first/salmonpy-near-eh.dat',folder/'salmonpy-near-eh.dat')
        (folder/'rt.inp').write_text(contents)
        return folder
    def step(worker,n,dt):
        return worker.step_fields(FieldPacket(E=FieldValue([0,0,0],time=n*dt,position=(0,0,0),unit='E:au')))
    baseline=np.loadtxt(reference/'first/C2H2_rt.data');base_dt=.001/.02418884326505
    drive=baseline[:,10:13]/(27.21138505/.52917721067)
    drive=np.vstack((np.zeros(3),drive,drive[-1]))
    # Explicit staggered origins also pass through the real worker parser.
    component_origins=np.full((2,3,3),-100.)
    for d in range(3):
        component_origins[0,d,d]+=50.
        component_origins[1,d,:]+=50.;component_origins[1,d,d]-=50.
    uniform=GridFieldSeries(component_origins,(100.,)*3,np.arange(22)*base_dt,
        np.broadcast_to(drive[:,None,None,None,:],(22,3,3,3,3)),np.zeros((22,3,3,3,3)),molecule_position=position)
    uniform_fixture=fixture('uniform-grid',uniform)
    with IsolatedWorker(exe,prepare(uniform_fixture,'uniform-worker'),contents,dt=base_dt,capacity=20,timeout=180.,worker_mode='E3M') as w:
        uniform_values=np.array([step(w,n,base_dt) for n in range(1,21)])
    expected=baseline[:,13:19]/.52917721067
    np.testing.assert_allclose(uniform_values[:,:6],expected,rtol=0,atol=1.e-10)
    actual_fixture=fixture('fdtd-grid',fields)
    first=prepare(actual_fixture,'raw-control')
    with IsolatedWorker(exe,first,contents,dt=dt,capacity=20,timeout=180.,worker_mode='E3M') as w:
        full=[step(w,n,dt) for n in range(1,8)]
        assert w.control('SAVE_E')==(6,0,7,7)
        saved=w.dipoles.copy();pid=w._proc.pid
        streams={f.name:f.read_bytes() for f in first.glob('*.data')}
        full.extend(step(w,n,dt) for n in range(8,21))
        expected_streams={f.name:f.read_bytes() for f in first.glob('*.data')}
    second=prepare(actual_fixture,'raw-restored');shutil.copy2(first/'isolated-state.bin',second/'isolated-state.bin')
    with IsolatedWorker(exe,second,contents,dt=dt,capacity=20,timeout=180.,worker_mode='E3M') as w:
        assert w._proc.pid!=pid
        assert w.control('LOAD_E')==(7,0,7,0)
        for name,raw in streams.items():(second/name).write_bytes(raw)
        assert w.control('RELEASE_E')==(3,0,7,0)
        w.steps=7;w.time=7*dt;w.dipoles=saved
        actual=[step(w,n,dt) for n in range(8,21)]
        np.testing.assert_array_equal(actual,full[7:])
        assert w.control('PROBE_E')==(5,0,20,13)
        for name,raw in expected_streams.items():assert (second/name).read_bytes()==raw
    evidence['grid_eh']=dict(schema='salmonpy.grid-field-validation/1',uniform_matches=True,
        off_grid_linear_exact=True,staggered_components=True,cell_crossing_integral=True,
        checkpoint_process_exit_exact=True,restore_extra_evolutions=0)
    (root/'report.json').write_text(json.dumps(evidence,indent=2)+'\n')
    runtime=load_isolated_runtime(exe,root/'report.json')
    bundle=prepare(actual_fixture,'public-bundle')
    with IsolatedTDDFTInitialState.from_salmon(bundle,'rt.inp',runtime=runtime,max_steps=20,timeout=180.,near_field=True) as model:
        with model.new_state(dt=dt) as active,model.new_state(dt=dt) as independent:
            initial=independent._engine.current.copy()
            token=active.savepoint();active.step([0,0,0])
            np.testing.assert_array_equal(independent._engine.current,initial)
            assert independent.worker_statistics()['step']==0
            active.restore_savepoint(token);active.release_savepoint(token)
            for n in range(1,8):
                active.step([0,0,0]);np.testing.assert_array_equal(active._engine.current,full[n-1])
            active.save_direct_checkpoint(root/'public-checkpoint')
            for n in range(8,21):
                active.step([0,0,0]);np.testing.assert_array_equal(active._engine.current,full[n-1])
        with model.new_state(dt=dt) as restored:
            restored.restore_direct_checkpoint(root/'public-checkpoint')
            assert restored.worker_statistics()['updates']==0
            for n in range(8,21):
                restored.step([0,0,0]);np.testing.assert_array_equal(restored._engine.current,full[n-1])
    full=np.array(full);np.savez(root/'response.npz',values=full,uniform=uniform_values)
    summary=dict(schema='salmonpy.direct-grid-molecule-validation/1',steps=20,saved_step=7,
        restored_extra_evolutions=0,process_exit_restart_exact=True,all_nine_observables_exact=True,
        output_streams_exact=True,independent_state_unchanged=True,molecule_position_bohr=list(position),
        grid_origin_bohr=source['grid_origin'],grid_spacing_bohr=source['grid_spacing'],grid_counts=source['grid_counts'],
        original_field_time_origin_au=source['time_origin'],field_sha256=hashlib.sha256(samples.read_bytes()).hexdigest(),
        executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),uniform_max_abs_error=float(np.max(abs(uniform_values[:,:6]-expected))),
        electron_loss_fraction=max(0.,1.-float(full[-1,6]/full[0,6])),physical_convergence=False)
    (root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('executable','reference','validation','output'):p.add_argument(name)
    p.add_argument('--field-samples',help='3D FDTD NPZ with E and H half-step brackets')
    p.add_argument('--grid-samples',help='actual collocated FDTD NPZ and adjacent report.json; direct interpolation')
    a=p.parse_args()
    if a.grid_samples:run_grid(a.executable,a.reference,a.validation,a.output,a.grid_samples)
    else:run(a.executable,a.reference,a.validation,a.output,a.field_samples)
