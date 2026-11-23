"""Synthetic native Γ records→GS/2→TDDFT parity. No SCF/actual GS collection.

python tools/test_gs_gamma_export_runtime.py CONVERTER_BUILD NEW_OUTPUT
Uses existing SHA-validated runtime/probe, builds nothing. Fixtures deliberately
include synthetic complete-result metadata only to exercise the API contract.
"""
import argparse,hashlib,json,os,shutil,subprocess
from pathlib import Path
import numpy as np
from salmonpy.gs_native import GammaRestartConverter
from salmonpy.dft import DFT,Pseudopotential,SCFSettings,GroundStateResult
from salmonpy.fdtd import SalmonExecutable


def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('build',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    record=json.loads(a.build.read_text());converter=GammaRestartConverter.from_build(a.build)
    engine=SalmonExecutable(record['target_executable'],record['target_sha256']);probe=Path(record['probe'])
    if digest(probe)!=record['probe_sha256']:raise ValueError('probe hash changed')
    root=a.output.resolve();root.mkdir(parents=True,exist_ok=False);env=dict(os.environ)
    package=Path(__import__('salmonpy').__file__).resolve().parent
    pp=Pseudopotential.from_file(package/'examples/si_intro_data/Si_rps.dat',atomic_number=14,local_channel=2)
    summaries=[]
    for kind in ('real','complex'):
        case=root/kind;native=case/'source/native';restart=native/'data_for_restart';restart.mkdir(parents=True)
        mode='write-gs-norm-real' if kind=='real' else 'write-gs-norm'
        with (case/'probe.log').open('w') as log:
            subprocess.run([str(probe),mode,str(restart)+os.sep],stdout=log,stderr=subprocess.STDOUT,env=env,check=True,timeout=60)
        model=DFT([14],np.eye(3)*8,[[.3,.4,.45]],[pp],4,2,(8,8,8),(1,1,1),'pz');settings=SCFSettings()
        (native/'gs.inp').write_text(model.input_text(settings));(native/'pseudo1_rps.dat').write_bytes(pp.data)
        (native/'gs_status.dat').write_text('SALMONPY_GS_V1\n1 1\n1e-10 1e-9 -1.0\nrho_dne\n')
        (native/'gs_eigen.data').write_text('#esp: single-particle energies (eigen energies)\n#occ: occupation numbers, io: orbital index\n# 1:io, 2:esp[a.u.], 3:occ\nk= 1,  spin= 1\n1 -.5 2\n2 .2 2\n')
        gs=dict(model._restart_signature(),schema='salmonpy.gs/1',status='complete',fractional_positions=model.fractional_positions.tolist(),
                executable_sha256=engine.sha256,energy_hartree=-1.,iterations=1,tolerance=settings.tolerance,max_iterations=settings.max_iterations,
                input_sha256=digest(native/'gs.inp'),files={x.relative_to(case/'source').as_posix():digest(x) for x in native.rglob('*') if x.is_file()})
        (case/'source/ground-state.json').write_text(json.dumps(gs,indent=2)+'\n')
        source=GroundStateResult.open(case/'source');archive=source.export_archive(case/'archive',converter=converter,orbital_block=1)
        converted=converter.write_bundle(archive,case/'reimported',executable=engine,pseudopotentials=[pp])
        check='read-gs-norm-cos' if kind=='real' else 'read-gs-norm'
        with (case/'probe.log').open('a') as log:
            subprocess.run([str(probe),check,str(converted/'data_for_restart')+os.sep],stdout=log,stderr=subprocess.STDOUT,env=env,check=True,timeout=60)
        values=np.array([[0,0,.0001*n] for n in range(1,7)])
        with source.tddft_initial_state(executable=engine,max_steps=12,timeout=60) as original,archive.tddft_initial_state(
             converter=converter,executable=engine,pseudopotentials=[pp],max_steps=12,timeout=60) as imported:
            with original.new_state(dt=.001) as one,imported.new_state(dt=.001) as two:
                first=[];second=[]
                for index,A in enumerate(values):
                    first.append(one.step(A));second.append(two.step(A))
                    if index==1:token=two.snapshot()
                left=np.array(first);right=np.array(second);assert left.tobytes()==right.tobytes()
                imported.close();two.restore(token)
                replay=np.array([two.step(A) for A in values[2:]])
                assert replay.tobytes()==right[2:].tobytes()
                for label,state in [('original',one),('imported',two)]:shutil.copy2(Path(state._work.name)/'run/worker.log',case/(label+'-worker.log'))
        density_electrons=float(archive.array('density').sum()*512/512)
        if not np.isclose(density_electrons,4.,atol=1e-12,rtol=0):raise ValueError('reconstructed density integral differs')
        summaries.append(dict(native_kind=kind,unchanged_native_reader_wave_parity=True,bitwise_native_vs_exported_currents=True,
            bitwise_replay=True,max_current_difference=float(np.max(abs(left-right))),density_integral=density_electrons,
            occupations=archive.array('occupations').tolist(),currents=right.tolist(),provenance=dict(archive.metadata['provenance'])))
    # MappingProxy/tuple nested provenance values are serialized by explicit conversion.
    def plain(x):
        from collections.abc import Mapping
        if isinstance(x,Mapping):return {k:plain(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [plain(v) for v in x]
        return x
    result=dict(scope='synthetic native records and complete-result metadata; real unchanged runtime TDDFT, NOT SCF/physical convergence',
                target_sha256=engine.sha256,exporter_sha256=converter.exporter_sha256,cases=summaries)
    (root/'result.json').write_text(json.dumps(plain(result),indent=2)+'\n');print(json.dumps(plain(result),indent=2))

if __name__=='__main__':main()
