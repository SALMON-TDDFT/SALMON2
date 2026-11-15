"""Local single-rank Si regression for relocated catalog-selected GS/1.

Run with an installed wheel using python -I; requires ASE/f90nml and an explicit
compatible DFT/TDDFT runtime. Uses at most two workers, one thread each. Blocks
Python socket calls, not OS networking. Retains GS files for inspection.
"""
import argparse
import hashlib
import time
import json
from pathlib import Path
import sys
import numpy as np
from salmonpy.examples.configured_si import run as prepare
from salmonpy.examples.catalog_tddft import run
from salmonpy.config import load_materials
from salmonpy.builtin_materials import builtin_material_registry
from salmonpy.preparation import prepare_inputs, verify_inputs
from salmonpy.initial_states import tddft_from_catalog
from salmonpy.runtime import load_runtime



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runtime')
    parser.add_argument('directory')
    args = parser.parse_args()
    root = Path(args.directory).resolve()
    root.mkdir(parents=True, exist_ok=False)
    runtime = args.runtime
    started = time.monotonic()
    def guard(event,args):
        if event.startswith('socket.'):
            raise RuntimeError('Python network socket prohibited in this regression')
    sys.addaudithook(guard)
    report=prepare(runtime,root/'original')
    source=root/'original'
    recipe=json.loads((source/'inputs/materials.json').read_text())
    manifest=source/'gs/ground-state.json'
    record=json.loads(manifest.read_text())
    digest=hashlib.sha256(manifest.read_bytes()).hexdigest()
    recipe['initial_states']={'Si-ground':{'material_id':'Si-starter','source':{
        'format':'salmonpy.gs/1','path':'saved','manifest_sha256':digest}}}
    files={name:source/'inputs'/name for name in ('si.xyz','Si_rps.dat','LICENSE.txt')}
    files['saved/ground-state.json']=manifest
    files.update({'saved/'+name:source/'gs'/name for name in record['files']})
    registry=builtin_material_registry()
    bundle=prepare_inputs(recipe,root/'transfer',registry=registry,base_directory=source/'inputs',files=files,
                          provenance={'saved/ground-state.json':{'source':'fresh local Si connection regression'}})
    anchor=hashlib.sha256((bundle/'manifest.json').read_bytes()).hexdigest()
    relocated=root/'relocated';bundle.rename(relocated)
    (source/'gs').rename(source/'gs-unavailable-at-original-path')
    verify_inputs(relocated,manifest_sha256=anchor)
    with np.load(source/'current.npz') as a:
        A=a['A'].copy(); reference=a['current'].copy(); sample_time=a['time'].copy()
    np.savez(root/'waveform.npz',A=A,dt=float(sample_time[1]-sample_time[0]))
    run(relocated/'materials.json','Si-ground',runtime,root/'waveform.npz',root/'selected')
    with np.load(root/'selected/current.npz') as a:
        np.testing.assert_array_equal(a['current'],reference)
    catalog=load_materials(relocated/'materials.json',registry=registry)
    engine=load_runtime(runtime,required=('dft','tddft','pz','single-rank'))
    with tddft_from_catalog(catalog,'Si-ground',base_directory=relocated,executable=engine,max_steps=len(A)-1) as model:
        with model.new_state(dt=float(sample_time[1]-sample_time[0])) as first, model.new_state(dt=float(sample_time[1]-sample_time[0])) as second:
            values=[first.current.copy()]
            for field in A[1:]:
                old=second.current.copy();first.step(field)
                np.testing.assert_array_equal(old,second.current)
                second.step(-field);values.append(first.current.copy())
            np.testing.assert_array_equal(values,reference)
            assert not np.array_equal(first.current,second.current)
            assert first.initial_state_id==second.initial_state_id=='Si-ground'
        metadata=model.record()
    report.update(relocated_gs_selection=True,selected_current_max_abs_error=0.,independent_states_exact=True,
                  identity=metadata,python_sockets_blocked=True,bundle_manifest_sha256=anchor)
    report['resources']['wall_seconds'] = time.monotonic()-started
    (root/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
