"""Explicit compiled C-plugin conformance and handle-lifetime regression.

Usage: python tools/test_c_material_plugin.py BUILD_JSON NEW_OUTPUT
No native build or GS collection. Input binary hash is rechecked by CMaterial.
"""
import ctypes as c
import json
from pathlib import Path
import struct
import sys
import numpy as np
from salmonpy.c_material import CMaterial,CABIError
from salmonpy.errors import InputError,RuntimeFailure
from salmonpy.materials import bind_material,check_material
from salmonpy.crystal import respond
from salmonpy.examples.custom_material import ExampleMaterial,waveform

record=json.loads(Path(sys.argv[1]).read_text());output=Path(sys.argv[2]);output.mkdir(parents=True,exist_ok=False)
if record.get('status')!='built-audited-and-host-tested':raise InputError('audited tested example build required')
with CMaterial(record['library'],record['sha256'],struct.pack('=dd',2.,.8)) as model:
    report=check_material(model,waveform())
    native=respond(model,waveform()).current
    reference=respond(ExampleMaterial(),waveform()).current
    np.testing.assert_allclose(native,reference,rtol=1e-12,atol=1e-14)
    binding=bind_material(model,key='c-plugin',identity=model.identity)
    adapter=binding.backend.create_state(binding.provider.prepare(None),1,.02,None)
    try:
        adapted=np.vstack([adapter.current,*[adapter.step(v) for v in waveform().values[1:]]])
        np.testing.assert_array_equal(adapted,native)
    finally:adapter.close()
    one=model.new_state(dt=.02);two=model.new_state(dt=.03)
    one.step([.1,0,0]);np.testing.assert_array_equal(two.current,np.zeros(3))
    two.step([0,.2,0]);np.testing.assert_array_equal(one.current,[.2,0,0])
    # Model close prevents construction but existing states remain alive.
    model.close()
    with np.testing.assert_raises(RuntimeFailure):model.new_state(dt=.02)
    one.step([.2,0,0]);two.step([0,.3,0]);one.close();two.close()
    assert model._handle.value==0

with CMaterial(record['library'],record['sha256'],struct.pack('=dd',2.,.8)) as model:
    with model.new_state(dt=.02) as state:
        before=state.current.copy();clock=state.time
        with np.testing.assert_raises(InputError):state.step([float('nan'),0,0])
        np.testing.assert_array_equal(state.current,before);assert state.time==clock
        # Finite input, overflow in the plugin: status4, no native/cache mutation.
        try:state.step([1e308,0,0])
        except CABIError as exc:assert exc.status==4
        else:raise AssertionError('native overflow accepted')
        np.testing.assert_array_equal(state.current,before);assert state.time==clock
        np.testing.assert_array_equal(state.step([.1,0,0]),[.2,0,0])

np.savez(output/'currents.npz',native=native,reference=reference,adapted=adapted)
summary={'library_sha256':record['sha256'],'maximum_native_reference_difference':float(np.max(abs(native-reference))),
         'maximum_repeat_difference':report.maximum_repeat_difference,'adapter_bitwise_equal':True,
         'independent_states':True,'deferred_model_destruction':True,'atomic_nan_and_overflow_rejection':True,
         'scope':'synthetic ABI conformance; not native SALMON TDDFT reentrancy or physical validation'}
(output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
