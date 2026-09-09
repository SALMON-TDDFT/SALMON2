#!/usr/bin/env python3
from pathlib import Path
import os, shutil, subprocess, tempfile

root=Path(__file__).resolve().parents[2]
identity=root/'src/rt/dg/rt_dg_hybrid_system_identity.f90'
sha=root/'src/common/dg_portable_sha256.f90'
gs=(root/'src/gs/main_dft.f90').read_text().lower()
rt=(root/'src/rt/main_tddft.f90').read_text().lower()
for text,label in ((gs,'GS publisher'),(rt,'RT initializer')):
    assert 'fingerprint_rt_dg_hybrid_system' in text,f'RED: {label} does not bind current system identity'
assert 'canonical_pp_fingerprint(pp)' in rt,'RED: RT does not compare current canonical PP identity'
identity_text=identity.read_text().lower()
for token in ('system%no','system%nk','system%vec_k','system%wtk','system%rocc','dg_sha256_mix_int64'):
    assert token in identity_text,f'RED: system identity omits {token}'
for source in (identity,root/'src/gs/dc/dg_canonical_pp_fingerprint.f90',root/'src/rt/dg/rt_dg_hybrid_checkpoint_v4.f90'):
    assert 'ieor(ishftc(hash,7)' not in source.read_text().lower(),f'RED: rotate/XOR digest remains in {source}'
with tempfile.TemporaryDirectory(prefix='hybrid-system-identity-') as name:
    build=Path(name);(build/'config.h').write_text('')
    exe=build/'identity'
    subprocess.run([shutil.which('mpifort'),'-cpp','-DUSE_MPI','-I',str(build),'-J',str(build),
      str(root/'src/common/structures.f90'),str(sha),str(identity),
      str(root/'tests/dg/test_rt_dg_hybrid_system_identity.f90'),'-o',str(exe)],check=True)
    env={**os.environ,'OMP_NUM_THREADS':'1','OMPI_MCA_rmaps_base_oversubscribe':'1'}
    for nrank in (1,2,4,8):
      run=subprocess.run([shutil.which('mpiexec'),'-n',str(nrank),str(exe)],capture_output=True,text=True,env=env,timeout=30)
      assert run.returncode==0,(nrank,run.stdout,run.stderr)
      assert f'PASS authoritative Hybrid system identity mutations ranks={nrank}' in run.stdout
print('PASS Hybrid GS/RT authoritative physical-system identity on 1/2/4/8 ranks')
