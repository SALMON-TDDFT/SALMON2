"""End-to-end Si: ASE -> GS -> reopened GS -> bands and independent TDDFT.

Usage: python tools/test_ase_ground_state.py EXE SI_RPS_DAT NEW_DIRECTORY
"""
from pathlib import Path
import sys,os,hashlib,json,shutil,subprocess
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ase.build import bulk
from salmonpy.ase import dft_from_atoms,BOHR_ANGSTROM,atoms_from_initial
from salmonpy.dft import Pseudopotential,SCFSettings,GroundStateResult,_read_status
from salmonpy.bands import BandPath
from salmonpy.fdtd import SalmonExecutable
exe=Path(sys.argv[1]).resolve();pp=Path(sys.argv[2]).resolve();out=Path(sys.argv[3]).resolve()
out.mkdir(parents=True,exist_ok=False)
engine=SalmonExecutable(exe,hashlib.sha256(exe.read_bytes()).hexdigest())
atoms=bulk('Si','diamond',a=10.26*BOHR_ANGSTROM,cubic=True)
model=dft_from_atoms(atoms,pseudopotentials=[Pseudopotential.from_file(pp,atomic_number=14,local_channel=2)],
                     electrons=32,states=32,grid=(12,12,12),kmesh=(2,2,2),xc='pz')
result=model.run(engine,out/'gs',threads=2)
reopened=GroundStateResult.open(result.directory)
assert reopened.energy_hartree==result.energy_hartree
# Direct fresh run uses exactly the same input/potential, no restart copied.
direct=out/'direct';direct.mkdir()
text=(out/'gs/native/gs.inp').read_text()
for p in (out/'gs/native').glob('*_rps.dat'):shutil.copy2(p,direct)
with (direct/'stdout.log').open('w') as log:
    subprocess.run([str(exe)],input=text,text=True,cwd=direct,stdout=log,stderr=subprocess.STDOUT,
                   env={**os.environ,'OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'1'},timeout=600,check=True)
status=_read_status(direct/'gs_status.dat')
assert status['energy_hartree']==result.energy_hartree
assert status['iterations']==result.iterations
# Reuse without manually constructing an input bundle.
with reopened.band_initial_state(executable=engine) as initial:
    np.testing.assert_allclose(atoms_from_initial(initial).cell.array,atoms.cell.array,rtol=0,atol=1e-12)
    bands=initial.run(BandPath([[0,0,0],[.5,0,0]],[3],['Γ','X']),out/'bands',bands=24,threads=2)
with reopened.tddft_initial_state(executable=engine,max_steps=4,threads=2) as initial:
    with initial.new_state(dt=.02) as first,initial.new_state(dt=.02) as second:
        np.testing.assert_array_equal(first.current,second.current)
        observed=[]
        for n in range(1,5):
            A=np.array([1e-5*n,0,0]);observed.append(first.step(A));second.step(-A)
    with initial.new_state(dt=.02) as replay:
        for n,J in enumerate(observed,1):np.testing.assert_array_equal(replay.step([1e-5*n,0,0]),J)
GroundStateResult.open(result.directory) # Every output hash unchanged after consumers.
try:model.run(engine,out/'failed',scf=SCFSettings(max_iterations=1,tolerance=1e-15),threads=2)
except RuntimeError as exc:assert 'converg' in str(exc)
else:raise AssertionError('unconverged GS accepted')
assert json.loads((out/'failed/ground-state.json').read_text())['status']=='failed'
# Relocation and corruption detection on a copy, never on the source result.
shutil.copytree(out/'gs',out/'relocated')
GroundStateResult.open(out/'relocated')
with (out/'relocated/native/data_for_restart/occupation.bin').open('ab') as stream:stream.write(b'tamper')
try:GroundStateResult.open(out/'relocated')
except ValueError:pass
else:raise AssertionError('tampered GS accepted')
summary={'energy_hartree':result.energy_hartree,'iterations':result.iterations,'direct_energy_difference_hartree':0.,
         'bands_points':len(bands.bands.energies),'independent_tddft_replay_exact':True,'gs_unchanged_after_consumers':True,
         'nonconvergence_rejected':True,'relocation_passed':True,'tampering_rejected':True,'binary_sha256':engine.sha256}
(out/'validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
