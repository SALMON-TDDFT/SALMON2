"""Real Si GS integration of exact disconnected ASE band path.

Usage: python tools/test_ase_bands.py EXE GS_BUNDLE GS_INPUT NEW_DIRECTORY
Requires optional ASE + f90nml and plotting dependencies.
"""
from pathlib import Path
import sys,os,shutil,hashlib,json,subprocess
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from salmonpy.ase import atoms_from_initial,run_bands
from salmonpy.bands import BandInitialState,BandStructure
from salmonpy.fdtd import SalmonExecutable
from ase.io.jsonio import read_json

def hashes(root):
    return {p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob('*') if p.is_file()}
exe,source=map(lambda v:Path(v).resolve(),sys.argv[1:3]);inp=sys.argv[3];out=Path(sys.argv[4]).resolve()
out.mkdir(parents=True,exist_ok=False)
original=hashes(source);engine=SalmonExecutable(exe,hashlib.sha256(exe.read_bytes()).hexdigest())
with BandInitialState.from_salmon(source,inp,executable=engine) as gs:
    atoms=atoms_from_initial(gs)
    path=atoms.cell.bandpath('GX,MR',npoints=17)
    before=path.kpts.copy();before_atoms=atoms.positions.copy()
    result=run_bands(gs,atoms,path,out/'ase',reference_ev=0.,bands=24,threads=2)
    np.testing.assert_array_equal(path.kpts,before)
    np.testing.assert_array_equal(atoms.positions,before_atoms)
    np.testing.assert_array_equal(result.band_structure.path.kpts,before)
    assert result.band_structure.path.path=='GX,MR'
    assert hashes(source)==original
    assert dict(gs.manifest)==hashes(gs._snapshot)
# Independent native execution from original GS and exact emitted path/input.
direct=out/'direct';shutil.copytree(source,direct)
shutil.copy2(out/'ase/native/salmonpy-kpoints.data',direct)
text=(out/'ase/native/band-input.inp').read_text()
with (direct/'stdout.log').open('w') as log:
    subprocess.run([str(exe)],input=text,text=True,cwd=direct,stdout=log,stderr=subprocess.STDOUT,
                   env={**os.environ,'OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'1'},timeout=300,check=True)
native=BandStructure.from_salmon(direct/'band.dat')
nk=len(before)
np.testing.assert_array_equal(result.native.bands.energies,native.energies[:nk,:24])
np.testing.assert_array_equal(result.native.bands.k_cartesian,native.k_cartesian[:nk])
loaded=read_json(out/'ase/ase-bands.json')
np.testing.assert_array_equal(loaded.energies,result.band_structure.energies)
assert loaded.path.path=='GX,MR'
import matplotlib
matplotlib.use('Agg')
result.plot(filename=str(out/'bands-ase.png'),show=False,emin=-7,emax=14)
summary={'ase_path':path.path,'physical_points':nk,'source_gs_unchanged':True,
         'caller_atoms_and_path_unchanged':True,'ase_json_roundtrip':True,
         'direct_native_max_difference_Ha':0.,'binary_sha256':engine.sha256}
(out/'validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
