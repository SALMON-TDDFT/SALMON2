"""Real GS integration: EXE GS_BUNDLE GS_INPUT NEW_OUTPUT_DIRECTORY."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from salmonpy.bands import BandInitialState,BandPath,BandStructure
from salmonpy.fdtd import SalmonExecutable

def hashes(root):
    return {p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob('*') if p.is_file()}
exe,source=map(lambda x:Path(x).resolve(),sys.argv[1:3]);inp=sys.argv[3];out=Path(sys.argv[4]).resolve()
out.mkdir(exist_ok=False,parents=True)
original=hashes(source);engine=SalmonExecutable(exe,hashlib.sha256(exe.read_bytes()).hexdigest())
path=BandPath([[0,0,0],[.5,0,0],[.5,.5,0],[.5,.5,.5],[0,0,0]],[4]*4,['Γ','X','M','R','Γ'])
with BandInitialState.from_salmon(source,inp,executable=engine) as initial:
    first=initial.run(path,out/'path-a',bands=24,threads=2)
    second=initial.run(BandPath([[0,0,0],[.5,0,0]],[7],['Γ','X']),out/'path-b',bands=24,threads=2)
    try:initial.run(path,out/'unconverged',bands=24,max_iterations=1,tolerance=1e-15,threads=2)
    except RuntimeError as exc:assert 'converg' in str(exc),str(exc)
    else:raise AssertionError('unconverged run accepted')
    assert json.loads((out/'unconverged/band-run.json').read_text())['status']=='failed'
    assert hashes(source)==original,'source GS mutated'
    assert dict(initial.manifest)==hashes(initial._snapshot),'snapshot mutated'
# Compare exact same native input invoked directly, from an independent original copy.
direct=out/'direct';shutil.copytree(source,direct)
text=(out/'path-a/native/band-input.inp').read_text()
with (direct/'stdout.log').open('w') as log:
    subprocess.run([str(exe)],cwd=direct,input=text,text=True,stdout=log,stderr=subprocess.STDOUT,
                   env={**os.environ,'OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'1'},check=True,timeout=300)
native=BandStructure.from_salmon(direct/'band.dat')
np.testing.assert_array_equal(first.bands.energies,native.energies[:17,:24])
np.testing.assert_array_equal(first.bands.k_cartesian,native.k_cartesian[:17])
import matplotlib
matplotlib.use('Agg')
reference=float(first.bands.energies[:,:16].max()*27.211386245988)
ax=first.plot(reference_ev=reference,title='Si: folded 8-atom cell bands; PZ, coarse grid')
ax.figure.savefig(out/'bands.png',dpi=160)
summary={'direct_native_max_difference_Ha':float(np.max(abs(first.bands.energies-native.energies[:17,:24]))),
         'first_points':len(first.bands.energies),'second_points':len(second.bands.energies),
         'nonconvergence_rejected':True,'original_gs_unchanged':True,'snapshot_unchanged':True,
         'executable_sha256':engine.sha256,'reference_ev':reference}
(out/'validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
