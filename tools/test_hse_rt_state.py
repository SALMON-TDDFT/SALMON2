"""Explicit native lifecycle probe from the pinned fast-PR HSE implementation.

python tools/test_hse_rt_state.py SOURCE REFERENCE_SOURCE NEW_OUTPUT --blas LIB
No GS or solver run; actual ACE/history algorithms, extracted unchanged legacy
Taylor lifecycle, synthetic occupied factors. Source/object audits precede links.
"""
import argparse,hashlib,json,shutil,subprocess,sys
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','reference','output']:p.add_argument(name,type=Path)
    p.add_argument('--fc',default='gfortran');p.add_argument('--blas',type=Path,required=True)
    a=p.parse_args();a.source=a.source.resolve();a.reference=a.reference.resolve();a.output=a.output.resolve();a.blas=a.blas.resolve()
    if 'native_kernel' in (a.source/'src/xc/exx_rt_state.f90').read_text():
        raise ValueError('this standalone probe is parent-profile-only; use test_hse_explicit_state.py BUILD NEW_OUTPUT --native-resources for the integrated profile')
    if not a.blas.is_file():raise ValueError('explicit existing BLAS/LAPACK library required')
    root=Path(__file__).resolve().parents[1];checker=Path(__file__).with_name('check_forbidden_eigensolvers.py')
    for tree in [a.source,a.reference]:subprocess.run([sys.executable,str(checker),'--source-root',str(tree/'src')],check=True)
    a.output.mkdir(parents=True,exist_ok=False);src=a.output/'src';src.mkdir()
    for name in ['exx_ace','exx_factor_history']:
        shutil.copy2(a.source/f'src/xc/{name}.f90',src/f'{name}.f90')
    # Missing state module is intentionally observed by the first red probe.
    state=a.source/'src/xc/exx_rt_state.f90'
    if state.exists():shutil.copy2(state,src/'exx_rt_state.f90')
    text=(a.reference/'src/xc/exx_native.f90').read_text();start=text.index('  subroutine exx_taylor_stage(');end=text.index('  end subroutine',start)+len('  end subroutine')
    (src/'legacy.f90').write_text('''module legacy_lifecycle
use exx_ace,only:s_exx_ace,exx_ace_clear
implicit none
type(s_exx_ace)::ace,initial_ace,midpoint_ace
complex(8),allocatable::full_source(:,:,:),initial_source(:,:,:)
logical::taylor_active=.false.,taylor_midpoint=.false.
character(32)::propagator
contains
'''+text[start:end]+'\nend module\n')
    shutil.copy2(root/'tests/native/test_hse_rt_state.f90',src/'probe.f90')
    commands=[]
    with (a.output/'build.log').open('w') as log:
        def run(cmd):
            commands.append([str(x) for x in cmd]);subprocess.run(cmd,cwd=a.output,stdout=log,stderr=subprocess.STDOUT,check=True)
        run([sys.executable,checker,'--source-root',src])
        names=['exx_ace','exx_factor_history']+(['exx_rt_state'] if state.exists() else [])+['legacy','probe']
        for name in names:run([a.fc,'-O2','-fopenmp','-ffree-line-length-none','-c',src/(name+'.f90'),'-o',a.output/(name+'.o')])
        run([sys.executable,checker,'--object-tree',a.output])
        run([a.fc,'-fopenmp',*[a.output/(name+'.o') for name in names],a.blas,'-o',a.output/'probe'])
        run([a.output/'probe'])
    (a.output/'result.json').write_text(json.dumps(dict(status='passed',scope='synthetic fast-PR ACE/factor-history lifecycle, not full native multi-material or speed benchmark',commands=commands,
        source_hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in src.iterdir()},blas_sha256=hashlib.sha256(a.blas.read_bytes()).hexdigest()),indent=2)+'\n')
    print('HSE RT state native lifecycle probe passed')


if __name__=='__main__':main()
