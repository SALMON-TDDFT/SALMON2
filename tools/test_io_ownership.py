"""Synthetic Fortran I/O ownership checks and optional old-output byte controls.

No SCF, GS archive, physical density or wavefunction data is generated.
"""
import argparse,hashlib,json,os,subprocess,sys
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--reference-build',type=Path)
    p.add_argument('--directory-scope',action='store_true')
    p.add_argument('--fc',default='mpifort');p.add_argument('--ranks',type=int,default=2)
    p.add_argument('--checker',type=Path,default=Path(__file__).with_name('check_forbidden_eigensolvers.py'))
    args=p.parse_args();root=Path(__file__).resolve().parents[1]
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    checker=args.checker.resolve();fixture=root/'tests/native/test_io_ownership.f90'
    subprocess.run([sys.executable,str(checker),'--source-root',str(fixture.parent)],check=True)
    def compile_fixture(receipt,label,fixture=fixture):
        record=json.loads(receipt.read_text());library=Path(record['library'])
        if hashlib.sha256(library.read_bytes()).hexdigest()!=record['sha256']:raise ValueError('runtime hash changed')
        command=record.get('configure_command') or record['commands'][0]
        build=Path(command[command.index('-B')+1])
        directory=out/label;directory.mkdir();obj=directory/'fixture.o';exe=directory/'fixture'
        subprocess.run([args.fc,'-I'+str(build),'-c',str(fixture),'-o',str(obj)],check=True)
        for objects in (directory,build.parent):
            subprocess.run([sys.executable,str(checker),'--object-tree',str(objects)],check=True)
        subprocess.run([args.fc,str(obj),str(library),'-o',str(exe)],check=True)
        return exe
    new=compile_fixture(args.build,'new')
    old=compile_fixture(args.reference_build,'old') if args.reference_build else None
    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
    def run(exe,mode,variant,label):
        work=out/(label+'-'+mode);work.mkdir()
        with (work/'run.log').open('w') as log:
            subprocess.run(['mpiexec','-n',str(args.ranks),str(exe),mode,variant],cwd=work,
                           env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
        if 'IO ownership PASS' not in (work/'run.log').read_text():raise ValueError('fixture incomplete')
        return work
    modes=('units','dos','pdos','band','avs','avs-unsplit')
    for mode in modes:
        occupied=run(new,mode,'occupied','new')
        if mode=='units':continue
        for f in occupied.glob('host-*'):
            if f.read_text()!='before\nafter\n':raise ValueError('host file modified')
        reference=run(old or new,mode,'reference','reference')
        def outputs(path):
            return {f.name:f.read_bytes() for f in path.iterdir() if f.is_file() and f.name not in ('run.log',) and not f.name.startswith('host-')}
        if outputs(occupied)!=outputs(reference):raise ValueError('output bytes changed: '+mode)
    if args.directory_scope:
        exe=compile_fixture(args.build,'directory-fixture',root/'tests/native/test_directory_context.f90')
        work=out/'directory-run';work.mkdir()
        with (work/'run.log').open('w') as log:
            subprocess.run(['mpiexec','-n',str(args.ranks),str(exe),str(work.resolve())],cwd=work,
                           env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
        if 'directory scope/host/absolute/relative/failure/sibling PASS' not in (work/'run.log').read_text():
            raise ValueError('directory fixture incomplete')
    (out/'result.json').write_text(json.dumps(dict(status='pass',ranks=args.ranks,
        synthetic_only=True,GS_data_collected=False,units_A_B_A=True,directory_scopes_checked=args.directory_scope,host_units_1_77_100_101_preserved=True,
        DOS_PDOS_band_AVS_outputs_exact=True,previous_runtime_outputs_exact=bool(old)),indent=2)+'\n')
    print('Synthetic IO ownership and output parity PASS')


if __name__=='__main__':main()
