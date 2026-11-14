"""Link a native explicit-state exchange probe to all audited HSE application objects.

python tools/test_hse_explicit_state.py BUILD_RECORD NEW_OUTPUT [--threads 1]
No full solver/source is modified. Audits all linked application objects before
linking; uses an existing CMake GNU serial HSE build with its compiler profile.
"""
import argparse,hashlib,json,re,shlex,shutil,subprocess,sys,os
from pathlib import Path


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('build',type=Path);p.add_argument('output',type=Path);p.add_argument('--threads',type=int,default=1)
    a=p.parse_args();record=json.loads(a.build.read_text());engine=Path(record['executable']).resolve();build=engine.parent;source=Path(record['source_root']).resolve()
    if record.get('schema')!='salmonpy.hse-rt-state-build/1' or record.get('status')!='built-and-audited' or digest(engine)!=record['sha256']:raise ValueError('unchanged audited HSE build required')
    if a.threads<1:raise ValueError('positive threads required')
    config=(build/'config.h').read_text()
    if not re.search(r'^#define USE_HSE\b',config,re.M) or any(re.search(r'^#define '+feature+r'\b',config,re.M) for feature in ['USE_MPI','USE_OPENACC','USE_EXX_CUFFT']):
        raise ValueError('serial CPU HSE build required for this probe')
    root=Path(__file__).resolve().parents[1];checker=Path(__file__).with_name('check_forbidden_eigensolvers.py');out=a.output.resolve()
    if out.is_relative_to(source) or out.is_relative_to(build) or out.is_relative_to(root):raise ValueError('new external output required')
    subprocess.run([sys.executable,str(checker),'--source-root',str(source/'src')],check=True)
    out.mkdir(parents=True,exist_ok=False);src=out/'src';src.mkdir();shutil.copy2(root/'tests/native/test_hse_explicit_state.f90',src/'probe.f90')
    cache=(build/'CMakeCache.txt').read_text();fc=re.search(r'^CMAKE_Fortran_COMPILER:[^=]*=(.*)$',cache,re.M)[1].strip()
    if 'GNU Fortran' not in subprocess.check_output([fc,'--version'],text=True):raise ValueError('GNU compiler profile required')
    flags=(build/'src/CMakeFiles/salmon.dir/flags.make').read_text();options=shlex.split(re.search(r'^Fortran_FLAGS\s*=\s*(.*)$',flags,re.M)[1]);options=[v for v in options if not v.startswith('-J')]
    includes=shlex.split(re.search(r'^Fortran_INCLUDES\s*=\s*(.*)$',flags,re.M)[1]);commands=[]
    with (out/'build.log').open('w') as log:
        def run(cmd,**kw):
            commands.append([str(x) for x in cmd]);subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,**kw)
        run([sys.executable,checker,'--source-root',src])
        run([fc,*options,*includes,'-I'+str(build),'-c',src/'probe.f90','-o',out/'probe.o'],cwd=out)
        run([sys.executable,checker,'--object-tree',build]);run([sys.executable,checker,'--object-tree',out])
        link=shlex.split((build/'src/CMakeFiles/salmon.dir/link.txt').read_text());link[link.index('-o')+1]=str(out/'probe')
        linked=[]
        for i,v in enumerate(link):
            if v=='CMakeFiles/salmon.dir/main.f90.o':link[i]=str(out/'probe.o')
            elif v.endswith('.o'):
                path=(build/'src'/v).resolve()
                if not path.is_relative_to(build):raise ValueError('application object outside audited tree')
                linked.append(path)
        if str(out/'probe.o') not in link:raise ValueError('unsupported executable link recipe')
        run(link,cwd=build/'src')
        run([out/'probe'],cwd=out,env={**os.environ,'OMP_NUM_THREADS':str(a.threads),'OPENBLAS_NUM_THREADS':'1'})
    (out/'result.json').write_text(json.dumps(dict(status='passed',scope='real dense ACE adapter explicit RT routing at common material configuration; not full native multi-material',threads=a.threads,commands=commands,probe_sha256=digest(out/'probe'),linked_application_objects={str(f):digest(f) for f in linked+[out/'probe.o']}),indent=2)+'\n');print('Explicit-state HSE adapter probe passed')


if __name__=='__main__':main()
