"""Audit/link prescribed-density native tests using an existing runtime build.

All linked application objects (including the new test object) are audited.
No full runtime rebuild, prohibited driver, GS calculation or state export.
"""
import argparse,hashlib,json,os,shlex,subprocess
from pathlib import Path

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('build',type=Path);p.add_argument('output',type=Path)
p.add_argument('--test',default='test_nlcc_xc_variation')
p.add_argument('--expect-failure',action='store_true')
p.add_argument('--ranks',type=int,nargs='+',default=[1])
a=p.parse_args();root=Path(__file__).resolve().parents[1];build=a.build.resolve()
source=dict(line.split('=',1) for line in (build/'CMakeCache.txt').read_text().splitlines()
            if '=' in line and not line.startswith(('#','//')))['CMAKE_HOME_DIRECTORY:INTERNAL']
checker=root/'tools/check_forbidden_eigensolvers.py'
def check(*args):subprocess.run(['python3',str(checker),*map(str,args)],check=True)
check('--source-root',Path(source)/'src')
check('--source-root',root/'tests/runtime')
a.output.mkdir(parents=True,exist_ok=False);out=a.output.resolve()
obj=out/(a.test+'.o');exe=out/a.test
link=shlex.split((build/'src/CMakeFiles/salmon.dir/link.txt').read_text())
compiler=link[0]
subprocess.run([compiler,'-fopenmp','-ffree-line-length-none','-I'+str(build),'-c',
                str(root/'tests/runtime'/(a.test+'.f90')),'-o',str(obj)],check=True)
check('--object-tree',build/'src')
check('--object-tree',out)
cmd=[];skip=False
for arg in link:
    if skip:skip=False;continue
    if arg=='-o':skip=True;continue
    if arg=='CMakeFiles/salmon.dir/main.f90.o':continue
    if arg.endswith('.o') and not (build/'src'/arg).resolve().is_relative_to(build/'src'):
        raise ValueError('unaudited object outside runtime object tree')
    cmd.append(arg)
cmd += [str(obj),'-o',str(exe)]
subprocess.run(cmd,cwd=build/'src',check=True)
record={'test':a.test,'source_root':source,'executable_sha256':hashlib.sha256(exe.read_bytes()).hexdigest()}
results=[]
for ranks in a.ranks:
 for threads in (1,2):
    command=[str(exe)] if ranks==1 else ['/opt/homebrew/bin/mpiexec','--bind-to','none','-n',str(ranks),str(exe)]
    run=subprocess.run(command,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
        env={**os.environ,'OMP_NUM_THREADS':str(threads),'OPENBLAS_NUM_THREADS':'1'},timeout=120)
    print(run.stdout,flush=True)
    results.append(dict(ranks=ranks,threads=threads,returncode=run.returncode,stdout=run.stdout))
record['runs']=results
(out/'test-record.json').write_text(json.dumps(record,indent=2)+'\n')
if any((r['returncode']==0)==a.expect_failure for r in results):raise SystemExit('unexpected native test outcome')
print(json.dumps({k:v for k,v in record.items() if k!='runs'}))
