from pathlib import Path
import subprocess,json,resource,sys,os,shlex,tempfile,argparse
root=Path(__file__).resolve().parents[2]
tests=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--source',type=Path,default=root/'src/gs/dc/lcfo_diag_chefsi_complex.f90')
s=parser.parse_args().source.read_text()
a=s.index('  interface salmon_all_finite')
b=s.index('  end interface',a)+len('  end interface')
helpers=s[s.index('  pure logical function finite_real_1d'):s.index('end module lcfo_diag_chefsi_complex')]
a2=s.index('    flag=0')
b2=s.index('    call sync_status(flag)',a2)
body=s[a2:b2].replace('flag=0','ok=.true.').replace('flag=1','ok=.false.')
folder=tempfile.TemporaryDirectory(prefix='lcfo-finite-')
out=Path(folder.name)
probe='module finite_probe\nimplicit none\n'+s[a:b]+'\ncontains\nsubroutine probe(a,ok)\nimplicit none\ncomplex(8),intent(in)::a(:,:,:,:)\nlogical,intent(out)::ok\ncomplex(8),allocatable::hdiag_sym(:,:,:)\nassociate(hrow=>a)\nallocate(hdiag_sym(2,2,1))\nhdiag_sym=(1d0,2d0)\n'+body+'end associate\nend subroutine\n'+helpers+'end module\n'
wrappers=''
for rank in (2,3):
 dims=','.join([':']*rank)
 expression='salmon_all_finite(a)' if 'finite_complex_'+str(rank)+'d' in s else \
    'salmon_all_finite(real(a,8)).and.salmon_all_finite(aimag(a))'
 wrappers+='subroutine probe'+str(rank)+'(a,ok)\nimplicit none\ncomplex(8),intent(in)::a('+dims+')\nlogical,intent(out)::ok\nok='+expression+'\nend subroutine\n'
probe=probe.replace('end module\n',wrappers+'end module\n')
(out/'extracted.f90').write_text(probe)
cmd=shlex.split(os.environ.get('FC','gfortran'))
cmd+=shlex.split(os.environ.get('FFLAGS','-O0 -fcheck=all -fstack-arrays'))
cmd+=['extracted.f90',str(tests/'probe.f90'),'-o','probe']
subprocess.run(cmd,cwd=out,check=True)
def limit(): resource.setrlimit(resource.RLIMIT_STACK,(8*1024*1024,8*1024*1024))
r=subprocess.run(['./probe'],cwd=out,preexec_fn=limit,capture_output=True,text=True)
record={'exit':r.returncode,'stdout':r.stdout,'stderr':r.stderr,'stack_bytes':8*1024*1024}
print(json.dumps(record,indent=2))
sys.exit(0 if r.returncode==0 else 1)
