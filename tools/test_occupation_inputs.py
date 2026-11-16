"""Explicit local occupation/original comparison: BUILD_RECORD NEW_OUTPUT."""
import argparse, hashlib, json, os, re, shlex, subprocess, sys
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build','output'):p.add_argument(name,type=Path)
    p.add_argument('--dip',action='store_true',help='test dipole electron-count diagnostics and strict bounds')
    a=p.parse_args();root=Path(__file__).resolve().parents[1]
    record=json.loads(a.build.read_text());exe=Path(record['executable']);build=exe.parent
    if record.get('status')!='built-and-audited' or hashlib.sha256(exe.read_bytes()).hexdigest()!=record['sha256']:raise ValueError('unchanged audited build required')
    out=a.output.resolve()
    if any(out.is_relative_to(p.resolve()) for p in (root,build,Path(record['source_root']))):raise ValueError('new external output required')
    if re.search(r'^#define (USE_MPI|USE_OPENACC)\b',(build/'config.h').read_text(),re.M):raise ValueError('serial CPU build required')
    target='rt/dip.f90' if a.dip else 'gs/occupation.f90'
    original=root/'runtime/sources/hse-fast-reference/src'/target
    manifest=json.loads((root/'runtime/source-trees/hse-explicit-routing-replacements.json').read_text())
    original_hash=hashlib.sha256(original.read_bytes()).hexdigest()
    if original_hash!=manifest['files']['src/'+target]['before_sha256']:raise ValueError('original source hash mismatch')
    out.mkdir(parents=True,exist_ok=False)
    module='dip' if a.dip else 'occupation'
    reference=re.sub(r'(?im)(\bmodule\s+)'+module+r'\b',r'\1reference_'+module,original.read_text())
    (out/'reference.f90').write_text(reference)
    (out/'probe.f90').write_bytes((root/('tests/native/test_dip_occupation_inputs.f90' if a.dip else 'tests/native/test_occupation_inputs.f90')).read_bytes())
    flags=(build/'src/CMakeFiles/salmon.dir/flags.make').read_text()
    options=[v for v in shlex.split(re.search(r'^Fortran_FLAGS\s*=\s*(.*)$',flags,re.M)[1]) if not v.startswith('-J')]
    includes=shlex.split(re.search(r'^Fortran_INCLUDES\s*=\s*(.*)$',flags,re.M)[1])
    compiler=re.search(r'^CMAKE_Fortran_COMPILER:[^=]*=(.*)$',(build/'CMakeCache.txt').read_text(),re.M)[1].strip()
    commands=[]
    with (out/'build.log').open('w') as log:
        def run(command,**kwargs):
            commands.append([str(v) for v in command]);subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT,timeout=120,**kwargs)
        checker=root/'tools/check_forbidden_eigensolvers.py'
        run([sys.executable,checker,'--source-root',out])
        for name in ('reference','probe'):run([compiler,*options,*includes,'-I'+str(out),'-J'+str(out),'-c',out/f'{name}.f90','-o',out/f'{name}.o'],cwd=out)
        run([sys.executable,checker,'--object-tree',build]);run([sys.executable,checker,'--object-tree',out])
        link=shlex.split((build/'src/CMakeFiles/salmon.dir/link.txt').read_text());link[link.index('-o')+1]=str(out/'probe')
        link=[str(out/'probe.o') if v=='CMakeFiles/salmon.dir/main.f90.o' else v for v in link];link.append(str(out/'reference.o'))
        run(link,cwd=build/'src')
    result=subprocess.run([str(out/'probe')],cwd=out,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),capture_output=True,text=True,timeout=30)
    (out/'run.log').write_text(result.stdout+result.stderr);result.check_returncode()
    if ('Dip occupation isolation PASS' if a.dip else 'Occupation input isolation PASS') not in result.stdout:raise AssertionError('probe did not finish')
    if a.dip:
        for mode in ('lower','upper'):
            for case_id in range(1,5):
                check=subprocess.run([str(out/'probe'),mode,str(case_id)],cwd=out,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),capture_output=True,text=True,timeout=30)
                (out/(mode+str(case_id)+'.log')).write_text(check.stdout+check.stderr)
                if check.returncode!=0 or 'Ne=' not in check.stdout or 'PASS' in check.stdout:raise AssertionError('strict boundary was not rejected by the legacy STOP')
    (out/'result.json').write_text(json.dumps({'status':'passed','scope':'dip diagnostic' if a.dip else 'occupation solver','build':record,'original_sha256':original_hash,'probe_sha256':hashlib.sha256((out/'probe.f90').read_bytes()).hexdigest(),'commands':commands},indent=2)+'\n')
    print('Occupation input isolation / original numerical equality PASS')

if __name__=='__main__':main()
