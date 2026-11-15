"""Original/current GS occupation comparison: BUILD_RECORD REFERENCE_RECORD NEW_OUTPUT."""
import argparse, hashlib, json, os, re, shutil, subprocess
from pathlib import Path
from test_gs_lifecycle import INPUT


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build','reference','output'):p.add_argument(name,type=Path)
    a=p.parse_args();root=Path(__file__).resolve().parents[1]
    records=[json.loads(v.read_text()) for v in (a.build,a.reference)]
    for r in records:
        if r.get('status')!='built-and-audited' or hashlib.sha256(Path(r['executable']).read_bytes()).hexdigest()!=r['sha256']:raise ValueError('unchanged audited builds required')
        if re.search(r'^#define (USE_MPI|USE_OPENACC)\b',(Path(r['executable']).parent/'config.h').read_text(),re.M):raise ValueError('serial CPU builds required')
    out=a.output.resolve()
    protected=[root]+[Path(r['source_root']) for r in records]+[Path(r['executable']).parent for r in records]
    if any(out.is_relative_to(p.resolve()) for p in protected):raise ValueError('new external output required')
    pseudo=root/'salmonpy/examples/si_intro_data/Si_rps.dat'
    if not pseudo.exists():pseudo=root.parent/'examples/si_intro_data/Si_rps.dat'
    out.mkdir(parents=True,exist_ok=False);cases=[]
    for xc,ne,temp in [('pz',3,-1),('pz',4,-1),('pz',4,0),('pz',3,.05),('hse06',4,-1),('hse06',4,.05)]:
        case=out/f'{xc}-ne{ne}-temp{temp}';case.mkdir()
        text=INPUT.replace("xc='hse06'",f"xc='{xc}'").replace('exx_pre_scf_threshold=1d8','exx_pre_scf_threshold=0')
        text=text.replace('nelec=4, nstate=2',f'nelec={ne}, nstate=4, temperature={temp}')
        text=text.replace('&atomic_red_coor','&parallel\n nproc_k=1, nproc_ob=1, nproc_rgrid=1,1,1\n/\n&atomic_red_coor').replace('&scf\n','&scf\n nscf_init_redistribution=0,\n')
        hashes={};residuals=[]
        try:
            for label,record in zip(('current','reference'),records):
                directory=case/label;directory.mkdir();shutil.copy2(pseudo,directory/'Si_rps.dat');(directory/'input.inp').write_text(text)
                with (directory/'run.log').open('w') as log:
                    result=subprocess.run([record['executable']],input=text,text=True,cwd=directory,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),stdout=log,stderr=subprocess.STDOUT,timeout=120)
                if result.returncode:raise RuntimeError(f'{case.name}/{label}: exit {result.returncode}; see {directory / "run.log"}')
                log=(directory/'run.log').read_text();residuals.append([l.split() for l in log.splitlines() if l.strip().startswith('iter and ')])
                if not residuals[-1] or temp>=0 and 'Fractional Occupation Numbers' not in log:raise AssertionError('required occupation path not reached')
                restart=directory/'data_for_restart'
                files=['gs_info.data','gs_eigen.data','gs_k.data']+[str(f.relative_to(directory)) for f in sorted(restart.rglob('*')) if f.is_file()]
                if not (restart/'occupation.bin').exists():raise AssertionError('missing occupation checkpoint')
                hashes[label]={f:hashlib.sha256((directory/f).read_bytes()).hexdigest() for f in files}
            if hashes['current']!=hashes['reference'] or residuals[0]!=residuals[1]:raise AssertionError(f'original/current mismatch: {case.name}')
            cases.append({'xc':xc,'nelec':ne,'temperature':temp,'sha256':hashes,'residuals':residuals[0]})
        finally:
            for label in ('current','reference'):
                shutil.rmtree(case/label/'data_for_restart',ignore_errors=True);(case/label/'Si_rps.dat').unlink(missing_ok=True)
    (out/'result.json').write_text(json.dumps({'status':'passed','cases':cases,'build':records[0],'reference':records[1]},indent=2)+'\n')
    print(f'GS occupation original comparison PASS ({len(cases)} cases)')

if __name__=='__main__':main()
