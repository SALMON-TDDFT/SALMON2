"""Interleaved full-support/adaptive RT measurements using archived GS seeds."""
import argparse, hashlib, json, math, os, platform, re, shlex, shutil
import statistics, subprocess, sys, time
from pathlib import Path
from run_rt import ROOT, HERE, SHARED, SHAPES, rt_input, parse_rt

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--seeds',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--repeat',type=int,default=3)
    p.add_argument('--mpiexec',default='mpiexec --bind-to none')
    a=p.parse_args()
    if a.repeat<1:p.error('repeat must be positive')
    a.binary=a.binary.resolve();a.seeds=a.seeds.resolve();a.output=a.output.resolve()
    a.output.mkdir(parents=True,exist_ok=False)
    previous=json.loads((a.seeds/'results.json').read_text())
    if not previous['complete']:raise RuntimeError('seed preparation incomplete')
    def verify_seeds():
        for prep in previous['preparations']:
            for name,sha in prep['payload_sha256'].items():
                if digest(a.seeds/prep['folder']/name)!=sha:raise RuntimeError('seed hash mismatch')
    verify_seeds()
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    sources=[HERE/'run_adaptive.py',HERE/'run_rt.py',HERE/'inputs.py',SHARED/'run.py',SHARED/'rank_probe.py',SHARED/'template.inp']
    data=dict(complete=False,production_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        platform=platform.platform(),binary=str(a.binary),binary_sha256=digest(a.binary),
        sources={str(f.relative_to(ROOT)):f.read_text() for f in sources},
        source_sha256={str(f.relative_to(ROOT)):digest(f) for f in sources},
        seed_results_sha256=digest(a.seeds/'results.json'),preparations=previous['preparations'],
        conditions=previous['conditions'],timing_note=previous['timing_note'],memory_note=previous['memory_note'],
        launcher=a.mpiexec,mpi_version=subprocess.check_output(shlex.split(a.mpiexec)+['--version'],text=True),
        threads={k:env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
        pseudo_sha256=digest(ROOT/'testsuites/pseudo/H_rps.dat'),
        build_cache=(a.binary.parent/'CMakeCache.txt').read_text(),runs=[],summary=[])
    data['conditions'].update(repeats=a.repeat,norm_fractions=[1.0,.999],full_support='paired reference',local_fft='auto')
    data['accuracy_note']='Energy width is recorded, not used to stop approximate-support timing runs. Original 1e-7 Ha threshold is reported as a flag.'
    def save():(a.output/'results.json').write_text(json.dumps(data,indent=2)+'\n')
    cases=[dict(shape=s,ranks=math.prod(s),suites=['weak']) for s in SHAPES]
    cases[6]['suites'].append('strong')
    cases.extend(dict(shape=(4,4,1),ranks=n,suites=['strong']) for n in (1,2,4,8))
    save()
    for case in cases:
        shape,ranks=case['shape'],case['ranks'];tag='x'.join(map(str,shape));seed=a.seeds/('gs-'+tag)
        for rep in range(1,a.repeat+1):
            # Reverse alternating repetitions to reduce systematic ordering bias.
            for mode,fraction in ([('full',1.),('adaptive',.999)] if rep%2 else [('adaptive',.999),('full',1.)]):
                name=f'{mode}-{tag}-mpi{ranks}-r{rep}';folder=a.output/name;folder.mkdir()
                inp=rt_input(shape,ranks).replace('&functional',f'&functional\n exx_mlwf_norm_fraction={fraction:.3f}d0\n exx_local_fft=\'auto\'')
                (folder/'inputfile').write_text(inp)
                shutil.copy(ROOT/'testsuites/pseudo/H_rps.dat',folder)
                (folder/'data_dcdft').symlink_to(os.path.relpath(seed/'data_dcdft',folder),target_is_directory=True)
                cmd=shlex.split(a.mpiexec)+['-n',str(ranks),sys.executable,str(SHARED/'rank_probe.py'),str(a.binary)]
                print('START',name,flush=True);start=time.perf_counter()
                with (folder/'inputfile').open() as fin,(folder/'output').open('w') as fout:
                    proc=subprocess.run(cmd,stdin=fin,stdout=fout,stderr=subprocess.STDOUT,cwd=folder,env=env,timeout=3600)
                if proc.returncode:raise RuntimeError('MPI failure '+name)
                wall=time.perf_counter()-start
                result=parse_rt(folder,ranks,max_energy_width=None)
                text=(folder/'output').read_text()
                pairs=[list(map(int,m)) for m in re.findall(r'EXX_ADAPTIVE local/global pairs/local FFT points \(orbital group 0\):\s+(\d+)\s+(\d+)\s+(\d+)',text)]
                supports=[list(map(float,m)) for m in re.findall(r'EXX_ADAPTIVE fraction/max radius/max norm loss:\s+(\S+)\s+(\S+)\s+(\S+)',text)]
                result.update(case,mode=mode,norm_fraction=fraction,repeat=rep,folder=name,launcher_wall_seconds=wall,
                    pair_counters=pairs,support_counters=supports,energy_width_within_1e_7=result['post_impulse_energy_width_ha']<=1e-7)
                data['runs'].append(result);save()
                print('DONE',name,'seconds',result['rt_max_seconds'],'MiB',result['peak_rank_bytes']/2**20,'energy width',result['post_impulse_energy_width_ha'],flush=True)
        for mode in ('full','adaptive'):
            rows=[r for r in data['runs'] if r['shape']==shape and r['ranks']==ranks and r['mode']==mode]
            summary=dict(case,mode=mode)
            for key in ('rt_max_seconds','propagation_max_seconds','peak_rank_bytes','launcher_wall_seconds','post_impulse_energy_width_ha'):
                vals=[r[key] for r in rows];summary[key]=dict(min=min(vals),median=statistics.median(vals),max=max(vals))
            data['summary'].append(summary)
        save()
    verify_seeds();data['complete']=True;save()

if __name__=='__main__':main()
