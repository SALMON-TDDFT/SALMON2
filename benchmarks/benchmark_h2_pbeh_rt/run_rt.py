"""Converged full-cell DC seed -> native mesh PBEh impulse, 16 RT steps."""
import argparse,hashlib,json,math,os,platform,re,shlex,shutil,statistics,subprocess,sys,time
from pathlib import Path
from inputs import SHAPES,LAYOUT,gs_input,rt_input
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SHARED=HERE.parent/'benchmark_h2_pbeh'

def table(path):
    rows=[[float(x) for x in line.split()] for line in path.read_text().splitlines() if line.strip() and not line.startswith('#')]
    if not rows or not all(math.isfinite(x) for row in rows for x in row):raise RuntimeError('nonfinite or empty '+str(path))
    return rows

def timing(text,label):
    m=re.search(r'^\s*'+re.escape(label)+r'\s+(\S+)\s+(\S+)',text,re.M)
    if m is None:raise RuntimeError('missing timer '+label)
    value=float(m[2])
    if not math.isfinite(value) or value<0:raise RuntimeError('invalid timer '+label)
    return value

def parse_rt(folder,ranks,max_energy_width=1e-7):
    text=(folder/'output').read_text()
    if 'end SALMON' not in text or 'end complex DC-LCFO wavefunction reconstruction' not in text or 'Native LCFO RT active' in text:
        raise RuntimeError('not a completed mesh RT job: '+str(folder))
    obs=table(folder/'h2_rt.data');energy=table(folder/'h2_rt_energy.data')
    if len(obs)!=16 or len(energy)!=17:raise RuntimeError('incorrect RT step count')
    if any(abs(row[0]-.02*(i+1))>1e-9 for i,row in enumerate(obs)):raise RuntimeError('incorrect time grid')
    if any(abs(row[1]-1e-4)>1e-12 for row in obs):raise RuntimeError('incorrect impulse')
    if any(abs(row[0]-.02*i)>1e-9 for i,row in enumerate(energy)):raise RuntimeError('incorrect energy time grid')
    rank_rows=[json.loads((folder/f'rank-{i}.json').read_text()) for i in range(ranks)]
    if sorted(r['rank'] for r in rank_rows)!=list(range(ranks)) or any(r['returncode']!=0 or r['peak_rss_bytes']<=0 for r in rank_rows):raise RuntimeError('rank failure')
    post=[row[1] for row in energy[1:]];drift=max(post)-min(post)
    if max_energy_width is not None and drift>max_energy_width:raise RuntimeError(f'post-impulse energy drift exceeds {max_energy_width} Ha')
    result=dict(rt_max_seconds=timing(text,'rt iterations'),propagation_max_seconds=timing(text,'time propagation'),
        total_root_seconds=float(re.search(r'total calculation time,(\S+)',text)[1]),
        peak_rank_bytes=max(r['peak_rss_bytes'] for r in rank_rows),per_rank=rank_rows,
        observables=obs,energies=energy,post_impulse_energy_width_ha=drift,
        localization_lines=[line for line in text.splitlines() if 'refresh/iterations/status/spread' in line])
    if not math.isfinite(result['total_root_seconds']) or result['total_root_seconds']<0:raise RuntimeError('invalid total timer')
    result['rt_seconds_per_step']=result['rt_max_seconds']/16
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--repeat',type=int,default=3);p.add_argument('--mpiexec',default='mpiexec --bind-to none')
    p.add_argument('--timeout',type=float,default=3600)
    a=p.parse_args();a.binary=a.binary.resolve();a.output=a.output.resolve()
    if a.repeat<1 or a.timeout<=0:p.error('repeat and timeout must be positive')
    a.output.mkdir(parents=True,exist_ok=True)
    if (a.output/'results.json').exists():p.error('use a fresh output directory')
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    sources={str(f.relative_to(ROOT)):f.read_text() for f in [HERE/'run_rt.py',HERE/'inputs.py',SHARED/'run.py',SHARED/'template.inp',SHARED/'rank_probe.py']}
    data=dict(complete=False,platform=platform.platform(),logical_cpus=os.cpu_count(),binary=str(a.binary),
        binary_sha256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),
        production_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        sources=sources,source_sha256={k:hashlib.sha256(v.encode()).hexdigest() for k,v in sources.items()},
        pseudo_sha256=hashlib.sha256((ROOT/'testsuites/pseudo/H_rps.dat').read_bytes()).hexdigest(),
        launcher=a.mpiexec,mpi_version=subprocess.check_output(shlex.split(a.mpiexec)+['--version'],text=True),
        threads={k:env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
        conditions=dict(functional='pbeh40',base_cell_bohr=[8,8,8],bond_bohr=1.4,spacing_bohr=.5,coulomb_radius_bohr=4,
            full_support=True,mlwf_interval=5,mlwf_maxiter=100,mlwf_tolerance=1e-7,kgrid=[1,1,1],
            dt_au=.02,steps=16,impulse_au=1e-4,direction=[1,0,0],fixed_ions=True,propagator='hse_taylor4',
            seed='one full-cell DC fragment with zero buffer; occupied-only; zero-temperature GS',
            repeats=a.repeat,orbital_groups=1,spatial_layouts=LAYOUT),
        timing_note='Minimum max-rank native RT-iteration timer; fixed 16 steps. Includes per-step energy/current output. Separate propagation timer and total-root timer. Launcher wall excludes output parsing.',
        memory_note='Median max-rank lifetime RSS of RT SALMON processes only, including seed loading/reconstruction; GS preparation is a separate process and excluded. Not isolated loop-only RSS or simultaneous node RSS.',
        preparations=[],runs=[],summary=[])
    cache=a.binary.parent/'CMakeCache.txt'
    if cache.exists():data['build_cache']=cache.read_text()
    def save():(a.output/'results.json').write_text(json.dumps(data,indent=2)+'\n')
    def execute(folder,inp,ranks,seed=None):
        folder.mkdir();(folder/'inputfile').write_text(inp);shutil.copy(ROOT/'testsuites/pseudo/H_rps.dat',folder)
        if seed:(folder/'data_dcdft').symlink_to(os.path.relpath(seed/'data_dcdft',folder),target_is_directory=True)
        cmd=shlex.split(a.mpiexec)+['-n',str(ranks),sys.executable,str(SHARED/'rank_probe.py'),str(a.binary)]
        print('START',folder.name,flush=True);start=time.perf_counter()
        with (folder/'inputfile').open() as inp,(folder/'output').open('w') as out:
            proc=subprocess.run(cmd,stdin=inp,stdout=out,stderr=subprocess.STDOUT,cwd=folder,env=env,timeout=a.timeout)
        if proc.returncode:raise RuntimeError('MPI failure: '+str(folder))
        return time.perf_counter()-start
    save();seeds={}
    for shape in SHAPES:
        ranks=min(math.prod(shape),8);name='gs-'+'x'.join(map(str,shape));folder=a.output/name
        wall=execute(folder,gs_input(shape,ranks),ranks)
        text=(folder/'output').read_text()
        scf=re.findall(r'DC #SCF =\s+(\d+)\s+Total Energy =\s+(\S+)\s+diff =\s+(\S+)',text)
        if 'end SALMON' not in text or 'end DC-LCFO complex' not in text or not scf or not math.isfinite(float(scf[-1][2])) or float(scf[-1][2])>=1e-10:
            raise RuntimeError('unconverged or incomplete seed: '+name)
        hashes={str(f.relative_to(folder)):hashlib.sha256(f.read_bytes()).hexdigest() for f in (folder/'data_dcdft').rglob('*') if f.is_file()}
        if not any(k.endswith('wavefunctions.bin') for k in hashes):raise RuntimeError('seed wavefunctions missing')
        data['preparations'].append(dict(shape=shape,ranks=ranks,folder=name,wall_seconds=wall,iterations=int(scf[-1][0]),
            final_dc_residual=float(scf[-1][2]),energy_ev_printed=float(scf[-1][1]),payload_sha256=hashes))
        seeds[shape]=folder;save();print('PREPARED',name,'iterations',scf[-1][0],flush=True)
    cases=[dict(shape=s,ranks=math.prod(s),suites=['weak']) for s in SHAPES];cases[6]['suites'].append('strong')
    cases.extend(dict(shape=(4,4,1),ranks=n,suites=['strong']) for n in (1,2,4,8))
    for case in cases:
        shape,ranks=case['shape'],case['ranks'];subset=[]
        for rep in range(a.repeat):
            name='rt-'+'x'.join(map(str,shape))+f'-mpi{ranks}-r{rep+1}';folder=a.output/name
            wall=execute(folder,rt_input(shape,ranks),ranks,seeds[shape]);result=parse_rt(folder,ranks)
            result.update(case,repeat=rep+1,folder=name,seed=seeds[shape].name,launcher_wall_seconds=wall)
            data['runs'].append(result);subset.append(result);save()
            print('RT DONE',name,'16-step seconds',result['rt_max_seconds'],'peak MiB',result['peak_rank_bytes']/2**20,flush=True)
        summary=dict(case)
        for key in ('rt_max_seconds','rt_seconds_per_step','propagation_max_seconds','total_root_seconds','peak_rank_bytes','launcher_wall_seconds'):
            v=[r[key] for r in subset];summary[key]=dict(min=min(v),median=statistics.median(v),max=max(v))
        data['summary'].append(summary);save()
    strong=[r for r in data['runs'] if 'strong' in r['suites']];ref=next(r for r in strong if r['ranks']==1)
    current=max(abs(a-b) for r in strong for row,base in zip(r['observables'],ref['observables']) for a,b in zip(row[13:16],base[13:16]))
    energy=max(abs(row[1]-base[1]) for r in strong for row,base in zip(r['energies'],ref['energies']))
    data['strong_current_max_error_au']=current;data['strong_energy_max_error_ha']=energy
    if current>1e-10 or energy>1e-8:raise RuntimeError('strong-layout RT trajectories disagree')
    # Confirm the common seed payloads were not modified by RT readers.
    for prep in data['preparations']:
        for f,h in prep['payload_sha256'].items():
            if hashlib.sha256((a.output/prep['folder']/f).read_bytes()).hexdigest()!=h:raise RuntimeError('seed payload changed')
    data['complete']=True;save()

if __name__=='__main__':main()
