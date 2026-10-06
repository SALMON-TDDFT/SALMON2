"""H2 PBEh40 supercell SCF weak and strong scaling, sequential fresh MPI jobs."""
import argparse,hashlib,json,math,os,platform,re,shlex,shutil,statistics,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SHAPES=[(1,1,1),(2,1,1),(4,1,1),(8,1,1),(16,1,1),(2,2,1),(4,4,1),(2,2,2)]
LAYOUT={1:(1,1,1),2:(1,2,1),4:(1,2,2),8:(1,4,2),16:(1,4,4)}

def input_text(shape,ranks):
    molecules=math.prod(shape)
    s=(HERE/'template.inp').read_text()
    s=s.replace('nproc_rgrid=1,1,1','nproc_rgrid='+','.join(map(str,LAYOUT[ranks])))
    s=s.replace('al=8d0,8d0,8d0','al='+','.join(str(8*n)+'d0' for n in shape))
    s=s.replace('num_rgrid=16,16,16','num_rgrid='+','.join(str(16*n) for n in shape))
    s=s.replace('nstate=1',f'nstate={molecules}').replace('nelec=2',f'nelec={2*molecules}').replace('natom=2',f'natom={2*molecules}')
    s=s[:s.index('&atomic_coor')]+ '&atomic_coor\n'
    for z in range(shape[2]):
        for y in range(shape[1]):
            for x in range(shape[0]):
                for delta in (-.7,.7):s+=f" 'H' {8*x+4+delta:.8f} {8*y+4:.8f} {8*z+4:.8f} 1\n"
    return s+'/\n'

def parse(folder,ranks):
    text=(folder/'output').read_text()
    match=re.search(r'#GS converged at\s+(\d+)\s+:\s+(\S+)',text)
    if not match or 'end SALMON' not in text:raise RuntimeError('SCF not converged: '+str(folder))
    info=(folder/'h2_info.data').read_text()
    iterations=int(re.search(r'Total number of iteration\s*=\s*(\d+)',info)[1])
    scf=re.search(r'^\s*scf iterations\s+(\S+)\s+(\S+)',text,re.M)
    if scf is None:raise RuntimeError('missing SCF timer: '+str(folder))
    rows=[json.loads((folder/f'rank-{i}.json').read_text()) for i in range(ranks)]
    if any(r['returncode']!=0 or r['peak_rss_bytes']<=0 for r in rows):raise RuntimeError('rank failure')
    localizations=[]
    for m in re.finditer(r'(EXX_SPATIAL|EXX_WANNIER) refresh/iterations/status/spread/gradient/overlap:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)',text):
        localizations.append(dict(backend=m[1],update=int(m[2]),iterations=int(m[3]),status=int(m[4]),spread=float(m[5]),gradient=float(m[6])))
    result=dict(iterations=iterations,residual=float(match[2]),energy_ev=float(re.search(r'Total energy \(eV\) =\s*(\S+)',info)[1]),
        total_root_seconds=float(re.search(r'total calculation time,(\S+)',text)[1]),scf_max_seconds=float(scf[2]),
        per_rank=rows,localizations=localizations)
    for k in ('residual','energy_ev','total_root_seconds','scf_max_seconds'):
        if not math.isfinite(result[k]):raise RuntimeError('nonfinite '+k)
    if result['residual']>=1e-10:raise RuntimeError('SCF tolerance failed')
    result['scf_seconds_per_iteration']=result['scf_max_seconds']/iterations
    result['peak_rank_bytes']=max(r['peak_rss_bytes'] for r in rows)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--repeat',type=int,default=3);p.add_argument('--mpiexec',default='mpiexec --bind-to none')
    p.add_argument('--case',action='append',dest='selected',help='Subset case, e.g. 4x4x1:8 (repeatable)')
    p.add_argument('--timeout',type=float,default=3600);p.add_argument('--pilot',action='store_true')
    a=p.parse_args();a.binary=a.binary.resolve();a.output=a.output.resolve()
    a.output.mkdir(parents=True,exist_ok=True)
    if (a.output/'results.json').exists():p.error('use a new output directory to preserve measurements')
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    cases=[dict(shape=s,ranks=math.prod(s),suites=['weak']) for s in SHAPES]
    for n in (1,2,4,8):cases.append(dict(shape=(4,4,1),ranks=n,suites=['strong']))
    cases[6]['suites'].append('strong')
    if a.selected:
        cases=[c for c in cases if 'x'.join(map(str,c['shape']))+':'+str(c['ranks']) in a.selected]
        if len(cases)!=len(set(a.selected)):p.error('unknown or duplicate case selection')
    if a.pilot:cases=[dict(shape=(4,4,1),ranks=16,suites=['pilot'])];a.repeat=1
    sources={str(f.relative_to(ROOT)):f.read_text() for f in HERE.iterdir() if f.suffix in ('.py','.inp')}
    data=dict(complete=False,platform=platform.platform(),logical_cpus=os.cpu_count(),binary=str(a.binary),
        binary_sha256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),
        production_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        pseudo_sha256=hashlib.sha256((ROOT/'testsuites/pseudo/H_rps.dat').read_bytes()).hexdigest(),
        sources=sources,launcher=a.mpiexec,mpi_version=subprocess.check_output(shlex.split(a.mpiexec)+['--version'],text=True),
        threads={k:env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
        conditions=dict(functional='pbeh40',base_cell_bohr=[8,8,8],bond_bohr=1.4,spacing_bohr=.5,coulomb_radius_bohr=4,
            mlwf_radius=0,mlwf_interval=5,mlwf_maxiter=100,mlwf_tolerance=1e-7,kgrid=[1,1,1],scf_threshold=1e-10,alpha_mb=.1,method_init_wf='gauss10',nscf=2000,
            repeats=a.repeat,orbital_groups=1,spatial_layouts=LAYOUT),runs=[],summary=[])
    cache=a.binary.parent/'CMakeCache.txt'
    if cache.exists():data['build_cache']=cache.read_text()
    def save():(a.output/'results.json').write_text(json.dumps(data,indent=2)+'\n')
    save()
    for case in cases:
        shape,ranks=case['shape'],case['ranks'];name='x'.join(map(str,shape))+f'-mpi{ranks}'
        subset=[]
        for rep in range(a.repeat):
            folder=a.output/(name+f'-r{rep+1}');folder.mkdir()
            (folder/'inputfile').write_text(input_text(shape,ranks))
            shutil.copy(ROOT/'testsuites/pseudo/H_rps.dat',folder)
            cmd=shlex.split(a.mpiexec)+['-n',str(ranks),sys.executable,str(HERE/'rank_probe.py'),str(a.binary)]
            print('START',folder.name,flush=True);start=time.perf_counter()
            with (folder/'inputfile').open() as inp,(folder/'output').open('w') as out:
                proc=subprocess.run(cmd,stdin=inp,stdout=out,stderr=subprocess.STDOUT,cwd=folder,env=env,timeout=a.timeout)
            if proc.returncode:raise RuntimeError('MPI failure: '+str(folder))
            result=parse(folder,ranks)
            result.update(case,repeat=rep+1,folder=folder.name,launcher_wall_seconds=time.perf_counter()-start)
            data['runs'].append(result);subset.append(result);save()
            print('DONE',folder.name,'iters',result['iterations'],'SCF s',result['scf_max_seconds'],flush=True)
        summary=dict(case)
        for key in ('total_root_seconds','scf_max_seconds','scf_seconds_per_iteration','peak_rank_bytes','launcher_wall_seconds','energy_ev','iterations'):
            vals=[r[key] for r in subset];summary[key]=dict(median=statistics.median(vals),min=min(vals),max=max(vals))
        data['summary'].append(summary);save()
    strong=[r for r in data['runs'] if 'strong' in r['suites']]
    if strong:
        delta=max(r['energy_ev'] for r in strong)-min(r['energy_ev'] for r in strong)
        data['strong_energy_range_ev']=delta
        if delta>1e-5:raise RuntimeError('strong-scaling energies disagree: '+str(delta))
    data['complete']=True;save()

if __name__=='__main__':main()
