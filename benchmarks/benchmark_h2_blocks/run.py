"""Buffered DC preparation and paired native RT scaling: eight H2 per core."""
import copy
import argparse, hashlib, json, math, os, platform, re, shlex, shutil, statistics, subprocess, sys, time
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE.parent/'benchmark_h2_pbeh_rt'))
from run_rt import parse_rt
from inputs import rt_input
from run import SHAPES, LAYOUT, input_text
SHARED=HERE.parent/'benchmark_h2_pbeh'

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def geometry(shape):return tuple(2*n for n in shape)
def fragment_states(shape):return math.prod(3 if n>1 else 2 for n in shape)+4

def gs_input(shape):
    ranks=math.prod(shape)
    s=input_text(geometry(shape),1).replace("yn_dc='n'","yn_dc='y'")
    s=s.replace("method_init_wf='gauss10'","method_init_wf='random'")
    s=s.replace('&system','&system\n temperature_k=300d0')
    s=s.replace('&functional',"&functional\n yn_exx_dc_mlwf='n'\n exx_pre_scf_threshold=1d-4\n exx_mlwf_norm_fraction=1d0")
    s+='\n&dc\n num_fragment='+','.join(map(str,shape))
    s+='\n num_rgrid_buffer='+','.join('8' if n>1 else '0' for n in shape)
    s+='\n nproc_rgrid_tot='+','.join(map(str,shape))
    s+=f"\n nstate_frag={fragment_states(shape)}\n yn_dc_lcfo='y'\n lcfo_eigensolver='scalapack'\n energy_cut=100d0\n lambda_cut=1d-7\n/\n"
    assert s.count(" 'H'")==16*ranks
    return s

def rt_block(shape,ranks,fraction,pair_tolerance=None,ace_support='occupied'):
    if ace_support not in ('occupied','source'):raise ValueError('invalid ACE support')
    if ace_support=='source' and pair_tolerance is not None:raise ValueError('source ACE requires exact support pairs')
    if pair_tolerance is not None and (not math.isfinite(pair_tolerance) or pair_tolerance<0):
        raise ValueError('pair tolerance must be finite and nonnegative')
    s=rt_input(geometry(shape),ranks).replace('exx_mlwf_tolerance=1d-7','exx_mlwf_tolerance=1d-6').replace('exx_mlwf_maxiter=100','exx_mlwf_maxiter=1000')
    s=s.replace('num_fragment=1,1,1','num_fragment='+','.join(map(str,shape)))
    s=s.replace('num_rgrid_buffer=0,0,0','num_rgrid_buffer='+','.join('8' if n>1 else '0' for n in shape))
    s=s.replace(f'nstate_frag={8*math.prod(shape)}',f'nstate_frag={fragment_states(shape)}')
    if fraction<1 and ace_support=='source':
        s=s.replace('&functional',"&functional\n exx_ace_support='source'")
    if fraction<1 and pair_tolerance is not None:
        value=f'{pair_tolerance:.17e}'.replace('e','d')
        s=s.replace('&functional',f"&functional\n exx_pair_screening='on'\n exx_pair_tolerance={value}")
    return s.replace('&functional',f"&functional\n exx_mlwf_norm_fraction={fraction}d0\n exx_local_fft='auto'")

def reference_record(row,origin,results_path,results_hash):
    """Preserve provenance when explicitly reusing a historical full reference."""
    if row['mode']!='full':raise ValueError('cross-binary reuse is only for full references')
    result=copy.deepcopy(row)
    result['binary_sha256']=row.get('binary_sha256',origin['binary_sha256'])
    result['production_commit']=row.get('production_commit',origin['production_commit'])
    result['reference_source']=dict(results_path=str(results_path),results_sha256=results_hash)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--full-source',type=Path,help='explicitly reuse completed historical full references, preserving their executable provenance')
    p.add_argument('--ace-support',choices=('occupied','source'),default='occupied',help='ACE construction vectors for adaptive RT only')
    p.add_argument('--pair-tolerance',type=float,help='screen adaptive RT only; default is off')
    p.add_argument('--large-repeat',type=int);p.add_argument('--rt-source',type=Path);p.add_argument('--prepare-only',action='store_true');p.add_argument('--repeat',type=int,default=3);p.add_argument('--pilot',action='store_true')
    p.add_argument('--mpiexec',default='/opt/homebrew/bin/mpiexec --bind-to none')
    p.add_argument('--allow-seed-binary-change',action='store_true');p.add_argument('--seed-source',type=Path);p.add_argument('--resume',action='store_true');p.add_argument('--generate-only',action='store_true')
    a=p.parse_args();a.binary=a.binary.resolve();a.output=a.output.resolve()
    if a.repeat<1 or (a.large_repeat is not None and a.large_repeat<1):p.error('repeat must be positive')
    if a.pair_tolerance is not None and (not math.isfinite(a.pair_tolerance) or a.pair_tolerance<0):p.error('invalid pair tolerance')
    if a.ace_support=='source' and a.pair_tolerance is not None:p.error('source ACE requires exact support pairs')
    shapes=SHAPES[:2] if a.pilot else SHAPES
    cases=[dict(shape=list(s),ranks=math.prod(s),suites=['weak']) for s in shapes]
    if not a.pilot:
        cases[6]['suites'].append('strong')
        cases.extend(dict(shape=[4,4,1],ranks=n,suites=['strong']) for n in (8,4,2,1))
    for c in cases:c['repeats']=a.large_repeat if a.large_repeat is not None and 8*math.prod(c['shape'])>=64 else a.repeat
    if a.generate_only:
        for shape in shapes:
            print('geometry',shape,'H2',8*math.prod(shape),'states/fragment',fragment_states(shape),'input bytes',len(gs_input(shape)))
        return
    a.output.mkdir(parents=True,exist_ok=True)
    result_file=a.output/'results.json'
    sources=[Path(__file__).resolve(),SHARED/'run.py',SHARED/'rank_probe.py',SHARED/'template.inp',HERE.parent/'benchmark_h2_pbeh_rt'/'inputs.py',HERE.parent/'benchmark_h2_pbeh_rt'/'run_rt.py']
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    binary_hash=digest(a.binary)
    if result_file.exists():
        if not a.resume:p.error('existing output requires --resume')
        data=json.loads(result_file.read_text())
        if data['conditions'].get('ace_support','occupied')!=a.ace_support:raise RuntimeError('ACE support changed')
        if data['conditions'].get('pair_tolerance')!=a.pair_tolerance:raise RuntimeError('screening configuration changed')
        if data['binary_sha256']!=binary_hash or data['repeats']!=a.repeat or data['cases']!=cases:raise RuntimeError('resume configuration mismatch')
        if data['sources']!={str(f.relative_to(ROOT)):f.read_text() for f in sources}:raise RuntimeError('source changed during experiment')
        if data['launcher']!=a.mpiexec or data['mpi_version']!=subprocess.check_output(shlex.split(a.mpiexec)+['--version'],text=True):raise RuntimeError('MPI configuration changed')
        if data['pseudo_sha256']!=digest(ROOT/'testsuites/pseudo/H_rps.dat'):raise RuntimeError('pseudopotential changed')
        if any(data['conditions']['threads'][k]!=env[k] for k in data['conditions']['threads']):raise RuntimeError('thread configuration changed')
    else:
        data=dict(complete=False,binary_sha256=binary_hash,production_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
          sources={str(f.relative_to(ROOT)):f.read_text() for f in sources},build_cache=(a.binary.parent/'CMakeCache.txt').read_text(),
          platform=platform.platform(),launcher=a.mpiexec,mpi_version=subprocess.check_output(shlex.split(a.mpiexec)+['--version'],text=True),
          pseudo_sha256=digest(ROOT/'testsuites/pseudo/H_rps.dat'),repeats=a.repeat,large_repeats=a.large_repeat,cases=cases,
          conditions=dict(core_bohr=[16]*3,spacing_bohr=.5,H2_per_core=8,buffer_split_bohr=4,coulomb_radius_bohr=4,functional='pbeh40',
            gs_temperature_k=300,dc_mlwf=False,pre_scf_threshold=1e-4,gs_threshold=1e-10,dt=.02,steps=16,impulse=1e-4,fractions=[1.,.999],
            mlwf_interval=5,mlwf_maxiter=1000,mlwf_tolerance=1e-6,pair_tolerance=a.pair_tolerance,ace_support=a.ace_support,threads={k:env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')}),
          preparations=[],runs=[])
    def save():
        tmp=result_file.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(result_file)
    def execute(name,inp,ranks,seed=None):
        folder=a.output/name
        if folder.exists():raise RuntimeError('incomplete prior job preserved; use a fresh output: '+str(folder))
        folder.mkdir();(folder/'inputfile').write_text(inp);shutil.copy(ROOT/'testsuites/pseudo/H_rps.dat',folder)
        if seed:(folder/'data_dcdft').symlink_to(os.path.relpath(seed/'data_dcdft',folder),target_is_directory=True)
        cmd=shlex.split(a.mpiexec)+['-n',str(ranks),sys.executable,str(SHARED/'rank_probe.py'),str(a.binary)]
        print('START',name,flush=True);start=time.perf_counter()
        with (folder/'inputfile').open() as fin,(folder/'output').open('w') as fout:
            proc=subprocess.run(cmd,stdin=fin,stdout=fout,stderr=subprocess.STDOUT,cwd=folder,env=env)
        if proc.returncode:raise RuntimeError('MPI failed '+name)
        rows=[json.loads((folder/f'rank-{i}.json').read_text()) for i in range(ranks)]
        if any(r['returncode'] or r['peak_rss_bytes']<=0 for r in rows):raise RuntimeError('rank failed '+name)
        return folder,time.perf_counter()-start,rows
    def verify(prep):
        for f,h in prep['payload_sha256'].items():
            if digest(a.output/prep['folder']/f)!=h:raise RuntimeError('seed changed')
    if a.seed_source and not data['preparations']:
        origin=a.seed_source.resolve();old=json.loads((origin/'results.json').read_text())
        if old['pseudo_sha256']!=data['pseudo_sha256']:raise RuntimeError('imported GS pseudopotential mismatch')
        if old['binary_sha256']!=binary_hash and not a.allow_seed_binary_change:raise RuntimeError('imported GS executable differs; explicit override required')
        for prep in old['preparations']:
            if prep['shape'] not in [list(s) for s in shapes]:continue
            folder=origin/prep['folder']
            if (folder/'inputfile').read_text()!=gs_input(prep['shape']):raise RuntimeError('imported GS input mismatch')
            for name,sha in prep['payload_sha256'].items():
                if digest(folder/name)!=sha:raise RuntimeError('imported GS hash mismatch')
            (a.output/prep['folder']).symlink_to(folder,target_is_directory=True)
            prep.setdefault('binary_sha256',old['binary_sha256'])
            data['preparations'].append(prep)
        data['imported_seed_source']=str(origin);data['imported_seed_results_sha256']=digest(origin/'results.json')
    if a.rt_source and not data['runs']:
        origin=a.rt_source.resolve();old=json.loads((origin/'results.json').read_text())
        if old['binary_sha256']!=binary_hash or old['pseudo_sha256']!=data['pseudo_sha256']:raise RuntimeError('imported RT executable/pseudopotential mismatch')
        if old['launcher']!=data['launcher'] or old['mpi_version']!=data['mpi_version'] or old['conditions']!=data['conditions']:raise RuntimeError('imported RT conditions mismatch')
        for r in old['runs']:
            if a.full_source and r['mode']=='full':continue
            c=next((c for c in cases if c['shape']==r['shape'] and c['ranks']==r['ranks']),None)
            if c is None or r['repeat']>c['repeats']:continue
            prep=next(v for v in data['preparations'] if v['shape']==r['shape'])
            oldprep=next(v for v in old['preparations'] if v['shape']==r['shape'])
            if prep['payload_sha256']!=oldprep['payload_sha256']:raise RuntimeError('imported RT seed differs')
            folder=origin/r['folder'];fraction=1. if r['mode']=='full' else .999
            if (folder/'inputfile').read_text()!=rt_block(r['shape'],r['ranks'],fraction,a.pair_tolerance,a.ace_support):raise RuntimeError('imported RT input mismatch')
            check=parse_rt(folder,r['ranks'],max_energy_width=None)
            for key in ('rt_max_seconds','peak_rank_bytes','observables','energies'):
                if check[key]!=r[key]:raise RuntimeError('imported RT payload differs')
            (a.output/r['folder']).symlink_to(folder,target_is_directory=True)
            data['runs'].append(r)
        data['imported_rt_source']=str(origin);data['imported_rt_results_sha256']=digest(origin/'results.json')
    if a.full_source:
        origin=a.full_source.resolve();old=json.loads((origin/'results.json').read_text())
        if old['pseudo_sha256']!=data['pseudo_sha256']:raise RuntimeError('full reference pseudopotential mismatch')
        if old['launcher']!=data['launcher'] or old['mpi_version']!=data['mpi_version']:raise RuntimeError('full reference MPI mismatch')
        if old['conditions']['threads']!=data['conditions']['threads']:raise RuntimeError('full reference thread mismatch')
        provenance_hash=digest(origin/'results.json')
        if data.get('full_reference_results_sha256',provenance_hash)!=provenance_hash:raise RuntimeError('full reference source changed')
        for r in old['runs']:
            if r['mode']!='full':continue
            c=next((c for c in cases if c['shape']==r['shape'] and c['ranks']==r['ranks']),None)
            if c is None or r['repeat']>c['repeats']:continue
            if any(v['folder']==r['folder'] for v in data['runs']):continue
            prep=next(v for v in data['preparations'] if v['shape']==r['shape'])
            oldprep=next(v for v in old['preparations'] if v['shape']==r['shape'])
            if prep['payload_sha256']!=oldprep['payload_sha256']:raise RuntimeError('full reference seed differs')
            folder=origin/r['folder']
            if (folder/'inputfile').read_text()!=rt_block(r['shape'],r['ranks'],1.,a.pair_tolerance,a.ace_support):raise RuntimeError('full reference input mismatch')
            check=parse_rt(folder,r['ranks'],max_energy_width=None)
            for key in ('rt_max_seconds','peak_rank_bytes','observables','energies'):
                if check[key]!=r[key]:raise RuntimeError('full reference payload differs')
            (a.output/r['folder']).symlink_to(folder,target_is_directory=True)
            record=reference_record(r,old,origin/'results.json',provenance_hash)
            record.update(c)
            data['runs'].append(record)
        data['full_reference_source']=str(origin);data['full_reference_results_sha256']=provenance_hash
        data['full_reference_note']='User-authorized historical full references; binaries and measurement times differ. Only missing full cases are newly measured.'
    save()
    seen_shapes=set()
    for shape in shapes+([(4,4,1)] if not a.pilot and not a.prepare_only else []):
        strong_phase=tuple(shape) in seen_shapes;seen_shapes.add(tuple(shape))
        tag='x'.join(map(str,shape));name='gs-'+tag;ranks=math.prod(shape)
        prep=next((r for r in data['preparations'] if r['folder']==name),None)
        if prep:verify(prep)
        else:
            folder,wall,rows=execute(name,gs_input(shape),ranks)
            text=(folder/'output').read_text()
            scf=re.findall(r'DC #SCF =\s+(\d+)\s+Total Energy =\s+(\S+)\s+diff =\s+(\S+)',text)
            if 'end SALMON' not in text or 'end DC-LCFO complex' not in text or not scf or not math.isfinite(float(scf[-1][2])) or float(scf[-1][2])>=1e-10:raise RuntimeError('GS not converged '+name)
            charge=float(re.findall(r'integral\(rho_tot\)=\s*(\S+)',text)[-1])
            if abs(charge-16*ranks)>1e-7:raise RuntimeError('electron count mismatch')
            hashes={str(f.relative_to(folder)):digest(f) for f in (folder/'data_dcdft').rglob('*') if f.is_file()}
            if not any(k.endswith('wavefunctions.bin') for k in hashes):raise RuntimeError('missing payload')
            prep=dict(shape=list(shape),ranks=ranks,folder=name,iterations=int(scf[-1][0]),energy_ev=float(scf[-1][1]),residual=float(scf[-1][2]),
              charge=charge,binary_sha256=binary_hash,wall_seconds=wall,per_rank=rows,peak_rank_bytes=max(r['peak_rss_bytes'] for r in rows),payload_sha256=hashes,
              timers=[l for l in text.splitlines() if re.match(r'\s*(scf iterations|total calculation time|DC|lcfo)',l)])
            data['preparations'].append(prep);save();print('PREPARED',name,wall,'seconds',flush=True)
        if a.prepare_only:continue
        for case in (c for c in cases if c['shape']==list(shape) and ('strong' if strong_phase else 'weak') in c['suites']):
            ranks=case['ranks']
            for rep in range(1,case['repeats']+1):
                for mode,fraction in ([('full',1.),('adaptive',.999)] if rep%2 else [('adaptive',.999),('full',1.)]):
                    name=f'{mode}-{tag}-mpi{ranks}-r{rep}'
                    if any(r['folder']==name for r in data['runs']):continue
                    folder,wall,_=execute(name,rt_block(shape,ranks,fraction,a.pair_tolerance,a.ace_support),ranks,a.output/prep['folder'])
                    r=parse_rt(folder,ranks,max_energy_width=None);text=(folder/'output').read_text()
                    r.update(case,folder=name,mode=mode,repeat=rep,launcher_wall_seconds=wall,binary_sha256=binary_hash,production_commit=data['production_commit'],
                      pair_counters=[list(map(int,m)) for m in re.findall(r'EXX_ADAPTIVE local/global pairs/local FFT points \(orbital group 0\):\s+(\d+)\s+(\d+)\s+(\d+)',text)],
                      support_counters=[list(map(float,m)) for m in re.findall(r'EXX_ADAPTIVE fraction/max radius/max norm loss:\s+(\S+)\s+(\S+)\s+(\S+)',text)],
                      retained_gauge_updates=text.count('retained accepted transported gauge'),
                      generated_pair_counters=[list(map(int,m)) for m in re.findall(r'EXX_PAIR generated grid products/catalogue entries:\s+(\d+)\s+(\d+)',text)],
                      pair_product_points=[int(m) for m in re.findall(r'EXX_PAIR evaluated product points:\s+(\d+)',text)],
                      source_ace_accepted=text.count('EXX_SUPPORT_ACE accepted: T'),source_ace_fallbacks=text.count('EXX_SUPPORT_ACE accepted: F'),
                      source_ace_counts=[list(map(int,m)) for m in re.findall(r'EXX_SUPPORT_ACE products/skipped/catalogue/product points:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)',text)],
                      pair_fallbacks=text.count('EXX_PAIR unscreened ACE fallback: T'),
                      accepted_pair_bounds=[float(m) for m in re.findall(r'EXX_PAIR unscreened ACE fallback: [TF] accepted action bound:\s+(\S+)',text)])
                    if mode=='adaptive' and a.ace_support=='source':
                        if r['source_ace_accepted']+r['source_ace_fallbacks']!=33 or len(r['source_ace_counts'])!=33:raise RuntimeError('missing source ACE diagnostics')
                    if mode=='adaptive' and a.pair_tolerance is not None:
                        if len(r['accepted_pair_bounds'])!=33 or len(r['generated_pair_counters'])!=33:raise RuntimeError('missing pair diagnostics')
                        if any(not math.isfinite(v) or v<0 or v>a.pair_tolerance*(1+1e-10) for v in r['accepted_pair_bounds']):raise RuntimeError('invalid accepted pair bound')
                    data['runs'].append(r);save();print('DONE',name,r['rt_max_seconds'],'seconds',r['peak_rank_bytes']/2**20,'MiB',flush=True)
        verify(prep)
    data['preparation_complete']=True;data['complete']=not a.prepare_only;save()

if __name__=='__main__':main()
