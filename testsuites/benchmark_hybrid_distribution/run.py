"""Fresh-process synthetic kernel scaling; never an end-to-end water benchmark."""
import argparse,hashlib,json,math,os,platform,shlex,statistics,subprocess,tempfile
from pathlib import Path


def run():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--grid',type=int,default=32)
    p.add_argument('--states',type=int,default=16)
    p.add_argument('--dense-dim',type=int,default=1024)
    p.add_argument('--repeat',type=int,default=3)
    p.add_argument('--timeout',type=int,default=300)
    p.add_argument('--mpiexec',default='mpiexec')
    p.add_argument('--link-flags',help='Override FFTW, ScaLAPACK and BLAS/LAPACK link flags')
    p.add_argument('--phase',choices=['all','exchange','dense'],default='all')
    p.add_argument('--smoke',action='store_true',help='Small two-layout validation, not a performance result')
    a=p.parse_args()
    if min(a.grid,a.states,a.dense_dim,a.repeat,a.timeout)<1:p.error('sizes, repetitions and timeout must be positive')
    b=a.build.resolve();root=Path(__file__).resolve().parents[2]
    cache={}
    for line in (b/'CMakeCache.txt').read_text().splitlines():
        if '=' in line and ':' in line and not line.startswith(('//','#')):
            key,value=line.split('=',1);cache[key.split(':',1)[0]]=value
    if cache.get('USE_SCALAPACK')!='ON' or cache.get('USE_HSE')!='ON':p.error('MPI/HSE/ScaLAPACK build required')
    if a.link_flags:
        libs=shlex.split(a.link_flags)
    else:
        vendor=cache.get('ScaLAPACK_VENDOR_FLAGS')
        fftw=cache.get('HSE_FFTW_LIBRARY')
        if not vendor or not fftw:p.error('provide --link-flags for this build; vendor/FFTW flags not found in cache')
        libs=shlex.split(vendor)+[fftw]
    fc=cache.get('CMAKE_Fortran_COMPILER','mpifort')
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    layouts=[(1,1),(4,1),(4,2),(4,4),(8,2),(8,4)]
    dense_ranks=[1,4,8]
    if a.smoke:
        a.grid=8;a.states=8;a.dense_dim=19;a.repeat=1
        layouts=[(1,1),(4,2)];dense_ranks=[1,4]
    if a.phase!='dense':
        if a.states<max(o for _,o in layouts):p.error('each orbital group needs a state')
        if a.states%2:p.error('pair rotations require an even state count')
        k1=int(a.states**(1/3)+0.5)
        while a.states%k1:k1-=1
        k2=math.isqrt(a.states//k1)
        while (a.states//k1)%k2:k2-=1
        kd=(k1,k2,a.states//k1//k2)
        if min(kd)<2 or max(kd)>a.grid:p.error('wave-packet factorization needs 2..grid modes per axis')
        if a.grid%2:p.error('these spatial layouts require an even grid size')
    result={'scope':'synthetic production-kernel phases; not SCF/RT or water',
            'platform':platform.platform(),'machine':platform.machine(),
            'production_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
            'compiler':subprocess.check_output([fc,'--version'],text=True).splitlines()[0],
            'build':str(b),'compiler_flags':cache.get('CMAKE_Fortran_FLAGS'),
            'build_type':cache.get('CMAKE_BUILD_TYPE'),'link_flags':libs,
            'build_type_flags':cache.get('CMAKE_Fortran_FLAGS_'+cache.get('CMAKE_BUILD_TYPE','').upper()),
            'probe_flags':['-O2','-fopenmp','-ffree-line-length-none'],
            'exchange_workload':{'source':'tensor Fourier wave packets with pair rotations, not water orbitals','h_bohr':1.,
                'occupations':2.,'mlwf_maxiter':3,'mlwf_tolerance':1e-7,
                'checksum':'unscaled Tr(Psi^dagger K Psi); no hybrid mixing, no rVV10'},
            'phase_labels':{'exchange':['refresh','exact_exchange','ace_build','ace_apply'],
                'dense':['replicated_block_input_and_layout','full_solve_and_diagnostics','unused','unused']},
            'threads':{k:env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
            'memory_note':'Per-process lifetime high-water RSS. Sum of rank peaks is not simultaneous node RSS. Baseline is after MPI init; deltas are not isolated allocation accounting.',
            'timing_note':'Cold first-call phases, including FFT plan initialization and validation collectives. Reported phases use max rank time; fresh MPI processes for every repeat.',
            'runs':[],'summary':[]}
    obj=b/'src/CMakeFiles/salmon.dir'
    names=['xc/exx_pair_candidates.f90.o','xc/exx_spatial.f90.o','xc/exx_spatial_local.f90.o','xc/exx_local_fft.f90.o','xc/exx_orbitals.f90.o','xc/exx_ace.f90.o',
           'xc/exx_wannier_gauge.f90.o','xc/fftw_pencils.f90.o','gs/dc/lcfo_scalapack.f90.o',
           'parallel/communication.f90.o','misc/nvtx_wrapper.f90.o']
    result['object_sha256']={name:hashlib.sha256((obj/name).read_bytes()).hexdigest() for name in names}
    result['probe_sha256']={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in
        (Path(__file__).resolve(),Path(__file__).with_name('probe.f90').resolve(),root/'testsuites/unit_lcfo_rt/peak_rss.c')}
    result['probe_source_snapshot']={path:(root/path).read_text() for path in result['probe_sha256']}
    result['mpi_version']=subprocess.run(shlex.split(a.mpiexec)+['--version'],text=True,capture_output=True).stdout.strip()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    def save():a.output.write_text(json.dumps(result,indent=2)+'\n')
    with tempfile.TemporaryDirectory(prefix='hybrid-scaling-') as tmp:
        tmp=Path(tmp);exe=tmp/'probe'
        subprocess.run(['cc','-O2','-c',str(root/'testsuites/unit_lcfo_rt/peak_rss.c'),'-o',str(tmp/'rss.o')],check=True)
        subprocess.run([fc,'-O2','-fopenmp','-ffree-line-length-none','-I'+str(b),str(Path(__file__).with_name('probe.f90')),
                        str(tmp/'rss.o'),*[str(obj/name) for name in names],*libs,'-o',str(exe)],check=True)
        tasks=[]
        for omega in (.11,0.):
            for ranks,orbitals in layouts:
                tasks.append(('exchange',ranks,orbitals,omega,[1,a.grid,a.states,orbitals,omega]))
        for ranks in dense_ranks:tasks.append(('dense',ranks,0,None,[2,a.dense_dim]))
        tasks=[task for task in tasks if a.phase=='all' or task[0]==a.phase]
        references={}
        for kind,ranks,orbitals,omega,args in tasks:
            subset=[]
            for repeat in range(a.repeat):
                cmd=shlex.split(a.mpiexec)+['-n',str(ranks),str(exe),*map(str,args)]
                proc=subprocess.run(cmd,cwd=tmp,env=env,text=True,capture_output=True,timeout=a.timeout)
                if proc.returncode:raise RuntimeError(proc.stdout+proc.stderr)
                localizations={}
                for line in proc.stdout.splitlines():
                    if line.startswith('LOCALIZATION '):
                        f=line.split()
                        localizations[int(f[1])]=dict(iterations=int(f[2]),status=int(f[3]),spread=float(f[4]),gradient=float(f[5]))
                if kind=='exchange' and sorted(localizations)!=list(range(ranks)):
                    raise RuntimeError('missing localization diagnostics')
                rows=[]
                for line in proc.stdout.splitlines():
                    if not line.startswith('MEASURE '):continue
                    f=line.split()
                    if len(f)!=10:raise RuntimeError('malformed measurement: '+line)
                    rows.append(dict(rank=int(f[1]),baseline_bytes=int(f[2]),peak_bytes=int(f[3]),
                                     phase_seconds=list(map(float,f[4:8])),checksum=float(f[8]),error=float(f[9])))
                if sorted(r['rank'] for r in rows)!=list(range(ranks)):raise RuntimeError('missing or duplicate rank reports')
                for row in rows:
                    if kind=='exchange':
                        row['localization']=localizations[row['rank']]
                        loc=row['localization']
                        if (loc['status']!=0 and loc['iterations']!=3) or not all(math.isfinite(loc[x]) for x in ('spread','gradient')):
                            raise RuntimeError('localization failed: '+str(loc))
                    if row['baseline_bytes']<=0 or row['peak_bytes']<row['baseline_bytes']:
                        raise RuntimeError('invalid RSS measurement')
                    if any(not math.isfinite(x) for x in row['phase_seconds']+[row['checksum'],row['error']]):
                        raise RuntimeError('nonfinite measurement')
                    if min(row['phase_seconds'])<0:raise RuntimeError('negative duration')
                key=(kind,omega)
                reference=references.setdefault(key,rows[0]['checksum'])
                if any(abs(r['checksum']-reference)>1e-9*max(1.,abs(reference)) or r['error']>1e-10 for r in rows):
                    raise RuntimeError('decomposition parity failed')
                record=dict(kind=kind,ranks=ranks,orbitals=orbitals,omega=omega,repeat=repeat,
                            grid=a.grid,states=a.states,dense_dim=a.dense_dim,per_rank=rows)
                result['runs'].append(record);subset.append(record);save()
            peak=[max(r['peak_bytes'] for r in x['per_rank']) for x in subset]
            phases=[[max(r['phase_seconds'][i] for r in x['per_rank']) for x in subset] for i in range(4)]
            summary=dict(kind=kind,ranks=ranks,orbitals=orbitals,omega=omega,
                         peak_rank_bytes_median=statistics.median(peak),peak_rank_bytes_min=min(peak),peak_rank_bytes_max=max(peak),
                         sum_rank_peaks_bytes_median=statistics.median(sum(r['peak_bytes'] for r in x['per_rank']) for x in subset),
                         phase_seconds_median=[statistics.median(x) for x in phases],
                         phase_seconds_min=[min(x) for x in phases],phase_seconds_max=[max(x) for x in phases])
            result['summary'].append(summary);save()
            print(json.dumps(summary),flush=True)
    result['complete']=True;save()


if __name__=='__main__':run()
