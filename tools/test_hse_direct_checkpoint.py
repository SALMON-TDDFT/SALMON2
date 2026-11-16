"""Direct fresh-process fast-PR HSE context restart, no elapsed propagation replay.

BUILD_JSON FIXTURES NEW_OUTPUT [--threads 1]. Uses only audited fast-PR HSE
objects, never a conventional HSE or comparison baseline execution engine.
"""
import argparse,hashlib,json,os,re,shlex,shutil,subprocess,sys
from pathlib import Path
import numpy as np


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build','fixtures','output'):p.add_argument(name,type=Path)
    p.add_argument('--threads',type=int,default=1)
    args=p.parse_args();root=Path(__file__).resolve().parents[1]
    record=json.loads(args.build.read_text());exe=Path(record['executable']);build=exe.parent
    source=Path(record['source_root']);marker=json.loads((source/'salmonpy-hse-rt-state.json').read_text())
    if (record.get('schema')!='salmonpy.hse-rt-state-build/1' or record.get('status')!='built-and-audited'
        or digest(exe)!=record['sha256'] or not marker.get('explicit_rt_routing')):
        raise ValueError('audited explicit fast-PR HSE build required')
    if args.threads<1:raise ValueError('positive threads required')
    out=args.output.resolve()
    if any(out.is_relative_to(p.resolve()) for p in (root,build,source,args.fixtures)):
        raise ValueError('new external output required')
    out.mkdir(parents=True,exist_ok=False);src=out/'src';src.mkdir()
    config=(build/'config.h').read_text()
    if not re.search(r'^#define USE_HSE\b',config,re.M) or re.search(r'^#define (USE_MPI|USE_OPENACC|USE_FFTW)\b',config,re.M):
        raise ValueError('serial CPU fast-PR HSE/FFTE build required')
    adapters=root/'salmonpy/native_adapters'
    if not adapters.exists():adapters=root.parent/'native_adapters'
    for original,name in [(adapters/'hse_checkpoint.f90','hse_checkpoint'),
                          (adapters/'tddft_context.f90','context'),
                          (root/'tests/native/test_hse_direct_checkpoint.f90','probe')]:
        shutil.copy2(original,src/(name+'.f90'))
    flags=(build/'src/CMakeFiles/salmon.dir/flags.make').read_text()
    options=[v for v in shlex.split(re.search(r'^Fortran_FLAGS\s*=\s*(.*)$',flags,re.M)[1]) if not v.startswith('-J')]
    includes=shlex.split(re.search(r'^Fortran_INCLUDES\s*=\s*(.*)$',flags,re.M)[1])
    fc=re.search(r'^CMAKE_Fortran_COMPILER:[^=]*=(.*)$',(build/'CMakeCache.txt').read_text(),re.M)[1].strip()
    checker=root/'tools/check_forbidden_eigensolvers.py';commands=[]
    with (out/'build.log').open('w') as log:
        def run(cmd,**kw):
            commands.append([str(x) for x in cmd]);subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,**kw)
        run([sys.executable,checker,'--source-root',source/'src'])
        run([sys.executable,checker,'--source-root',src])
        for name in ('hse_checkpoint','context','probe'):
            run([fc,*options,'-I'+str(out),*includes,'-J'+str(out),'-c',src/(name+'.f90'),'-o',out/(name+'.o')],cwd=out)
        run([sys.executable,checker,'--object-tree',build]);run([sys.executable,checker,'--object-tree',out])
        link=shlex.split((build/'src/CMakeFiles/salmon.dir/link.txt').read_text())
        link[link.index('-o')+1]=str(out/'probe')
        link=[str(out/'probe.o') if v=='CMakeFiles/salmon.dir/main.f90.o' else v for v in link]
        link[1:1]=[str(out/'context.o'),str(out/'hse_checkpoint.o')]
        run(link,cwd=build/'src')
    from salmonpy.bands import _input_groups
    results=[];rejections=[]
    env={**os.environ,'OMP_NUM_THREADS':str(args.threads),'OPENBLAS_NUM_THREADS':'1'}
    for nk,mesh in [(2,(2,1,1)),(3,(3,1,1))]:
        fixture=args.fixtures/f'k{nk}'/'reimported'
        for prop,factor in [('hse_taylor4','off'),('hse_taylor4','learned'),
                            ('hse_taylor4_full','off'),('hse_ptcn','learned')]:
            case=out/f'k{nk}-{prop}-{factor}';case.mkdir()
            groups=_input_groups((fixture/'rt.inp').read_text())
            groups['calculation']="\n theory='tddft_response'\n"
            groups['functional']=f"\n xc='hse06', exx_factor_mode='{factor}', exx_factor_exact_interval={1 if prop=='hse_ptcn' else 2}, exx_factor_n_teach=1\n"
            groups['kgrid']='\n num_kgrid='+','.join(map(str,mesh))+'\n'
            pc='n' if prop=='hse_ptcn' else 'y'
            groups['propagation']=f"\n propagator='{prop}', yn_predictor_corrector='{pc}', n_hamil=4\n"
            groups['tgrid']='\n dt=.001, nt=12\n'
            groups['emfield']="\n ae_shape1='impulse',e_impulse=0,trans_longi='tr'\n"
            groups['analysis']='\n out_rt_energy_step=1\n'
            groups['parallel']='\n nproc_k=1,nproc_ob=1,nproc_rgrid=1,1,1\n'
            text='\n'.join('&'+k+v+'/' for k,v in groups.items())+'\n'
            def execute(label,mode,split,path,input_text=None):
                d=case/label;d.mkdir();(d/'data_for_restart').mkdir()
                for n in ('info.bin','occupation.bin','wfn.bin'):shutil.copy2(fixture/'data_for_restart'/n,d/'data_for_restart'/n)
                shutil.copy2(fixture/'pseudo1_rps.dat',d/'pseudo1_rps.dat');(d/'inputfile').write_text(text if input_text is None else input_text)
                with (d/'run.log').open('w') as log,(d/'inputfile').open() as inp:
                    r=subprocess.run([str(out/'probe'),mode,str(split),str(path)],cwd=d,stdin=inp,
                        stdout=log,stderr=subprocess.STDOUT,env=env,timeout=180)
                if r.returncode:raise RuntimeError(f'{d}: exit {r.returncode}; see run.log')
                return d
            baseline=execute('continuous','continuous',0,case/'unused')
            for split in (6,7):
                archive=case/f'checkpoint-{split}.bin'
                saved=execute(f'save-{split}','save',split,archive)
                sha=digest(archive)
                receipt={'schema':'salmonpy.fast-hse-checkpoint-receipt/1',
                    'checkpoint_sha256':sha,'runtime_sha256':record['sha256'],
                    'probe_sha256':digest(out/'probe'),'input_sha256':hashlib.sha256(text.encode()).hexdigest(),
                    'initial_files':{n:digest(fixture/'data_for_restart'/n) for n in ('info.bin','occupation.bin','wfn.bin')},
                    'pseudo_sha256':digest(fixture/'pseudo1_rps.dat')}
                archive.with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n')
                assert digest(archive)==receipt['checkpoint_sha256']
                assert digest(exe)==receipt['runtime_sha256'] and digest(out/'probe')==receipt['probe_sha256']
                assert all(digest(fixture/'data_for_restart'/n)==sha for n,sha in receipt['initial_files'].items())
                assert digest(fixture/'pseudo1_rps.dat')==receipt['pseudo_sha256']
                assert hashlib.sha256(text.encode()).hexdigest()==receipt['input_sha256']
                resumed=execute(f'resume-{split}','resume',split,archive)
                assert digest(archive)==sha
                calls=np.loadtxt(resumed/'calls.dat',dtype=int)
                np.testing.assert_array_equal(calls,[split,12-split,12])
                before=np.loadtxt(baseline/'api.dat');after=np.vstack((np.loadtxt(saved/'api.dat'),np.loadtxt(resumed/'api.dat')))
                assert np.isfinite(before).all() and np.isfinite(after).all()
                np.testing.assert_array_equal(before,after)
                assert (baseline/'trace.bin').read_bytes()==(saved/'trace.bin').read_bytes()+(resumed/'trace.bin').read_bytes(), 'internal state trace mismatch'
                results.append({'nk':nk,'propagator':prop,'factor':factor,'split':split,
                    'resumed_calls':int(calls[1]),'api_max_error':float(abs(before-after).max()),
                    'internal_trace_equal':True,'checkpoint_sha256':sha,'checkpoint_bytes':archive.stat().st_size})
            broken=case/'truncated.bin'
            content=(case/'checkpoint-7.bin').read_bytes();broken.write_bytes(content[:len(content)//2])
            trailing=case/'trailing.bin';trailing.write_bytes(content+b'x')
            for label,path,input_text in [('reject-truncated',broken,None),
                    ('reject-settings',case/'checkpoint-7.bin',text.replace('dt=.001','dt=.002')),
                    ('reject-late',trailing,None)]:
                rejected=execute(label,'reject-late' if label=='reject-late' else 'reject',7,path,input_text)
                calls=np.loadtxt(rejected/'calls.dat',dtype=int)
                assert calls[1]==0
                if label=='reject-late':assert calls[2]==-1
                rejections.append({'nk':nk,'propagator':prop,'factor':factor,'case':label,'propagation_calls':0})
    (out/'result.json').write_text(json.dumps({'status':'passed','scope':'serial fast-PR HSE direct context restart; synthetic input, no coupled MPI claim',
        'threads':args.threads,'engine_sha256':record['sha256'],'source_hashes':{p.name:digest(p) for p in src.iterdir()},
        'cases':results,'rejections':rejections,'commands':commands},indent=2)+'\n')
    print('Fast-PR HSE direct checkpoint PASS',len(results),flush=True)


if __name__=='__main__':main()
