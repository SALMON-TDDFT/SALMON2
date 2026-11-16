"""Local full-step two-context regression against isolated native processes.

BUILD_RECORD FIXTURES NEW_OUTPUT [--threads N]. No CI or network is used.
Only synthetic restart fixtures; outputs and all linked objects are audited.
"""
import argparse, hashlib, json, os, re, shlex, shutil, subprocess, sys
from pathlib import Path
import numpy as np


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['build','fixtures','output']:p.add_argument(n,type=Path)
    p.add_argument('--different-capacity',action='store_true',help='A accepts 12 steps; B accepts 6')
    p.add_argument('--different-openmp',action='store_true',help='A uses orbital, B uses real-space OpenMP selection')
    p.add_argument('--different-dt',action='store_true',help='A and B use independent fixed timesteps; poison legacy global dt')
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--reference',type=Path,required=True,help='immutable native fast-HSE full-driver build record')
    a=p.parse_args();root=Path(__file__).resolve().parents[1]
    record=json.loads(a.build.read_text());exe=Path(record['executable']);build=exe.parent
    if record.get('schema')!='salmonpy.hse-rt-state-build/1' or record.get('status')!='built-and-audited' or digest(exe)!=record['sha256']:raise ValueError('unchanged audited build required')
    reference=json.loads(a.reference.read_text())
    if reference.get('schema')!='salmonpy.hse-rt-state-build/1' or reference.get('status')!='built-and-audited' or digest(reference['executable'])!=reference['sha256']:raise ValueError('unchanged audited reference required')
    if a.threads<1:raise ValueError('positive threads required')
    out=a.output.resolve()
    if any(out.is_relative_to(p.resolve()) for p in [root,build,Path(record['source_root']),a.fixtures]):
        raise ValueError('new external output directory required')
    out.mkdir(parents=True,exist_ok=False);src=out/'src';src.mkdir()
    checker=Path(__file__).with_name('check_forbidden_eigensolvers.py')
    adapter=root/'salmonpy/native_adapters/tddft_context.f90'
    if not adapter.exists():adapter=root.parent/'native_adapters/tddft_context.f90'
    shutil.copy2(adapter,src/'context.f90');shutil.copy2(adapter.with_name('hse_checkpoint.f90'),src/'hse_checkpoint.f90');shutil.copy2(root/'tests/native/test_tddft_contexts.f90',src/'probe.f90')
    config=(build/'config.h').read_text()
    if not re.search(r'^#define USE_HSE\b',config,re.M) or re.search(r'^#define (USE_MPI|USE_OPENACC|USE_FFTW)\b',config,re.M):
        raise ValueError('serial CPU HSE/FFTE build required')
    flags=(build/'src/CMakeFiles/salmon.dir/flags.make').read_text()
    options=[v for v in shlex.split(re.search(r'^Fortran_FLAGS\s*=\s*(.*)$',flags,re.M)[1]) if not v.startswith('-J')]
    includes=shlex.split(re.search(r'^Fortran_INCLUDES\s*=\s*(.*)$',flags,re.M)[1])
    fc=re.search(r'^CMAKE_Fortran_COMPILER:[^=]*=(.*)$',(build/'CMakeCache.txt').read_text(),re.M)[1].strip()
    commands=[]
    with (out/'build.log').open('w') as log:
        def run(cmd,**kw):
            commands.append([str(v) for v in cmd]);subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,**kw)
        run([sys.executable,checker,'--source-root',record['source_root']+'/src'])
        run([sys.executable,checker,'--source-root',src])
        for name in ['hse_checkpoint','context','probe']:
            run([fc,*options,'-I'+str(out),*includes,'-J'+str(out),'-c',src/(name+'.f90'),'-o',out/(name+'.o')],cwd=out)
        run([sys.executable,checker,'--object-tree',build]);run([sys.executable,checker,'--object-tree',out])
        link=shlex.split((build/'src/CMakeFiles/salmon.dir/link.txt').read_text());link[link.index('-o')+1]=str(out/'probe')
        link=[str(out/'probe.o') if v=='CMakeFiles/salmon.dir/main.f90.o' else v for v in link];link.insert(1,str(out/'context.o'));link.insert(1,str(out/'hse_checkpoint.o'))
        run(link,cwd=build/'src')
    from salmonpy.bands import _input_groups
    results=[]
    for nk,mesh in [(2,(2,1,1)),(3,(3,1,1))]:
        fixture=a.fixtures/f'k{nk}'/'reimported'
        for propagator,learned in [('hse_taylor4','off'),('hse_taylor4','learned'),('hse_taylor4_full','off'),('hse_ptcn','learned')]:
            case=out/f'k{nk}-{propagator}-{learned}';case.mkdir()
            for mode in ['solo_a','solo_b','interleaved','reverse','reopen_a','reopen_b','native']:
                work=case/mode;work.mkdir();(work/'data_for_restart').mkdir()
                for name in ['info.bin','occupation.bin','wfn.bin']:shutil.copy2(fixture/'data_for_restart'/name,work/'data_for_restart'/name)
                shutil.copy2(fixture/'pseudo1_rps.dat',work/'pseudo1_rps.dat')
                groups=_input_groups((fixture/'rt.inp').read_text())
                groups['calculation']="\n theory='tddft_response'\n"
                groups['functional']=f"\n xc='hse06', exx_factor_mode='{learned}', exx_factor_exact_interval={1 if propagator=='hse_ptcn' else 2}, exx_factor_n_teach=1\n"
                groups['kgrid']='\n num_kgrid='+','.join(map(str,mesh))+'\n'
                groups['propagation']=f"\n propagator='{propagator}', yn_predictor_corrector='{"n" if propagator=='hse_ptcn' else "y"}', n_hamil=4\n"
                groups['tgrid']='\n dt='+('.0015' if a.different_dt and mode=='solo_b' else '.001')+', nt='+('6' if a.different_capacity and mode=='solo_b' else '12')+'\n'
                groups['emfield']="\n ae_shape1='impulse', e_impulse=0, trans_longi='tr'\n"
                groups['analysis']='\n out_rt_energy_step=1\n'
                groups['parallel']='\n nproc_k=1, nproc_ob=1, nproc_rgrid=1,1,1\n'
                (work/'inputfile').write_text('\n'.join('&'+k+v+'/' for k,v in groups.items())+'\n')
                with (work/'run.log').open('w') as log,(work/'inputfile').open() as inp:
                    r=subprocess.run([out/'probe',mode,'different_dt' if a.different_dt else 'same_dt','different_omp' if a.different_openmp else 'same_omp','different_cap' if a.different_capacity else 'same_cap'],cwd=work,stdin=inp,stdout=log,stderr=subprocess.STDOUT,timeout=120,
                                     env={**os.environ,'OMP_NUM_THREADS':str(a.threads),'OPENBLAS_NUM_THREADS':'1'})
                if r.returncode:raise RuntimeError(f'{work}: native exit {r.returncode}; see run.log')
                shutil.rmtree(work/'data_for_restart');(work/'pseudo1_rps.dat').unlink()
            for label in ['a','b']:
                expected_dt=.0015 if a.different_dt and label=='b' else .001
                api=np.loadtxt(case/('solo_'+label)/(label+'_api.dat'))
                np.testing.assert_allclose(api[:,0],np.arange(1,7 if a.different_capacity and label=="b" else 13)*expected_dt,rtol=0,atol=1e-17)
                for mode in ['interleaved','reverse']:
                    solo=case/('solo_'+label);both=case/mode
                    for suffix in ['_rt.data','_rt_energy.data','_api.dat']:
                        x=np.loadtxt(solo/(label+suffix));y=np.loadtxt(both/(label+suffix))
                        if not np.isfinite(x).all() or x.tobytes()!=y.tobytes():raise ValueError('native output differs: '+label+suffix)
                    if (solo/(label+'_trace.bin')).read_bytes()!=(both/(label+'_trace.bin')).read_bytes():raise ValueError('full-state trace differs: '+label)
            for reopened,unaffected in [('a','b'),('b','a')]:
                work=case/('reopen_'+reopened)
                solo=case/('solo_'+unaffected)
                if (work/(unaffected+'_trace.bin')).read_bytes()!=(solo/(unaffected+'_trace.bin')).read_bytes():
                    raise ValueError('close/reopen disturbed surviving state')
                fresh=(work/(reopened+'_reopened_trace.bin')).read_bytes()
                expected=(case/('solo_'+reopened)/(reopened+'_trace.bin')).read_bytes()
                if not fresh or expected[:len(fresh)]!=fresh:raise ValueError('reopened state differs from fresh process')
            # The immutable reference rejects sampled-input fields. Use its supported
            # Acos2 pulse for a separate full-driver comparison, without modifying it.
            work=case/'native_reference';work.mkdir();(work/'data_for_restart').mkdir()
            for name in ['info.bin','occupation.bin','wfn.bin']:shutil.copy2(fixture/'data_for_restart'/name,work/'data_for_restart'/name)
            shutil.copy2(fixture/'pseudo1_rps.dat',work/'pseudo1_rps.dat')
            native_groups=dict(groups)
            native_groups['emfield']="\n ae_shape1='Acos2',E_amplitude1=.01,omega1=20,tw1=.024,epdir_re1=0,0,1,trans_longi='tr'\n"
            (work/'inputfile').write_text('\n'.join('&'+k+v+'/' for k,v in native_groups.items())+'\n')
            with (work/'run.log').open('w') as log,(work/'inputfile').open() as inp:
                r=subprocess.run([reference['executable']],cwd=work,stdin=inp,stdout=log,stderr=subprocess.STDOUT,timeout=120,
                                 env={**os.environ,'OMP_NUM_THREADS':str(a.threads),'OPENBLAS_NUM_THREADS':'1','SALMONPY_HSE_CONTEXT':'0'})
            if r.returncode:raise RuntimeError(f'{work}: native reference exit {r.returncode}; see run.log')
            native=np.loadtxt(work/'rt_rt.data');isolated=np.loadtxt(case/'native'/'a_rt.data')
            np.testing.assert_allclose(native[:,13:16],isolated[:,13:16],rtol=1e-12,atol=1e-14)
            np.testing.assert_allclose(np.loadtxt(work/'rt_rt_energy.data'),np.loadtxt(case/'native'/'a_rt_energy.data'),rtol=1e-12,atol=1e-14)
            # The development entry returns charge current, opposite to native matter current.
            api=np.loadtxt(case/'native'/'a_api.dat')
            np.testing.assert_allclose(api[:,1:4],-native[:,13:16],rtol=1e-12,atol=1e-14)
            np.testing.assert_allclose(api[:,4],np.loadtxt(work/'rt_rt_energy.data')[1:,1],rtol=1e-12,atol=1e-14)
            shutil.rmtree(work/'data_for_restart');(work/'pseudo1_rps.dat').unlink()
            results.append(dict(nk=nk,propagator=propagator,factor=learned,steps=12,dt_a=.001,dt_b=.0015 if a.different_dt else .001,identical=True));print(results[-1],flush=True)
    (out/'result.json').write_text(json.dumps(dict(status='passed',scope='full-step synthetic serial HSE contexts; owned fixed dt and otherwise common immutable configuration',different_dt=a.different_dt,different_openmp=a.different_openmp,different_capacity=a.different_capacity,threads=a.threads,source_hashes={f.name:digest(f) for f in src.iterdir()},engine_sha256=record['sha256'],reference_sha256=reference['sha256'],probe_sha256=digest(out/'probe'),commands=commands,cases=results),indent=2)+'\n')


if __name__=='__main__':main()
