"""Compare whole fast-PR HSE solvers on synthetic native orbitals, without SCF.

python tools/test_hse_rt_state_runtime.py NEW_BUILD_JSON REF_BUILD_JSON FIXTURES NEW_OUTPUT
Requires recorded source/object audits. Existing binaries are checked by hash.
The builtin regular k mesh is used: fast-PR HSE rejects file_kw custom meshes.
No original fixture is modified, no GS output is collected, no physical accuracy
or full native multi-material/speed conclusion is implied by this regression.
"""
import argparse,hashlib,json,os,shutil,subprocess,time
from pathlib import Path
import numpy as np


def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()



def _check_scientific_trace(reference,changed,propagator,factor_mode):
    if reference!=changed:raise ValueError('factor/solver scientific diagnostics differ')
    predictions=sum(line.startswith('FACTOR_HISTORY step/predicted:') and line.endswith('T') for line in changed)
    if propagator=='hse_taylor4' and factor_mode=='learned' and predictions==0:
        raise ValueError('learned Taylor path installed no prediction')
    if propagator=='hse_ptcn' and not any(line.startswith('HSE_PT_CN residual=') for line in changed):
        raise ValueError('PT-CN solver residual trace missing')
    return predictions


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['build','reference','fixtures','output']:p.add_argument(name,type=Path)
    p.add_argument('--threads',type=int,default=1)
    a=p.parse_args()
    if a.threads<1:raise ValueError('positive thread count required')
    a.output.mkdir(parents=True,exist_ok=False)
    builds=[json.loads(x.read_text()) for x in [a.reference,a.build]]
    for r in builds:
        if r.get('schema')!='salmonpy.hse-rt-state-build/1' or r.get('status')!='built-and-audited':raise ValueError('audited HSE build required')
        if digest(r['executable'])!=r['sha256']:raise ValueError('HSE runtime hash changed')
    from salmonpy.bands import _input_groups
    results=[]
    for nk,mesh in [(2,(2,1,1)),(3,(3,1,1)),(6,(3,2,1)),(8,(2,2,2))]:
        fixture=a.fixtures/('k'+str(nk))/'reimported'
        template=_input_groups((fixture/'rt.inp').read_text())
        for dt in [.001,.002]:
            for mode,learned in [('hse_taylor4','off'),('hse_taylor4','learned'),('hse_taylor4_full','off'),('hse_ptcn','learned')]:
                case=a.output/(f'k{nk}-dt{dt}-{mode}-{learned}');case.mkdir()
                tables=[];elapsed=[];scientific=[]
                for label,build in zip(['reference','changed'],builds):
                    work=case/label;work.mkdir();(work/'data_for_restart').mkdir()
                    for name in ['info.bin','occupation.bin','wfn.bin']:shutil.copy2(fixture/'data_for_restart'/name,work/'data_for_restart'/name)
                    shutil.copy2(fixture/'pseudo1_rps.dat',work/'pseudo1_rps.dat')
                    groups=dict(template)
                    groups['calculation']="\n theory='tddft_response'\n"
                    interval=1 if mode=='hse_ptcn' else 2
                    groups['functional']=f"\n xc='hse06', exx_factor_mode='{learned}', exx_factor_exact_interval={interval}, exx_factor_n_teach=1\n"
                    groups['kgrid']='\n num_kgrid='+','.join(str(x) for x in mesh)+'\n'
                    pc='n' if mode=='hse_ptcn' else 'y'
                    groups['propagation']=f"\n propagator='{mode}', yn_predictor_corrector='{pc}', n_hamil=4\n"
                    groups['tgrid']=f'\n dt={dt:.17e}, nt=12\n'
                    groups['emfield']="\n ae_shape1='impulse', e_impulse=.0001, epdir_re1=0,0,1, trans_longi='tr'\n"
                    groups['analysis']='\n out_rt_energy_step=1\n'
                    groups['parallel']='\n nproc_k=1, nproc_ob=1, nproc_rgrid=1,1,1\n'
                    (work/'inputfile').write_text('\n'.join('&'+key+v+'/' for key,v in groups.items())+'\n')
                    start=time.perf_counter()
                    with (work/'run.log').open('w') as log,(work/'inputfile').open() as inp:
                        run=subprocess.run([build['executable']],cwd=work,stdin=inp,stdout=log,stderr=subprocess.STDOUT,
                            env={**os.environ,'OMP_NUM_THREADS':str(a.threads),'OPENBLAS_NUM_THREADS':'1'},timeout=120)
                    elapsed.append(time.perf_counter()-start)
                    if run.returncode:raise RuntimeError(f'{work}: native exit {run.returncode}; see run.log')
                    rt=np.loadtxt(work/'rt_rt.data');energy=np.loadtxt(work/'rt_rt_energy.data')
                    if rt.shape!=(12,16) or energy.shape!=(13,3):raise ValueError('unexpected native sample count')
                    if not np.isfinite(rt).all() or not np.isfinite(energy).all():raise ValueError('nonfinite native result')
                    np.testing.assert_allclose(rt[:,0],dt*np.arange(1,13),rtol=0,atol=1e-10)
                    tables.append((rt,energy))
                    lines=(work/'run.log').read_text().splitlines()
                    scientific.append([line.strip() for line in lines if line.strip().startswith(('FACTOR_HISTORY step/predicted:', 'FACTOR_HISTORY k/teacher/step/count/coeff/ready:', 'HSE_PT_CN residual=', 'HSE_PT_CN_ORTH pre_gram='))])
                for one,two in zip(tables[0],tables[1]):np.testing.assert_allclose(one,two,rtol=1e-12,atol=1e-14)
                predictions=_check_scientific_trace(scientific[0],scientific[1],mode,learned)
                result=dict(nk=nk,mesh=mesh,dt=dt,propagator=mode,factor_mode=learned,steps=12,threads=a.threads,
                    current_max_abs_difference=float(np.max(abs(tables[0][0][:,13:16]-tables[1][0][:,13:16]))),
                    energy_max_abs_difference=float(np.max(abs(tables[0][1][:,1:]-tables[1][1][:,1:]))),
                    identical_parsed_current=tables[0][0].tobytes()==tables[1][0].tobytes(),identical_parsed_energy=tables[0][1].tobytes()==tables[1][1].tobytes(),
                    identical_current_output=(case/'reference/rt_rt.data').read_bytes()==(case/'changed/rt_rt.data').read_bytes(),
                    identical_energy_output=(case/'reference/rt_rt_energy.data').read_bytes()==(case/'changed/rt_rt_energy.data').read_bytes(),
                    identical_factor_solver_diagnostics=True,
                    installed_taylor_predictions=predictions,
                    end_to_end_seconds=dict(zip(['reference','changed'],elapsed)))
                results.append(result);print(nk,dt,mode,learned,result['current_max_abs_difference'],flush=True)
                # Keep logs/input/scalar tables; remove only the owned synthetic
                # native initial copies to avoid collecting a GS payload.
                for label in ['reference','changed']:
                    shutil.rmtree(case/label/'data_for_restart')
                    (case/label/'pseudo1_rps.dat').unlink()
    (a.output/'result.json').write_text(json.dumps(dict(scope='synthetic fast-PR HSE formatted-output numerical regression; not internal-bit, converged GS, full contexts or speed benchmark',
        reference_sha256=builds[0]['sha256'],changed_sha256=builds[1]['sha256'],cases=results),indent=2)+'\n')


if __name__=='__main__':main()
