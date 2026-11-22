"""Verify launched symmetry/shared-space Si diagnostics without GS collection."""
from pathlib import Path
import argparse,hashlib,json,re
import numpy as np
p=argparse.ArgumentParser();p.add_argument('work',type=Path);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();results={}
for name in ['off384','on384','on256']:
    d=a.work/name;text=(d/'run.log').read_text()
    assert 'end SALMON' in text and 'HSEX_INTERNAL_REDUCED_DIAGNOSTIC_PASS no W publication' in text
    assert not (d/'hsex-w.bin').exists()
    assert re.search(r'number of threads\s*:\s*12',text)
    def values(label):
        return next(line[len(label):].split() for line in text.splitlines() if line.startswith(label))
    sizes=list(map(int,values('HSEX_REDUCED_SIZE grid/RHS/maxrank ')))
    assert sizes[:2]==[1728,1560]
    counts=list(map(int,values('HSEX_REDUCED_COUNTS CG_RHS/CG_H/shared_rank/shared_H/validation_H ')))
    assert counts[0]==1560 and counts[3]==counts[2]+counts[4]
    errors=[float(line.split()[-2]) for line in text.splitlines() if line.startswith('HSEX_REDUCED_TARGET_CHI_ERROR ')]
    assert len(errors)==6 and np.isfinite(errors).all()
    times=list(map(float,values('HSEX_REDUCED_SECONDS setup/SVD/CG/shared_build/validation/total ')))
    assert np.isfinite(times).all() and min(times)>=0
    err,residual,status=map(float,values('HSEX_REDUCED_ERROR solution/true_residual/status '))
    assert int(status) in [0,2] and np.isfinite([err,residual]).all()
    candidate=times[0]+times[3]+times[4]
    row={'grid_dimension':sizes[0],'full_RHS':sizes[1],'maxrank':sizes[2],'actual_rank':counts[2],
         'full_CG_Hpsi':counts[1],'shared_Hpsi':counts[3],'representative_validation_Hpsi':counts[4],
         'max_target_chi_relative_Frobenius_error':max(errors),'representative_solution_error':err,
         'representative_true_residual':residual,'solver_status':int(status),'solver_converged':int(status)==0,
         'setup_seconds':times[0],'diagnostic_SVD_seconds':times[1],'full_CG_seconds':times[2],
         'shared_build_and_chi_validation_seconds':times[3],'representative_H_validation_seconds':times[4],
         'candidate_seconds_including_diagnostics':candidate,'CG_to_candidate_ratio':times[2]/candidate,
         'input_sha256':hashlib.sha256((d/'input').read_bytes()).hexdigest()}
    if name.startswith('on'):
        operations,nq,full_rhs,nrep,rep_rhs=map(int,values('HSEX_REDUCED_SYMMETRY_COUNTS operations/q/full_RHS/rep_q/rep_RHS '))
        assert (operations,nq,full_rhs,nrep,rep_rhs)==(6,6,1560,3,780)
        assert counts[4]==rep_rhs
        exact=[float(line.split()[-1]) for line in text.splitlines() if line.startswith('HSEX_REDUCED_SYMMETRY_CHI_ERROR ')]
        assert len(exact)==6 and max(exact)<1e-6 and np.isfinite(exact).all()
        validation=list(map(float,values('HSEX_REDUCED_SYMMETRY_VALIDATION H_probes/H_cov/projector/energy/chi ')))
        h,sec=values('HSEX_REDUCED_SYMMETRY_CG retained_H/seconds ')
        row.update({'operations':operations,'selected_q':nq,'representative_q':nrep,'representative_RHS':rep_rhs,
                    'max_symmetry_only_chi_error':max(exact),'symmetry_H_probes':int(validation[0]),
                    'max_H_covariance_error':validation[1],'max_projector_error':validation[2],
                    'max_energy_intertwining_error_Ha':validation[3],
                    'representative_CG_Hpsi':int(h),'representative_CG_solve_only_seconds':float(sec)})
    else:
        assert counts[4]==1560
        row.update({'selected_q':6,'representative_q':6,'representative_RHS':1560,'symmetry_H_probes':0})
    exe='salmon-final' if name=='on256' else 'salmon'
    row['executable_sha256']=hashlib.sha256((a.work/exe).read_bytes()).hexdigest()
    row['candidate_Hpsi_including_symmetry_probes']=row['shared_Hpsi']+row['symmetry_H_probes']
    results[name]=row
results['comparison']={'same_rank_symmetry_speed_ratio':results['off384']['candidate_seconds_including_diagnostics']/results['on384']['candidate_seconds_including_diagnostics']}
a.output.write_text(json.dumps(results,indent=2)+'\n')
print('SI_RESPONSE_SYMMETRY_VERIFIED',len(results)-1,'cases; no W publication')
for name,row in results.items():
    print(name,row)
