"""Validate a launched fixed-target Si diagnostic; collect no GS payloads."""
from pathlib import Path
import argparse,hashlib,json,re
import numpy as np
parser=argparse.ArgumentParser()
parser.add_argument('work',type=Path)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
results={}
for name in ['rank384','rank768','final256','final384','target19-rank384']:
    case=args.work/name
    text=(case/'run.log').read_text()
    assert 'end SALMON' in text and 'HSEX_INTERNAL_REDUCED_DIAGNOSTIC_PASS no W publication' in text
    assert re.search(r'number of threads\s*:\s*12',text)
    assert not (case/'hsex-w.bin').exists(),f'diagnostic published W: {name}'
    assert re.search(r'HSEX_REDUCED_DIAGNOSTIC_PASS owner=\d+ converged=[01]',text)
    def fields(label):
        line=next(line for line in text.splitlines() if line.startswith(label))
        return line[len(label):].split()
    ng,nrhs,limit=map(int,fields('HSEX_REDUCED_SIZE grid/RHS/maxrank '))
    counts=list(map(int,fields('HSEX_REDUCED_COUNTS CG_RHS/CG_H/shared_rank/shared_H/validation_H ')))
    cg_rhs,cg_h,dimension,shared_h,validate_h=counts
    assert (ng,nrhs)==(1728,1040) and cg_rhs==nrhs
    assert dimension<=limit and shared_h==dimension+validate_h and validate_h==nrhs
    error,residual,status=map(float,fields('HSEX_REDUCED_ERROR solution/true_residual/status '))
    assert int(status) in [0,2]
    times=np.array(fields('HSEX_REDUCED_SECONDS setup/SVD/CG/shared_build/validation/total '),float)
    assert np.isfinite(times).all() and (times>=0).all()
    rows=[line.split()[-3:] for line in text.splitlines()
          if line.startswith('HSEX_REDUCED_TARGET_CHI_ERROR q/relative/norm ')]
    assert len(rows)==4
    chi_errors=np.array([float(row[1]) for row in rows])
    assert np.isfinite(chi_errors).all() and (chi_errors>=0).all()
    assert np.isfinite([error,residual]).all()
    source_lines=[line for line in text.splitlines() if line.startswith('HSEX_REDUCED_Q target/')]
    target=int(source_lines[0].split()[-7])
    total_shared=float(times[0]+times[3]+times[4])
    # Both comparison costs include shared geometry/vertex setup. Only the
    # candidate also pays the diagnostic RHS SVD, conservatively included.
    reference_with_setup=float(times[2]+times[0]-times[1])
    executable='salmon-initial-validated' if name.startswith('rank') else 'salmon-reviewed'
    results[name]={
        'target_k':target,'grid_dimension':ng,'RHS':nrhs,'rank_limit':limit,'dimension':dimension,
        'CG_Hpsi':cg_h,'shared_Hpsi_including_validation':shared_h,
        'validation_Hpsi':validate_h,'CG_seconds':float(times[2]),'setup_seconds':float(times[0]),
        'diagnostic_SVD_seconds':float(times[1]),'shared_build_seconds':float(times[3]),
        'validation_seconds':float(times[4]),'shared_seconds_including_setup_SVD_validation':total_shared,
        'CG_seconds_including_common_setup':reference_with_setup,
        'fixed_target_speed_ratio':reference_with_setup/total_shared,
        'target_chi_max_relative_Frobenius_error':float(chi_errors.max()),
        'solution_relative_error':error,'max_normalized_true_residual':residual,
        'solver_status':int(status),'solver_converged':status==0,
        'q_records':source_lines,'input_sha256':hashlib.sha256((case/'input').read_bytes()).hexdigest(),
        'executable_sha256':hashlib.sha256((args.work/executable).read_bytes()).hexdigest()}
args.output.write_text(json.dumps(results,indent=2)+'\n')
print('FIXED_TARGET_SI_DIAGNOSTIC_VERIFIED',len(results),'cases; no W publication')
for name,row in results.items():
    print(name,'chi_error=',row['target_chi_max_relative_Frobenius_error'],
          'seconds=',row['shared_seconds_including_setup_SVD_validation'],
          'ratio=',row['fixed_target_speed_ratio'],'converged=',row['solver_converged'])
