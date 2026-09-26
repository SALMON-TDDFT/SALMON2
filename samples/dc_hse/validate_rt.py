"""Small occupied-only GS/RT smoke test and U-refresh interval comparison."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import numpy as np
from validate import ROOT


def launch(exe, path, text):
    path.mkdir(parents=True, exist_ok=False)
    (path/'inputfile').write_text(text)
    shutil.copy(ROOT/'testsuites/pseudo/H_rps.dat',path)
    with (path/'inputfile').open('rb') as inp, (path/'outputfile').open('wb') as out:
        p=subprocess.run(['mpiexec','-n','2',str(exe)],cwd=path,stdin=inp,stdout=out,
                         stderr=subprocess.STDOUT,timeout=180,
                         env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
    output=(path/'outputfile').read_text()
    assert p.returncode==0 and 'end SALMON' in output, str(path/'outputfile')
    return output


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('executable',type=Path)
    parser.add_argument('work',type=Path)
    args=parser.parse_args(); exe=args.executable.resolve(); work=args.work.resolve()
    base=(ROOT/'testsuites/422_H_dcdft_hse/inputfile').read_text()
    base=re.sub(r'&dc\n.*?\n/\n','',base,flags=re.S)
    base=base.replace("yn_dc='y'","yn_dc='n'").replace('nstate=4','nstate=2')
    base=base.replace('temperature_k=300d0','').replace("xc='hse06'","xc='hse06'\n yn_hse_wannier='y'")
    base=base.replace('threshold=1d-8','threshold=1d-10').replace("write_gs_restart_data='no'","write_gs_restart_data='yes'")
    output=launch(exe,work/'gs',base)
    assert '#GS converged at' in output
    rt=base.replace("theory='dft'","theory='tddft_response'")
    rt=rt.replace("write_gs_restart_data='yes'","directory_read_data='../gs/data_for_restart/'")
    rt=rt.replace("sysname='H_dc_hse'","sysname='H_dc_hse'\n checkpoint_interval=2")
    rt+='''
&tgrid
 nt=4
 dt=0.02d0
/
&emfield
 trans_longi='tr'
 ae_shape1='impulse'
 e_impulse=1d-4
 epdir_re1=0d0,1d0,0d0
/
&analysis
 out_rt_energy_step=1
/
'''
    traces=[]
    for interval in [1,10]:
        text=rt.replace("yn_hse_wannier='y'",f"yn_hse_wannier='y'\n hse_mlwf_interval={interval}")
        path=work/f'rt_interval{interval}'
        output=launch(exe,path,text)
        assert 'HSE_WANNIER' in output
        energy=np.loadtxt(path/'H_dc_hse_rt_energy.data')
        current=np.loadtxt(path/'H_dc_hse_rt.data')
        assert np.isfinite(energy).all() and np.isfinite(current).all()
        assert len(energy)>=4 and len(current)>=4
        traces.append((energy,current))
    np.testing.assert_allclose(traces[0][0],traces[1][0],rtol=0,atol=1e-9)
    np.testing.assert_allclose(traces[0][1],traces[1][1],rtol=0,atol=1e-9)
    restart=rt.replace("directory_read_data='../gs/data_for_restart/'",
                       "directory_read_data='../rt_interval10/checkpoint_rt_000002/'\n yn_restart='y'")
    launch(exe,work/'restart',restart)
    restart_energy=np.atleast_2d(np.loadtxt(work/'restart/H_dc_hse_rt_energy.data'))
    restart_current=np.atleast_2d(np.loadtxt(work/'restart/H_dc_hse_rt.data'))
    np.testing.assert_allclose(restart_energy[-1],traces[1][0][-1],rtol=0,atol=1e-9)
    np.testing.assert_allclose(restart_current[-1],traces[1][1][-1],rtol=0,atol=1e-9)
    result=dict(steps=4,restart_after_step=2,
                max_restart_energy_difference=float(np.max(abs(restart_energy[-1]-traces[1][0][-1]))),intervals=[1,10],max_energy_trace_difference=float(np.max(abs(traces[0][0]-traces[1][0]))),
                max_current_trace_difference=float(np.max(abs(traces[0][1]-traces[1][1]))))
    (work/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
