"""Bounded DC-HSE comparisons. Each case uses a fresh directory; no stored references."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[2]


def run(exe,work,name,text,nproc):
    path=work/name
    path.mkdir(parents=True,exist_ok=False)
    (path/'inputfile').write_text(text)
    shutil.copy(ROOT/'testsuites/pseudo/H_rps.dat',path)
    start=time.perf_counter()
    with (path/'inputfile').open('rb') as inp,(path/'outputfile').open('wb') as out:
        p=subprocess.run(['mpiexec','-n',str(nproc),str(exe)],cwd=path,stdin=inp,stdout=out,stderr=subprocess.STDOUT,
                         env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),timeout=180)
    output=(path/'outputfile').read_text()
    if p.returncode or 'end SALMON' not in output:
        raise AssertionError(f'{name}: calculation failed; see {path}/outputfile')
    if "yn_dc='y'" in text:
        records=[line.split() for line in output.splitlines() if line.split()[:3]==['DC','#SCF','=']]
        assert records and float(records[-1][10])<1e-8, name+' did not converge'
        info=path/'data_dcdft/fragments/000001/H_dc_hse_info.data'
    else:
        assert '#GS converged at' in output,name+' did not converge'
        info=path/'H_dc_hse_info.data'
    energy=float(re.search(r'Total energy \(eV\)\s*=\s*(\S+)',info.read_text())[1])
    result=dict(name=name,nproc=nproc,energy_eV=energy,seconds=time.perf_counter()-start)
    eigen=path/'data_dcdft/total/H_dc_hse_lcfo_complex_eigen.data'
    if eigen.exists():
        result['eigenvalues_Ha']=np.loadtxt(eigen)[:,3].tolist()
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('executable',type=Path)
    parser.add_argument('work',type=Path)
    args=parser.parse_args();exe=args.executable.resolve();work=args.work.resolve()
    work.mkdir(parents=True,exist_ok=True)
    base=(ROOT/'testsuites/422_H_dcdft_hse/inputfile').read_text().replace('threshold=1d-8','threshold=1d-10')
    results=[]
    results.append(run(exe,work,'dc_k2',base,4))
    results.append(run(exe,work,'dc_k1',base.replace('nproc_k=2','nproc_k=1').replace('nproc_rgrid_tot=4,1,1','nproc_rgrid_tot=2,1,1'),2))
    assert abs(results[0]['energy_eV']-results[1]['energy_eV'])<1e-7
    np.testing.assert_allclose(results[0]['eigenvalues_Ha'],results[1]['eigenvalues_Ha'],atol=1e-7,rtol=0)
    one=base.replace('num_fragment=2,1,1','num_fragment=1,1,1').replace('num_rgrid_buffer=4,0,0','num_rgrid_buffer=0,0,0')
    one=one.replace('nproc_rgrid_tot=4,1,1','nproc_rgrid_tot=2,1,1')
    results.append(run(exe,work,'one_fragment',one,2))
    conventional=one.replace("yn_dc='y'","yn_dc='n'").replace("xc='hse06'","xc='hse06'\n yn_hse_wannier='y'")
    results.append(run(exe,work,'conventional',conventional,2))
    assert abs(results[-1]['energy_eV']-results[-2]['energy_eV'])<1e-5, 'one-fragment energy parity failed'
    gamma=base.replace('num_kgrid=1,2,1','num_kgrid=1,1,1').replace('nproc_k=2','nproc_k=1')
    gamma=gamma.replace('nproc_rgrid_tot=4,1,1','nproc_rgrid_tot=2,1,1')
    results.append(run(exe,work,'dc_gamma',gamma,2))
    smaller=base.replace('al=16d0,8d0,8d0','al=32d0,8d0,8d0').replace('num_rgrid=16,8,8','num_rgrid=32,8,8')
    smaller=smaller.replace('11.3d0','19.3d0').replace('12.7d0','20.7d0')
    results.append(run(exe,work,'smaller_fragment',smaller,4))
    (work/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))

if __name__=='__main__':
    main()
