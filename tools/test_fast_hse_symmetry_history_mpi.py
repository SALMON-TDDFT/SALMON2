"""Verify fast-PR symmetry histories and PT-CN: BUILD FULL REDUCED NEW_OUTPUT.

Explicit local MPI tests, no compile/download. Synthetic fixtures, not convergence.
Requires matching history capability and read-only rank-local history diagnostics.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build','full','reduced','output'):p.add_argument(name,type=Path)
    args=p.parse_args();root=Path(__file__).resolve().parents[1];out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1');commands=[];records=[]
    def run(tool,arguments,ranks,label):
        cmd=['mpiexec','-n',str(ranks),sys.executable,str(root/'tools'/tool),*map(str,arguments)];commands.append(cmd)
        with (out/(label+'.log')).open('w') as log:subprocess.run(cmd,env=env,cwd=root,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
    drive=np.zeros((32,2,3));drive[:,:,2]=(.0001*np.sin(.3*np.arange(1,33)))[:,None]
    path=out/'drive.npz';np.savez(path,A_over_c=drive)
    for label,propagator,factor,history in [('ptcn','hse_ptcn','off',None),('taylor-strict','hse_taylor4','strict',None),('taylor-learned','hse_taylor4','learned','taylor'),('ptcn-learned','hse_ptcn','learned','ptcn')]:
        bundles=[]
        for kind,source in [('full',args.full),('reduced',args.reduced)]:
            target=out/'fixtures'/label/kind;shutil.copytree(source,target);q=target/'rt.inp';text=q.read_text()
            interval=1 if propagator=='hse_ptcn' else 2
            old="exx_factor_mode='off'";new=f"exx_factor_mode='{factor}',exx_factor_exact_interval={interval},exx_factor_n_teach=2,exx_factor_warmup_steps=4"
            if old not in text:raise ValueError('canonical factor-off fixture required')
            text=text.replace(old,new).replace("propagator='hse_taylor4'","propagator='"+propagator+"'")
            if propagator=='hse_ptcn':text=text.replace("yn_predictor_corrector='y'","yn_predictor_corrector='n'")
            q.write_text(text);bundles.append(target)
        extra=['--require-history',history] if history else []
        for size in (1,2):
            name=label+'-k'+str(size)
            run('test_fast_hse_low_gradient_mpi.py',[args.build.resolve(),*bundles,out/name,'--group-size',size,'--drive',path,*extra],2*size,name)
            records.append(json.loads((out/name/'result.json').read_text())|{'mode':label})
        for split in (3,7,13):
            checkpoint=out/(label+'-checkpoint-'+str(split))
            for mode in ('save','resume'):
                name=label+'-'+str(split)+'-'+mode
                run('test_fast_hse_symmetry_mpi.py',[args.build.resolve(),bundles[1],out/name,'--mode',mode,'--checkpoint',checkpoint,'--group-size',2,'--split',split,*extra],4,name)
            traces=[]
            for rank in range(4):
                with np.load(checkpoint/f'rank-{rank}/trace.npz') as a:traces.append(a['trace'].copy())
            np.testing.assert_array_equal(traces[0],traces[1]);np.testing.assert_array_equal(traces[2],traces[3]);assert not np.array_equal(traces[0],traces[2])
    (out/'report.json').write_text(json.dumps({'status':'pass','cases':records,'commands':commands,'replayed_steps':0,'caveat':'Synthetic numerical/history validation, not material convergence.'},indent=2)+'\n')
    print('Symmetry PT-CN, learned histories, MPI and direct resume PASS',flush=True)


if __name__=='__main__':main()
