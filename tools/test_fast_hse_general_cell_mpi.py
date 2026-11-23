"""Synthetic general-cell regression for explicitly built fast-PR HSE.

BUILD FULL_BUNDLE REDUCED_BUNDLE NEW_OUTPUT [--old-build OLD_BUILD].
Runs local MPI groups1/2, full/reduced comparisons and direct checkpoint resume.
Existing fixtures must use the canonical 8-cube grid,2x2x2 mesh and cos(z) orbitals.
No compilation/download, physical convergence claim or input axis normalization.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np


def check_independent_groups(checkpoint):
    traces=[]
    for rank in range(4):
        with np.load(checkpoint/f'rank-{rank}/trace.npz') as data:traces.append(data['trace'].copy())
    np.testing.assert_array_equal(traces[0],traces[1])
    np.testing.assert_array_equal(traces[2],traces[3])
    if np.array_equal(traces[0],traces[2]):raise AssertionError('different screening parameters produced identical state traces')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('build','full','reduced','output'):parser.add_argument(name,type=Path)
    parser.add_argument('--old-build',type=Path)
    parser.add_argument('--roundoff-reference-build',type=Path,help='comparison-only old kernels with admission tolerance fix; the uncorrected old core rejects this case')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    build=args.build.resolve();env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    commands=[];records=[]
    def run(tool,arguments,ranks,label):
        cmd=['mpiexec','-n',str(ranks),sys.executable,str(root/'tools'/tool),*map(str,arguments)]
        commands.append(cmd)
        with (out/(label+'.log')).open('w') as log:
            subprocess.run(cmd,env=env,cwd=root,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
    cases={
        'rotated': 'al_vec1=0,8,0,al_vec2=-8,0,0,al_vec3=0,0,8',
        'left-handed': 'al_vec1=0,8,0,al_vec2=8,0,0,al_vec3=0,0,8',
        'shear': 'al_vec1=8,0,0,al_vec2=1,8,0,al_vec3=0,0,8',
        'permuted-z': 'al_vec1=0,0,8,al_vec2=8,0,0,al_vec3=0,8,0',
    }
    for case,vectors in cases.items():
        bundles=[]
        for kind,source in [('full',args.full),('reduced',args.reduced)]:
            target=out/'fixtures'/case/kind;shutil.copytree(source,target)
            p=target/'rt.inp';original=p.read_text();text=original.replace('al=8,8,8','al=8,8,8,'+vectors)
            if text==original:raise ValueError('canonical al=8,8,8 fixture required')
            if case=='permuted-z' and kind=='reduced':
                text=text.replace("yn_symmetry='yyn'","yn_symmetry='nyy'")
                ops=np.zeros((2,3,4));ops[:,:,:3]=np.eye(3);ops[1,:,:3]=np.diag([1,-1,-1])
                np.savetxt(target/'sym.dat',ops.reshape(6,4),fmt='%.1f')
            p.write_text(text);bundles.append(target)
        for size in (1,2):
            label=case+'-k'+str(size)
            run('test_fast_hse_low_gradient_mpi.py',[build,*bundles,out/label,'--group-size',size],2*size,label)
            records.append(json.loads((out/label/'result.json').read_text())|{'case':case})
        # Each independent material group uses its own omega, demonstrating no
        # cross-group communication or checkpoint/state contamination.
        checkpoint=out/(case+'-checkpoint')
        common=[build,bundles[1]]
        run('test_fast_hse_symmetry_mpi.py',[*common,out/(case+'-save'),'--mode','save','--checkpoint',checkpoint,'--group-size',2],4,case+'-save')
        run('test_fast_hse_symmetry_mpi.py',[*common,out/(case+'-resume'),'--mode','resume','--checkpoint',checkpoint,'--group-size',2],4,case+'-resume')
        check_independent_groups(checkpoint)
    # Actual three-component A inputs, with matching retained little groups.
    # Circular xy uses mirror-z, which preserves the entire polarization plane.
    for polarization,diagonal,flags in [('x',[1,-1,-1],'nyy'),('y',[-1,1,-1],'yny'),('circular-xy',[1,1,-1],'yyy')]:
        bundles=[]
        for kind,source in [('full',args.full),('reduced',args.reduced)]:
            target=out/'fixtures'/polarization/kind;shutil.copytree(source,target)
            if kind=='reduced':
                p=target/'rt.inp';p.write_text(p.read_text().replace("yn_symmetry='yyn'","yn_symmetry='"+flags+"'"))
                ops=np.zeros((2,3,4));ops[:,:,:3]=np.eye(3);ops[1,:,:3]=np.diag(diagonal)
                np.savetxt(target/'sym.dat',ops.reshape(6,4),fmt='%.1f')
            bundles.append(target)
        step=np.arange(1,33);drive=np.zeros((32,2,3))
        if polarization=='circular-xy':
            drive[:,:,0]=(.0001*(np.cos(.3*step)-1))[:,None]
            drive[:,:,1]=(.0001*np.sin(.3*step))[:,None]
        else:drive[:,:,'xyz'.index(polarization)]=(.0001*np.sin(.3*step))[:,None]
        path=out/(polarization+'-drive.npz');np.savez(path,A_over_c=drive)
        for size in (1,2):
            label=polarization+'-k'+str(size)
            run('test_fast_hse_low_gradient_mpi.py',[build,*bundles,out/label,'--group-size',size,'--drive',path],2*size,label)
            records.append(json.loads((out/label/'result.json').read_text())|{'case':polarization})
        checkpoint=out/(polarization+'-checkpoint')
        for mode in ('save','resume'):
            run('test_fast_hse_symmetry_mpi.py',[build,bundles[1],out/(polarization+'-'+mode),'--mode',mode,'--checkpoint',checkpoint,'--group-size',2,'--polarization',polarization],4,polarization+'-'+mode)
        check_independent_groups(checkpoint)
    # Ordinary diagonal cells retain the prior optimized arithmetic, including
    # nonbinary-exact lengths for which n*h does not exactly reconstruct a.
    for label,counts,lengths in [('diagonal',(8,8,8),(8.,8.,8.)),('diagonal-roundoff',(12,12,12),(7.3,7.3,7.3))]:
        target=out/'fixtures'/label;shutil.copytree(args.reduced,target)
        if counts[0]!=8:
            p=target/'rt.inp';p.write_text(p.read_text().replace('al=8,8,8','al='+','.join(map(str,lengths))).replace('num_rgrid=8,8,8','num_rgrid='+','.join(map(str,counts))))
            z=np.arange(counts[2])[None,None,:];psi=np.empty((*counts,2),dtype='=c16',order='F')
            volume=float(np.prod(lengths));psi[...,0]=1/np.sqrt(volume);psi[...,1]=np.sqrt(2/volume)*np.cos(2*np.pi*z/counts[2])
            (target/'data_for_restart/wfn.bin').write_bytes(np.repeat(psi[...,None],4,axis=4).tobytes(order='F'))
        reference=args.old_build if label=='diagonal' else args.roundoff_reference_build
        for which,engine in [('new',build),('old',reference.resolve() if reference else None)]:
            if engine is None:continue
            run('test_fast_hse_symmetry_mpi.py',[engine,target,out/(label+'-'+which),'--mode','save','--checkpoint',out/(label+'-'+which+'-checkpoint'),'--group-size',2],4,label+'-'+which)
        if reference:
            for rank in range(4):
                with np.load(out/(label+'-new-checkpoint')/f'rank-{rank}/trace.npz') as a,np.load(out/(label+'-old-checkpoint')/f'rank-{rank}/trace.npz') as b:
                    np.testing.assert_array_equal(a['trace'],b['trace'])
    (out/'report.json').write_text(json.dumps({'status':'pass','cases':records,'diagonal_bit_equal':bool(args.old_build),'roundoff_old_kernel_bit_equal':bool(args.roundoff_reference_build),'commands':commands,'caveat':'Synthetic initialization/numerical regression, not converged material physics.'},indent=2)+'\n')
    print('General-cell full/reduced MPI, independent states and direct resume PASS',flush=True)


if __name__=='__main__':main()
