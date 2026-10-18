#!/usr/bin/env python3
"""Prepare independent fixed-ion W2Se4 HSE06/SOC/NLCC convergence inputs (no submission)."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil

BOHR = 0.529177210903
AUTIME_FS = 0.024188843265857
ENGINE = '7eca90d9eafc36b22c1d16dbef03d3f5e294c0f7'
PSEUDOS = {
 'W': ('9baa14856929efb1187962fdc332993f1ab846b802cb0729e1407e9f73fe0959', 12),
 'Se': ('594f5442860fabf86e0379ea823b6bfeea88cbb3b420a3170e046efd2f15cb48', 6),
}
REMOTE = '/vol0001/hp260497/si-3d-weak-scaling/WSe2_HSE_SOC_NLCC_e1e1a89_20261018'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def f(value):
    return ('%.16e' % value).replace('e', 'd')

def common(case):
    a, b, c = 3.29, math.sqrt(3)*3.29, case['lz_A']
    z = 1.67/c
    atoms = [('W', 0., 0., .5, 1), ('W', .5, .5, .5, 1),
             ('Se', 0., 1/3, .5-z, 2), ('Se', 0., 1/3, .5+z, 2),
             ('Se', .5, 5/6, .5-z, 2), ('Se', .5, 5/6, .5+z, 2)]
    return """&units unit_system='au' /
&system yn_periodic='y', spin='noncollinear', yn_spinorbit='y',yn_symmetry='yyy',
 al=%s,%s,%s, nelem=2,natom=6,nelec=48,nstate=48,temperature=-1d0 /
&pseudo izatom=74,34,
 file_pseudo(1)='../../pseudo/W_PBE19.vps',
 file_pseudo(2)='../../pseudo/Se_PBE19.vps' /
&functional xc='hse06',hse_omega=.11d0,exx_factor_mode='off',
 hse_block_rows=1,hse_q_reduction=1,1,1,
 yn_exx_fft_orbits='n',yn_exx_k_pair_symmetry='n',
 hse_sr_group_fft='off',exx_pre_scf_threshold=%s,exx_pre_scf_steps=3 /
&rgrid num_rgrid=%s /
&kgrid num_kgrid=%s,dk_shift=.5d0,.5d0,0d0,yn_gamma_centered='n' /
&parallel nproc_k=%d,nproc_ob=1,nproc_rgrid=1,1,1 /
&atomic_red_coor
%s
/
""" % (f(a/BOHR), f(b/BOHR), f(c/BOHR), f(case['pre_scf']),
         ','.join(map(str, case['grid'])), ','.join(map(str,case['k'])), case['ranks'],
         '\n'.join(" '%s' %s %s %s %d" % (s,f(x),f(y),f(z),el) for s,x,y,z,el in atoms))

def job(case, phase):
    ranks = case['ranks']
    return """#!/bin/bash
#PJM -N "WSe2%s"
#PJM -g hp260497
#PJM -L "rscgrp=%s"
#PJM -L "node=%d"
#PJM -L "elapse=24:00:00"
#PJM --mpi proc=%d
#PJM --mpi max-proc-per-node=1
#PJM -j
set -eu
cd %s/cases/%s/%s
python3 ../../preflight.py --case .. --binary ../../../build/salmon
%s
mkdir run.lock
date -Is > started.txt
export OMP_NUM_THREADS=48 PARALLEL=48 OMP_SCHEDULE=static
export FLIB_FASTOMP=FALSE FLIB_NOHARDBARRIER=TRUE FLIB_CNTL_BARRIER_ERR=FALSE
for v in $(compgen -v SALMON_); do unset "$v"; done
set +e
/usr/bin/time -v -o launcher-time.txt mpiexec -n %d -stdout launcher-stdout -stderr launcher-stderr ../../rank-run.sh
rc=$?
set -e
printf '%%s\\n' "$rc" > launcher-exit.txt
date -Is > ended.txt
[ "$rc" = 0 ] || exit "$rc"
python3 ../../check_run.py . --phase %s
""" % (phase.replace('-','')[:6], 'small' if ranks<=384 else 'large',ranks,ranks,REMOTE,case['name'],phase,
          ('python3 ../../check_run.py ../gs-x --phase gs' if phase=='rt-x' else 'python3 ../../check_run.py ../gs --phase gs') if phase.startswith('rt-') else ':', ranks,
          'rt' if phase.startswith('rt-') else 'gs')

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pseudo',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--grid',type=int,nargs=3,default=[18,32,144])
    args=ap.parse_args()
    dest=args.output.resolve(); dest.mkdir(parents=True,exist_ok=False)
    (dest/'pseudo').mkdir()
    potentials={}
    for symbol,(digest,valence) in PSEUDOS.items():
        src=args.pseudo/(symbol+'_PBE19.vps')
        if sha(src)!=digest: raise ValueError('Incorrect original potential: '+str(src))
        shutil.copy2(src,dest/'pseudo'/src.name)
        potentials[symbol]={'sha256':digest,'valence':valence,
          'url':'https://www.openmx-square.org/vps_pao2019/%s/%s_PBE19.vps'%(symbol,symbol)}
    here=Path(__file__).resolve().parent
    for name in ('preflight.py','check_run.py','compare_spectra.py','README_JA.md','LITERATURE_JA.md','rank-run.sh','prepare.py','validate_symmetry.py'):
        shutil.copy2(here/name,dest/name)
    cases=[]
    def add(name,k=(24,18,1),grid=None,lz=25.2,dt=.00025,duration=60.,kick=1e-4,pre=1e-4):
        nk=math.prod(k)
        target=432 if nk>=1440 or math.prod(grid or args.grid)>=300000 else 144
        case=dict(name=name,k=list(k),grid=list(grid or args.grid),lz_A=lz,
                  dt_fs=dt,duration_fs=duration,kick_au=kick,pre_scf=pre)
        ng=math.prod(case['grid'])
        case['k_metrics']={'points':nk,'primitive_equal_density_N':math.sqrt(2*nk),
          'dx_Ainv':2*math.pi/3.29/k[0], 'dy_Ainv':2*math.pi/(math.sqrt(3)*3.29)/k[1],
          'exact_hex_folding':False,'contains_Gamma_K_Kprime':True}
        case['spacing_A']=[v/n for v,n in zip((3.29,math.sqrt(3)*3.29,lz),case['grid'])]
        root=dest/name; root.mkdir()
        case['stage_ranks']={}
        case['irreducible_kpoints']={}
        for phase in ('gs','gs-x','rt-x','rt-y'):
            xgroup=phase in ('gs-x','rt-x')
            nrep=nk if xgroup else (k[0]//2+1)*k[1]
            stage=dict(case)
            def memory_for(nrank):
                nfull=nrep//nrank*(1 if xgroup else 2)
                tile=3*ng*nk*16+min(ng,nrank)*ng*16
                matrix=2*ng*48*nfull*16
                return {'density_tiles':tile/2**30,'one_spinor_matrix':matrix/2**30,
                  'planning_estimate_24_matrices_plus_tiles':(24*matrix+tile)/2**30,
                  'note':'Heuristic with expanded little-group sources; not a proven RSS upper bound.'}
            choices=[r for r in range(1,nrep+1) if nrep%r==0 and memory_for(r)['planning_estimate_24_matrices_plus_tiles']<27]
            if not choices:raise ValueError('No memory planning layout for '+name+' '+phase)
            below=[r for r in choices if r<=target]
            stage['ranks']=max(below) if below else min(choices)
            case.setdefault('stage_memory_model_GiB',{})[phase]=memory_for(stage['ranks'])
            case['stage_ranks'][phase]=stage['ranks']
            case['irreducible_kpoints'][phase]=nrep
            mask='nyy' if xgroup else ('yny' if phase=='rt-y' else 'yyy')
            com=common(stage).replace("yn_symmetry='yyy'","yn_symmetry='%s'"%mask)
            folder=root/phase;folder.mkdir()
            text="! Prepared convergence candidate; numerical convergence is not established.\n"
            text+="&calculation theory='%s',yn_dc='n' /\n" % ('tddft_response' if phase.startswith('rt-') else 'dft')
            text+="&control sysname='WSe2',yn_restart='n',checkpoint_interval=-1,time_shutdown=0d0,\n"
            text+=(" write_gs_restart_data='wfn' /\n" if not phase.startswith('rt-') else " directory_read_data='../%s/data_for_restart/' /\n"%('gs-x' if phase=='rt-x' else 'gs'))
            text+=com
            if not phase.startswith('rt-'):
                text+="&scf method_init_wf='gauss',method_init_density='wf',nscf_init_no_diagonal=0,\n nscf=600,ncg=12,threshold=1d-8,method_mixing='simple_dm',mixrate=.15d0 /\n"
            else:
                text=text.replace('exx_pre_scf_threshold='+f(pre),'exx_pre_scf_threshold=0d0')
                text+="&propagation propagator='hse_taylor4' /\n"
                text+='&tgrid dt=%s,nt=%d /\n'%(f(dt/AUTIME_FS),round(duration/dt))
                text+="&emfield trans_longi='tr',ae_shape1='impulse',e_impulse=%s,epdir_re1=%s /\n"%(f(kick),'1d0,0d0,0d0' if phase=='rt-x' else '0d0,1d0,0d0')
                text+="&analysis de=%s,nenergy=800,out_rt_energy_step=100,yn_out_rt_energy_components='y' /\n"%f(.01/27.211386245988)
            (folder/'inputfile').write_text(text)
            (folder/'job.sh').write_text(job(stage,phase))
            ops=[]
            for tx,ty in ((0.,0.),(.5,.5)):
                for sx in ((1,) if xgroup else (1,-1)):
                    for sz in (1,-1):
                        ops.extend(['%d 0 0 %.1f'%(sx,tx), '0 1 0 %.1f'%ty, '0 0 %d 0'%sz])
            (folder/'sym.dat').write_text('\n'.join(ops)+'\n')
        (root/'case.json').write_text(json.dumps(case,indent=2)+'\n');cases.append(case)
    for nx,ny in ((12,8),(18,12),(24,18),(36,24),(48,30),(60,36)):
        add('k%dx%d'%(nx,ny),k=(nx,ny,1))
    add('grid24x42x192',grid=(24,42,192))
    add('grid28x48x216',grid=(28,48,216))
    add('vacuum31p5',lz=31.5,grid=(args.grid[0],args.grid[1],round(args.grid[2]*1.25)))
    add('vacuum37p8',lz=37.8,grid=(args.grid[0],args.grid[1],round(args.grid[2]*1.5)))
    add('dt-half',dt=.000125)
    add('kick-half',kick=5e-5)
    add('time120fs',duration=120.)
    add('gs-no-pre',pre=0.)
    record={'engine_commit':ENGINE,'remote_root':REMOTE,'status':'prepared_not_submitted',
       'cell':'orthogonal single-layer 1H W2Se4, not two-layer bulk 2H',
       'a_A':3.29,'Se_height_A':1.67,'potentials':potentials,'cases':cases,
       'baseline':'k24x18','source_generator':str(here),'binary_sha256':None}
    coords=[('W',0,0,12.6),('W',3.29/2,math.sqrt(3)*3.29/2,12.6),
       ('Se',0,math.sqrt(3)*3.29/3,12.6-1.67),('Se',0,math.sqrt(3)*3.29/3,12.6+1.67),
       ('Se',3.29/2,5*math.sqrt(3)*3.29/6,12.6-1.67),('Se',3.29/2,5*math.sqrt(3)*3.29/6,12.6+1.67)]
    (dest/'structure.xyz').write_text('6\nW2Se4 single layer; Angstrom; cell 3.29 sqrt(3)*3.29 25.2\n'+
      '\n'.join('%s %.12f %.12f %.12f'%r for r in coords)+'\n')
    (dest/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    sums=[]
    for p in sorted(dest.rglob('*')):
        if p.is_file(): sums.append('%s  %s'%(sha(p),p.relative_to(dest)))
    (dest/'SHA256SUMS').write_text('\n'.join(sums)+'\n')
    print('Prepared %d independent GS + x/y RT cases in %s'%(len(cases),dest))

if __name__=='__main__':main()
