#!/usr/bin/env python3
"""Small rectangular Si2 developer fixture: SOC/NLCC symmetry on/off GS and x/y RT."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--salmon',type=Path,required=True);ap.add_argument('--pseudo',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--threads',type=int,default=2)
    ap.add_argument('--pair-fft', choices=('y','n'), help='Legacy comparison switch')
    ap.add_argument('--k-algorithm', choices=('auto','reference','pair','orbits'), default='auto')
    ap.add_argument('--resume',action='store_true',help='Reuse only completed directories with exactly identical input')
    ap.add_argument('--steps',type=int,default=100)
    ap.add_argument('--y-steps',type=int)
    args=ap.parse_args();root=args.output.resolve();root.mkdir(parents=True,exist_ok=args.resume)
    args.y_steps=args.y_steps or args.steps
    assert args.steps>0 and args.y_steps>0
    binary=args.salmon.resolve();pp=args.pseudo.resolve()
    assert hashlib.sha256(pp.read_bytes()).hexdigest()=='9a274b3f50cd247d53e735b8ad1669d491286934b61c295bcd6e8b6ba555d80a'
    digest=hashlib.sha256(binary.read_bytes()).hexdigest()
    pin=root/'execution-binary.sha256'
    if args.resume:
        assert pin.read_text().strip()==digest,'Cannot reuse results from a different binary'
    else:pin.write_text(digest+'\n')
    report={'binary_sha256':digest,'threads':args.threads,'yn_exx_spinor_pair_fft':args.pair_fft,'exx_k_algorithm':args.k_algorithm,
            'fixture':'coarse rectangular bcc-like Si2 developer cell; no material accuracy claim'}
    base="""&units unit_system='au' /
&system yn_periodic='y',al=8d0,10d0,12d0,nelem=1,natom=2,nelec=8,nstate=8,
 spin='noncollinear',yn_spinorbit='y',yn_symmetry='MASK' /
&pseudo izatom(1)=14,file_pseudo(1)='PP' /
&functional xc='hse06',exx_pre_scf_threshold=0d0,hse_block_rows=1, SELECTOR /
&rgrid num_rgrid=8,8,8 /
&kgrid num_kgrid=3,1,1 /
&parallel nproc_k=1,nproc_ob=1,nproc_rgrid=1,1,1 /
&atomic_red_coor
 'Si' 0d0 0d0 0d0 1
 'Si' .5d0 .5d0 .5d0 1
/
""".replace('PP',str(pp)).replace('SELECTOR', ','.join(x for x in [
        "exx_k_algorithm='%s'"%args.k_algorithm if args.k_algorithm!='auto' else '',
        "yn_exx_spinor_pair_fft='%s'"%args.pair_fft if args.pair_fft else ''] if x))
    def run(name,text,sxlist):
        d=root/name
        if args.resume and d.exists():
            assert (d/'inputfile').read_text()==text,'Cannot reuse mismatched input '+name
            log=(d/'run.log').read_text()
            assert 'end SALMON' in log and 'ERROR STOP' not in log,'Cannot reuse incomplete '+name
            return d,log
        d.mkdir();(d/'inputfile').write_text(text)
        ops=[]
        for tx,ty,tz in ((0,0,0),(.5,.5,.5)):
            for sx in sxlist:
                for sz in (1,-1):
                    ops.extend(['%d 0 0 %.1f'%(sx,tx),'0 1 0 %.1f'%ty,'0 0 %d %.1f'%(sz,tz)])
        (d/'sym.dat').write_text('\n'.join(ops)+'\n')
        with (d/'inputfile').open() as inp,(d/'run.log').open('w') as log:
            rc=subprocess.run([str(binary)],cwd=d,stdin=inp,stdout=log,stderr=subprocess.STDOUT,
                env=dict(os.environ,OMP_NUM_THREADS=str(args.threads),OPENBLAS_NUM_THREADS='1'),timeout=1800).returncode
        log=(d/'run.log').read_text()
        assert rc==0 and 'end SALMON' in log and 'ERROR STOP' not in log,name+' failed'
        expected='spinor_pair' if args.k_algorithm=='pair' or args.pair_fft=='y' or (args.k_algorithm=='auto' and args.pair_fft!='n') else 'reference'
        assert 'selected='+expected in log, name+' algorithm selection not logged'
        return d,log
    energies={}
    for mode,mask,rot in [('off','n',(1,)),('on','yyy',(1,-1)),('x','nyy',(1,))]:
        text="&calculation theory='dft' /\n&control sysname='Si',write_gs_restart_data='wfn' /\n"+base.replace('MASK',mask)
        text+="&scf nscf=600,ncg=12,threshold=1d-8,method_mixing='simple_dm',mixrate=.2d0,method_init_density='wf' /\n"
        d,log=run('gs-'+mode,text,rot)
        residual=float(re.findall(r'iter and int_x.*?=\s*\d+\s+([0-9.Ee+-]+)',log)[-1])
        assert residual<1e-8 and '#GS converged at' in log,'GS nonconvergence'
        electrons=float(re.findall(r'Ne=\s*([0-9.Ee+-]+)',log)[-1]);assert abs(electrons-8)<1e-10
        energies[mode]=float(re.findall(r'Total Energy=\s*([0-9.Ee+-]+)',log)[-1])
        report['gs-'+mode]={'residual':residual,'printed_energy_eV':energies[mode]}
    report['gs_on_off_energy_difference_eV']=abs(energies['on']-energies['off'])
    report['gs_x_off_energy_difference_eV']=abs(energies['x']-energies['off'])
    assert max(abs(energies[k]-energies['off']) for k in ('x','on'))<1e-4,'GS energy differs'
    currents={}
    for axis in ('x','y'):
      for mode in ('off','on','direct'):
        gs='off' if mode=='off' else ('x' if axis=='x' else 'on')
        mask='n' if mode=='off' else ('nyy' if axis=='x' else 'yny')
        field='1d0,0d0,0d0' if axis=='x' else '0d0,1d0,0d0'
        text="&calculation theory='tddft_response' /\n&control sysname='Si',directory_read_data='../gs-%s/data_for_restart/' /\n"%gs
        text+=base.replace('MASK',mask)
        text+="&propagation propagator='%s' /\n"%('hse_taylor4_full' if mode=='direct' else 'hse_taylor4')
        steps=args.steps if axis=='x' else args.y_steps
        text+="&tgrid dt=.02d0,nt=STEPS /\n&emfield trans_longi='tr',ae_shape1='impulse',e_impulse=1d-4,epdir_re1=%s /\n"%field
        text=text.replace('nt=STEPS','nt=%d'%steps)
        text+="&analysis out_rt_energy_step=1,yn_out_rt_energy_components='y',de=.001d0,nenergy=200 /\n"
        d,log=run('rt-%s-%s'%(axis,mode),text,(1,) if axis=='x' or mode=='off' else (1,-1))
        rows=[[float(x) for x in s.split()] for s in (d/'Si_rt.data').read_text().splitlines() if s.strip() and not s.lstrip().startswith('#')]
        assert len(rows)==steps and all(math.isfinite(x) for r in rows for x in r),'RT output incomplete/nonfinite'
        currents[axis,mode]=rows
      # Current is columns 14..16; other non-time fields must also match for the identical impulse.
      delta=max(abs(a-b) for r,s in zip(currents[axis,'off'],currents[axis,'on']) for a,b in zip(r[1:],s[1:]))
      exact=max(abs(a-b) for r,s in zip(currents[axis,'direct'],currents[axis,'on']) for a,b in zip(r[1:],s[1:]))
      report[axis+'_steps']=steps
      report[axis+'_rt_off_on_max_difference']=delta;report[axis+'_rt_ace_direct_max_difference']=exact
      assert delta<1e-7,'symmetry RT differs'
      assert exact<1e-8,'symmetry ACE/direct differs'
    report['result']='PASS';(root/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
