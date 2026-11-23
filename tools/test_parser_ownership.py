"""Compile/read actual parser cases; compare owned A→B→A with fresh legacy parses."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

A="""&calculation theory='dft' /
&control sysname='a' /
&units unit_system='A_eV_fs' /
&system yn_periodic='y', al=8.,8.,8., nstate=4, nelec=2, nelem=1, natom=2 /
&pseudo file_pseudo(1)='Si-a.dat', izatom(1)=14 /
&functional xc='libxc_pbe', cval=.7 /
&rgrid num_rgrid=8,8,8 /
&kgrid num_kgrid=2,2,2 /
&scf mixrate=.23, nmemory_mb=5 /
&emfield E_amplitude1=.07 /
&atomic_coor
Si 0. 0. 0. 1
Si 1. 1. 1. 1
/
"""
B="""&calculation theory='dft' /
&control sysname='b' /
&system yn_periodic='y', al=9.,9.,9., nstate=4, nelec=2, nelem=1, natom=1, file_atom_red_coor='b.coords' /
&pseudo file_pseudo(1)='Si-b.dat', izatom(1)=14 /
&rgrid num_rgrid=8,8,8 /
"""
BLANK="""&calculation theory='dft' /
&system yn_periodic='y', al=8.,8.,8., nelem=0, natom=0 /
&rgrid num_rgrid=8,8,8 /
"""


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--fc',default='mpifort');p.add_argument('--mpiexec',default='mpiexec')
    p.add_argument('--ranks',type=int,default=2)
    p.add_argument('--checker',type=Path,default=Path(__file__).with_name('check_forbidden_eigensolvers.py'))
    p.add_argument('--reference-build',type=Path,help='unchanged audited library for fresh legacy parse controls')
    p.add_argument('--fast-hse',action='store_true',help='exercise fast-PR HSE-specific defaults')
    a=p.parse_args();record=json.loads(a.build.read_text());library=Path(record['library'])
    if hashlib.sha256(library.read_bytes()).hexdigest()!=record['sha256']:raise ValueError('library checksum mismatch')
    configure=record.get('configure_command') or record['commands'][0]
    build=Path(configure[configure.index('-B')+1])
    root=Path(__file__).resolve().parents[1];out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    obj=out/'test.o';exe=out/'test-parser'
    subprocess.run([a.fc,'-cpp',*(['-DSALMONPY_FAST_PARSER_TEST'] if a.fast_hse else []),'-I'+str(build),'-c',str(root/'tests/native/test_parser_ownership.f90'),'-o',str(obj)],check=True)
    subprocess.run([sys.executable,str(a.checker),'--object-tree',str(out)],check=True)
    subprocess.run([a.fc,str(obj),str(library),'-o',str(exe)],check=True)
    for mode in ('owned','legacy-a','legacy-b'):
        work=out/mode;work.mkdir()
        for name,text in (('a.inp',A.replace("cval=.7","cval=.7, hse_omega=.17, exx_gpu_batch_size=7") if a.fast_hse else A),('b.inp',B),('blank.inp',BLANK),('b.coords','Si .125 .25 .375 1\n')):
            (work/name).write_text(text)
        sentinels=('.namelist.tmp','.atomic_coor.tmp','.atomic_red_coor.tmp','.reenetrance.tmp')
        if mode=='owned':
            for name in sentinels:(work/name).write_text('previous work; preserve\n')
        command=[a.mpiexec,'-n',str(a.ranks),str(exe),mode]
        with (out/(mode+'.log')).open('w') as log:
            subprocess.run(command,cwd=work,env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'},
                           check=True,stdout=log,stderr=subprocess.STDOUT,timeout=120)
        marker='model parser ownership PASS' if mode=='owned' else 'legacy parse PASS'
        if marker not in (out/(mode+'.log')).read_text():raise ValueError('parser test did not reach completion')
        if mode=='owned':
            for name in sentinels:assert (work/name).read_text()=='previous work; preserve\n'
    for first,second in (('a.variables.log','again.variables.log'),):
        assert (out/'owned'/first).read_bytes()==(out/'owned'/second).read_bytes()
    for mode,name in (('legacy-a','a.variables.log'),('legacy-b','b.variables.log')):
        assert (out/mode/'variables.log').read_bytes()==(out/'owned'/name).read_bytes()
    if a.reference_build:
        reference=json.loads(a.reference_build.read_text());ref_lib=Path(reference['library'])
        if hashlib.sha256(ref_lib.read_bytes()).hexdigest()!=reference['sha256']:
            raise ValueError('reference library checksum mismatch')
        configure=reference.get('configure_command') or reference['commands'][0]
        ref_build=Path(configure[configure.index('-B')+1])
        control=out/'legacy-control.f90'
        control.write_text("program legacy_parser_control\n"
            "use inputoutput,only:read_input\n"
            "use math_constants,only:set_math_constants\n"
            "use parallelization,only:setup_parallel,end_parallel\n"
            "implicit none\ncharacter(4096)::filename\n"
            "call set_math_constants\ncall setup_parallel\n"
            "call get_command_argument(1,filename)\ncall read_input(trim(filename))\n"
            "call end_parallel\nend program\n")
        subprocess.run([sys.executable,str(a.checker),'--source-root',str(out)],check=True)
        ref_obj=out/'legacy-control.o';ref_exe=out/'legacy-control'
        subprocess.run([a.fc,'-I'+str(ref_build),'-c',str(control),'-o',str(ref_obj)],check=True)
        subprocess.run([sys.executable,str(a.checker),'--object-tree',str(out)],check=True)
        subprocess.run([a.fc,str(ref_obj),str(ref_lib),'-o',str(ref_exe)],check=True)
        for name in ('a','b'):
            work=out/('reference-'+name);work.mkdir()
            for f in ('a.inp','b.inp','blank.inp','b.coords'):shutil.copy2(out/'owned'/f,work/f)
            with (out/('reference-'+name+'.log')).open('w') as log:
                subprocess.run([a.mpiexec,'-n',str(a.ranks),str(ref_exe),name+'.inp'],cwd=work,
                    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'},
                    check=True,stdout=log,stderr=subprocess.STDOUT,timeout=120)
            assert (work/'variables.log').read_bytes()==(out/'owned'/(name+'.variables.log')).read_bytes()
    (out/'result.json').write_text(json.dumps({'status':'pass','ranks':a.ranks,'A_B_A':True,
        'omitted_defaults_independent':True,'scratch_sources_interleaved':True,'legacy_logs_exact':True,
        'preexisting_tempfiles_preserved':True,'units_and_atom_sources_independent':True,'fast_hse_settings':a.fast_hse,
        'old_legacy_logs_exact':bool(a.reference_build)},indent=2)+'\n')
    print('Parser ownership and legacy parity PASS')


if __name__=='__main__':main()
