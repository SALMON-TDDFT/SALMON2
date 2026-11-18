"""Synthetic SCF ownership and extracted convergence parity, without SCF runs.

BUILD_RECORD ORIGINAL_ITERATION_SOURCE NEW_OUTPUT. Requires an explicitly built
standard MPI runtime and the pinned pre-refactor source as comparison evidence.
No archive, GS, executable download, implicit rebuild or solver call is made.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

REFERENCE_SHA256='70137c4a3fb995f35815ebd3df4674325714af6b115815be8cdaa97e1cf1c5c8'


def wrapper(name,text,owned):
    def routine(label):
        found=re.search(r'(?ims)^  subroutine '+label+r'\(.*?^  end subroutine '+label+r'\s*$',text)
        if found is None:raise ValueError('missing extracted convergence routine: '+label)
        return found[0]
    declarations='' if owned else 'type(s_scalar) :: rho_old,Vlocal_old\nreal(8) :: rNe\n'
    return f'''module {name}
use structures
use mpi
implicit none
contains
subroutine measure(oldrho,newrho,oldpot,newpot,mode,spin,electrons,spin_electrons,result)
implicit none
real(8),intent(in) :: oldrho(:,:,:),newrho(:,:,:),oldpot(:,:,:),newpot(:,:,:)
character(*),intent(in) :: mode
integer,intent(in) :: spin,electrons,spin_electrons(2)
real(8),intent(out) :: result
type(s_rgrid) :: lg,mg
type(s_dft_system) :: system
type(s_parallel_info) :: info
type(s_scalar) :: rho,V_local(spin)
type(s_mixing) :: mixing
type(s_dcdft) :: dc
integer :: ix,iy,iz,nelec,nelec_spin(2)
character(12) :: convergence
character :: yn_dc
real(8) :: sum1
{declarations}
mg%is=[-1,2,0];mg%ie=mg%is+shape(oldrho)-1
lg%num=shape(oldrho);system%nspin=spin;system%Hvol=.25d0
info%icomm_r=MPI_COMM_WORLD
nelec=electrons;nelec_spin=spin_electrons;yn_dc='n';convergence=mode
call allocate_scalar(mg,rho)
call allocate_scalar(mg,V_local(1))
rho%f=oldrho;V_local(1)%f=oldpot
call init_convergence_check()
rho%f=newrho;V_local(1)%f=newpot
call convergence_check(lg,mg,system,info,rho,V_local)
result=sum1
contains
{routine('init_convergence_check')}
{routine('convergence_check')}
subroutine comm_summation(local,total,comm)
real(8),intent(in) :: local
real(8),intent(out) :: total
integer,intent(in) :: comm
integer :: ierr
call MPI_Allreduce(local,total,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
if(ierr/=MPI_SUCCESS)stop 90
end subroutine
end subroutine
end module
'''


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build',type=Path)
    parser.add_argument('reference',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args();root=Path(__file__).resolve().parents[1]
    record=json.loads(args.build.read_text());library=Path(record['library'])
    if record.get('schema')!='salmonpy.embedded-tddft/1' or hashlib.sha256(library.read_bytes()).hexdigest()!=record['sha256']:
        raise ValueError('unchanged explicit standard build required')
    if hashlib.sha256(args.reference.read_bytes()).hexdigest()!=REFERENCE_SHA256:
        raise ValueError('pinned pre-refactor iteration source required')
    build=args.build.resolve().parent/'build'
    source=Path(record['configure_command'][record['configure_command'].index('-S')+1])
    current=source/'src/gs/scf_iteration_dft.f90'
    if hashlib.sha256(current.read_bytes()).hexdigest()!=record['source_hashes']['src/gs/scf_iteration_dft.f90']:
        raise ValueError('built source changed')
    # Compare the entire executable iteration body after reversing only the
    # ownership substitution and removing the new end-of-call scalar capture.
    before=args.reference.read_text();after=current.read_text()
    for owned,local in [('rho_old','rho_old'),('Vlocal_old','Vlocal_old'),
                        ('electron_normalization','rNe'),('converged','flag_conv'),
                        ('no_mix_converged_count','icnt_conv_nomix')]:
        after=after.replace('mixing%scf_progress%'+owned,local)
    def body(text):
        text=text[text.index('call init_convergence_check'):text.index('end subroutine scf_iteration_dft')]
        lines=[line for line in text.splitlines() if not line.lstrip().startswith('!')
               and not line.startswith('mixing%scf_progress%')]
        return ''.join(''.join(lines).split()).lower()
    if body(before)!=body(after):raise AssertionError('iteration arithmetic/order changed beyond ownership substitution')
    out=args.output.resolve()
    if any(out.is_relative_to(p.resolve()) for p in (root,source,build)):
        raise ValueError('new external output required')
    out.mkdir(parents=True,exist_ok=False);probe=out/'src';probe.mkdir()
    shutil.copy2(root/'tests/native/test_scf_progress.f90',probe/'ownership.f90')
    (probe/'baseline.f90').write_text(wrapper('baseline_convergence',args.reference.read_text(),False))
    (probe/'owned.f90').write_text(wrapper('owned_convergence',current.read_text(),True))
    shutil.copy2(root/'tests/native/test_scf_convergence.f90',probe/'convergence.f90')
    flags=(build/'src/CMakeFiles/salmon.dir/flags.make').read_text()
    options=[v for v in shlex.split(re.search(r'^Fortran_FLAGS\s*=\s*(.*)$',flags,re.M)[1]) if not v.startswith('-J')]
    includes=shlex.split(re.search(r'^Fortran_INCLUDES\s*=\s*(.*)$',flags,re.M)[1])
    compiler=re.search(r'^CMAKE_Fortran_COMPILER:[^=]*=(.*)$',(build/'CMakeCache.txt').read_text(),re.M)[1].strip()
    from build_embedded_tddft import mpi_link_command
    checker=root/'tools/check_forbidden_eigensolvers.py';commands=[]
    with (out/'build.log').open('w') as log:
        def run(command,**kwargs):
            command=list(map(str,command));commands.append(command)
            subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120,**kwargs)
        run([sys.executable,checker,'--source-root',source/'src'])
        run([sys.executable,checker,'--source-root',probe])
        for name in ('ownership','baseline','owned','convergence'):
            run([compiler,*options,*includes,'-I'+str(out),'-I'+str(build),'-J'+str(out),
                 '-c',probe/(name+'.f90'),'-o',out/(name+'.o')],cwd=out)
        # Every linked application object is covered, including structures.o.
        run([sys.executable,checker,'--object-tree',build])
        run([sys.executable,checker,'--object-tree',out])
        structures=build/'src/CMakeFiles/salmon.dir/common/structures.f90.o'
        for name,objects in [('ownership',['ownership']),('convergence',['baseline','owned','convergence'])]:
            run(mpi_link_command([compiler,'-fopenmp',*[str(out/(n+'.o')) for n in objects],str(structures),'-o',str(out/name)]),cwd=out)
        for threads in (1,2):
            for ranks in (1,2):
                env=dict(os.environ,OMP_NUM_THREADS=str(threads),OPENBLAS_NUM_THREADS='1')
                for name in ('ownership','convergence'):
                    run(['mpiexec','-n',ranks,out/name],cwd=out,env=env)
    result={'status':'pass','SCF_calls':0,'MPI_ranks':[1,2],'OMP_threads':[1,2],
            'reference_source_sha256':REFERENCE_SHA256,'runtime_sha256':record['sha256'],
            'ownership_independent':True,'convergence_bit_exact':True,'iteration_body_equivalent':True,
            'scope':'synthetic ownership and extracted convergence routines; not SCF continuation',
            'commands':commands}
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(out/'result.json')


if __name__=='__main__':main()
