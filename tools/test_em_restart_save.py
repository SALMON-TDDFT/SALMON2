"""Explicit GS-free test of maintained EM restart writers (not run on import).

Build exact maintained procedures in a small MPI harness. Optional reference
source compares every supported payload's bytes. This is per-file publication,
not an atomic multi-file archive, crash durability, or a cross-platform claim.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('output', type=Path)
    p.add_argument('--source-root', type=Path)
    p.add_argument('--reference-source', type=Path)
    p.add_argument('--compiler', default='mpifort')
    p.add_argument('--launcher', default='mpiexec')
    a = p.parse_args()
    repo = Path(__file__).resolve().parents[1]
    source = a.source_root or repo/'runtime/sources/solver'
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    records = []

    def run(cmd, cwd, name):
        env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
        with (cwd/(name+'.log')).open('w') as log:
            subprocess.run([str(x) for x in cmd], cwd=cwd, env=env,
                           stdout=log, stderr=subprocess.STDOUT, check=True, timeout=120)

    for label, tree in [('current', source), ('reference', a.reference_source)]:
        if tree is None:
            continue
        tree = tree.resolve()
        run(['python3', repo/'tools/check_forbidden_eigensolvers.py', '--source-root', tree/'src'], out, label+'-source-audit')
        text = (tree/'src/maxwell/common_maxwell.f90').read_text()
        begin = text.index('  subroutine output_r_txt_em(')
        last = 'restart_save_status' if label=='current' else 'output_r_bin_em'
        end = text.index('  end subroutine '+last, begin)+len('  end subroutine '+last)
        procedures = text[begin:end]
        calls = []
        for mode in ('single', 'all'):
            for key, value in [('i0d','7'), ('i1d','ints'), ('r0d','2d0'), ('r1d','x(:,1,1,1,1)')]:
                calls.append(f"call output_r_txt_em(600,1,2,'{mode}','{mode}_{key}',{key}={value},io_status=ios)")
            for kind, array in [('r','x'), ('c','z')]:
                for dims in range(1,6):
                    section = ','.join([':']*dims+['1']*(5-dims))
                    calls.append(f"call output_r_bin_em(600,lo,hi,'{mode}','{mode}_{kind}{dims}', &\n {kind}{dims}d={array}({section}),io_status=ios)")
        calls.append("call output_r_bin_em(600,lo,hi,'all','empty',r5d=x,ipe=0,io_status=ios)")
        driver = '''module parallelization
 integer::nproc_id_global=0
end module
module communication
contains
 logical function comm_is_root(i)
 integer,intent(in)::i
 comm_is_root=(i==0)
 end function
end module
module restart_test
contains
 function resolve_path(path) result(resolved)
 character(*),intent(in)::path
 character(:),allocatable::resolved
 resolved=path
 end function
'''+procedures+'''
end module
program probe
 use mpi
 use parallelization
 use restart_test
 implicit none
 integer::err,rank,unit,ios=0,bad,total,i
 integer::lo(5)=1,hi(5)=2,ints(2)=[2,3]
 real(8)::x(2,2,2,2,2)
 complex(8)::z(2,2,2,2,2)
 character(32)::mode,name
 logical::opened
 call MPI_Init(err)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,err)
 nproc_id_global=rank
 call get_command_argument(1,mode)
 write(name,'(A,I0)')'host-',rank
 open(unit=600,file=trim(name),status='new')
 write(600,*)'before'
 x=reshape([(dble(i),i=1,32)],shape(x));z=cmplx(x,-x,kind=8)
 if(trim(mode)=='normal')then
'''+ '\n'.join(calls)+'''
 else
   call output_r_bin_em(600,lo,hi,'all','probe',r5d=x,io_status=ios)
 endif
 inquire(unit=600,opened=opened)
 if(.not.opened)error stop 'caller unit closed'
 write(600,*)'after'
 close(600)
 write(name,'(A,I0)')'status-',rank
 open(newunit=unit,file=trim(name),status='new')
 write(unit,*)ios
 close(unit)
 bad=0
 if(ios/=0)bad=1
 call MPI_Allreduce(bad,total,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,err)
 if(trim(mode)=='normal'.and.total/=0)error stop 'normal save failed'
 call MPI_Finalize(err)
end program
'''
        if label=='reference':
            driver = driver.replace(',io_status=ios','')
        build = out/(label+'-build');build.mkdir()
        (build/'probe.f90').write_text(driver)
        run([a.compiler,'-c','probe.f90'],build,'compile')
        run(['python3',repo/'tools/check_forbidden_eigensolvers.py','--object-tree',build],build,'object-audit')
        run([a.compiler,'probe.o','-o','probe'],build,'link')
        for ranks in (1,2):
            normal = out/f'{label}-normal-{ranks}';normal.mkdir();(normal/'data_for_restart').mkdir()
            run([a.launcher,'-n',ranks,build/'probe','normal'],normal,'run')
            if label=='reference':
                control = out/f'current-normal-{ranks}'/'data_for_restart'
                files = list((normal/'data_for_restart').iterdir())
                assert {f.name for f in files}=={f.name for f in control.iterdir()}
                assert all(f.read_bytes()==(control/f.name).read_bytes() for f in files)
                records.append({'reference_ranks':ranks,'files_byte_equal':len(files)})
                continue
            for fault in ('pending-file','pending-directory','target-directory','info-directory','missing-parent'):
                d=out/f'{fault}-{ranks}';d.mkdir();data=d/'data_for_restart'
                if fault!='missing-parent':data.mkdir()
                for rank in range(ranks):
                    if not data.exists():continue
                    target=data/f'probe_{rank}.bin';target.write_bytes(b'previous checkpoint')
                    if rank!=ranks-1:continue
                    if fault=='pending-file':Path(str(target)+'.salmonpy-pending').write_bytes(b'previous failed attempt')
                    elif fault=='pending-directory':Path(str(target)+'.salmonpy-pending').mkdir()
                    elif fault=='target-directory':target.unlink();target.mkdir();(target/'sentinel').write_bytes(b'old directory')
                    elif fault=='info-directory':(data/f'probe_{rank}_inf.bin').mkdir()
                before={str(f.relative_to(d)):f.read_bytes() for f in d.rglob('*') if f.is_file()}
                run([a.launcher,'-n',ranks,build/'probe','failure'],d,'failed-attempt')
                codes=[int((d/f'status-{rank}').read_text()) for rank in range(ranks)]
                assert codes[-1]!=0
                for name,value in before.items():
                    # Other healthy ranks may publish their files. Only failed owner is checked.
                    if f'probe_{ranks-1}' in name:assert (d/name).read_bytes()==value
                evidence=d/'failed-attempt';evidence.mkdir()
                for f in list(d.iterdir()):
                    if f.name in ('failed-attempt','failed-attempt.log'):continue
                    f.rename(evidence/f.name)
                shutil.copytree(evidence/'data_for_restart',data) if (evidence/'data_for_restart').exists() else data.mkdir()
                for f in list(data.iterdir()):
                    if f.name.endswith('.salmonpy-pending') or f.is_dir():
                        # Move repair obstacles aside; keep all failed evidence.
                        f.rename(d/('repair-obstacle-'+f.name))
                run([a.launcher,'-n',ranks,build/'probe','normal'],d,'repair-retry')
                # Retry same probe target, separately from the comprehensive normal payload set.
                for rank in range(ranks):
                    for name in (f'host-{rank}',f'status-{rank}'):(d/name).rename(d/('normal-'+name))
                run([a.launcher,'-n',ranks,build/'probe','failure'],d,'probe-retry')
                assert all(int((d/f'status-{rank}').read_text())==0 for rank in range(ranks))
                records.append({'fault':fault,'ranks':ranks,'status':codes,'repair_retry':True,'failed_owner_existing_bytes_preserved':True})
        records.append({'build':label,'executable_sha256':hashlib.sha256((build/'probe').read_bytes()).hexdigest()})
    (out/'summary.json').write_text(json.dumps({'cases':records,'scope':'per-file stage/close/rename; no multi-file transaction or crash durability'},indent=2)+'\n')
    print('PASS',len(records),'records')


if __name__=='__main__':
    main()
