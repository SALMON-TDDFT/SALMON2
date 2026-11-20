"""Synthetic Fortran I/O ownership checks and optional old-output byte controls.

No SCF, GS archive, physical density or wavefunction data is generated.
"""
import argparse,hashlib,json,os,re,subprocess,sys
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--reference-build',type=Path)
    p.add_argument('--directory-scope',action='store_true')
    p.add_argument('--host-unit-library',action='store_true')
    p.add_argument('--host-timer-library',action='store_true')
    p.add_argument('--host-mpi-library',action='store_true')
    p.add_argument('--checkpoint-scratch',action='store_true')
    p.add_argument('--scalar-inquiries',action='store_true')
    p.add_argument('--finite-source',type=Path)
    p.add_argument('--field-close',action='store_true')
    p.add_argument('--serializer-reference',type=Path)
    p.add_argument('--fc',default='mpifort');p.add_argument('--ranks',type=int,default=2)
    p.add_argument('--checker',type=Path,default=Path(__file__).with_name('check_forbidden_eigensolvers.py'))
    args=p.parse_args();root=Path(__file__).resolve().parents[1]
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    checker=args.checker.resolve();fixture=root/'tests/native/test_io_ownership.f90'
    subprocess.run([sys.executable,str(checker),'--source-root',str(fixture.parent)],check=True)
    if args.field_close:
        field=args.finite_source or root/'runtime/sources/solver/src/maxwell/embedded_eh.f90'
        text=field.read_text();start=text.index('  integer(c_int) function sp_eh_close()')
        end=text.index('  end function',start)+len('  end function')
        source=out/'field_close_probe.f90'
        source.write_text("""module field_close_probe
use iso_c_binding
implicit none
type buffer_type
 real(8),allocatable::buffer(:)
end type
type fs_type
 integer::srg_ng=0
end type
type field_type
 type(buffer_type)::fe,ld
 type(fs_type)::fs
 integer::step=0,capacity=1,mpi_child_owner=7
end type
type(field_type),allocatable::state
type(buffer_type),allocatable::pending,step_saved
logical::step_restored=.false.
integer::finalized=0,nproc_group_global=0
contains
subroutine comm_resume_children(owner,status)
 integer,intent(in)::owner
 integer,intent(out)::status
 status=0
end subroutine
subroutine comm_get_max(status,group)
 integer,intent(inout)::status
 integer,intent(in)::group
end subroutine
subroutine comm_end_children(owner)
 integer,intent(in)::owner
end subroutine
subroutine comm_release_children(owner)
 integer,intent(in)::owner
 if(allocated(state))then
  if(state%fs%srg_ng/=0)stop 31
 endif
end subroutine
subroutine eh_finalize(fs,fe,ld)
 type(fs_type),intent(inout)::fs
 type(buffer_type),intent(inout)::fe,ld
 finalized=finalized+1
end subroutine
subroutine ld_finalize(ld)
 type(buffer_type),intent(inout)::ld
 if(allocated(ld%buffer))deallocate(ld%buffer)
end subroutine
subroutine dealloc_cache(value)
 integer,intent(inout)::value
 value=0
end subroutine
"""+text[start:end]+'\nend module\n')
        subprocess.run([sys.executable,str(checker),'--source-root',str(out)],check=True)
        for name,file,defines in [('helpers',source,[]),('probe',fixture,['-DSALMONPY_FIELD_CLOSE_PROBE'])]:
            subprocess.run([args.fc,'-cpp',*defines,'-I'+str(out),'-J'+str(out),'-c',str(file),'-o',str(out/(name+'.o'))],check=True)
        subprocess.run([sys.executable,str(checker),'--object-tree',str(out)],check=True)
        exe=out/'field-close-probe'
        subprocess.run([args.fc,str(out/'helpers.o'),str(out/'probe.o'),'-o',str(exe)],check=True)
        subprocess.run([str(exe)],cwd=out,check=True)
        (out/'result.json').write_text(json.dumps({'status':'passed','lifecycle_cases':8,
            'actual_close_body':True,'synthetic_allocations':True,'GS_data_collected':False})+'\n')
        return
    if args.scalar_inquiries:
        field=args.finite_source or root/'runtime/sources/solver/src/maxwell/embedded_eh.f90'
        adapter=root/'salmonpy/native_adapters/hse_checkpoint.f90'
        if not adapter.exists():adapter=root.parent/'native_adapters/hse_checkpoint.f90'
        field_text=field.read_text();history_text=adapter.read_text()
        for text in (field_text,history_text):
            if re.search(r'(?i)\b(?:all|any)\s*\(\s*ieee_is_(?:finite|nan)\s*\(',text):
                raise ValueError('array IEEE inquiry violates scalar policy')
        helpers=[]
        for rank in range(1,5):
            start=field_text.index('  subroutine accumulate_finite_'+str(rank)+'d(')
            end=field_text.index('  end subroutine',start)+len('  end subroutine')
            helpers.append(field_text[start:end])
        start=history_text.index(' logical function valid_history(')
        end=history_text.index(' end function',start)+len(' end function')
        generated=out/'scalar_ieee_probe.f90'
        generated.write_text('module scalar_ieee_probe\nuse exx_factor_history,only:s_factor_history\n'
            'use,intrinsic::ieee_arithmetic,only:ieee_is_finite\nimplicit none\ncontains\n'+
            '\n'.join(helpers)+'\n'+history_text[start:end]+'\nend module\n')
        subprocess.run([sys.executable,str(checker),'--source-root',str(out)],check=True)
        record=json.loads(args.build.read_text());lib=Path(record['library'])
        if hashlib.sha256(lib.read_bytes()).hexdigest()!=record['sha256']:raise ValueError('runtime hash changed')
        command=record.get('configure_command') or record['commands'][0]
        build=Path(command[command.index('-B')+1])
        for name,file,defines in [('helpers',generated,[]),('probe',fixture,['-DSALMONPY_SCALAR_IEEE_PROBE'])]:
            subprocess.run([args.fc,'-cpp',*defines,'-I'+str(build),'-I'+str(out),'-J'+str(out),
                '-c',str(file),'-o',str(out/(name+'.o'))],check=True)
        for objects in (out,build.parent):
            subprocess.run([sys.executable,str(checker),'--object-tree',str(objects)],check=True)
        exe=out/'scalar-probe'
        subprocess.run([args.fc,str(out/'helpers.o'),str(out/'probe.o'),str(lib),'-o',str(exe)],check=True)
        subprocess.run([str(exe)],cwd=out,check=True)
        (out/'result.json').write_text(json.dumps({'status':'passed','finite_and_nonfinite':True,
            'all_positions_and_ranks_1_to_4':True,'history_metadata_preserved':True,
            'Fujitsu_compiled':False,'GS_data_collected':False})+'\n')
        return
    if args.checkpoint_scratch:
        record=json.loads(args.build.read_text());lib=Path(record['library'])
        if hashlib.sha256(lib.read_bytes()).hexdigest()!=record['sha256']:raise ValueError('runtime hash changed')
        command=record.get('configure_command') or record['commands'][0]
        build=Path(command[command.index('-B')+1])
        adapter=root/'salmonpy/native_adapters/tddft_context.f90'
        if not adapter.exists():adapter=root.parent/'native_adapters/tddft_context.f90'
        def probe(text,directory,reference=False):
            directory.mkdir(exist_ok=False)
            payload=text[text.index('type :: checkpoint_payload'):text.index('end type',text.index('type :: checkpoint_payload'))+len('end type')]
            helpers=text[text.index('subroutine checkpoint_io('):text.rindex('end module')]
            generated=directory/'checkpoint_scope_probe.f90'
            legacy='integer(8)::checkpoint_array_limit=0_8\n' if reference else ''
            generated.write_text('module checkpoint_scope_probe\nuse structures\nimplicit none\n'+legacy+payload+'\ncontains\n'+helpers+'end module\n')
            test=directory/'probe.f90';test_text=fixture.read_text()
            if reference:
                test_text=re.sub(r'call (checkpoint_io\w*)\(([^\n]+),ios,(\d+)_8\)',
                    r'checkpoint_array_limit=\3_8;call \1(\2,ios)',test_text)
            test.write_text(test_text)
            subprocess.run([sys.executable,str(checker),'--source-root',str(directory)],check=True)
            for name,file,defines in [('helpers',generated,[]),('probe',test,['-DSALMONPY_CHECKPOINT_SCOPE_PROBE'])]:
                subprocess.run([args.fc,'-cpp',*defines,'-I'+str(build),'-I'+str(directory),'-J'+str(directory),
                    '-c',str(file),'-o',str(directory/(name+'.o'))],check=True)
            for objects in (directory,build.parent):
                subprocess.run([sys.executable,str(checker),'--object-tree',str(objects)],check=True)
            exe=directory/'checkpoint-probe'
            subprocess.run([args.fc,str(directory/'helpers.o'),str(directory/'probe.o'),str(lib),'-o',str(exe)],check=True)
            subprocess.run([str(exe)],cwd=directory,check=True)
        probe(adapter.read_text(),out/'current')
        if args.serializer_reference:
            probe(args.serializer_reference.read_text(),out/'reference',True)
            for name in ['small.bin','large.bin','vector.bin']:
                if (out/'current'/name).read_bytes()!=(out/'reference'/name).read_bytes():
                    raise ValueError('serializer bytes changed: '+name)
        (out/'result.json').write_text(json.dumps({'status':'passed','synthetic_only':True,
            'GS_data_collected':False,'small_large_small_exact':True,'oversize_rejected':True,
            'old_writer_bytes_exact':bool(args.serializer_reference)})+'\n')
        return
    if args.host_unit_library or args.host_timer_library or args.host_mpi_library:
        obj=out/'host-unit.o';lib=out/('host-unit'+('.dll' if sys.platform=='win32' else '.dylib' if sys.platform=='darwin' else '.so'))
        extra=[];libraries=[]
        if args.host_timer_library or args.host_mpi_library:
            record=json.loads(args.build.read_text());runtime=Path(record['library'])
            if hashlib.sha256(runtime.read_bytes()).hexdigest()!=record['sha256']:raise ValueError('runtime hash changed')
            command=record.get('configure_command') or record['commands'][0]
            build=Path(command[command.index('-B')+1])
            extra=['-DSALMONPY_HOST_MPI_PROBE' if args.host_mpi_library else '-DSALMONPY_HOST_TIMER_PROBE','-I'+str(build)]
            libraries=[str(runtime)]
            subprocess.run([sys.executable,str(checker),'--object-tree',str(build.parent)],check=True)
        subprocess.run([args.fc,'-cpp','-DSALMONPY_HOST_UNIT_PROBE',*extra,'-fPIC','-J'+str(out),
                        '-c',str(fixture),'-o',str(obj)],check=True)
        subprocess.run([sys.executable,str(checker),'--object-tree',str(out)],check=True)
        subprocess.run([args.fc,'-dynamiclib' if sys.platform=='darwin' else '-shared',str(obj),*libraries,'-o',str(lib)],check=True)
        (out/'result.json').write_text(json.dumps({'library':str(lib),
            'sha256':hashlib.sha256(lib.read_bytes()).hexdigest(),'runtime_sha256':record['sha256'] if args.host_timer_library or args.host_mpi_library else None,'status':'built-and-audited'})+'\n')
        return
    def compile_fixture(receipt,label,fixture=fixture):
        record=json.loads(receipt.read_text());library=Path(record['library'])
        if hashlib.sha256(library.read_bytes()).hexdigest()!=record['sha256']:raise ValueError('runtime hash changed')
        command=record.get('configure_command') or record['commands'][0]
        build=Path(command[command.index('-B')+1])
        directory=out/label;directory.mkdir();obj=directory/'fixture.o';exe=directory/'fixture'
        subprocess.run([args.fc,'-cpp','-I'+str(build),'-c',str(fixture),'-o',str(obj)],check=True)
        for objects in (directory,build.parent):
            subprocess.run([sys.executable,str(checker),'--object-tree',str(objects)],check=True)
        subprocess.run([args.fc,str(obj),str(library),'-o',str(exe)],check=True)
        return exe
    new=compile_fixture(args.build,'new')
    old=compile_fixture(args.reference_build,'old') if args.reference_build else None
    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
    def run(exe,mode,variant,label):
        work=out/(label+'-'+mode);work.mkdir()
        with (work/'run.log').open('w') as log:
            subprocess.run(['mpiexec','-n',str(args.ranks),str(exe),mode,variant],cwd=work,
                           env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
        if 'IO ownership PASS' not in (work/'run.log').read_text():raise ValueError('fixture incomplete')
        return work
    modes=('units','dos','pdos','band','avs','avs-unsplit')
    for mode in modes:
        occupied=run(new,mode,'occupied','new')
        if mode=='units':continue
        for f in occupied.glob('host-*'):
            if f.read_text()!='before\nafter\n':raise ValueError('host file modified')
        reference=run(old or new,mode,'reference','reference')
        def outputs(path):
            return {f.name:f.read_bytes() for f in path.iterdir() if f.is_file() and f.name not in ('run.log',) and not f.name.startswith('host-')}
        if outputs(occupied)!=outputs(reference):raise ValueError('output bytes changed: '+mode)
    if args.directory_scope:
        exe=compile_fixture(args.build,'directory-fixture',root/'tests/native/test_directory_context.f90')
        work=out/'directory-run';work.mkdir()
        with (work/'run.log').open('w') as log:
            subprocess.run(['mpiexec','-n',str(args.ranks),str(exe),str(work.resolve())],cwd=work,
                           env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
        if 'directory scope/host/absolute/relative/failure/sibling PASS' not in (work/'run.log').read_text():
            raise ValueError('directory fixture incomplete')
        if 'MPI resumed owners/late children/sibling/host PASS' not in (work/'run.log').read_text():
            raise ValueError('MPI owner fixture incomplete')
    (out/'result.json').write_text(json.dumps(dict(status='pass',ranks=args.ranks,
        synthetic_only=True,GS_data_collected=False,units_A_B_A=True,directory_scopes_checked=args.directory_scope,mpi_resumed_owners_checked=args.directory_scope,host_units_1_77_100_101_preserved=True,
        DOS_PDOS_band_AVS_outputs_exact=True,previous_runtime_outputs_exact=bool(old)),indent=2)+'\n')
    print('Synthetic IO ownership and output parity PASS')


if __name__=='__main__':main()
