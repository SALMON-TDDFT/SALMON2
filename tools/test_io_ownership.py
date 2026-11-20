"""Synthetic Fortran I/O ownership checks and optional old-output byte controls.

No SCF/GS export is performed. Optional ABI-state mode reads existing GS via
references and stores only ordinary RT currents/energies, never disk checkpoints.
"""
import argparse,hashlib,json,os,re,subprocess,sys
from pathlib import Path
from build_embedded_tddft import mpi_link_command


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
    p.add_argument('--abi-contract',action='store_true')
    p.add_argument('--legacy-settings-warning',action='store_true',help='check one-time HSE compatibility warning without GS creation')
    p.add_argument('--tbmbj-input',action='store_true',help='check public TB-mBJ input on synthetic density, without GS')
    p.add_argument('--parser-host',action='store_true',help='check status-only versus explicit native-host parser selection')
    p.add_argument('--abi-state',action='store_true')
    p.add_argument('--fixture',type=Path)
    p.add_argument('--state-output',type=Path)
    p.add_argument('--mpiexec',default='mpiexec')
    p.add_argument('--fast-hse',action='store_true')
    p.add_argument('--state-symmetry',action='store_true',help='exercise reduced-k HSE host lifecycle with symmetry-compatible z drive')
    p.add_argument('--cc',default='mpicc')
    p.add_argument('--serializer-reference',type=Path)
    p.add_argument('--fc',default='mpifort');p.add_argument('--ranks',type=int,default=2)
    p.add_argument('--checker',type=Path,default=Path(__file__).with_name('check_forbidden_eigensolvers.py'))
    args=p.parse_args()
    if args.state_symmetry and not (args.abi_state and args.fast_hse):
        p.error('state-symmetry is a fast-HSE ABI-state fixture option')
    if args.build.is_dir():args.build=args.build/'build.json'
    if not args.build.is_file():p.error(f'expected a build.json record or a directory containing it: {args.build}')
    root=Path(__file__).resolve().parents[1]
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    checker=args.checker.resolve();fixture=root/'tests/native/test_io_ownership.f90'
    subprocess.run([sys.executable,str(checker),'--source-root',str(fixture.parent)],check=True)
    if args.tbmbj_input:
        import numpy as np
        from importlib.resources import files
        from salmonpy.dft import DFT,Pseudopotential
        from dataclasses import replace
        record=json.loads(args.build.read_text());lib=Path(record['library'])
        digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
        if digest(lib)!=record['sha256']:raise ValueError('runtime hash changed')
        command=record.get('configure_command') or record['commands'][0]
        build=Path(command[command.index('-B')+1]);obj=out/'tbmbj.o';exe=out/'tbmbj-host'
        subprocess.run([args.fc,'-cpp','-DSALMONPY_TBMBJ_INPUT_PROBE','-I'+str(build),'-c',str(fixture),'-o',str(obj)],check=True)
        for objects in (out,build.parent):
            subprocess.run([sys.executable,str(checker),'--object-tree',str(objects)],check=True)
        subprocess.run(mpi_link_command([args.fc,str(obj),str(lib),'-o',str(exe)]),check=True)
        exe_sha=digest(exe)
        pseudo=Pseudopotential.from_file(Path(str(files('salmonpy.examples')/'si_intro_data/Si_rps.dat')),atomic_number=14,local_channel=2)
        model=DFT((14,),np.eye(3)*8,[[0,0,0]],(pseudo,),4,2,(8,8,8),(2,1,1),xc='tbmbj')
        for mode,coefficient in [('auto',None),('fixed',1.2)]:
            work=out/mode;work.mkdir()
            (work/'inputfile').write_text(replace(model,cval=coefficient).input_text())
            for name,data in model.input_assets:(work/name).write_bytes(data)
            if digest(exe)!=exe_sha or digest(lib)!=record['sha256']:raise ValueError('validated executable/runtime changed')
            with (work/'stdout.log').open('w') as log:
                subprocess.run([args.mpiexec,'-n','1',str(exe),mode],cwd=work,check=True,timeout=120,
                    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'},stdout=log,stderr=subprocess.STDOUT)
            assert 'TB-mBJ input and independent XC evaluation PASS' in (work/'stdout.log').read_text()
        (out/'result.json').write_text(json.dumps({'status':'passed','runtime_sha256':record['sha256'],
            'executable_sha256':exe_sha,'libxc_enabled':'#define USE_LIBXC' in (build/'config.h').read_text(),
            'automatic_and_fixed_input':True,'owned_XC_unchanged_by_global_settings':True,
            'GS_data_collected':False,'scope':'synthetic XC/parser contract; no SCF or physical-material validation'},indent=2)+'\n')
        return
    if args.legacy_settings_warning:
        record=json.loads(args.build.read_text());lib=Path(record['library'])
        if record.get('schema')!='salmonpy.fast-hse-mpi/1':p.error('legacy-settings-warning needs a fast-HSE build')
        digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
        if digest(lib)!=record['sha256']:raise ValueError('runtime hash changed')
        command=record.get('configure_command') or record['commands'][0]
        build=Path(command[command.index('-B')+1]);obj=out/'warning.o';exe=out/'warning-host'
        subprocess.run([args.fc,'-cpp','-DSALMONPY_HSE_LEGACY_WARNING_PROBE','-I'+str(lib.parent),'-I'+str(build),'-c',str(fixture),'-o',str(obj)],check=True)
        for objects in (out,build.parent):
            subprocess.run([sys.executable,str(checker),'--object-tree',str(objects)],check=True)
        link=mpi_link_command([args.fc,str(obj),str(lib),'-o',str(exe)])
        subprocess.run(link,check=True);exe_sha=digest(exe)
        for ranks in sorted({1,args.ranks}):
            if digest(exe)!=exe_sha or digest(lib)!=record['sha256']:raise ValueError('validated executable/runtime changed')
            result=subprocess.run([args.mpiexec,'-n',str(ranks),str(exe)],capture_output=True,text=True,timeout=120,
                env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
            (out/f'stdout-{ranks}.log').write_text(result.stdout)
            (out/f'stderr-{ranks}.log').write_text(result.stderr)
            result.check_returncode()
            assert result.stderr.count('legacy_settings fallback')==1
            assert result.stderr.index('explicit-owner-complete')<result.stderr.index('legacy_settings fallback')
            assert 'explicit owner attempted before any fallback' in result.stdout
        (out/'result.json').write_text(json.dumps({'status':'passed','runtime_sha256':record['sha256'],
            'executable_sha256':exe_sha,'ranks':sorted({1,args.ranks}),'warning_count_per_MPI_job':1,
            'solver_states_created':False},indent=2)+'\n')
        return
    if args.parser_host:
        record=json.loads(args.build.read_text());lib=Path(record['library'])
        digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
        if digest(lib)!=record['sha256']:raise ValueError('runtime hash changed')
        command=record.get('configure_command') or record['commands'][0]
        build=Path(command[command.index('-B')+1]);obj=out/'parser-host.o';exe=out/'parser-host'
        subprocess.run([args.fc,'-cpp','-DSALMONPY_PARSER_HOST_PROBE','-I'+str(build),'-c',str(fixture),'-o',str(obj)],check=True)
        for objects in (out,build.parent):
            subprocess.run([sys.executable,str(checker),'--object-tree',str(objects)],check=True)
        link=mpi_link_command([args.fc,str(obj),str(lib),'-o',str(exe)])
        subprocess.run(link,check=True);exe_sha=digest(exe)
        for ranks in sorted({1,args.ranks}):
            work=out/f'ranks-{ranks}';work.mkdir()
            (work/'inputfile').write_text("&calculation theory='maxwell' /\n&system yn_periodic='n' /\n&propagation propagator='user_extension' /\n")
            if digest(exe)!=exe_sha or digest(lib)!=record['sha256']:raise ValueError('validated executable/runtime changed')
            with (work/'stdout.log').open('w') as log:
                subprocess.run([args.mpiexec,'-n',str(ranks),str(exe)],cwd=work,check=True,timeout=120,
                    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'},stdout=log,stderr=subprocess.STDOUT)
            assert 'parser explicit native-host contract PASS' in (work/'stdout.log').read_text()
        (out/'result.json').write_text(json.dumps({'status':'passed','runtime_sha256':record['sha256'],
            'executable_sha256':exe_sha,'link_command':link,'ranks':sorted({1,args.ranks}),
            'status_only_keeps_electronic_checks':True,'explicit_false_keeps_electronic_checks':True,
            'native_host_true_selects_field_checks':True,'solver_states_created':False},indent=2)+'\n')
        return
    if args.abi_state:
        if args.fixture is None or args.state_output is None:
            p.error('abi-state needs an existing fixture and state-output')
        if args.ranks not in (1,2):p.error('this host fixture covers ranks 1/2')
        fixture_root=args.fixture.resolve();gs=(fixture_root/'data_for_restart').resolve()
        pp=(fixture_root/'pseudo1_rps.dat').resolve()
        if not gs.is_dir() or not pp.is_file():raise ValueError('existing GS and pseudopotential references required')
        text=(fixture_root/'rt.inp').read_text()
        record=json.loads(args.build.read_text());lib=Path(record['library'])
        expected='salmonpy.fast-hse-mpi/1' if args.fast_hse else 'salmonpy.embedded-tddft/1'
        if record.get('schema')!=expected:raise ValueError('build record does not match selected host fixture')
        def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
        if digest(lib)!=record['sha256']:raise ValueError('runtime hash changed')
        before={str(path.relative_to(gs)):digest(path) for path in gs.rglob('*') if path.is_file()}
        pp_before=digest(pp)
        symmetry=(fixture_root/'sym.dat').resolve() if args.state_symmetry else None
        symmetry_before=digest(symmetry) if symmetry is not None else None
        command=record.get('configure_command') or record['commands'][0]
        build=Path(command[command.index('-B')+1]);headers=root/'salmonpy/native_adapters'
        if not headers.exists():headers=root.parent/'native_adapters'
        executables={}
        for language,compiler,file,flags in [('c',args.cc,'test_native_state.c',['-I'+str(headers)]),
                                           ('fortran',args.fc,'test_native_state.f90',['-cpp','-J'+str(out)])]:
            if args.fast_hse:flags=[*flags,'-DSALMONPY_FAST_HSE']
            if args.state_symmetry:flags=[*flags,'-DSALMONPY_STATE_SYMMETRY=1']
            obj=out/(language+'.o');exe=out/('state-'+language)
            subprocess.run([compiler,*flags,'-c',str(fixture.parent/file),'-o',str(obj)],check=True)
            for objects in (out,build.parent):
                subprocess.run([sys.executable,str(checker),'--object-tree',str(objects)],check=True)
            link=mpi_link_command([compiler,str(obj),str(lib),'-Wl,-rpath,'+str(lib.parent),'-o',str(exe)])
            subprocess.run(link,check=True)
            executables[language]={'path':str(exe),'sha256':digest(exe),'link_command':link}
        (out/'executables.json').write_text(json.dumps(executables,indent=2)+'\n')
        results=args.state_output.resolve();results.mkdir(parents=True,exist_ok=False)
        import numpy as np
        env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'}
        cases=[];traces={};symmetry_counts=[]
        for ranks in sorted({1,args.ranks}):
            for language,point in [('c',0),('fortran',0),('c',1),('c',2)]:
                tag=f'{language}-n{ranks}-point{point}';work=results/tag;work.mkdir()
                (work/'data_for_restart').symlink_to(gs,target_is_directory=True)
                (work/'pseudo1_rps.dat').symlink_to(pp)
                if symmetry is not None:(work/'sym.dat').symlink_to(symmetry)
                input_text=text
                if args.fast_hse:
                    # This fixture exercises k-only groups; do not change supplied GS.
                    input_text,count=re.subn(r'(?i)\bnproc_k\s*=\s*\d+',f'nproc_k={ranks}',text)
                    if count!=1:raise ValueError('HSE host fixture needs one explicit nproc_k setting')
                (work/'rt.inp').write_text(input_text)
                exe=Path(executables[language]['path'])
                if digest(exe)!=executables[language]['sha256'] or digest(lib)!=record['sha256']:
                    raise ValueError('validated executable/runtime changed')
                with (work/'stdout.log').open('w') as log:
                    subprocess.run([args.mpiexec,'-n',str(ranks),str(exe),str(work),str(point)],cwd=work,
                                   env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=240)
                if args.state_symmetry:
                    observed=(work/'stdout.log').read_text()
                    blocks=re.findall(r'init_sym_kvector\(start\)(.*?)init_sym_kvector\(end\)',observed,re.S)
                    counts=[len(re.findall(r'^\s*\d+\s+[-+0-9.eEdD]+\s+',block,re.M)) for block in blocks]
                    if not counts or not all(counts):raise ValueError('reduced-k host did not report symmetry representatives')
                    symmetry_counts.extend(counts)
                data=np.loadtxt(work/'trace.tsv');assert data.shape==(24 if point==0 else 12,8)
                if not np.isfinite(data).all():raise ValueError('nonfinite host output')
                traces[(language,ranks,point)]=data;cases.append(tag);print(tag,'PASS',flush=True)
            paired=traces[('c',ranks,0)]
            np.testing.assert_allclose(paired,traces[('fortran',ranks,0)],rtol=2e-10,atol=2e-12)
            for point in (1,2):
                np.testing.assert_allclose(paired[paired[:,1]==point],traces[('c',ranks,point)],rtol=2e-10,atol=2e-12)
            if ranks!=1:
                np.testing.assert_allclose(paired,traces[('c',1,0)],rtol=2e-10,atol=2e-12)
        after={str(path.relative_to(gs)):digest(path) for path in gs.rglob('*') if path.is_file()}
        if before!=after or pp_before!=digest(pp):raise ValueError('referenced input was modified')
        if symmetry is not None and symmetry_before!=digest(symmetry):raise ValueError('referenced symmetry was modified')
        profile=('fast-PR HSE reduced-k' if args.state_symmetry else 'fast-PR HSE full-k') if args.fast_hse else 'standard PZ'
        summary={'status':'passed','profile':profile+' fixture','symmetry_enabled':args.state_symmetry,
                 'drive':'z (fixture retained symmetry)' if args.state_symmetry else 'xyz',
                 'observed_symmetry_k_counts':sorted(set(symmetry_counts)),'cases':cases,
                 'C_Fortran_currents_match':True,'independent_controls_match':True,
                 'bad_old_A_preserves_state':True,'stale_generation_rejected':True,
                 'fresh_state_matches_original_start':True,'MPI_caller_survives_detach':True,
                 'GS_bytes_unchanged':True,'GS_data_collected':False,'electronic_disk_checkpoints_created':False,
                 'runtime_sha256':record['sha256'],'executables':executables,'rtol':2e-10,'atol':2e-12,
                 'scope':profile+' state-lifecycle ABI; not full profile/physics/platform conformance'}
        (results/'result.json').write_text(json.dumps(summary,indent=2)+'\n')
        (out/'result.json').write_text(json.dumps(summary,indent=2)+'\n')
        return
    if args.abi_contract:
        record=json.loads(args.build.read_text());lib=Path(record['library'])
        if hashlib.sha256(lib.read_bytes()).hexdigest()!=record['sha256']:raise ValueError('runtime hash changed')
        command=record.get('configure_command') or record['commands'][0]
        build=Path(command[command.index('-B')+1])
        headers=root/'salmonpy/native_adapters'
        if not headers.exists():headers=root.parent/'native_adapters'
        defines=['-DSALMONPY_ABI_FAST_HSE'] if args.fast_hse else []
        outputs={}
        for language,compiler,source,flags in [
                ('c',args.cc,fixture.with_name('test_native_protocol.c'),['-I'+str(headers)]),
                ('fortran',args.fc,fixture,['-cpp','-DSALMONPY_ABI_PROTOCOL_PROBE','-J'+str(out)])]:
            obj=out/(language+'.o');exe=out/('probe-'+language)
            subprocess.run([compiler,*defines,*flags,'-c',str(source),'-o',str(obj)],check=True)
            for objects in (out,build.parent):
                subprocess.run([sys.executable,str(checker),'--object-tree',str(objects)],check=True)
            link=mpi_link_command([compiler,str(obj),str(lib),'-Wl,-rpath,'+str(lib.parent),'-o',str(exe)])
            subprocess.run(link,check=True)
            text=subprocess.check_output([str(exe)],cwd=out,text=True)
            outputs[language]=[int(value) for value in text.split()]
            if outputs[language]!=[1,1,1,1]:raise ValueError('protocol query mismatch')
        (out/'result.json').write_text(json.dumps({'status':'passed','profile':'fast_hse' if args.fast_hse else 'standard',
            'protocols':outputs,'runtime_sha256':record['sha256'],'MPI_initialized':False,
            'numerical_states_created':False,'GS_data_collected':False,'scope':'query-only C/Fortran ABI; not full state/step conformance'},indent=2)+'\n')
        return
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
        subprocess.run(mpi_link_command([args.fc,str(out/'helpers.o'),str(out/'probe.o'),'-o',str(exe)]),check=True)
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
        subprocess.run(mpi_link_command([args.fc,str(out/'helpers.o'),str(out/'probe.o'),str(lib),'-o',str(exe)]),check=True)
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
            subprocess.run(mpi_link_command([args.fc,str(directory/'helpers.o'),str(directory/'probe.o'),str(lib),'-o',str(exe)]),check=True)
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
        subprocess.run(mpi_link_command([args.fc,str(obj),*libraries,'-o',str(lib)],shared=True),check=True)
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
        subprocess.run(mpi_link_command([args.fc,str(obj),str(library),'-o',str(exe)]),check=True)
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
