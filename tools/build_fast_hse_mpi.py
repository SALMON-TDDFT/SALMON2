"""Explicit host-MPI fast-PR HSE library build from a new prepared source.

PREPARED_FAST_SOURCE NEW_OUTPUT --jobs 2 --cmake-arg=... . Never runs on import.
Only the pinned explicit-routing fast-PR engine and host MPI are accepted.
"""
import argparse,hashlib,json,os,re,shlex,shutil,subprocess,sys
from pathlib import Path

FAST_COMMIT='c9e8ce431f43eefd4a88998134adb4796f8ec837'

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--jobs',type=int,default=2);p.add_argument('--cmake-arg',action='append',default=[])
    a=p.parse_args();root=Path(__file__).resolve().parents[1];source=a.source.resolve();out=a.output.resolve()
    marker_path=source/'salmonpy-hse-rt-state.json';marker=json.loads(marker_path.read_text())
    if (marker.get('upstream_commit')!=FAST_COMMIT or not marker.get('explicit_rt_routing')
        or marker.get('status')!='prepared-and-source-audited'):
        raise ValueError('new prepared explicit-routing fast-PR source required')
    if os.name!='posix' or a.jobs<1:raise ValueError('POSIX host and positive build workers required')
    if out.is_relative_to(source) or out.is_relative_to(root):raise ValueError('new external output directory required')
    checker=root/'tools/check_forbidden_eigensolvers.py'
    subprocess.run([sys.executable,str(checker),'--source-root',str(source/'src')],check=True)
    out.mkdir(parents=True,exist_ok=False);build=out/'build';staged=out/'src';staged.mkdir()
    pic=out/'pic.cmake';pic.write_text('add_compile_options(-fPIC)\n')
    commands=[]
    with (out/'build.log').open('w') as log:
        def run(cmd,**kw):
            cmd=list(map(str,cmd));commands.append(cmd)
            subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,**kw)
        run(['cmake','-G','Unix Makefiles','-S',source,'-B',build,*a.cmake_arg,
             '-DCMAKE_PROJECT_SALMON_INCLUDE='+str(pic),'-DCMAKE_POSITION_INDEPENDENT_CODE=ON',
             '-DUSE_HSE=ON','-DUSE_MPI=ON'])
        run(['cmake','--build',build,'-j',str(a.jobs)])
        config=(build/'config.h').read_text()
        if not all(re.search(r'^#define '+f+r'\b',config,re.M) for f in ('USE_HSE','USE_MPI')):
            raise ValueError('MPI fast HSE build required')
        for f in ('USE_OPENACC','USE_FFTW','USE_SCALAPACK','USE_EIGENEXA'):
            if re.search(r'^#define '+f+r'\b',config,re.M):raise ValueError('unsupported profile: '+f)
        adapters=root/'salmonpy/native_adapters'
        if not adapters.exists():adapters=root.parent/'native_adapters'
        names=['hse_model_inputs','hse_checkpoint','tddft_context','fast_hse_mpi','parallel_bridge']
        if (build/'src/CMakeFiles/salmon.dir/salmonpy/parallel_bridge.f90.o').is_file():
            names.remove('parallel_bridge')
        if not (build/'src/CMakeFiles/salmon.dir/io/step_output.f90.o').is_file():
            names.insert(1,'step_output')
        for name in names:
            path=(root/'runtime/sources/standard/src/salmonpy/parallel_bridge.f90' if name=='parallel_bridge'
                  else root/'runtime/sources/standard/src/io/step_output.f90' if name=='step_output'
                  else adapters/(name+'.f90'))
            shutil.copy2(path,staged/(name+'.f90'))
        shutil.copy2(adapters/'fast_hse_handles.h',staged/'fast_hse_handles.h')
        flags=(build/'src/CMakeFiles/salmon.dir/flags.make').read_text()
        options=[x for x in shlex.split(re.search(r'^Fortran_FLAGS\s*=\s*(.*)$',flags,re.M)[1]) if not x.startswith('-J')]
        includes=shlex.split(re.search(r'^Fortran_INCLUDES\s*=\s*(.*)$',flags,re.M)[1])
        fc=re.search(r'^CMAKE_Fortran_COMPILER:[^=]*=(.*)$',(build/'CMakeCache.txt').read_text(),re.M)[1].strip()
        run([sys.executable,checker,'--source-root',staged])
        for name in names:
            run([fc,*options,'-fPIC','-I'+str(out),*includes,'-J'+str(out),'-c',staged/(name+'.f90'),'-o',out/(name+'.o')],cwd=out)
        run([sys.executable,checker,'--object-tree',build]);run([sys.executable,checker,'--object-tree',out])
        link=shlex.split((build/'src/CMakeFiles/salmon.dir/link.txt').read_text())
        library=out/('libsalmonpy_fast_hse.dylib' if sys.platform=='darwin' else 'libsalmonpy_fast_hse.so')
        link[link.index('-o')+1]=str(library)
        link=[x for x in link if x!='CMakeFiles/salmon.dir/main.f90.o']
        link[1:1]=['-shared',*[str(out/(n+'.o')) for n in names]]
        run(link,cwd=build/'src')
    r={'schema':'salmonpy.fast-hse-mpi/1','status':'built-and-audited','library':str(library),
       'sha256':digest(library),'origin_sha256':digest(build/'salmon'),'source':marker,
       'source_marker_sha256':digest(marker_path),'source_hashes':{str(x.relative_to(source)):digest(x) for x in sorted((source/'src').rglob('*')) if x.is_file()},
       'adapter_hashes':{x.name:digest(x) for x in staged.iterdir()},'commands':commands,
       'scope':'host MPI fast-PR k-only HSE; legacy Python one lifetime; shared-model generation-checked C contexts; same-topology rank-local checkpoint',
       'handle_abi':{'protocol':1,'header':'fast_hse_handles.h','shared_immutable_model':True,'main_thread_collective':True}}
    # One build can serve field and material rank groups. Keep the executable
    # fingerprint separate from the host-library fingerprint.
    if (build/'src/CMakeFiles/salmon.dir/maxwell/embedded_eh.f90.o').is_file():
        import ctypes as c
        probe=c.CDLL(str(library))
        probe.sp_eh_checkpoint_protocol.argtypes=[]
        probe.sp_eh_checkpoint_protocol.restype=c.c_int
        if probe.sp_eh_checkpoint_protocol()!=1:raise ValueError('unsupported compiled field protocol')
        r['field_abi']={'protocol':1,'scheme':'salmonpy.eh-centered-current/1'}
    if (build/'src/CMakeFiles/salmon.dir/gs/main_dft.f90.o').is_file():
        import ctypes as c
        probe=c.CDLL(str(library))
        protocol=getattr(probe,'sp_gs_protocol',None)
        if protocol is not None:
            protocol.argtypes=[];protocol.restype=c.c_int
            if protocol()!=1:raise ValueError('unsupported compiled GS protocol')
            r['gs_abi']={'protocol':1,'modes':['import'],'arrays':['density','occupations','orbitals']}
            memory=getattr(probe,'sp_fast_hse_memory_gs_protocol',None)
            if memory is not None:
                memory.argtypes=[];memory.restype=c.c_int
                if memory()!=1:raise ValueError('unsupported compiled memory GS protocol')
                r['gs_abi']['memory_rt']=1
    executable_record={'schema':'salmonpy.runtime-build/1','status':'built-and-audited',
        'executable':str(build/'salmon'),'sha256':r['origin_sha256'],
        'build_configuration':{'mpi_enabled':True,'hse_enabled':True},
        'source':marker,'source_marker_sha256':r['source_marker_sha256']}
    (out/'executable-build.json').write_text(json.dumps(executable_record,indent=2)+'\n')
    (out/'build.json').write_text(json.dumps(r,indent=2)+'\n');print(out/'build.json')


if __name__=='__main__':main()
