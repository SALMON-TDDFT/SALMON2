"""Explicit host-MPI TDDFT library build; never run on import.

python tools/build_embedded_tddft.py PREPARED_SOURCE NEW_DIRECTORY --jobs 2 \
  --cmake-arg=-DCMAKE_Fortran_COMPILER=mpifort --cmake-arg=-DCMAKE_C_COMPILER=mpicc
The caller supplies host BLAS/OpenMP/toolchain configuration as needed.
"""
import argparse,hashlib,json,os,re,shlex,subprocess,sys
from pathlib import Path


def link_library(build,output,source):
    checker=Path(__file__).with_name('check_forbidden_eigensolvers.py').resolve()
    subprocess.run([sys.executable,str(checker),'--source-root',str(source)],check=True)
    subprocess.run([sys.executable,str(checker),'--object-tree',str(build)],check=True)
    config=(build/'config.h').read_text()
    if not re.search(r'^#define USE_MPI\b',config,re.M):raise ValueError('MPI build required')
    for feature in ('USE_OPENACC','USE_SCALAPACK','USE_EIGENEXA','USE_FFTW','USE_HSE'):
        if re.search(r'^#define '+feature+r'\b',config,re.M):raise ValueError('unsupported embedded build: '+feature)
    command=shlex.split((build/'src/CMakeFiles/salmon.dir/link.txt').read_text())
    out_index=command.index('-o')
    library=output/('libsalmonpy_tddft.dylib' if sys.platform=='darwin' else 'libsalmonpy_tddft.so')
    command[out_index+1]=str(library)
    command=[v for v in command if v!='CMakeFiles/salmon.dir/main.f90.o']
    command.insert(1,'-shared')
    subprocess.run(command,cwd=build/'src',check=True)
    source_hashes={str(p.relative_to(source)):hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in sorted((source/'src').rglob('*')) if p.is_file()}
    record={'schema':'salmonpy.embedded-tddft/1','library':str(library),
            'sha256':hashlib.sha256(library.read_bytes()).hexdigest(),
            'source':json.loads((source/'salmonpy-runtime-source.json').read_text()),
            'source_hashes':source_hashes,'link_command':command,
            'scope':{'one_point_per_group':True,'one_lifetime_per_process':True,
                     'host_mpi':True,'eh_coupling_scheme':'salmonpy.eh-centered-current/1'}}
    (output/'build.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source');p.add_argument('directory');p.add_argument('--jobs',type=int,default=2)
    p.add_argument('--cmake-arg',action='append',default=[])
    a=p.parse_args();source=Path(a.source).resolve();out=Path(a.directory).resolve()
    if os.name!='posix':p.error('host library build currently supports macOS/Linux')
    if a.jobs<1:p.error('positive jobs required')
    if out.is_relative_to(source):p.error('build outside the prepared source')
    if not all((source/f).is_file() for f in ('src/rt/embedded_tddft.f90','src/maxwell/embedded_eh.f90')):
        p.error('prepare a runtime with embedded TDDFT and E/H support')
    out.mkdir(parents=True,exist_ok=False);build=out/'build'
    # CMake's executable PIC property emits -fPIE, which is insufficient when
    # these objects are subsequently linked into a Linux shared library.
    # Directory compile options follow that generated flag, so require -fPIC
    # for both languages without replacing caller optimization/toolchain flags.
    pic_options=out/'embedded-pic.cmake'
    pic_options.write_text('add_compile_options(-fPIC)\n')
    command=['cmake','-G','Unix Makefiles','-S',str(source),'-B',str(build),*a.cmake_arg,
             '-DCMAKE_PROJECT_SALMON_INCLUDE='+str(pic_options),
             '-DUSE_MPI=ON','-DCMAKE_POSITION_INDEPENDENT_CODE=ON']
    with (out/'build.log').open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
        subprocess.run(['cmake','--build',str(build),'-j',str(a.jobs)],stdout=log,stderr=subprocess.STDOUT,check=True)
    record=link_library(build,out,source)
    record['configure_command']=command
    (out/'build.json').write_text(json.dumps(record,indent=2)+'\n')
    print(out/'build.json')


if __name__=='__main__':main()
