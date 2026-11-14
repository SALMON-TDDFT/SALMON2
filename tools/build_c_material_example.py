"""Explicit audited build of the manufactured C plugin and C++/Fortran hosts.

No SALMON solver source/build or GS data. Existing destinations are never reused.
Usage: python tools/build_c_material_example.py NEW_DIRECTORY
"""
import argparse,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
import time


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path)
    p.add_argument('--cc',default='cc');p.add_argument('--cxx',default='c++');p.add_argument('--fc',default='gfortran')
    a=p.parse_args();root=Path(__file__).resolve().parents[1]
    include=root/'salmonpy/include';source=root/'salmonpy/examples/c_material'
    # Installed tools live under salmonpy/_runtime_source/tools.
    if not source.is_dir():
        package=Path(__file__).resolve().parents[2];include=package/'include';source=package/'examples/c_material'
    checker=root/'tools/check_forbidden_eigensolvers.py'
    if not checker.is_file():checker=Path(__file__).with_name('check_forbidden_eigensolvers.py')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    # Preserve exactly the sources used by this build; never rewrite old builds.
    snapshot=out/'sources';snapshot.mkdir()
    shutil.copytree(include,snapshot/'include')
    shutil.copytree(source,snapshot/'example')
    include=snapshot/'include';source=snapshot/'example'
    commands=[]
    with (out/'build.log').open('w') as log:
        def run(args):
            commands.append(args);log.write(json.dumps(args)+'\n');log.flush()
            subprocess.run(args,cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
        run([sys.executable,str(checker),'--source-root',str(source)])
        run([sys.executable,str(checker),'--source-root',str(include)])
        run([a.cc,'-std=c11','-O2','-ffp-contract=off','-fPIC','-I',str(include),'-c',str(source/'memory.c'),'-o','memory.o'])
        run([a.cxx,'-std=c++11','-O2','-ffp-contract=off','-I',str(include),'-c',str(source/'host.cpp'),'-o','host_cpp.o'])
        run([a.fc,'-O2','-ffp-contract=off','-J',str(out),'-I',str(out),'-c',str(include/'salmonpy_material_iso.f90'),'-o','iso.o'])
        run([a.fc,'-O2','-ffp-contract=off','-I',str(out),'-c',str(source/'host.f90'),'-o','host_fortran.o'])
        run([sys.executable,str(checker),'--object-tree',str(out)])
        suffix='.dylib' if sys.platform=='darwin' else '.so'
        library=out/('libexample_material'+suffix)
        run([a.cc,'-dynamiclib' if sys.platform=='darwin' else '-shared',str(out/'memory.o'),'-lm','-o',str(library)])
        rpath=['-Wl,-rpath,'+str(out)]
        run([a.cxx,str(out/'host_cpp.o'),str(library),*rpath,'-o',str(out/'host_cpp')])
        run([a.fc,str(out/'iso.o'),str(out/'host_fortran.o'),str(library),*rpath,'-o',str(out/'host_fortran')])
        run([str(out/'host_cpp')]);run([str(out/'host_fortran')])
    record={'schema':'salmonpy.c-material-example/1','status':'built-audited-and-host-tested',
            'library':str(library),'sha256':hashlib.sha256(library.read_bytes()).hexdigest(),
            'commands':commands,'objects':{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in out.glob('*.o')},
            'sources':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [source/'memory.c',source/'host.cpp',source/'host.f90',include/'salmonpy_material.h',include/'salmonpy_material_iso.f90']}}
    (out/'build.json').write_text(json.dumps(record,indent=2)+'\n');print(out/'build.json')

if __name__=='__main__':main()
