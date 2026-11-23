"""Explicit GNU converter build tied to a trusted audited serial runtime record.

python tools/build_gs_restart_converter.py RUNTIME_BUILD_JSON NEW_OUTPUT
Uses the recorded runtime source and actual target flags; no solver rebuild/GS.
"""
import argparse,hashlib,json,os,re,shlex,shutil,subprocess,sys
from pathlib import Path
import numpy as np


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('runtime_build',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();rp=a.runtime_build.resolve();record=json.loads(rp.read_text())
    environment=dict(os.environ)
    if environment.get('GFORTRAN_CONVERT_UNIT','').lower() not in ('','native','swap','big_endian','little_endian'):
        raise ValueError('unit-specific native I/O conversion settings unsupported')
    if record.get('schema')!='salmonpy.runtime-build/1' or record.get('status')!='built-and-audited':raise ValueError('audited runtime build record required')
    engine=Path(record['executable']).resolve()
    if digest(engine)!=record['sha256']:raise ValueError('runtime hash changed')
    if record.get('build_configuration',{}).get('mpi_enabled') is not False:raise ValueError('initial converter requires an explicitly serial runtime')
    build=engine.parent;cache=build/'CMakeCache.txt'
    def cmake(key):
        m=re.search(r'^'+key+r':[^=]*=(.*)$',cache.read_text(),re.M)
        if not m:raise ValueError('missing compiler profile: '+key)
        return m.group(1).strip()
    fc=Path(cmake('CMAKE_Fortran_COMPILER')).resolve()
    version=subprocess.check_output([str(fc),'--version'],text=True)
    if 'GNU Fortran' not in version:raise ValueError('only GNU compiler-bound native record conversion validated')
    flags_path=build/'src/CMakeFiles/salmon.dir/flags.make'
    m=re.search(r'^Fortran_FLAGS\s*=\s*(.*)$',flags_path.read_text(),re.M)
    if not m:raise ValueError('actual runtime target flags are required (Makefiles generator)')
    flags=[];skip=False
    for token in shlex.split(m.group(1)):
        if skip:skip=False;continue
        if token in ('-J','-I'):skip=True;continue
        if token.startswith(('-J','-I')):continue
        flags.append(token)
    source=Path(record.get('source_root',rp.parent/'source')).resolve()
    root=Path(__file__).resolve().parents[1];package=root/'salmonpy';checker=root/'tools/check_forbidden_eigensolvers.py'
    if not (package/'native_adapters').is_dir():
        package=Path(__file__).resolve().parents[2];checker=Path(__file__).with_name('check_forbidden_eigensolvers.py')
    out=a.output.resolve()
    if any(out.is_relative_to(parent) for parent in (source,build,root,package)):
        raise ValueError('converter output must be outside source, package and fixed runtime build trees')
    out.mkdir(parents=True,exist_ok=False)
    commands=[];state={'schema':'salmonpy.gs-native-build/1','status':'started','target_executable':str(engine),'target_sha256':record['sha256']}
    (out/'build.json').write_text(json.dumps(state,indent=2)+'\n')
    with (out/'build.log').open('w') as log:
        def run(args,**kw):
            args=list(map(str,args));commands.append(args);log.write(json.dumps(args)+'\n');log.flush()
            subprocess.run(args,cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120,env=environment,**kw)
        try:
            # Audit full referenced SALMON source BEFORE extracting/compiling its procedures.
            run([sys.executable,checker,'--source-root',source/'src'])
            snapshot=out/'sources';snapshot.mkdir()
            for name in ('gs_gamma_converter.f90','gs_gamma_exporter.f90','gs_gamma_probe.f90'):shutil.copy2(package/'native_adapters'/name,snapshot/name)
            shutil.copy2(package/'include/salmonpy_gs_stream.f90',snapshot/'salmonpy_gs_stream.f90')
            native=source/'src/io/checkpoint_restart.f90';text=native.read_text()
            for name in ('read_wavefunction','write_wavefunction'):
                match=re.search(r'^subroutine '+name+r'\(.*?^end subroutine '+name+r'\s*$',text,re.M|re.S)
                if not match:raise ValueError('native wavefunction procedure not found')
                (snapshot/(name+'.f90')).write_text('! Extracted unchanged from audited external SALMON source.\n'+match.group(0)+'\n')
            run([sys.executable,checker,'--source-root',snapshot])
            opts=[*flags,'-cpp','-UUSE_MPI','-J'+str(out),'-I'+str(out)]
            for name in ('salmonpy_gs_stream','gs_gamma_converter','gs_gamma_exporter','gs_gamma_probe','read_wavefunction','write_wavefunction'):
                run([fc,*opts,'-c',snapshot/(name+'.f90'),'-o',out/(name+'.o')])
            run([sys.executable,checker,'--object-tree',out])
            converter=out/'gs_gamma_converter';exporter=out/'gs_gamma_exporter';probe=out/'native_probe'
            run([fc,*flags,out/'salmonpy_gs_stream.o',out/'gs_gamma_converter.o','-o',converter])
            run([fc,*flags,out/'salmonpy_gs_stream.o',out/'gs_gamma_exporter.o','-o',exporter])
            run([fc,*flags,*[out/(n+'.o') for n in ('gs_gamma_probe','read_wavefunction','write_wavefunction')],'-o',probe])
            converted=out/'synthetic-converted';reference=out/'synthetic-reference'
            (converted/'data_for_restart').mkdir(parents=True);reference.mkdir()
            (converted/'request.txt').write_text('8 8 8 2\n2\n')
            i=np.arange(1,513,dtype=float)[None,:];b=np.arange(1,3)[:,None]
            (converted/'portable_real.bin').write_bytes(np.asarray((i+1000*b)/1024,dtype='<f8').tobytes())
            (converted/'portable_imag.bin').write_bytes(np.asarray((i-1000*b)/2048,dtype='<f8').tobytes())
            (converted/'portable_occupations.bin').write_bytes(np.array([2.,0.],dtype='<f8').tobytes())
            # Retain complete commands/logs; native probe uses its own isolated working dirs.
            commands.append([str(converter)])
            subprocess.run([str(converter)],cwd=converted,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60,env=environment)
            run([probe,'write',str(reference)+os.sep])
            run([probe,'read',str(converted/'data_for_restart')+os.sep])
            for name in ('info.bin','occupation.bin','wfn.bin'):
                if (converted/'data_for_restart'/name).read_bytes()!=(reference/name).read_bytes():raise ValueError('native writer bytes differ: '+name)
            for kind in ('complex','real'):
                native_case=out/('native-export-'+kind);native_case.mkdir()
                (native_case/'native_input').mkdir()
                run([probe,'write' if kind=='complex' else 'write-real',str(native_case/'native_input')+os.sep])
                (native_case/'request.txt').write_text('8 8 8 2\n2\n0\n')
                commands.append([str(exporter)])
                subprocess.run([str(exporter)],cwd=native_case,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60,env=environment)
                expected_real=np.asarray((i+1000*b)/1024,dtype='<f8').tobytes()
                expected_imag=np.asarray((i-1000*b)/2048,dtype='<f8').tobytes() if kind=='complex' else bytes(8*512*2)
                if (native_case/'portable_real.bin').read_bytes()!=expected_real or (native_case/'portable_imag.bin').read_bytes()!=expected_imag:
                    raise ValueError('native export bits differ: '+kind)
                if (native_case/'portable_occupations.bin').read_bytes()!=np.array([2.,0.],dtype='<f8').tobytes():
                    raise ValueError('native occupation export differs')
            for nk in (2,8):
                case=out/('multi-k-'+str(nk));case.mkdir();(case/'native_input').mkdir()
                run([probe,'write',str(case/'native_input')+os.sep,str(nk)])
                (case/'request.txt').write_text(f'8 8 8 2\n2\n0\n{nk}\n')
                (case/'portable_weights.bin').write_bytes(np.full(nk,1./nk,dtype='<f8').tobytes())
                subprocess.run([str(exporter)],cwd=case,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60,env=environment)
                kk=np.arange(nk)[:,None,None];ii=np.arange(1,513)[None,None,:];bb=np.arange(1,3)[None,:,None]
                er=np.asarray((ii+1000*bb+10000*kk)/1024,dtype='<f8').tobytes()
                ei=np.asarray((ii-1000*bb-5000*kk)/2048,dtype='<f8').tobytes()
                if (case/'portable_real.bin').read_bytes()!=er or (case/'portable_imag.bin').read_bytes()!=ei:raise ValueError('multi-k export order differs')
                if (case/'portable_occupations.bin').read_bytes()!=np.tile(np.array([2.,0.],dtype='<f8'),nk).tobytes():raise ValueError('multi-k occupations differ')
                (case/'request.txt').write_text(f'8 8 8 2\n2\n{nk}\n');(case/'data_for_restart').mkdir()
                subprocess.run([str(converter)],cwd=case,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60,env=environment)
                run([probe,'read',str(case/'data_for_restart')+os.sep,str(nk)])
                for name in ('info.bin','occupation.bin','wfn.bin'):
                    if (case/'native_input'/name).read_bytes()!=(case/'data_for_restart'/name).read_bytes():raise ValueError('multi-k native round trip differs')
            weighted=out/'multi-k-weighted';weighted.mkdir();(weighted/'data_for_restart').mkdir()
            shutil.copy2(out/'multi-k-2/portable_real.bin',weighted/'portable_real.bin')
            shutil.copy2(out/'multi-k-2/portable_imag.bin',weighted/'portable_imag.bin')
            (weighted/'portable_weights.bin').write_bytes(np.array([.25,.75],dtype='<f8').tobytes())
            (weighted/'portable_occupations.bin').write_bytes(np.array([[2.,0.],[2.,2.]],dtype='<f8').tobytes())
            (weighted/'request.txt').write_text('8 8 8 2\n3.5\n2\n')
            subprocess.run([str(converter)],cwd=weighted,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60,env=environment)
            run([probe,'read-weighted',str(weighted/'data_for_restart')+os.sep,'2'])
            state.update(status='built-audited-and-probed',converter=str(converter),sha256=digest(converter),
                         exporter=str(exporter),exporter_sha256=digest(exporter),native_export_parity=True,
                         probe=str(probe),probe_sha256=digest(probe),compiler=str(fc),compiler_sha256=digest(fc),compiler_version=version,
                         flags=flags,flags_make_sha256=digest(flags_path),runtime_build_sha256=digest(rp),
                         native_source_sha256=digest(native),sources={p.name:digest(p) for p in snapshot.glob('*.f90')},
                         objects={p.name:digest(p) for p in out.glob('*.o')},native_byte_parity=True,multi_k_byte_parity=True,commands=commands)
        except BaseException as exc:
            state.update(status='failed',error=str(exc),commands=commands);raise
        finally:(out/'build.json').write_text(json.dumps(state,indent=2)+'\n')
    print(out/'build.json')

if __name__=='__main__':main()
