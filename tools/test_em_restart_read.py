"""Explicit GS-free restart-reader failures/repair under MPI1/2.

Requires a built, source/object-audited runtime. Checks missing/truncated/invalid
text, missing/truncated real/complex payloads/metadata, and remap transfers.
No production numerical kernel is implemented in Python.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from build_embedded_tddft import mpi_link_command


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--batch',action='store_true',help='also run local LD+FDTD midpoint/resume and damaged archive tests')
    p.add_argument('--reference-executable',type=Path)
    p.add_argument('--reference-sha256')
    p.add_argument('--fc',default='mpifort')
    p.add_argument('--launcher',default='mpiexec')
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    record=json.loads(a.build.read_text());lib=Path(record['library'])
    assert hashlib.sha256(lib.read_bytes()).hexdigest()==record['sha256']
    command=record.get('configure_command') or record['commands'][0]
    build=Path(command[command.index('-B')+1])
    fixture=root/'tests/native/test_em_restart_read.f90'
    checker=root/'tools/check_forbidden_eigensolvers.py'
    subprocess.run([sys.executable,str(checker),'--source-root',str(fixture.parent)],check=True)
    obj=out/'probe.o';exe=out/'probe'
    subprocess.run([a.fc,'-I'+str(build),'-c',str(fixture),'-o',str(obj)],check=True)
    subprocess.run([sys.executable,str(checker),'--object-tree',str(out)],check=True)
    subprocess.run(mpi_link_command([a.fc,str(obj),str(lib),'-o',str(exe)]),check=True)
    cases=[]

    def run(work,ranks,mode,target='',old=None):
        with (work/'run.log').open('w') as log:
            subprocess.run([a.launcher,'-n',str(ranks),str(exe),mode,target,str(old or ranks)],
                cwd=work,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),
                stdout=log,stderr=subprocess.STDOUT,check=True,timeout=45)
        assert 'caller survival and retry PASS' in (work/'run.log').read_text()
        assert all((work/f'host-{rank}').read_text()=='before\nafter\n' for rank in range(ranks))

    for ranks in (1,2):
        seed=out/f'seed-{ranks}';seed.mkdir();(seed/'data_for_restart').mkdir()
        run(seed,ranks,'prepare')
        for mode in ('text','single-real','single-complex','local-real','local-complex','remap-real','remap-complex'):
            defects=['missing','truncated','directory']+(['invalid-text'] if mode=='text' else [])
            members=['payload','metadata'] if mode.startswith(('local','remap')) else ['payload']
            for member in members:
                for defect in defects:
                    name=f'{mode}-{member}-{defect}-{ranks}';d=out/name;d.mkdir()
                    shutil.copytree(seed/'data_for_restart',d/'restart')
                    kind='complex' if 'complex' in mode else 'real'
                    target=('value.data' if mode=='text' else kind+'.bin' if mode.startswith('single') else
                            f'{kind}_{ranks-1}'+('_inf.bin' if member=='metadata' else '.bin'))
                    path=d/'restart'/target;original=path.read_bytes()
                    Path(str(path)+'.repair').write_bytes(original)
                    if defect=='missing':path.unlink()
                    elif defect=='directory':path.unlink();path.mkdir()
                    elif defect=='invalid-text':path.write_text('not-an-integer\n')
                    else:path.write_bytes(original[:1])
                    run(d,ranks,mode,'restart/'+target)
                    assert path.read_bytes()==original
                    statuses=[int((d/f'failure-status-{rank}').read_text()) for rank in range(ranks)]
                    assert len(set(statuses))==1 and statuses[0]>0
                    cases.append({'case':name,'failure_status':statuses[0],'repaired_retry':True,'caller_MPI_and_unit_retained':True})
        readonly=out/f'readonly-archive-{ranks}';readonly.mkdir()
        shutil.copytree(seed/'data_for_restart',readonly/'restart')
        for file in (readonly/'restart').iterdir():file.chmod(0o444)
        (readonly/'restart').chmod(0o555)
        run(readonly,ranks,'remap-complex')
        cases.append({'case':readonly.name,'read_only_archive_accepted':True})
        # Changed-size normal remap: both directions use the same prepared values.
        changed=3-ranks;d=out/f'remap-{ranks}-to-{changed}';d.mkdir()
        shutil.copytree(seed/'data_for_restart',d/'restart')
        run(d,changed,'remap-real',old=ranks)
        cases.append({'case':d.name,'normal_remap':True})
    if a.batch:
        from salmonpy.fdtd import FDTD,PlaneWave
        from salmonpy.geometry import Grid,Geometry,Shape
        from salmonpy.fields import AcosPulse
        from salmonpy.ld import LorentzDrude
        grid=Grid((40,40,40),(2.,2.,2.))
        calc=FDTD(Geometry(grid,(Shape('box',(2.,4.,4.),(7.,1.,1.),material=1),)),
            (LorentzDrude(plasma_frequency=20.,strength=[1.],resonance=[0.],damping=[2.]),),
            .003,128,source=PlaneWave(AcosPulse(.001,20.,.15),0,-16.),
            observations=((1.,1.,1.),),observation_stride=8)
        base=calc.input_text().replace("sysname = 'fdtd'","sysname = 'fdtd', checkpoint_interval=64")
        def launch(work,ranks,mode='batch-full',binary=exe):
            with (work/'inputfile').open() as inp,(work/'run.log').open('w') as log:
                q=subprocess.run([a.launcher,'-n',str(ranks),str(binary)]+([mode] if binary==exe else []),
                    cwd=work,stdin=inp,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),
                    stdout=log,stderr=subprocess.STDOUT,timeout=90)
            return q.returncode
        def text(ranks,restart=None):
            result=base+f'\n&parallel\n nproc_k=1,nproc_ob=1,nproc_rgrid=1,{ranks},1\n/\n'
            if restart:result=result.replace('&control',"&control\n yn_restart='y', directory_read_data='restart'")
            return result
        def equal(work,reference):
            count=0;whitespace=[]
            for f in reference.iterdir():
                if f.suffix=='.data' or f.name=='shape.cube':
                    actual=(work/f.name).read_bytes();expected=f.read_bytes()
                    if actual!=expected and f.name.endswith('_at_point_rt.data'):
                        # Existing txtfile_copy_em trims old rows on restart.
                        assert [v.rstrip() for v in actual.splitlines()]==[v.rstrip() for v in expected.splitlines()],str(f)
                        whitespace.append(f.name)
                    else:
                        assert actual==expected,str(f);count+=1
            for f in (reference/'data_for_restart').iterdir():
                assert (work/'data_for_restart'/f.name).read_bytes()==f.read_bytes(),str(f);count+=1
            return {'byte_equal':count,'legacy_point_copy_trailing_spaces_only':whitespace}
        controls={};halves={}
        for ranks in (1,2):
            full=out/f'batch-full-{ranks}';full.mkdir();(full/'inputfile').write_text(text(ranks))
            assert launch(full,ranks)==0;controls[ranks]=full
            if a.reference_executable:
                assert a.reference_sha256 and hashlib.sha256(a.reference_executable.read_bytes()).hexdigest()==a.reference_sha256
                reference=out/f'batch-prior-{ranks}';reference.mkdir();(reference/'inputfile').write_text(text(ranks))
                assert launch(reference,ranks,binary=a.reference_executable)==0
                cases.append({'case':reference.name,'prior_normal_files_byte_equal':equal(full,reference)})
            half=out/f'batch-half-{ranks}';half.mkdir();(half/'inputfile').write_text(text(ranks))
            assert launch(half,ranks,mode='batch-half')==0;halves[ranks]=half
            assert int((half/'data_for_restart/iter_now.data').read_text())==64
        for old in (1,2):
            for ranks in (1,2):
                d=out/f'batch-resume-{old}-to-{ranks}';d.mkdir()
                shutil.copytree(halves[old]/'data_for_restart',d/'restart')
                before={str(f.relative_to(d/'restart')):hashlib.sha256(f.read_bytes()).hexdigest() for f in (d/'restart').iterdir()}
                (d/'inputfile').write_text(text(ranks,restart=True))
                assert launch(d,ranks)==0
                assert all(hashlib.sha256((d/'restart'/name).read_bytes()).hexdigest()==sha for name,sha in before.items())
                cases.append({'case':d.name,'continuous_files_byte_equal':equal(d,controls[ranks]),'input_archive_unchanged':True})
        for ranks in (1,2):
            for defect in ('missing-text','invalid-text','missing-info','truncated-info','missing-payload','truncated-payload'):
                d=out/f'batch-{defect}-{ranks}';d.mkdir();shutil.copytree(halves[ranks]/'data_for_restart',d/'restart')
                name='nproc_size_global.data' if 'text' in defect else f'ex_y_{ranks-1}'+('_inf.bin' if 'info' in defect else '.bin')
                target=d/'restart'/name
                if defect.startswith('missing'):target.rename(d/('retained-'+name))
                elif defect=='invalid-text':target.write_text('not-an-integer\n')
                else:target.write_bytes(target.read_bytes()[:1])
                before={str(f.relative_to(d/'restart')):hashlib.sha256(f.read_bytes()).hexdigest() for f in (d/'restart').iterdir()}
                (d/'inputfile').write_text(text(ranks,restart=True));code=launch(d,ranks)
                assert code!=0 and 'SALMONPY_RESTART_READ_ERROR' in (d/'run.log').read_text()
                assert all(hashlib.sha256((d/'restart'/name).read_bytes()).hexdigest()==sha for name,sha in before.items())
                cases.append({'case':d.name,'exit_code':code,'no_timeout':True,'input_archive_unchanged':True})
        # A failed repartition read must not overwrite metadata in its input set.
        for old,ranks in ((1,2),(2,1)):
            for broken in (False,True):
                d=out/f'batch-shared-path-{old}-to-{ranks}-{"broken" if broken else "normal"}'
                d.mkdir();shutil.copytree(halves[old]/'data_for_restart',d/'data_for_restart')
                archive=d/'data_for_restart'
                if broken:(archive/f'ex_y_{old-1}.bin').rename(d/f'retained-ex_y_{old-1}.bin')
                before={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in archive.iterdir()}
                inp=text(ranks,restart=True).replace("directory_read_data='restart'","directory_read_data='data_for_restart'")
                (d/'inputfile').write_text(inp);code=launch(d,ranks)
                if broken:
                    assert code!=0 and 'SALMONPY_RESTART_READ_ERROR' in (d/'run.log').read_text()
                    after={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in archive.iterdir()}
                    assert before==after,'failed shared-path read modified input archive'
                    cases.append({'case':d.name,'exit_code':code,'input_archive_unchanged':True})
                else:
                    assert code==0
                    cases.append({'case':d.name,'continuous_files_byte_equal':equal(d,controls[ranks])})
    (out/'summary.json').write_text(json.dumps({'cases':cases,'runtime_sha256':record['sha256'],
        'probe_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'GS_collected':False,
        'scope':'GNU MPI1/2 single-node helper I/O failure and explicit repaired retry; failed buffers are not published state'},indent=2)+'\n')
    print('PASS',len(cases),'restart-reader cases')


if __name__=='__main__':
    main()
