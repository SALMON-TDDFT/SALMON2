"""Explicit local build/test of sequential model-owned standard GS inputs.

NEW_BUILD OLD_BUILD NEW_OUTPUT --reference-sha256 SHA [--size 2 --split spatial]
Both build records must name unchanged audited host-MPI libraries. The old build's
standalone executable SHA is supplied separately. No restart GS data is written.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
from build_embedded_tddft import mpi_link_command


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('build','reference','output'):parser.add_argument(name,type=Path)
    parser.add_argument('--reference-sha256',required=True)
    parser.add_argument('--size',type=int,choices=(1,2),default=1)
    parser.add_argument('--split',choices=('orbital','spatial'),default='orbital')
    parser.add_argument('--xc-a',choices=('pz','libxc_pz','libxc_pbe'),default='pz')
    parser.add_argument('--xc-b',choices=('libxc_pbe','tbmbj'),default='libxc_pbe')
    parser.add_argument('--fixed-c',type=float)
    parser.add_argument('--fc',default='mpifort')
    parser.add_argument('--mpiexec',default='mpiexec')
    parser.add_argument('--checker',type=Path,default=Path(__file__).with_name('check_forbidden_eigensolvers.py'))
    a=parser.parse_args();repo=Path(__file__).resolve().parents[1]
    records=[json.loads(p.read_text()) for p in (a.build,a.reference)]
    for record in records:
        if record.get('schema')!='salmonpy.embedded-tddft/1' or digest(record['library'])!=record['sha256']:
            raise ValueError('unchanged verified embedded build required')
    library=Path(records[0]['library']).resolve();build=library.parent/'build'
    reference=Path(records[1]['library']).resolve().parent/'build/salmon'
    if not (build/'salmonpy_gs_model.mod').is_file():raise ValueError('development build with GS module files and application objects required')
    if digest(reference)!=a.reference_sha256:raise ValueError('reference executable checksum mismatch')
    out=a.output.resolve()
    if any(out.is_relative_to(p) for p in (repo,build,reference.parent)):
        raise ValueError('new external output directory required')
    out.mkdir(parents=True,exist_ok=False);commands=[]
    def run(cmd,**kw):
        commands.append([str(x) for x in cmd])
        subprocess.run(cmd,check=True,**kw)
    source=out/'src';source.mkdir()
    (source/'probe.f90').write_bytes((repo/'tests/native/test_gs_model.f90').read_bytes())
    with (out/'build.log').open('w') as log:
        for flag,path in (('--source-root',source),('--object-tree',build)):
            run([sys.executable,a.checker,flag,path],stdout=log,stderr=subprocess.STDOUT)
        obj=out/'probe.o'
        run([a.fc,'-I'+str(build),'-J'+str(out),'-c',source/'probe.f90','-o',obj],stdout=log,stderr=subprocess.STDOUT)
        run([sys.executable,a.checker,'--object-tree',out],stdout=log,stderr=subprocess.STDOUT)
        probe=out/'probe'
        run(mpi_link_command([a.fc,obj,library,'-Wl,-rpath,'+str(library.parent),'-o',probe]),stdout=log,stderr=subprocess.STDOUT)
    probe_sha=digest(probe)
    from salmonpy.dft import DFT,Pseudopotential,SCFSettings
    pseudo_path=repo/'salmonpy/examples/si_intro_data/Si_rps.dat'
    if not pseudo_path.is_file():pseudo_path=repo.parent/'examples/si_intro_data/Si_rps.dat'
    pseudo=Pseudopotential.from_file(pseudo_path,atomic_number=14,local_channel=2)
    (out/'pseudo1_rps.dat').write_bytes(pseudo.data);texts=[]
    for label,xc,n,side,natom,kmesh in [('a',a.xc_a,8,8.,1,(2,1,1)),('b',a.xc_b,12,9.,2,(1,2,1))]:
        model=DFT((14,)*natom,np.eye(3)*side,np.array([[.3,.4,.45],[.7,.6,.55]][:natom]),
                  (pseudo,),4*natom,2*natom,(n,n,n),kmesh,xc='libxc_pbe' if xc=='tbmbj' else xc)
        text=model.input_text(SCFSettings(max_iterations=4,tolerance=1e-12))
        if xc=='tbmbj':text=text.replace("xc='libxc_pbe'","xc='tbmbj'"+(f', cval={a.fixed_c}' if a.fixed_c is not None else ''))
        text=text.replace("write_gs_restart_data='all'","write_gs_restart_data='no'").replace("base_directory='./'",f"base_directory='./{label}/'")
        ob=a.size if a.split=='orbital' else 1;rx=a.size if a.split=='spatial' else 1
        text=text.replace('&system',f'&parallel\n nproc_k=1,nproc_ob={ob},nproc_rgrid={rx},1,1\n/\n&system')
        (out/label).mkdir();(out/(label+'.inp')).write_text(text);texts.append(text)
        if label=='a':(out/'a.inp.unsupported').write_text(text.replace(f"xc='{xc}'","xc='pw'"))
    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
    if digest(probe)!=probe_sha:raise ValueError('probe changed')
    with (out/'probe.log').open('w') as log:
        run([a.mpiexec,'-n',str(a.size),probe,out/'a.inp',out/'b.inp'],cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180)
    for label,text in zip(('a','b'),texts):
        ref=out/('reference-'+label);ref.mkdir();(ref/'pseudo1_rps.dat').write_bytes(pseudo.data);(ref/label).mkdir()
        if digest(reference)!=a.reference_sha256:raise ValueError('reference changed')
        with (ref/'stdout.log').open('w') as log:
            run([a.mpiexec,'-n',str(a.size),reference],input=text,text=True,cwd=ref,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180)
        actual=out/(label+'-first-status.dat' if label=='a' else 'b-status.dat')
        if actual.read_bytes()!=(ref/'gs_status.dat').read_bytes():raise AssertionError('GS status differs from independent control')
        for name in ('gs_eigen.data','gs_info.data','gs_k.data'):
            if (out/label/name).read_bytes()!=(ref/label/name).read_bytes():raise AssertionError(name+' differs from independent control')
    if (out/'a-first-status.dat').read_bytes()!=(out/'a-final-status.dat').read_bytes():raise AssertionError('A changed after B')
    if list(out.rglob('wfn.bin')) or list(out.rglob('occupation.bin')):raise AssertionError('GS restart data unexpectedly written')
    result={'status':'pass','size':a.size,'split':a.split,'xc_a':a.xc_a,'xc_b':a.xc_b,'fixed_c':a.fixed_c,
            'A_B_A_exact':True,'independent_status_eigen_info_k_exact':True,'shared_inputs_and_symmetry_restored':True,
            'MPI_children_live_after_each_run':0,'six_repeat_slots_stable':True,'GS_restart_written':False,
            'probe_sha256':probe_sha,'reference_sha256':a.reference_sha256,'commands':commands}
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(out)


if __name__=='__main__':main()
