#!/usr/bin/env python3
"""Small sequential full-mesh RT integration; requires NumPy, MPI, HSE binary and H pseudo."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--pseudo',type=Path,required=True)
p.add_argument('--mpirun',default='mpirun')
p.add_argument('--work-root',type=Path,default=Path(tempfile.gettempdir()))
a=p.parse_args()
root=Path(tempfile.mkdtemp(prefix='hse-grid-native-',dir=a.work_root))
repo=Path(__file__).resolve().parents[2]
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')

def run(name,text,rt=False,settings=None,orbital=1,spatial=2,reject=None):
    d=root/name;d.mkdir();shutil.copy(a.pseudo,d/'H_rps.dat')
    if rt:(d/'data_dcdft').symlink_to(root/'gs/data_dcdft',target_is_directory=True)
    controls=dict(settings or {})
    if rt:controls.setdefault('yn_hse_realspace_rt','y')
    for key,value in controls.items():
        value="'"+value+"'" if key.startswith('yn_') else str(value)
        pattern=r'(\b'+re.escape(key)+r'\s*=\s*)(?:\x27[^\x27]*\x27|[^,\s/]+)'
        if re.search(pattern,text):text=re.sub(pattern,lambda m:m[1]+value,text)
        else:text=text.replace('&functional','&functional\n '+key+'='+value,1)
    if rt:text=text.replace('nproc_rgrid=2,1,1','nproc_rgrid='+str(spatial)+',1,1')
    if orbital>1:text=text.replace('nproc_ob=1','nproc_ob='+str(orbital))
    (d/'inputfile').write_text(text)
    with (d/'inputfile').open('rb') as inp,(d/'run.log').open('wb') as log:
        result=subprocess.run([a.mpirun,'-np',str(spatial*orbital),str(a.binary.resolve())],cwd=d,env=env,stdin=inp,stdout=log,stderr=subprocess.STDOUT,timeout=240)
    log=(d/'run.log').read_text()
    if reject:
        assert reject.lower() in log.lower(),(d,log[-2000:])
        assert 'Real-space HSE RT active' not in log,d
    else:
        assert result.returncode==0 and 'end SALMON' in log,(d,log[-2000:])
        if rt:
            if controls.get('yn_hse_realspace_rt')=='y':
                assert 'Real-space HSE RT active' in log and 'Grid HSE ACE build' in log,d
            else:
                assert 'Real-space HSE RT active' not in log,d
                assert 'HSE_WANNIER refresh/' in log,d
            assert 'Native LCFO RT active' not in log,d
    return d

def rows(d):
    z=np.loadtxt(d/'H_dc_hse_rt.data');assert z.shape[0]==4 and np.isfinite(z).all(),d
    return z

def compare(x,y,tol=1e-9):
    assert x.shape==y.shape
    delta=float(np.max(np.abs(x-y)));assert delta<tol,delta
    return delta

source=(repo/'testsuites/422_H_dcdft_hse/inputfile').read_text()
source=source.replace("xc='hse06'","xc='hse06'\n hse_mlwf_maxiter=200")
source=source.replace('nproc_rgrid_tot=4,1,1','nproc_rgrid_tot=2,1,1').replace('nproc_k=2','nproc_k=1').replace('num_kgrid=1,2,1','num_kgrid=1,1,1')
run('gs',source)
rt=source.replace("theory='dft'","theory='tddft_response'").replace("yn_dc='y'","yn_dc='n'\n yn_conventional_from_dcdft='y'")
rt=rt.replace('nproc_rgrid=1,1,1','nproc_rgrid=2,1,1').replace(' nstate=4',' nstate=2').replace(' temperature_k=300d0\n','')
rt=rt.replace('&control',"&control\n write_rt_wfn_k='y'\n method_wf_distributor='single'")
rt+='''
&tgrid
 nt=4
 dt=0.02d0
/
&emfield
 ae_shape1='impulse'
 e_impulse=0.0001d0
 epdir_re1=1d0,0d0,0d0
/
&analysis
 yn_out_dns_rt='y'
 out_dns_rt_step=4
 out_rt_energy_step=1
 nenergy=20
 de=0.01d0
/
'''
base=run('full',rt,rt=True,settings={'yn_hse_wannier':'n'})
wf=run('wf_full',rt,rt=True,settings={'yn_hse_wannier':'y','hse_rt_wf_radius':0})
parity=compare(rows(base),rows(wf))
reference=run('ordinary_native_reference',rt,rt=True,spatial=1,
              settings={'yn_hse_realspace_rt':'n','yn_hse_wannier':'y'})
reference_current=compare(rows(base)[:,13:16],rows(reference)[:,13:16])
def density(d):
    lines=(d/'H_dc_hse_dns_000004.cube').read_text().splitlines()
    n=abs(int(lines[2].split()[0]))
    z=np.array([float(x) for line in lines[6+n:] for x in line.split()])
    assert np.isfinite(z).all() and z.size==16*8*8
    return z
reference_density=compare(density(base),density(reference))
# Compare the complete energy table, including all recorded physical steps.
reference_energy=compare(np.loadtxt(base/'H_dc_hse_rt_energy.data'),
                         np.loadtxt(reference/'H_dc_hse_rt_energy.data'),tol=1e-8)
orb=run('wf_orbital2',rt,rt=True,orbital=2,settings={'yn_hse_wannier':'y','hse_rt_wf_radius':0})
mpi_parity=compare(rows(wf),rows(orb))
planned=run('wf_distributed_measured',rt,rt=True,settings={
    'yn_hse_wannier':'y','yn_hse_rt_seed_distributed':'y',
    'hse_rt_fft_batch':2,'yn_hse_rt_fft_measure':'y'})
compare(rows(wf),rows(planned))

local=run('wf_radius3',rt,rt=True,settings={'yn_hse_wannier':'y','hse_rt_wf_radius':3})
rows(local)
reuse=run('ace2_u2',rt,rt=True,settings={'yn_hse_wannier':'y','hse_rt_ace_interval':2,'hse_rt_u_interval':2})
rows(reuse)
assert (reuse/'run.log').read_text().count('Grid HSE ACE build')<(wf/'run.log').read_text().count('Grid HSE ACE build')
run('reject_fixed_lcfo',rt,rt=True,settings={'yn_hse_lcfo_rt':'y'},reject='fixed')

# Parse the complete saved fragment bases, not only occupied orbitals. Their cores
# partition the mesh; least-squares projection is independent of basis normalization.
def outside_initial_basis(folder):
    files=[x for x in folder.rglob('wfn.bin') if 'data_dcdft' not in x.parts]
    assert len(files)==1,(folder,files)
    psi=np.fromfile(files[0],dtype=np.complex128).reshape((16,8,8,2),order='F')
    assert np.isfinite(psi).all()
    # Detect wrong final buffer or checkpoint layout before interpreting projection.
    rho=2*np.sum(abs(psi)**2,axis=3)
    compare(rho.reshape(-1,order='C'),density(folder),tol=1e-10)
    compare(np.sum(abs(psi)**2,axis=(0,1,2)),np.ones(2),tol=1e-8)
    imported=np.zeros_like(psi)
    coverage=np.zeros(psi.shape[:3],dtype=int)
    initial_norm=0.
    residual=0.;norm=0.;initial_residual=0.
    for bfile in sorted((root/'gs/data_dcdft').rglob('basis_functions.bin')):
        with bfile.open('rb') as f:
            assert f.read(16).rstrip()==b'SLCFO_COMPLEX_V1'
            control=np.fromfile(f,np.int32,6);assert control[1]==0x01020304
            header=int(np.fromfile(f,np.int64,1)[0]);f.read(96)
            meta=np.fromfile(f,np.int32,20);assert meta[17]==1 and meta[16]==1
            f.seek(header);assert np.fromfile(f,np.int32,1)[0]==1
            np.fromfile(f,np.int64,1)
            spin,nb=np.fromfile(f,np.int32,2);assert spin==1
            shape=tuple(meta[6:9]);ng=int(np.prod(shape))
            basis=np.fromfile(f,np.complex128,ng*nb).reshape((ng,nb),order='F')
        with (bfile.parent/'rgrid_index.bin').open('rb') as f:
            dims=np.fromfile(f,np.int32,3);np.fromfile(f,np.int32,3)
            indices=[np.fromfile(f,np.int32,int(n))[:int(k)]-1 for n,k in zip(dims,shape)]
        with (bfile.parent/'wavefunctions.bin').open('rb') as f:
            f.seek(header+12)
            spin_c,nb_c,nmat=np.fromfile(f,np.int32,3)
            counts=np.fromfile(f,np.int32,int(np.prod(meta[9:12])))
            row_ids=np.fromfile(f,np.int32,int(nb_c))
            assert spin_c==spin and nb_c==nb and sum(counts)==nmat
            coeff=np.fromfile(f,np.complex128,int(nb_c*meta[19])).reshape((nb_c,meta[19]),order='F')[:,:2]
        initial_values=basis@coeff
        imported[np.ix_(*indices,np.arange(2))]=initial_values.reshape(shape+(2,),order='F')
        coverage[np.ix_(*indices)]+=1
        values=psi[np.ix_(*indices,np.arange(2))].reshape((ng,2),order='F')
        q=np.linalg.qr(basis,mode='reduced')[0]
        residual+=float(np.linalg.norm(values-q@(q.conj().T@values))**2)
        norm+=float(np.linalg.norm(values)**2)
        initial_residual+=float(np.linalg.norm(initial_values-q@(q.conj().T@initial_values))**2)
        initial_norm+=float(np.linalg.norm(initial_values)**2)
    assert norm>0 and (initial_residual/initial_norm)**.5<1e-12
    assert np.all(coverage==1)
    compare(np.sum(abs(imported)**2,axis=(0,1,2)),np.ones(2),tol=1e-10)
    fraction=(residual/norm)**.5
    assert fraction>1e-8,('Propagation did not escape the full initial LCFO span',fraction)
    return fraction
escape=outside_initial_basis(base)
print(json.dumps({'work':str(root),'full_wf_parity':parity,'native_reference_current':reference_current,'native_reference_density':reference_density,'native_reference_energy':reference_energy,'orbital_mpi_parity':mpi_parity,'outside_full_initial_lcfo_span':escape},indent=2))
