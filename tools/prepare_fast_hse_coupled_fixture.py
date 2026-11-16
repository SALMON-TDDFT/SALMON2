"""Explicit zero-current synthetic three-k-point fixture for fast-HSE coupling regressions.
Not a converged GS or production material recipe. No engine run/download/build.
Usage: python -I tools/prepare_fast_hse_coupled_fixture.py PSEUDOPOTENTIAL NEW_DIR
"""
from pathlib import Path
import hashlib,json,struct,shutil,sys
import numpy as np

pseudo=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=False)
r=out/'data_for_restart';r.mkdir()
def record(stream,payload):
    stream.write(struct.pack('=i',len(payload)));stream.write(payload);stream.write(struct.pack('=i',len(payload)))
with (r/'info.bin').open('wb') as stream:
    for n in (3,2,0,1,0):record(stream,struct.pack('=i',n))
with (r/'occupation.bin').open('wb') as stream:record(stream,np.full((2,3),2.,dtype='=f8').tobytes(order='F'))
x=np.arange(8)[:,None,None]
psi=np.empty((8,8,8,2),dtype='=c16',order='F')
psi[:,:,:,0]=1/np.sqrt(512.)
psi[:,:,:,1]=np.sqrt(2/512.)*np.cos(2*np.pi*x/8)
np.testing.assert_allclose(psi.reshape(512,2,order='F').conj().T@psi.reshape(512,2,order='F'),np.eye(2),atol=1e-14)
(r/'wfn.bin').write_bytes(np.repeat(psi[...,None],3,axis=4).tobytes(order='F'));shutil.copy2(pseudo,out/'pseudo1_rps.dat')
for prop,factor,tag in [('hse_taylor4','learned','taylor'),('hse_ptcn','learned','ptcn')]:
    text=f"""&calculation
 theory='tddft_response'
/
&control
 sysname='rt',directory_read_data='./data_for_restart/',yn_restart='n'
/
&units
 unit_system='au'
/
&system
 yn_periodic='y',al=8,8,8,nelem=1,natom=1,nelec=4,nstate=2
/
&pseudo
 file_pseudo(1)='pseudo1_rps.dat',izatom(1)=14,lloc_ps(1)=2
/
&functional
 xc='hse06',exx_factor_mode='{factor}',exx_factor_exact_interval={1 if prop=='hse_ptcn' else 2},exx_factor_n_teach=1
/
&rgrid
 num_rgrid=8,8,8
/
&kgrid
 num_kgrid=3,1,1
/
&tgrid
 dt=.001,nt=256
/
&propagation
 propagator='{prop}',yn_predictor_corrector='{'n' if prop=='hse_ptcn' else 'y'}',n_hamil=4
/
&emfield
 ae_shape1='impulse',e_impulse=0,trans_longi='tr'
/
&analysis
 out_rt_energy_step=1
/
&parallel
 nproc_k=1,nproc_ob=1,nproc_rgrid=1,1,1
/
&atomic_red_coor
 'Si' .3 .4 .45 1
/
"""
    (out/f'{tag}.inp').write_text(text)
shutil.copy2(out/'taylor.inp',out/'rt.inp')
(out/'fixture.json').write_text(json.dumps({'scope':'zero-current real-orbital time-reversal synthetic continuation regression; not converged GS',
    'files':{str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in out.rglob('*') if p.is_file()}},indent=2)+'\n')
