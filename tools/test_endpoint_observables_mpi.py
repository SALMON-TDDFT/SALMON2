"""Manual frozen-wavefunction endpoint-A audit (not required physics CI).

--instrument PREPARED_SOURCE NEW_SOURCE copies a prepared source and adds the
Fortran probe below to that NEW diagnostic tree only. Build it with the usual
source/object audits. No probe is added to the distributed solver/API.

mpiexec -n N python -m tools.test_endpoint_observables_mpi BUILD FIXTURE NEW_OUT
Use --resume CHECKPOINT_DIRECTORY for a fresh-process direct restart at step4.
An existing periodic PZ Si rt.inp/pseudo1_rps.dat/data_for_restart is required.
GS is referenced in place, never calculated or copied. Only electronic RT
checkpoints, observable tables and wavefunction SHA256 values are written.
"""
import argparse
import ctypes as c
import hashlib
import json
import shutil
from pathlib import Path

PROBE = r'''
! Test-only frozen-state observation probe: never included in release sources.
integer(c_int) function sp_test_endpoint_probe(wave,capacity,values) bind(C) result(code)
  use total_energy,only:calc_eigen_energy
  use communication,only:comm_resume_children,comm_end_children
  complex(c_double_complex),intent(out)::wave(*)
  integer(c_int),value::capacity
  real(c_double),intent(out)::values(18)
  type(s_orbital)::psi,ht,tt
  type(s_pp_grid)::saved_pp
  type(s_dft_energy)::energy
  real(8)::saved_A(3),endpoint(3),curr(3,state%system%nspin,1),ep,em
  integer::n,i,ix,iy,iz,ispin,io,ik,im,axis,status
  code=-1;values=0d0
  if(.not.associated(state))return
  if(mod(state%itt,2)==1)then
    psi=state%spsi_out
  else
    psi=state%spsi_in
  endif
  n=product(state%mg%ie-state%mg%is+1)*state%system%nspin &
    *(state%info%io_e-state%info%io_s+1)*(state%info%ik_e-state%info%ik_s+1)
  if(capacity==0)then
    code=n;return
  endif
  if(capacity/=n)return
  call comm_resume_children(state%mpi_child_owner,status)
  if(status/=0)return
  call legacy_timer_enter()
  i=0
  do im=1,1
  do ik=state%info%ik_s,state%info%ik_e
  do io=state%info%io_s,state%info%io_e
  do ispin=1,state%system%nspin
  do iz=state%mg%is(3),state%mg%ie(3)
  do iy=state%mg%is(2),state%mg%ie(2)
  do ix=state%mg%is(1),state%mg%ie(1)
    i=i+1;wave(i)=psi%zwf(ix,iy,iz,ispin,io,ik,im)
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  saved_pp=state%ppg;saved_A=state%system%vec_Ac
  endpoint=state%rt%Ac_tot(:,state%itt)
  ht=state%tpsi;tt=state%tpsi
  values(1:3)=state%rt%curr(:,state%itt)
  call calc_current(state%system,state%mg,state%stencil,state%info,state%srg,psi,state%ppg,curr)
  values(4:6)=sum(curr(:,:,1),dim=2)
  values(10:11)=[state%energy%E_kin,state%energy%E_ion_nloc]
  values(14)=state%energy%E_tot
  state%system%vec_Ac=endpoint
  call update_kvector_nonlocalpt(state%info%ik_s,state%info%ik_e,state%system,state%ppg)
  call calc_current(state%system,state%mg,state%stencil,state%info,state%srg,psi,state%ppg,curr)
  values(7:9)=sum(curr(:,:,1),dim=2)
  call evaluate_energy()
  values(12:13)=[energy%E_kin,energy%E_ion_nloc]
  values(15)=values(14)+sum(values(12:13)-values(10:11))
  do axis=1,3
    state%system%vec_Ac=endpoint;state%system%vec_Ac(axis)=endpoint(axis)+1d-5
    call update_kvector_nonlocalpt(state%info%ik_s,state%info%ik_e,state%system,state%ppg)
    call evaluate_energy();ep=energy%E_kin+energy%E_ion_nloc
    state%system%vec_Ac(axis)=endpoint(axis)-1d-5
    call update_kvector_nonlocalpt(state%info%ik_s,state%info%ik_e,state%system,state%ppg)
    call evaluate_energy();em=energy%E_kin+energy%E_ion_nloc
    values(15+axis)=(ep-em)/(2d-5*state%system%det_a)
  enddo
  state%system%vec_Ac=saved_A;state%ppg=saved_pp
  call legacy_timer_leave()
  call comm_end_children(state%mpi_child_owner)
  code=0
contains
  subroutine evaluate_energy()
    energy=state%energy
    call calc_eigen_energy(energy,psi,ht,tt,state%system,state%info,state%mg,state%V_local, &
      state%stencil,state%srg,state%ppg)
  end subroutine
end function
'''


def instrument(source, output):
    source=source.resolve(strict=True);output=output.resolve()
    if output.exists():raise FileExistsError(output)
    shutil.copytree(source,output)
    path=output/'src/rt/embedded_tddft.f90'
    text=path.read_text();assert text.rstrip().endswith('end module')
    path.write_text(text[:text.rfind('end module')]+PROBE+'\nend module\n')
    (output/'endpoint-audit-instrumentation.json').write_text(json.dumps({
        'base_source':str(source),'diagnostic_only':True,
        'modified_file':'src/rt/embedded_tddft.f90',
        'sha256':hashlib.sha256(path.read_bytes()).hexdigest()},indent=2)+'\n')



def link_probe(base_record, source, output):
    """Compile only diagnostic module; audit every borrowed application object."""
    import shlex
    import subprocess
    import sys
    base_record=base_record.resolve();source=source.resolve();output=output.resolve()
    rec=json.loads(base_record.read_text());base=base_record.parent/'build'
    assert hashlib.sha256(Path(rec['library']).read_bytes()).hexdigest()==rec['sha256']
    assert (source/'endpoint-audit-instrumentation.json').is_file()
    changed={'src/rt/time_evolution_step.f90','src/rt/embedded_tddft.f90'}
    for name,digest in rec['source_hashes'].items():
        if name not in changed:
            assert hashlib.sha256((source/name).read_bytes()).hexdigest()==digest,name
    output.mkdir(exist_ok=False)
    checker=Path(__file__).with_name('check_forbidden_eigensolvers.py')
    def audit(*args):subprocess.run([sys.executable,str(checker),*map(str,args)],check=True)
    audit('--source-root',source/'src')
    flags=(base/'src/CMakeFiles/salmon.dir/flags.make').read_text()
    import re
    inc=shlex.split(re.search(r'^Fortran_INCLUDES = (.*)$',flags,re.M)[1])
    opts=shlex.split(re.search(r'^Fortran_FLAGS = (.*)$',flags,re.M)[1])
    opts=[x for x in opts if not x.startswith('-J')]
    replacements={}
    for name in ('src/rt/time_evolution_step.f90','src/rt/embedded_tddft.f90'):
        obj=output/(Path(name).name+'.o')
        subprocess.run(['mpifort',*inc,*opts,'-I'+str(base),'-J'+str(output),'-c',
                        str(source/name),'-o',str(obj)],check=True)
        replacements['CMakeFiles/salmon.dir/'+name.removeprefix('src/')+'.o']=obj
    command=list(rec['link_command']);idx=command.index('-o')
    library=output/Path(rec['library']).name;command[idx+1]=str(library)
    objects=[]
    for i,value in enumerate(command):
        if value.endswith('.o'):
            path=Path(value)
            if not path.is_absolute():path=base/'src'/path
            if value in replacements:path=replacements[value]
            command[i]=str(path);objects.append(str(path))
    audit('--object-tree',base)
    audit('--object-tree',output)
    subprocess.run(command,cwd=base/'src',check=True)
    rec.update(library=str(library),sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
               diagnostic_only=True,diagnostic_source=str(source),
               diagnostic_source_sha256={name:hashlib.sha256((source/name).read_bytes()).hexdigest()
                   for name in ('src/rt/embedded_tddft.f90','src/rt/time_evolution_step.f90')},
               link_command=command,audited_application_objects=objects)
    for name in changed:rec['source_hashes'][name]=hashlib.sha256((source/name).read_bytes()).hexdigest()
    (output/'build.json').write_text(json.dumps(rec,indent=2)+'\n')


def run(a):
    import numpy as np
    from mpi4py import MPI
    from types import SimpleNamespace
    from salmonpy.distributed import NativeParallelContext
    from salmonpy._native_group import _directory_call
    comm=MPI.COMM_WORLD
    try:
        rec=json.loads(a.build.read_text());library=Path(rec['library'])
        assert hashlib.sha256(library.read_bytes()).hexdigest()==rec['sha256']
        lib=c.CDLL(str(library));out=a.output.resolve();fixture=a.fixture.resolve()
        if comm.rank==0:
            out.mkdir(parents=True,exist_ok=False)
            text=(fixture/'rt.inp').read_text()
            import re
            text=re.sub(r"xc\s*=\s*'[^']+'", "xc='"+a.xc+"'", text)
            (out/'rt.inp').write_text(text)
            for name in ('data_for_restart','pseudo1_rps.dat'):
                (out/name).symlink_to(fixture/name,target_is_directory=name=='data_for_restart')
        comm.Barrier()
        owner=SimpleNamespace(solver_comm=comm,closed=False,_dependents=0)
        native=NativeParallelContext(owner,rec['library'],rec['sha256'])
        def call(fn,*args):return _directory_call(lib,comm,out,fn,*args)
        vec=np.ctypeslib.ndpointer(dtype=np.float64,shape=(3,),flags='C_CONTIGUOUS')
        dec_t=np.ctypeslib.ndpointer(dtype=np.int32,shape=(5,),flags='C_CONTIGUOUS')
        def setup(name,args):
            f=getattr(lib,'sp_tddft_'+name);f.argtypes=args;f.restype=c.c_int;return f
        init=setup('init',[c.c_char_p,c.c_double,c.c_int,vec,dec_t,vec])
        step=setup('step_performance' if a.performance else 'step',[vec,vec,vec]);close=setup('close',[])
        checkpoint=setup('checkpoint',[c.c_char_p,c.c_int])
        if a.handles:
            if not a.no_probe:raise ValueError('handle route uses the production no-probe check')
            token=c.c_int64()
            hi=setup('handle_init',[c.c_char_p,c.c_double,c.c_int,vec,dec_t,c.POINTER(c.c_int64),vec])
            hs=setup('handle_step_performance' if a.performance else 'handle_step',[c.c_int64,vec,vec,vec])
            hc=setup('handle_close',[c.c_int64]);hp=setup('handle_checkpoint',[c.c_int64,c.c_char_p,c.c_int])
            init=lambda path,dt,cap,A,dec,j:hi(path,dt,cap,A,dec,c.byref(token),j)
            step=lambda old,new,j:hs(token,old,new,j)
            checkpoint=lambda path,mode:hp(token,path,mode)
            close=lambda:hc(token)

        probe=None if a.no_probe else lib.sp_test_endpoint_probe
        if probe:probe.argtypes=[c.c_void_p,c.c_int,np.ctypeslib.ndpointer(dtype=np.float64,shape=(18,),flags='C_CONTIGUOUS')]
        if probe:probe.restype=c.c_int
        dec=np.array([comm.size,1,1,1,1] if a.split=='k' else
                     [1,comm.size,1,1,1] if a.split=='orbital' else [1,1,comm.size,1,1],dtype=np.int32)
        zero=np.zeros(3);j=zero.copy();old=zero.copy();rows=[];hashes=[]
        assert call(init,b'rt.inp',.02,8,zero,dec,j)==0
        start=1
        def drive(i):return .01*np.array([np.sin(.31*i),np.cos(.23*i)-1,np.sin(.17*i)])
        if a.resume:
            assert call(checkpoint,str(a.resume.resolve()/f'electronic-{comm.rank}.bin').encode(),1)==0
            old=drive(4);start=5
        try:
            for i in range(start,9):
                new=drive(i);assert call(step,old,new,j)==0
                v=np.zeros(18);v[:3]=j
                if probe:
                    n=call(probe,None,0,v);assert n>0
                    wave=np.empty(n,dtype=np.complex128)
                    assert call(probe,wave.ctypes.data,n,v)==0
                    # Legacy ABI returns electron current; handle API returns charge current.
                    np.testing.assert_array_equal(j,v[:3])
                    hashes.append(hashlib.sha256(wave.tobytes()).hexdigest())
                rows.append(v)
                if i==4:
                    assert call(checkpoint,str(out/f'electronic-{comm.rank}.bin').encode(),0)==0
                old=new
        finally:
            assert call(close)==0;native.close()
        assert not MPI.Is_finalized()
        rows=np.asarray(rows)
        np.savez(out/f'observations-{comm.rank}.npz',values=rows)
        (out/f'wave-hashes-{comm.rank}.json').write_text(json.dumps(hashes)+'\n')
        if not probe:
            if comm.rank==0:(out/'result.json').write_text(json.dumps({'production_library':True,
                'MPI_alive':True,'ranks':comm.size,'split':a.split,'xc':a.xc,
                'handle_api':a.handles,'current_sign':'charge' if a.handles else 'electron',
                'library_sha256':rec['sha256'],'restored_step':4 if a.resume else 0})+'\n')
            return
        errors={'current':float(np.max(np.abs(rows[:,:3]-rows[:,6:9]))),
                'energy':float(np.max(np.abs(rows[:,9:11]-rows[:,11:13]))),
                'dE_dA_over_volume':float(np.max(np.abs(rows[:,6:9]-rows[:,15:18])))}
        if comm.rank==0:
            (out/'result.json').write_text(json.dumps({'errors':errors,'ranks':comm.size,'split':a.split,'xc':a.xc,
                'library_sha256':rec['sha256'],'resumed_without_steps':4 if a.resume else 0,
                'GS_calculated':False,'GS_copied':False,'MPI_alive':True},indent=2)+'\n')
            print(json.dumps(errors),flush=True)
        np.testing.assert_allclose(rows[:,:3],rows[:,6:9],atol=2e-12,rtol=2e-10)
        np.testing.assert_allclose(rows[:,9:11],rows[:,11:13],atol=2e-11,rtol=2e-10)
        np.testing.assert_allclose(rows[:,6:9],rows[:,15:18],atol=2e-10,rtol=2e-6)
    except BaseException:
        import traceback
        traceback.print_exc();comm.Abort(1)


def main():
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='--instrument':
        instrument(Path(sys.argv[2]),Path(sys.argv[3]));return
    if len(sys.argv)>1 and sys.argv[1]=='--link-probe':
        link_probe(*map(Path,sys.argv[2:5]));return
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('build','fixture','output'):p.add_argument(key,type=Path)
    p.add_argument('--resume',type=Path)
    p.add_argument('--no-probe',action='store_true')
    p.add_argument('--handles',action='store_true')
    p.add_argument('--performance',action='store_true')
    p.add_argument('--xc',choices=('pz','libxc_pz','libxc_pbe','tbmbj'),default='pz')
    p.add_argument('--split',choices=('k','orbital','space'),default='k')
    run(p.parse_args())


if __name__=='__main__':main()
