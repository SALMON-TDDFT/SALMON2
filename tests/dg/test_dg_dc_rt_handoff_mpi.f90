program test_dg_dc_rt_handoff_mpi
  use mpi
  use dg_wpw_checkpoint,only:s_dg_wpw_checkpoint_state,write_dg_wpw_checkpoint,&
    write_dg_wpw_checkpoint_manifest,dg_wpw_checkpoint_checksum
  use rt_dg_wpw_checkpoint_handoff,only:s_rt_dg_wpw_checkpoint_handoff,load_rt_dg_wpw_checkpoint_handoff
  implicit none
  type(s_dg_wpw_checkpoint_state)::state
  type(s_rt_dg_wpw_checkpoint_handoff)::handoff
  integer::rank,ierr,info
  integer(8)::checksum,checksums(2)
  character(128)::path
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call make_state(state,rank)
  write(path,'(a,i6.6,a)')'tiny.dg_wpw.rank_',rank,'.chk'
  call write_dg_wpw_checkpoint(path,state,info);if(info/=0)error stop 1
  checksum=dg_wpw_checkpoint_checksum(state)
  call MPI_Gather(checksum,1,MPI_INTEGER8,checksums,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
  if(rank==0)call write_dg_wpw_checkpoint_manifest('tiny.dg_wpw.manifest',checksums,info)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call load_rt_dg_wpw_checkpoint_handoff('tiny.dg_wpw.manifest','tiny.dg_wpw.rank_',&
    MPI_COMM_WORLD,1d-12,handoff,info)
  if(info/=0.or..not.handoff%valid)error stop 2
  if(handoff%state%layout_fingerprint/=state%layout_fingerprint)error stop 3
  if(handoff%occupied_residual>1d-12.or.handoff%metric_orthonormality>1d-12)error stop 4
  if(rank==0)then
    open(unit=91,file='tiny.dg_wpw.manifest',access='stream',form='unformatted',status='old',position='append')
    write(91)1
    close(91)
  endif
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call load_rt_dg_wpw_checkpoint_handoff('tiny.dg_wpw.manifest','tiny.dg_wpw.rank_',&
    MPI_COMM_WORLD,1d-12,handoff,info)
  if(info==0)error stop 5
  if(rank==0)print '(a)','PASS checkpoint-backed DG-DC to RT H/S identity handoff'
  call MPI_Finalize(ierr)
contains
  subroutine make_state(s,r)
    type(s_dg_wpw_checkpoint_state),intent(out)::s
    integer,intent(in)::r
    integer::wid,pid
    wid=r+1;pid=3+r
    s%schema_version=2;s%operator_epoch=4;s%layout_fingerprint=111_8+r
    s%ownership_kind='fragment_root_v1';s%metric_convention='orthonormal_ww'
    s%operator_convention='windowed_kg_sipg_v1';s%n_occ=1
    s%peer_ranks=[integer::]
    s%owned_w_ids=[wid];s%owned_p_ids=[pid];s%support_w_ids=[wid];s%support_p_ids=[pid]
    s%eigenvalues=[1d0,3d0];s%occupations=[1d0,0d0]
    allocate(s%coeff_w(1,2),s%coeff_p(1,2));s%coeff_w=0;s%coeff_p=0
    if(r==0)s%coeff_w(1,1)=1d0
    s%potential=[0d0]
    s%ww_r=[wid];s%ww_c=[wid];s%wp_w=[integer::];s%wp_p=[integer::]
    s%pp_r=[pid];s%pp_c=[pid]
    s%ww_h=[cmplx(dble(wid),0d0,8)];s%ww_s=[(1d0,0d0)]
    s%wp_h=[complex(8)::];s%wp_s=[complex(8)::]
    s%pp_h=[cmplx(dble(pid),0d0,8)];s%pp_s=[(1d0,0d0)]
    s%ww_z=[cmplx(0.1d0*dble(wid),0d0,8)]
    s%wp_z=[complex(8)::];s%pp_z=[cmplx(0.1d0*dble(pid),0d0,8)]
  end subroutine
end program
