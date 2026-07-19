module rt_dg_wpw_checkpoint_handoff
  use mpi,only:MPI_Comm_rank,MPI_Comm_size,MPI_Allreduce,MPI_INTEGER,MPI_MAX,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_SUCCESS
  use dg_wpw_checkpoint,only:s_dg_wpw_checkpoint_state,read_dg_wpw_checkpoint,&
    read_dg_wpw_checkpoint_manifest,dg_wpw_checkpoint_checksum
  use dg_wpw_bounded_operator,only:s_dg_wpw_bounded_operator,initialize_dg_wpw_bounded_operator,&
    apply_h_dg_wpw_bounded,apply_s_dg_wpw_bounded,global_gram_dg_wpw_bounded
  implicit none
  private
  type,public::s_rt_dg_wpw_checkpoint_handoff
    logical::valid=.false.
    real(8)::occupied_residual=huge(1d0),metric_orthonormality=huge(1d0)
    real(8)::position_hermiticity=huge(1d0),commutator_norm=huge(1d0)
    complex(8),allocatable::position_reduced(:,:)
    type(s_dg_wpw_checkpoint_state)::state
    type(s_dg_wpw_bounded_operator)::operator
    type(s_dg_wpw_bounded_operator)::position_operator
  end type
  public::load_rt_dg_wpw_checkpoint_handoff
contains
  subroutine load_rt_dg_wpw_checkpoint_handoff(manifest_path,rank_prefix,comm,tolerance,handoff,info)
    character(*),intent(in)::manifest_path,rank_prefix
    integer,intent(in)::comm
    real(8),intent(in)::tolerance
    type(s_rt_dg_wpw_checkpoint_handoff),intent(out)::handoff
    integer,intent(out)::info
    integer::rank,nrank,ierr,local_bad,global_bad,load_info,nocc,nstate,i
    integer(8),allocatable::checksums(:)
    character(1024)::rank_path
    complex(8),allocatable::hw(:,:),hp(:,:),sw(:,:),sp(:,:),resw(:,:),resp(:,:),gram(:,:)
    complex(8),allocatable::zw(:,:),zp(:,:),hzw(:,:),hzp(:,:),zhw(:,:),zhp(:,:),zgram(:,:),allhw(:,:),allhp(:,:)
    real(8)::local_residual,global_residual,local_commutator,global_commutator
    info=1;local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)return
    call read_dg_wpw_checkpoint_manifest(manifest_path,checksums,load_info)
    if(load_info/=0.or.size(checksums)/=nrank)local_bad=1
    if(local_bad==0)then
      write(rank_path,'(a,i6.6,a)')trim(rank_prefix),rank,'.chk'
      call read_dg_wpw_checkpoint(trim(rank_path),handoff%state,load_info)
      if(load_info/=0)then
        local_bad=1
      else if(dg_wpw_checkpoint_checksum(handoff%state)/=checksums(rank+1))then
        local_bad=1
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call initialize_dg_wpw_bounded_operator(handoff%operator,comm,handoff%state%operator_epoch,&
      handoff%state%layout_fingerprint,handoff%state%ownership_kind,handoff%state%metric_convention,&
      handoff%state%operator_convention,handoff%state%peer_ranks,handoff%state%owned_w_ids,&
      handoff%state%owned_p_ids,handoff%state%support_w_ids,handoff%state%support_p_ids,&
      handoff%state%ww_r,handoff%state%ww_c,handoff%state%ww_h,handoff%state%ww_s,&
      handoff%state%wp_w,handoff%state%wp_p,handoff%state%wp_h,handoff%state%wp_s,&
      handoff%state%pp_r,handoff%state%pp_c,handoff%state%pp_h,handoff%state%pp_s,load_info)
    if(load_info/=0)return
    call initialize_dg_wpw_bounded_operator(handoff%position_operator,comm,handoff%state%operator_epoch,&
      handoff%state%layout_fingerprint,handoff%state%ownership_kind,handoff%state%metric_convention,&
      'length_gauge_volume_v1',handoff%state%peer_ranks,handoff%state%owned_w_ids,&
      handoff%state%owned_p_ids,handoff%state%support_w_ids,handoff%state%support_p_ids,&
      handoff%state%ww_r,handoff%state%ww_c,handoff%state%ww_z,handoff%state%ww_s,&
      handoff%state%wp_w,handoff%state%wp_p,handoff%state%wp_z,handoff%state%wp_s,&
      handoff%state%pp_r,handoff%state%pp_c,handoff%state%pp_z,handoff%state%pp_s,load_info)
    if(load_info/=0)return
    nocc=handoff%state%n_occ
    allocate(hw(size(handoff%state%coeff_w,1),nocc),hp(size(handoff%state%coeff_p,1),nocc),&
      sw(size(handoff%state%coeff_w,1),nocc),sp(size(handoff%state%coeff_p,1),nocc))
    allocate(resw(size(hw,1),nocc),resp(size(hp,1),nocc),gram(nocc,nocc))
    call apply_h_dg_wpw_bounded(handoff%operator,handoff%state%operator_epoch,&
      handoff%state%layout_fingerprint,handoff%state%coeff_w(:,1:nocc),handoff%state%coeff_p(:,1:nocc),hw,hp,load_info)
    if(load_info/=0)return
    call apply_s_dg_wpw_bounded(handoff%operator,handoff%state%operator_epoch,&
      handoff%state%layout_fingerprint,handoff%state%coeff_w(:,1:nocc),handoff%state%coeff_p(:,1:nocc),sw,sp,load_info)
    if(load_info/=0)return
    resw=hw;resp=hp
    do i=1,nocc
      resw(:,i)=resw(:,i)-handoff%state%eigenvalues(i)*sw(:,i)
      resp(:,i)=resp(:,i)-handoff%state%eigenvalues(i)*sp(:,i)
    enddo
    local_residual=sum(abs(resw)**2)+sum(abs(resp)**2)
    call MPI_Allreduce(local_residual,global_residual,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    handoff%occupied_residual=sqrt(max(0d0,global_residual))
    call global_gram_dg_wpw_bounded(handoff%operator,handoff%state%coeff_w(:,1:nocc),&
      handoff%state%coeff_p(:,1:nocc),sw,sp,gram,load_info)
    if(load_info/=0)return
    do i=1,nocc;gram(i,i)=gram(i,i)-1d0;enddo
    handoff%metric_orthonormality=maxval(abs(gram))
    nstate=size(handoff%state%eigenvalues)
    allocate(zw(size(hw,1),nstate),zp(size(hp,1),nstate),hzw(size(hw,1),nstate),hzp(size(hp,1),nstate),&
      zhw(size(hw,1),nstate),zhp(size(hp,1),nstate),zgram(nstate,nstate),&
      allhw(size(hw,1),nstate),allhp(size(hp,1),nstate))
    call apply_h_dg_wpw_bounded(handoff%position_operator,handoff%state%operator_epoch,&
      handoff%state%layout_fingerprint,handoff%state%coeff_w,handoff%state%coeff_p,zw,zp,load_info)
    if(load_info/=0)return
    call global_gram_dg_wpw_bounded(handoff%operator,handoff%state%coeff_w,&
      handoff%state%coeff_p,zw,zp,zgram,load_info)
    if(load_info/=0)return
    handoff%position_hermiticity=maxval(abs(zgram-conjg(transpose(zgram))))/&
      max(1d0,maxval(abs(zgram)))
    handoff%position_reduced=0.5d0*(zgram+conjg(transpose(zgram)))
    call apply_h_dg_wpw_bounded(handoff%operator,handoff%state%operator_epoch,&
      handoff%state%layout_fingerprint,zw,zp,hzw,hzp,load_info)
    if(load_info/=0)return
    call apply_h_dg_wpw_bounded(handoff%operator,handoff%state%operator_epoch,&
      handoff%state%layout_fingerprint,handoff%state%coeff_w,handoff%state%coeff_p,allhw,allhp,load_info)
    if(load_info/=0)return
    call apply_h_dg_wpw_bounded(handoff%position_operator,handoff%state%operator_epoch,&
      handoff%state%layout_fingerprint,allhw,allhp,zhw,zhp,load_info)
    if(load_info/=0)return
    local_commutator=sum(abs(hzw-zhw)**2)+sum(abs(hzp-zhp)**2)
    call MPI_Allreduce(local_commutator,global_commutator,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    handoff%commutator_norm=sqrt(max(0d0,global_commutator))
    local_bad=merge(0,1,handoff%occupied_residual<=tolerance.and.&
      handoff%metric_orthonormality<=tolerance.and.handoff%position_hermiticity<=tolerance.and.&
      handoff%commutator_norm<huge(1d0))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    handoff%valid=.true.;info=0
  end subroutine
end module
