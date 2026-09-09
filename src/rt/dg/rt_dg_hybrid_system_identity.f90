module rt_dg_hybrid_system_identity
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use structures,only:s_dft_system
  use dg_portable_sha256,only:s_dg_sha256_context,dg_sha256_schema,dg_sha256_init,&
    dg_sha256_update_int64,dg_sha256_update_real64,dg_sha256_update_logical,dg_sha256_final
  implicit none
  private
  public::fingerprint_rt_dg_hybrid_system
contains
  pure function fingerprint_rt_dg_hybrid_system(system,grid_num,periodic,projector_count,&
      pp_digest,xctype,spinorbit,plus_u,hse,fix_func,jm) result(digest)
    type(s_dft_system),intent(in)::system
    integer,intent(in)::grid_num(3),projector_count,xctype(:)
    logical,intent(in)::periodic,spinorbit,plus_u,hse,fix_func,jm
    integer(int64),intent(in)::pp_digest(4)
    integer(int64)::digest(4)
    type(s_dg_sha256_context)::hash
    integer::i,j,k,active_count
    real(real64),parameter::occupation_zero_tolerance=64d0*epsilon(1d0)
    real(real64)::electron_count
    digest=0_int64
    if(system%nion<1.or.system%nspin<1.or.system%no<1.or.system%nk<1.or.any(grid_num<1).or.&
      projector_count<0.or.all(pp_digest==0_int64))return
    if(.not.allocated(system%Rion).or..not.allocated(system%kion).or.&
      .not.allocated(system%vec_k).or..not.allocated(system%wtk).or..not.allocated(system%rocc))return
    if(size(system%Rion,1)/=3.or.size(system%Rion,2)/=system%nion.or.size(system%kion)/=system%nion)return
    if(size(system%vec_k,1)/=3.or.size(system%vec_k,2)/=system%nk.or.size(system%wtk)/=system%nk)return
    if(size(system%rocc,1)<system%no.or.size(system%rocc,2)/=system%nk.or.size(system%rocc,3)/=system%nspin)return
    if(.not.all(ieee_is_finite(system%Rion)).or..not.all(ieee_is_finite(system%vec_k)).or.&
      .not.all(ieee_is_finite(system%wtk)).or..not.all(ieee_is_finite(system%rocc)))return
    call dg_sha256_init(hash);call dg_sha256_update_int64(hash,int(z'53595354454D5635',int64))
    call dg_sha256_update_int64(hash,dg_sha256_schema)
    call dg_sha256_update_int64(hash,int(system%nion,int64));call dg_sha256_update_int64(hash,int(system%nspin,int64))
    call dg_sha256_update_int64(hash,int(system%nk,int64));call dg_sha256_update_int64(hash,int(system%ngrid,int64))
    call dg_sha256_update_int64(hash,int(projector_count,int64))
    do i=1,4;call dg_sha256_update_int64(hash,pp_digest(i));enddo
    call dg_sha256_update_logical(hash,periodic);call dg_sha256_update_logical(hash,spinorbit)
    call dg_sha256_update_logical(hash,plus_u);call dg_sha256_update_logical(hash,hse)
    call dg_sha256_update_logical(hash,fix_func);call dg_sha256_update_logical(hash,jm)
    do i=1,3
      call dg_sha256_update_int64(hash,int(grid_num(i),int64));call dg_sha256_update_real64(hash,system%hgs(i))
      do j=1,3;call dg_sha256_update_real64(hash,system%primitive_a(j,i));enddo
    enddo
    do i=1,system%nion
      call dg_sha256_update_int64(hash,int(system%kion(i),int64))
      do j=1,3;call dg_sha256_update_real64(hash,system%Rion(j,i));enddo
    enddo
    do i=1,system%nk
      do j=1,3;call dg_sha256_update_real64(hash,system%vec_k(j,i));enddo
      call dg_sha256_update_real64(hash,system%wtk(i))
    enddo
    ! Solver-only trailing zero states are excluded.  The exact bit patterns of
    ! every physically occupied value above the documented threshold are bound.
    electron_count=0d0
    do k=1,system%nspin;do j=1,system%nk
      active_count=count(abs(system%rocc(:,j,k))>occupation_zero_tolerance)
      call dg_sha256_update_int64(hash,int(active_count,int64))
      do i=1,size(system%rocc,1)
        if(abs(system%rocc(i,j,k))>occupation_zero_tolerance)then
          call dg_sha256_update_real64(hash,system%rocc(i,j,k))
          electron_count=electron_count+system%wtk(j)*system%rocc(i,j,k)
        endif
      enddo
    enddo;enddo
    ! Bind the physical charge explicitly as well as the canonical occupied list.
    call dg_sha256_update_real64(hash,electron_count)
    call dg_sha256_update_int64(hash,int(size(xctype),int64))
    do i=1,size(xctype);call dg_sha256_update_int64(hash,int(xctype(i),int64));enddo
    call dg_sha256_final(hash,digest)
  end function fingerprint_rt_dg_hybrid_system
end module rt_dg_hybrid_system_identity
