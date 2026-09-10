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
      pp_digest,xctype,spinorbit,plus_u,hse,fix_func,jm,electron_count,spin_electron_count) result(digest)
    type(s_dft_system),intent(in)::system
    integer,intent(in)::grid_num(3),projector_count,xctype(:)
    logical,intent(in)::periodic,spinorbit,plus_u,hse,fix_func,jm
    integer,intent(in)::electron_count,spin_electron_count(2)
    integer(int64),intent(in)::pp_digest(4)
    integer(int64)::digest(4)
    type(s_dg_sha256_context)::hash
    integer::i,j
    digest=0_int64
    if(system%nion<1.or.system%nspin<1.or.system%no<1.or.system%nk<1.or.any(grid_num<1).or.&
      projector_count<0.or.all(pp_digest==0_int64))return
    if(electron_count<1.or.any(spin_electron_count<0).or.&
      (sum(spin_electron_count)>0.and.(system%nspin/=2.or.sum(spin_electron_count)/=electron_count)))return
    if(.not.allocated(system%Rion).or..not.allocated(system%kion).or.&
      .not.allocated(system%vec_k).or..not.allocated(system%wtk))return
    if(size(system%Rion,1)/=3.or.size(system%Rion,2)/=system%nion.or.size(system%kion)/=system%nion)return
    if(size(system%vec_k,1)/=3.or.size(system%vec_k,2)/=system%nk.or.size(system%wtk)/=system%nk)return
    if(.not.all(ieee_is_finite(system%Rion)).or..not.all(ieee_is_finite(system%vec_k)).or.&
      .not.all(ieee_is_finite(system%wtk)))return
    call dg_sha256_init(hash);call dg_sha256_update_int64(hash,int(z'53595354454D5635',int64))
    call dg_sha256_update_int64(hash,dg_sha256_schema)
    call dg_sha256_update_int64(hash,int(system%nion,int64));call dg_sha256_update_int64(hash,int(system%nspin,int64))
    call dg_sha256_update_int64(hash,int(system%nk,int64));call dg_sha256_update_int64(hash,int(system%ngrid,int64))
    call dg_sha256_update_int64(hash,int(electron_count,int64))
    call dg_sha256_update_logical(hash,sum(spin_electron_count)>0)
    do i=1,2;call dg_sha256_update_int64(hash,int(spin_electron_count(i),int64));enddo
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
    ! Thermal and solver-specific occupations are intentionally excluded here.
    ! The LCFO occupations are authenticated separately in the checkpoint payload.
    call dg_sha256_update_int64(hash,int(size(xctype),int64))
    do i=1,size(xctype);call dg_sha256_update_int64(hash,int(xctype(i),int64));enddo
    call dg_sha256_final(hash,digest)
  end function fingerprint_rt_dg_hybrid_system
end module rt_dg_hybrid_system_identity
