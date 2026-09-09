module rt_dg_hybrid_system_identity
  use,intrinsic::iso_fortran_env,only:int64,real64
  use structures,only:s_dft_system
  use dg_portable_sha256,only:dg_sha256_schema,dg_sha256_mix_int64,dg_sha256_mix_real64,dg_sha256_mix_logical
  implicit none
  private
  public::fingerprint_rt_dg_hybrid_system
contains
  pure integer(int64) function fingerprint_rt_dg_hybrid_system(system,grid_num,periodic,&
      projector_count,pp_fingerprint,xctype,spinorbit,plus_u,hse,fix_func,jm) result(hash)
    type(s_dft_system),intent(in)::system
    integer,intent(in)::grid_num(3),projector_count,xctype(:)
    logical,intent(in)::periodic,spinorbit,plus_u,hse,fix_func,jm
    integer(int64),intent(in)::pp_fingerprint
    integer::i,j
    hash=int(z'243F6A8885A308D3',int64)
    if(system%nion<1.or.system%nspin<1.or.system%no<1.or.system%nk<1.or.any(grid_num<1).or.projector_count<0.or.&
      pp_fingerprint==0_int64.or..not.allocated(system%Rion).or..not.allocated(system%kion).or.&
      .not.allocated(system%vec_k).or..not.allocated(system%wtk).or..not.allocated(system%rocc).or.&
      size(system%Rion,1)/=3.or.size(system%Rion,2)/=system%nion.or.size(system%kion)/=system%nion.or.&
      size(system%vec_k,1)/=3.or.size(system%vec_k,2)/=system%nk.or.size(system%wtk)/=system%nk.or.&
      size(system%rocc,1)/=system%no.or.size(system%rocc,2)/=system%nk.or.&
      size(system%rocc,3)/=system%nspin)then
      hash=0_int64;return
    endif
    call mix_i(hash,int(z'53595354454D5634',int64))
    call mix_i(hash,dg_sha256_schema)
    call mix_i(hash,int(system%nion,int64));call mix_i(hash,int(system%nspin,int64))
    call mix_i(hash,int(system%no,int64));call mix_i(hash,int(system%nk,int64))
    call mix_i(hash,int(system%ngrid,int64));call mix_i(hash,int(projector_count,int64))
    call mix_i(hash,pp_fingerprint);call mix_l(hash,periodic);call mix_l(hash,spinorbit)
    call mix_l(hash,plus_u);call mix_l(hash,hse);call mix_l(hash,fix_func);call mix_l(hash,jm)
    do i=1,3
      call mix_i(hash,int(grid_num(i),int64));call mix_r(hash,system%hgs(i))
      do j=1,3;call mix_r(hash,system%primitive_a(j,i));enddo
    enddo
    do i=1,system%nion
      call mix_i(hash,int(system%kion(i),int64))
      do j=1,3;call mix_r(hash,system%Rion(j,i));enddo
    enddo
    do i=1,system%nk
      do j=1,3;call mix_r(hash,system%vec_k(j,i));enddo
      call mix_r(hash,system%wtk(i))
    enddo
    do j=1,system%nspin
      do i=1,system%nk
        call mix_occupations(hash,system%rocc(:,i,j))
      enddo
    enddo
    call mix_i(hash,int(size(xctype),int64))
    do i=1,size(xctype);call mix_i(hash,int(xctype(i),int64));enddo
    if(hash==0_int64)hash=1_int64
  contains
    pure subroutine mix_i(target,value)
      integer(int64),intent(inout)::target
      integer(int64),intent(in)::value
      call dg_sha256_mix_int64(target,value)
    end subroutine mix_i
    pure subroutine mix_r(target,value)
      integer(int64),intent(inout)::target
      real(real64),intent(in)::value
      call dg_sha256_mix_real64(target,value)
    end subroutine mix_r
    pure subroutine mix_l(target,value)
      integer(int64),intent(inout)::target
      logical,intent(in)::value
      call dg_sha256_mix_logical(target,value)
    end subroutine mix_l
    pure subroutine mix_occupations(target,values)
      integer(int64),intent(inout)::target
      real(real64),intent(in)::values(:)
      integer::orbital
      do orbital=1,size(values);call mix_r(target,values(orbital));enddo
    end subroutine mix_occupations
  end function fingerprint_rt_dg_hybrid_system
end module rt_dg_hybrid_system_identity
