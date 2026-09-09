module rt_dg_hybrid_system_identity
  use,intrinsic::iso_fortran_env,only:int64,real64
  use structures,only:s_dft_system
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
    if(system%nion<1.or.system%nspin<1.or.any(grid_num<1).or.projector_count<0.or.&
      pp_fingerprint==0_int64.or..not.allocated(system%Rion).or..not.allocated(system%kion).or.&
      size(system%Rion,1)/=3.or.size(system%Rion,2)/=system%nion.or.size(system%kion)/=system%nion)then
      hash=0_int64;return
    endif
    call mix_i(hash,int(z'53595354454D5634',int64))
    call mix_i(hash,int(system%nion,int64));call mix_i(hash,int(system%nspin,int64))
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
    call mix_i(hash,int(size(xctype),int64))
    do i=1,size(xctype);call mix_i(hash,int(xctype(i),int64));enddo
    if(hash==0_int64)hash=1_int64
  contains
    pure subroutine mix_i(target,value)
      integer(int64),intent(inout)::target
      integer(int64),intent(in)::value
      integer::byte
      do byte=0,7
        target=ieor(ishftc(target,7),int(ibits(value,8*byte,8),int64))
      enddo
    end subroutine mix_i
    pure subroutine mix_r(target,value)
      integer(int64),intent(inout)::target
      real(real64),intent(in)::value
      integer(int64)::bits
      bits=transfer(value,bits);call mix_i(target,bits)
    end subroutine mix_r
    pure subroutine mix_l(target,value)
      integer(int64),intent(inout)::target
      logical,intent(in)::value
      call mix_i(target,int(merge(1,0,value),int64))
    end subroutine mix_l
  end function fingerprint_rt_dg_hybrid_system
end module rt_dg_hybrid_system_identity
