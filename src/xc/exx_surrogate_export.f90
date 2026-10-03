! Small portable strict B snapshot. Caller must invoke on one designated rank.
! Complete marker makes truncated output invalid to the importer.
module exx_surrogate_export
  use iso_fortran_env,only:real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::surrogate_write_endpoint
contains
  logical function valid_sha(value)result(ok)
    character(*),intent(in)::value
    integer::i
    ok=.false.
    if(len_trim(value)/=64)return
    do i=1,64
      if(index('0123456789abcdef',value(i:i))==0)return
    enddo
    ok=.true.
  end function

  subroutine surrogate_write_endpoint(path,step,epoch,dt,b,strict,accepted,source_sha,q_sha,status)
    character(*),intent(in)::path,source_sha,q_sha
    integer,intent(in)::step,epoch
    real(real64),intent(in)::dt
    complex(real64),intent(in)::b(:,:)
    logical,intent(in)::strict,accepted
    integer,intent(out)::status
    integer::unit,ios,i,j,n,close_status
    real(real64)::scale
    status=1;n=size(b,1)
    if(.not.strict.or..not.accepted)return
    if(step<0.or.epoch<0.or.n<1.or.size(b,2)/=n)return
    if(.not.ieee_is_finite(dt).or.dt<=0)return
    if(.not.valid_sha(source_sha).or..not.valid_sha(q_sha))return
    do j=1,n
      do i=1,n
        if(.not.ieee_is_finite(real(b(i,j),real64)).or..not.ieee_is_finite(aimag(b(i,j))))return
      enddo
    enddo
    scale=max(maxval(abs(b)),1d0)
    if(maxval(abs(b-conjg(transpose(b))))>1d-12*scale)return
    open(newunit=unit,file=path,status='new',action='write',form='formatted',iostat=ios)
    if(ios/=0)return
    write(unit,'(a)',iostat=ios)'SALMON_STRICT_B_V1'
    if(ios==0)write(unit,*,iostat=ios)step,epoch,n,dt
    if(ios==0)write(unit,'(a)',iostat=ios)source_sha
    if(ios==0)write(unit,'(a)',iostat=ios)q_sha
    do j=1,n
      do i=1,n
        if(ios==0)write(unit,'(2(es26.17e3,1x))',iostat=ios)real(b(i,j),real64),aimag(b(i,j))
      enddo
    enddo
    if(ios==0)write(unit,'(a)',iostat=ios)'END_STRICT_B'
    close(unit,iostat=close_status)
    if(ios==0.and.close_status==0)status=0
    ! Failed output remains for diagnosis and must never be overwritten/reused.
  end subroutine
end module
