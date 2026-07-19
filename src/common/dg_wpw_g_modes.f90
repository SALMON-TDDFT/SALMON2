module dg_wpw_g_modes
  use,intrinsic::iso_fortran_env,only:int64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  public::select_dg_wpw_g_modes
contains
  subroutine select_dg_wpw_g_modes(box_length,energy_cutoff,target_count,indices,g_vectors,info)
    real(8),intent(in)::box_length(3),energy_cutoff
    integer,intent(in)::target_count
    integer,allocatable,intent(out)::indices(:,:)
    real(8),allocatable,intent(out)::g_vectors(:,:)
    integer,intent(out)::info
    integer::nk(3),ix,iy,iz,ncandidate,nselected,i,j,key,tmp_key,tmp_vec(3)
    integer::shell_start,shell_end,previous_total,next_total,nkeep
    integer,allocatable::candidate(:,:),keys(:)
    integer(int64)::candidate64
    real(8)::g(3),kcut
    real(8),parameter::pi=acos(-1d0)
    info=0
    if(any(box_length<=0d0).or..not.all(ieee_is_finite(box_length)).or.&
       energy_cutoff<=0d0.or..not.ieee_is_finite(energy_cutoff).or.target_count<=0)then
      info=1;return
    endif
    kcut=sqrt(2d0*energy_cutoff)
    nk=ceiling(kcut*box_length/(2d0*pi))+1
    candidate64=int(2*nk(1)+1,int64)*int(2*nk(2)+1,int64)*int(2*nk(3)+1,int64)
    if(candidate64>int(huge(ncandidate),int64))then;info=2;return;endif
    ncandidate=int(candidate64);allocate(candidate(3,ncandidate));nselected=0
    do iz=-nk(3),nk(3);do iy=-nk(2),nk(2);do ix=-nk(1),nk(1)
      if(ix==0.and.iy==0.and.iz==0)cycle
      g=2d0*pi*[dble(ix)/box_length(1),dble(iy)/box_length(2),dble(iz)/box_length(3)]
      if(0.5d0*sum(g*g)<=energy_cutoff*(1d0+32d0*epsilon(1d0)))then
        nselected=nselected+1;candidate(:,nselected)=[ix,iy,iz]
      endif
    enddo;enddo;enddo
    if(nselected==0)then;allocate(indices(3,0),g_vectors(3,0));return;endif
    allocate(keys(nselected))
    do i=1,nselected;keys(i)=sum(candidate(:,i)*candidate(:,i));enddo
    do i=2,nselected
      do j=i,2,-1
        if(keys(j)<keys(j-1).or.(keys(j)==keys(j-1).and.cubic_order_less(candidate(:,j),candidate(:,j-1))))then
          tmp_key=keys(j);keys(j)=keys(j-1);keys(j-1)=tmp_key
          tmp_vec=candidate(:,j);candidate(:,j)=candidate(:,j-1);candidate(:,j-1)=tmp_vec
        else;exit;endif
      enddo
    enddo
    nkeep=0;shell_start=1
    do while(shell_start<=nselected)
      shell_end=shell_start
      do while(shell_end<nselected.and.keys(shell_end+1)==keys(shell_start));shell_end=shell_end+1;enddo
      previous_total=nkeep;next_total=shell_end
      if(next_total<=target_count)then
        nkeep=next_total
      else
        if(previous_total==0.or.abs(next_total-target_count)<abs(previous_total-target_count))nkeep=next_total
        exit
      endif
      shell_start=shell_end+1
    enddo
    allocate(indices(3,nkeep),g_vectors(3,nkeep));indices=candidate(:,1:nkeep)
    do i=1,nkeep
      g_vectors(:,i)=2d0*pi*[dble(indices(1,i))/box_length(1),dble(indices(2,i))/box_length(2),&
        dble(indices(3,i))/box_length(3)]
    enddo
  end subroutine
  logical function cubic_order_less(a,b)result(less)
    integer,intent(in)::a(3),b(3);integer::sa,sb
    less=.false.;sa=sum(abs(a));sb=sum(abs(b))
    if(sa/=sb)then;less=sa<sb
    else if(abs(a(3))/=abs(b(3)))then;less=abs(a(3))<abs(b(3))
    else if(abs(a(2))/=abs(b(2)))then;less=abs(a(2))<abs(b(2))
    else if(abs(a(1))/=abs(b(1)))then;less=abs(a(1))<abs(b(1))
    else if(a(3)/=b(3))then;less=a(3)<b(3)
    else if(a(2)/=b(2))then;less=a(2)<b(2)
    else if(a(1)/=b(1))then;less=a(1)<b(1)
    endif
  end function
end module
