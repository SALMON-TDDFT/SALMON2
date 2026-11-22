! Gamma-cell closure from finite-q density response, never from raw current.
module hsex_density_gamma
  implicit none
  private
  public :: density_gamma_stencil,density_gamma_replace
contains
  subroutine density_gamma_stencil(mesh,rec,selected,status)
    use hsex_screen_math, only: screen_sites,screen_finite
    implicit none
    integer,intent(in) :: mesh(3)
    real(8),intent(in) :: rec(3,3)
    integer,allocatable,intent(out) :: selected(:)
    integer,intent(out) :: status
    integer,allocatable :: sites(:,:)
    real(8),allocatable :: radii(:)
    real(8) :: nearest(52),radius
    integer :: nq,i,j,nvalid,code
    status=1
    if(any(mesh<6).or..not.screen_finite(rec))return
    nq=product(mesh)
    allocate(sites(3,nq),radii(nq),stat=code)
    if(code/=0)return
    call screen_sites(mesh,sites)
    radii=huge(1d0);nearest=huge(1d0);nvalid=0
    do i=1,nq
      if(all(sites(:,i)==0).or.any(2*abs(sites(:,i))>=mesh))cycle
      radius=sum(matmul(rec,dble(sites(:,i))/mesh)**2)
      if(.not.screen_finite([radius]).or.radius<=0d0)return
      radii(i)=radius;nvalid=nvalid+1
      if(radius>=nearest(52))cycle
      j=52
      do while(j>1)
        if(nearest(j-1)<=radius)exit
        nearest(j)=nearest(j-1);j=j-1
      enddo
      nearest(j)=radius
    enddo
    if(nvalid<52)return
    radius=nearest(52)*(1d0+1d-10)
    nvalid=count(radii<=radius)
    allocate(selected(nvalid),stat=code)
    if(code/=0)return
    j=0
    do i=1,nq
      if(radii(i)>radius)cycle
      j=j+1;selected(j)=i
    enddo
    status=0
  end subroutine

  subroutine density_gamma_replace(coarse,first,na,mesh,rec,volume,tolerance,icomm,rank,status)
    use hsex_screen_math, only: screen_finite,screen_q2_head
    use communication, only: comm_get_max,comm_summation,comm_bcast,comm_is_root
    implicit none
    complex(8),intent(inout) :: coarse(:,:)
    integer,intent(in) :: first,na,mesh(3),icomm,rank
    real(8),intent(in) :: rec(3,3),volume,tolerance
    integer,intent(out) :: status
    integer,allocatable :: selected(:),owned(:),owners(:)
    complex(8),allocatable :: samples(:,:,:),local(:,:),total(:,:),point(:,:),average(:,:)
    real(8) :: qc
    integer :: code,j,s,a,b,nlocal,nq,stencil_status
    status=1;code=0;nlocal=size(coarse,1);nq=product(mesh)
    if(na<1.or.first<1.or.size(coarse,2)/=nq)code=1
    if(int(na,8)*na>huge(1))code=1
    if(.not.screen_finite([volume,tolerance]).or..not.screen_finite(rec))code=1
    if(volume<=0d0.or.tolerance<=0d0.or.tolerance>=1d0)code=1
    if(.not.screen_finite(coarse))code=1
    call comm_get_max(code,icomm)
    if(code/=0)return
    if(nlocal>0.and.first+nlocal-1>na*na)code=1
    call density_gamma_stencil(mesh,rec,selected,stencil_status)
    code=max(code,stencil_status);call comm_get_max(code,icomm)
    if(code/=0)return
    allocate(owned(na*na),owners(na*na),local(na,na),total(na,na),point(na,na),average(na,na),stat=code)
    call comm_get_max(code,icomm)
    if(code/=0)return
    owned=0
    do j=1,nlocal
      owned(first+j-1)=1
    enddo
    call comm_summation(owned,owners,size(owned),icomm)
    if(any(owners/=1))code=1
    call comm_get_max(code,icomm)
    if(code/=0)return
    if(comm_is_root(rank))allocate(samples(na,na,size(selected)),stat=code)
    call comm_get_max(code,icomm)
    if(code/=0)return
    do s=1,size(selected)
      local=0d0
      do j=1,nlocal
        a=mod(first+j-2,na)+1;b=(first+j-2)/na+1
        local(a,b)=coarse(j,selected(s))
      enddo
      call comm_summation(local,total,na*na,icomm,dest=0)
      if(comm_is_root(rank))samples(:,:,s)=total
    enddo
    status=0
    if(comm_is_root(rank))then
      qc=(6d0*acos(-1d0)**2/(volume*nq))**(1d0/3d0)
      call screen_q2_head(samples,mesh,rec,mesh,tolerance,qc,point,average,status,selected)
      write(*,'(a,2i8)')'HSEX_DENSITY_GAMMA_STENCIL samples/status ',size(selected),status
    endif
    call comm_bcast(status,icomm,0)
    if(status/=0)return
    call comm_bcast(average,icomm,0)
    ! Commit only after every check passes. No separate scalar head is added.
    do j=1,nlocal
      a=mod(first+j-2,na)+1;b=(first+j-2)/na+1
      coarse(j,1)=average(a,b)
    enddo
  end subroutine
end module
