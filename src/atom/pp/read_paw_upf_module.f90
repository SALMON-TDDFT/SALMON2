!--------------------------------------------
! Unified Pseudopotential Format (UPF)
! This potential is adopted in QUANTUM ESPRESSO 
! The unit of energy is in Rydberg, and converted to Hartree in this routine
!--------------------------------------------
module read_paw_upf_module

  use structures,only : s_pp_info
  use paw_variables,only : paw
  use radial_1d
  use builtin_pz, only: PZxc

  implicit none

  private
  public :: read_paw_upf
  public :: get_mr_norb_upf

contains

  subroutine read_paw_upf( g, pp, ik )
    implicit none
    integer,intent(in) :: g, ik
    type(s_pp_info),intent(inout) :: pp
    character(200) :: cbuf
    integer :: nrr, norb, nelem, i, j, l, n, ll,l0
    integer :: nwfc
    real(8),allocatable :: work(:),qij(:)
    real(8),allocatable :: vh(:),vxc(:)
    real(8) :: tmp,exc,trho,pi4,tmp2,tmp1

    nrr=maxval( pp%mr )
    norb=maxval( pp%num_orb )
    nelem=size( pp%mr )

    write(*,*) "nrr,norb,nelem",nrr,norb,nelem
    allocate( qij(nrr) ); qij=0.0d0

    allocate( paw%rad(nrr,nelem)       ); paw%rad=0.0d0
    allocate( paw%rab(nrr,nelem)       ); paw%rab=0.0d0
    allocate( paw%aewf(nrr,norb,nelem) ); paw%aewf=0.0d0
    allocate( paw%pswf(nrr,norb,nelem) ); paw%pswf=0.0d0
    allocate( paw%beta(nrr,norb,nelem) ); paw%beta=0.0d0
    allocate( paw%occ(norb,nelem)      ); paw%occ=0.0d0
    allocate( paw%rhops(nrr,nelem)     ); paw%rhops=0.0d0
    allocate( paw%rhoae(nrr,nelem)     ); paw%rhoae=0.0d0
    allocate( paw%rhocps(nrr,nelem)    ); paw%rhocps=0.0d0
    allocate( paw%rhocae(nrr,nelem)    ); paw%rhocae=0.0d0
    allocate( paw%vlocps(nrr,nelem)    ); paw%vlocps=0.0d0
    allocate( paw%vlocae(nrr,nelem)    ); paw%vlocae=0.0d0
    allocate( paw%orb_2_l(norb,nelem)  ); paw%orb_2_l=-1000
    allocate( paw%orb_2_n(norb,nelem)  ); paw%orb_2_n=0

    call find_tag( g, "z_valence", cbuf )
    call get_dnum_from_string( cbuf, "z_valence=", pp%zion )
    pp%zps(ik)=nint(pp%zion)
    write(*,*) cbuf(1:20)

    call find_tag( g, "l_max", cbuf )
    call get_num_from_string( cbuf, "l_max=", pp%mlps(ik) )
    write(*,*) cbuf(1:20)


    allocate( paw%DIJ(norb,norb,nelem) ); paw%DIJ=0.0d0
    allocate( paw%QIJL(nrr,norb,norb,0:4,nelem) ); paw%QIJL=0.0d0
    allocate( paw%Q(norb*norb,nelem) ); paw%Q=0.0d0
    allocate( paw%MULTIPOLES(norb*norb*5,nelem) ); paw%MULTIPOLES=0

    call find_tag( g, "mesh_size", cbuf )
    call get_num_from_string( cbuf, "mesh_size=", nrr )
    write(*,*) cbuf(1:20)

    call find_tag( g, "number_of_wfc", cbuf )
    call get_num_from_string( cbuf, "number_of_wfc=", nwfc )
    write(*,*) cbuf(1:20)

    call find_tag( g, "number_of_proj", cbuf )
    call get_num_from_string( cbuf, "number_of_proj=", norb )
    write(*,*) cbuf(1:20)

    call find_tag( g, "<PP_R ", cbuf )
    read(g,*) pp%rad(1:nrr,ik)
    write(*,*) cbuf(1:20)

    paw%rad(1:nrr,ik)=pp%rad(1:nrr,ik)

    call find_tag( g, "<PP_RAB ", cbuf )
    read(g,*) paw%rab(1:nrr,ik)
    write(*,*) cbuf(1:20)

    call find_tag( g, "<PP_NLCC ", cbuf )
    read(g,*) paw%rhocps(1:nrr,ik)
    write(*,*) cbuf(1:20)

    call find_tag( g, "<PP_LOCAL ", cbuf )
    read(g,*) paw%vlocps(1:nrr,ik)
    write(*,*) cbuf(1:20)

    call find_tag( g, "<PP_NONLOCAL>", cbuf )
    write(*,*) cbuf(1:20)
    do i=1,norb
       read(g,'(a)') cbuf
       write(*,*) cbuf
       call get_num_from_string( cbuf, "angular_momentum=", l )
       paw%orb_2_l(i,ik)=l
       paw%orb_2_n(i,ik)=count(paw%orb_2_l(:,ik)==l)
       write(*,*) i,paw%orb_2_l(i,ik),paw%orb_2_n(i,ik)
       do
          read(g,'(a)') cbuf
          if ( index(cbuf,">") /= 0 ) exit
       end do
       read(g,*) paw%beta(1:nrr,i,ik)
       read(g,'(a)') cbuf
    end do ! i

    call find_tag( g, "<PP_DIJ ", cbuf ); write(*,*) cbuf(1:20)
    read(g,*) paw%DIJ(1:norb,1:norb,ik)

    call find_tag( g, "<PP_Q ", cbuf ); write(*,*) cbuf(1:20)
    read(g,*) paw%Q(:,ik)

    call find_tag( g, "<PP_MULTIPOLES ", cbuf ); write(*,*) cbuf(1:20)
    read(g,*) paw%MULTIPOLES(1:norb*norb*(pp%mlps(ik)+1),ik)

    do
       call find_tag( g, "<PP_QIJL", cbuf ); write(*,*) cbuf(1:20)
       call get_num_from_string( cbuf, "first_index=", i )
       call get_num_from_string( cbuf, "second_index=", j )
       call get_num_from_string( cbuf, "angular_momentum=", l )
       read(g,*) paw%QIJL(1:nrr,i,j,l,ik)
       read(g,'(a)') cbuf ;write(*,*) "aaaa",cbuf(1:20)
       read(g,'(a)') cbuf ;write(*,*) "aaa",cbuf(1:20)
       if ( index(cbuf,"</PP_AUGMENTATION") /= 0 ) then
          exit
       else
          backspace(g)
       end if
    end do

    call find_tag( g, "<PP_PSWFC", cbuf ); write(*,*) cbuf(1:20)
    do i=1,nwfc
       call find_tag( g, "<PP_CHI", cbuf ); write(*,*) cbuf(1:20)
       do
          read(g,'(a)') cbuf
          if ( index(cbuf,">") /= 0 ) exit
       end do
       read(g,*) pp%upp_f(1:nrr,i,ik)   
       call find_tag( g, "</PP_CHI", cbuf )       
    end do

    call find_tag( g, "<PP_FULL_WFC ", cbuf ); write(*,*) cbuf(1:20)
    do i=1,norb
       call find_tag( g, "<PP_AEWFC", cbuf ); write(*,*) cbuf(1:20)
       read(g,*) paw%aewf(1:nrr,i,ik)
       call find_tag( g, "</PP_AEWFC", cbuf )
    end do
    do i=1,norb
       call find_tag( g, "<PP_PSWFC", cbuf ); write(*,*) cbuf(1:20)
       read(g,*) paw%pswf(1:nrr,i,ik)
       call find_tag( g, "</PP_PSWFC", cbuf )
    end do

    call find_tag( g, "<PP_RHOATOM ", cbuf ); write(*,*) cbuf(1:20)
    read(g,*) paw%rhops(1:nrr,ik)

    call find_tag( g, "<PP_OCCUPATIONS ", cbuf ); write(*,*) cbuf(1:20)
    read(g,*) paw%occ(1:norb,ik)

    call find_tag( g, "<PP_AE_NLCC ", cbuf ); write(*,*) cbuf(1:20)
    read(g,*) paw%rhocae(1:nrr,ik)

    call find_tag( g, "<PP_AE_VLOC ", cbuf ); write(*,*) cbuf(1:20)
    read(g,*) paw%vlocae(1:nrr,ik)

!---

    pp%lref(ik) = ubound(pp%vpp,2)

    if ( paw%rad(1,ik) == 0.0d0 ) then
      call stop_program('stop@read_paw_upf(1)')
    else
      pp%mr(ik) = nrr + 1
      do i=2,nrr+1
        pp%rad(i,ik) = paw%rad(i-1,ik)
      end do
      pp%rad(1,ik) = 0.0d0
      do j=1,norb
        do i=1,nrr
          pp%vpp(i,j-1) = sqrt(0.5d0)*paw%beta(i,j,ik)
        end do
        pp%vpp(0,j-1) = pp%vpp(1,j-1)
      end do
      do i=1,nrr
        pp%vpp(i,pp%lref(ik)) = 0.5d0*paw%vlocps(i,ik)
      end do
      pp%vpp(0,pp%lref(ik)) = pp%vpp(1,pp%lref(ik))
    end if

    if ( all(paw%rab(2:nrr,ik)==paw%rab(1,ik)) ) then
      tmp = pp%rad(2,ik) - pp%rad(1,ik)
      do i = pp%mr(ik)+1, pp%nrmax
        pp%rad(i,ik) = (i-1)*tmp
      end do
    else
      tmp1 = log(paw%rad(1,ik))
      tmp2 = log(paw%rad(2,ik)) - log(paw%rad(1,ik))
      do i = pp%mr(ik)+1,pp%nrmax
        tmp = tmp1 + (i-2)*tmp2
        pp%rad(i,ik) = exp(tmp)
      end do
    end if

    pp%nrps(ik)=0
    do j=1,norb
      do i=nrr,1,-1
        if ( paw%beta(i,j,ik) /= 0.0d0 ) then
          pp%nrps(ik) = max( pp%nrps(ik), i )
          exit
        end if
      end do
    end do

    pp%rps(ik) = pp%rad(pp%nrps(ik),ik)
    write(*,*) 'pp%nrps,pp%rps',pp%nrps(ik), pp%rps(ik)

    write(*,*) "norb=",norb

    pp%nproj(:,ik)=0
    do j=1,norb
      l=paw%orb_2_l(j,ik)
      pp%nproj(l,ik) = pp%nproj(l,ik) + 1
    end do

    pp%mlps(ik) = maxval( paw%orb_2_l(:,ik) )

    write(*,*) "pp%mlps=",pp%mlps(ik)

    i=0
    l0=0
    do ll = 0, pp%mlps(ik)
    do l = l0, l0+pp%nproj(ll,ik)-1
      i=i+1
      write(*,*) i,ll
    end do
    end do

!---

    l=maxval( paw%orb_2_l(:,ik) )
    n=maxval( paw%orb_2_n(:,ik) )
    allocate( paw%nl_2_orb(0:l,n,ik) ); paw%nl_2_orb=0
    do i=1,norb
       l=paw%orb_2_l(i,ik)
       n=paw%orb_2_n(i,ik)
       paw%nl_2_orb(l,n,ik)=i
    end do

!---

    paw%rhoae(:,ik)=0.0d0
    do j=1,norb
       paw%rhoae(:,ik)=paw%rhoae(:,ik)+paw%occ(j,ik)*paw%aewf(:,j,ik)**2
    end do

!---

    allocate( work(nrr) ); work=0.0d0

    rewind 10
    do i=1,nrr
       write(10,'(1x,5g20.10)') pp%rad(i,ik),(paw%beta(i,j,ik),j=1,norb)
    end do

    rewind 11
    do i=1,nrr
       write(11,'(1x,5g20.10)') pp%rad(i,ik),paw%vlocps(i,ik),paw%vlocae(i,ik)
    end do

    rewind 12
    do i=1,nrr
       write(12,'(1x,5g20.10)') pp%rad(i,ik),(paw%pswf(i,j,ik),j=1,norb)
    end do

    rewind 13
    do i=1,nrr
       write(13,'(1x,5g20.10)') pp%rad(i,ik),(paw%aewf(i,j,ik),j=1,norb)
    end do

    rewind 14
    do i=1,nrr
       write(14,'(1x,5g20.10)') pp%rad(i,ik),paw%rhops(i,ik),paw%rhoae(i,ik)
    end do

    work=paw%aewf(:,1,ik)**2-paw%pswf(:,1,ik)**2
    rewind 15
    do i=1,nrr
       write(15,'(1x,5g20.10)') pp%rad(i,ik),qij(i),work(i)
    end do
    
    pi4=4.0d0*acos(-1.0d0)
    do l=0,2*pp%mlps(ik)
    do j=1,norb
    do i=1,norb
    write(*,'(1x,"qijl",3i3,es20.10)') i,j,l,sum(paw%QIJL(:,i,j,l,ik)*paw%rad(:,ik)**l*paw%rab(:,ik))
    end do
    end do
    end do
    write(*,*) '|aewf|**2',sum(paw%aewf(:,1,ik)**2*paw%rab(:,ik))
    write(*,*) '|aewf|**2',sum(paw%aewf(:,3,ik)**2*paw%rab(:,ik))
    write(*,*) '|pswf|**2',sum(paw%pswf(:,3,ik)**2*paw%rab(:,ik))
    write(*,*) "rhops",sum(paw%rhops(:,ik)*paw%rab(:,ik))
    write(*,*) "rhoae",sum(paw%rhoae(:,ik)*paw%rab(:,ik))
    write(*,*) "rhocps",sum(paw%rhocps(:,ik)*paw%rab(:,ik)) &
         ,sum(paw%rhocps(:,ik)*pp%rad(:,ik)**2*paw%rab(:,ik))*pi4
    write(*,*) "rhocae",sum(paw%rhocae(:,ik)*paw%rab(:,ik)) &
         ,sum(paw%rhocae(:,ik)*pp%rad(:,ik)**2*paw%rab(:,ik))*pi4

    do j=1,norb
       work=work+paw%occ(j,ik)*paw%pswf(:,j,ik)**2
    end do
    write(*,*) 'sum(density_from_pswf)',sum(work*paw%rab(:,ik))

    do j=1,norb-1
       write(*,*) '|beta*pswf|',sum(paw%beta(:,j+1,ik)*paw%pswf(:,j,ik)*paw%rab(:,ik))
    end do

    allocate( vh(nrr)  ); vh=0.0d0
    allocate( vxc(nrr) ); vxc=0.0d0

    call calc_hartree_radial_1d &
         ( pp%rad(1:nrr,ik),paw%rab(1:nrr,ik),paw%rhocae(1:nrr,ik),vh,tmp )

    exc=0.0d0
    do i=1,nrr
!       trho=paw%rhops(i,ik)/pp%rad(i,ik)**2/pi4 +paw%rhocps(i,ik)
       trho=paw%rhocae(i,ik)
       call PZxc( trho, tmp, tmp2 )
       exc=exc+tmp*trho*paw%rab(i,ik)*pp%rad(i,ik)**2*pi4
       vxc(i)=tmp+trho*tmp2
    end do
    write(*,*) "exc=",exc

    rewind 40
    do i=1,nrr
       write(40,'(1x,5g20.10)') pp%rad(i,ik),vh(i),vxc(i),paw%vlocps(i,ik)
    end do
    rewind 41
    do i=1,nrr
       write(41,'(1x,5g20.10)') pp%rad(i,ik),vh(i),vxc(i) &
            ,-6.0d0/pp%rad(i,ik),paw%vlocae(i,ik)
    end do

  end subroutine read_paw_upf


  subroutine get_num_from_string( buf, str, n )
    implicit none
    character(*),intent(in) :: buf, str
    integer,intent(out) :: n
    integer :: i, j, l, len_str
    character(8) :: str_out
    l = index( buf, str )
    if ( l /= 0 ) then
       len_str = len_trim( adjustl(str) )
       i=l+len_str+1
       j=index( buf(i:), '"' ) + i - 2
       str_out = buf(i:j)
       read(str_out,*) n
    end if
  end subroutine get_num_from_string

  subroutine get_dnum_from_string( buf, str, d )
    implicit none
    character(*),intent(in) :: buf, str
    real(8),intent(out) :: d
    integer :: i, j, l, len_str
    character(8) :: str_out
    l = index( buf, str )
    if ( l /= 0 ) then
       len_str = len_trim( adjustl(str) )
       i=l+len_str+1
       j=index( buf(i:), '"' ) + i - 2
       str_out = buf(i:j)
       read(str_out,*) d
    end if
  end subroutine get_dnum_from_string


  subroutine get_mr_norb_upf( g, mr, norb )
    implicit none
    integer,intent(in) :: g
    integer,intent(out) :: mr, norb
    character(200) :: cbuf
    do
       call find_tag( g, "mesh_size", cbuf )
       call get_num_from_string( cbuf, "mesh_size=", mr )
       call find_tag( g, "number_of_proj", cbuf )
       call get_num_from_string( cbuf, "number_of_proj=", norb )
       exit
    end do
  end subroutine get_mr_norb_upf


  subroutine find_tag( g, indx, cbuf )
    implicit none
    integer,intent(in) :: g
    character(*),intent(in) :: indx
    character(*),intent(out) :: cbuf
    integer :: l
    l=len(indx)
    do
       read(g,'(a)',END=99) cbuf
       cbuf=adjustl(cbuf)
       if ( cbuf(1:l) == indx ) return
    end do
99  continue
  end subroutine find_tag


end module read_paw_upf_module
