!--------------------------------------------
! Unified Pseudopotential Format (UPF)
! This potential is adopted in QUANTUM ESPRESSO 
! The unit of energy is in Rydberg, and converted to Hartree in this routine
!--------------------------------------------
module read_ps_upf_module

  implicit none

  private
  public :: read_ps_upf

contains

  subroutine read_ps_upf( pp,rrc,rhor_nlcc,flag_nlcc_element,ik,ps_file )
    use structures, only: s_pp_info
    implicit none
    integer,intent(in) :: ik
    type(s_pp_info) :: pp
    real(8),intent(out) :: rrc(0:)
    real(8),intent(out) :: rhor_nlcc(0:,0:)
    logical,intent(inout) :: flag_nlcc_element(:)
    character(*),intent(in) :: ps_file
    character(30) :: cbuf
    integer,parameter :: g=4

    open(g,file=ps_file,status='old')

    read(g,'(a)',END=10) cbuf
    cbuf = adjustl( cbuf )
    !write(*,*) cbuf

    if ( cbuf(1:21) == '<UPF version="2.0.1">' ) then
       rewind g
       call read_ps_upf_ver201( g, pp, rrc, rhor_nlcc, ik )
    else
       write(*,*) "This may be an old UPF format"
       stop 'stop@read_ps_upf'
    end if

    if ( any(rhor_nlcc/=0.0d0) ) flag_nlcc_element(ik)=.true.

    return
10  stop 'Format is invalid (stop@read_ps_upf)'

    close(g)

  end subroutine read_ps_upf


  subroutine read_ps_upf_ver201( g, pp, rrc, rhor_nlcc, ik )
    use structures, only: s_pp_info
    implicit none
    integer,intent(in) :: g
    type(s_pp_info),intent(inout) :: pp
    real(8),intent(out) :: rrc(0:)
    real(8),intent(out) :: rhor_nlcc(0:,0:)
    integer,intent(in) :: ik
    integer,parameter :: max_loop=1000000, max_array_size = 16
    integer :: loop,i,j,k,l,ir,ic,nr,l0,ll
    integer,allocatable :: lo(:),no(:)
    character(100) :: cbuf, ckey
    integer :: norb,nrr,nsize,ltmp,columns
    real(8) :: Zps,r,dx,x1,x
    real(8),allocatable :: rr(:),rx(:),vql(:),cdc(:),cdd(:)
    real(8),allocatable :: viod(:,:),anorm(:),Dij(:,:),work(:)
    logical :: flag_spin_orb = .false.
    integer,allocatable :: i2l(:)

    write(*,'(a40," read_ps_upf_ver201")') repeat("-",40)

    nrr=0
    norb=0

    do loop=1,max_loop

       read(g,'(a)',END=10) cbuf
       ckey = adjustl( cbuf )

       if ( ckey(1:6) == "</UPF>" ) then
          write(*,*) ckey(1:6)
          exit
       end if

       if ( ckey(1:11) == "<PP_HEADER " ) then

          write(*,*) ckey(1:11)

          backspace(g)

          do i=1,max_loop

             read(g,'(a)') cbuf

             j = index( cbuf, "z_valence=" )
             if ( j /= 0 ) then
                j = j + len_trim( "z_valence=" )
                k = index( cbuf(j+1:), '"' )
                ckey=cbuf(j+1:j+k-1)
                read(ckey,*) Zps
             end if

             j = index( cbuf, "mesh_size=" )
             if ( j /= 0 ) then
                j = j + len_trim( "mesh_size=" )
                k = index( cbuf(j+1:), '"' )
                ckey=cbuf(j+1:j+k-1)
                read(ckey,*) nrr
             end if

             j = index( cbuf, "number_of_proj=" )
             if ( j /= 0 ) then
                j = j + len_trim( "number_of_proj=" )
                k = index( cbuf(j+1:), '"' )
                ckey=cbuf(j+1:j+k-1)
                read(ckey,*) norb
             end if

             j = index( cbuf, "/>" )
             if ( j /= 0 ) exit

          end do ! i

       end if ! </PP_HEADER>

       if ( nrr > 0 ) then
          if ( .not.allocated(rr) ) then
             allocate( rr(nrr)  ) ; rr =0.0d0
             allocate( rx(nrr)  ) ; rx =0.0d0
             allocate( cdc(nrr) ) ; cdc=0.0d0
             allocate( vql(nrr) ) ; vql=0.0d0
             allocate( cdd(nrr) ) ; cdd=0.0d0
          end if
       end if

       if ( ckey(1:8) == "<PP_MESH" ) then

          !write(*,*) ckey(1:8)

          do i=1,max_loop
             read(g,'(a)') cbuf
             ckey = adjustl( cbuf )
             if ( ckey(1:6) == "<PP_R " ) exit
          end do
          read(g,*) rr(1:nrr)

          do i=1,max_loop
             read(g,'(a)') cbuf
             ckey = adjustl( cbuf )
             if ( ckey(1:8) == "<PP_RAB " ) exit
          end do
          read(g,*) rx(1:nrr)

       end if ! </PP_MESH>

       if ( ckey(1:9) == "<PP_NLCC " ) then

          !write(*,*) ckey(1:9)

          read(g,*) cdc(1:nrr)

       end if

       if ( ckey(1:10) == "<PP_LOCAL " ) then

          !write(*,*) ckey(1:10)

          read(g,*) vql(1:nrr)

       end if ! </PP_LOCAL>

       if ( nrr > 0 ) then
          if ( .not.allocated(lo) ) then
             allocate( lo(max_array_size)       ) ; lo=0
             allocate( no(0:max_array_size)     ) ; no=0
             allocate( viod(nrr,max_array_size) ) ; viod=0.0d0
          end if
       end if

       if ( ckey(1:13) == "<PP_NONLOCAL>" ) then

          !write(*,*) ckey(1:13)

          j=0
          do i=1,max_loop

             read(g,'(a)') cbuf
             ckey = adjustl( cbuf )

             if ( ckey(1:9) == "<PP_BETA." ) then

                do k=1,max_loop

                   call get_num_from_string( cbuf, "angular_momentum=", ltmp )

                   l=index( cbuf, ">" )
                   if ( l /= 0 ) then
                      j=j+1
                      read(g,*) viod(1:nrr,j)
                      lo(j)=ltmp
                      no(lo(j)) = no(lo(j)) + 1
                      exit
                   end if

                   read(g,'(a)') cbuf

                end do ! k

                if ( j == norb ) exit

             else if ( ckey(1:8) == "<PP_DIJ " ) then
                backspace(g)
                !write(*,*) ckey,"exit"
                exit
             end if

          end do ! i

          if ( norb > 0 ) then
             if ( .not.allocated(Dij) ) then
                allocate( Dij(norb,norb) ) ; Dij=0.0d0
                allocate( anorm(norb)    ) ; anorm=0.0d0
             end if
          end if

          do i=1,max_loop

             read(g,'(a)') cbuf
             ckey = adjustl( cbuf )

             if ( ckey(1:8) == "<PP_DIJ " ) then
                write(*,*) ckey(1:8)
                call get_num_from_string( cbuf, "size=", nsize )
                call get_num_from_string( cbuf, "columns=", columns )
                write(*,*) "nsize,columns=",nsize,columns
                if ( nsize > 0 ) then
! ---
                allocate( work(nsize) ) ; work=0.0d0
                nr=nsize/columns
                do ir=1,nr
                   read(g,*) work((ir-1)*columns+1:ir*columns)
                end do
                if ( nr*columns < nsize ) read(g,*) work(nr*columns+1:nsize)
                do ic=1,norb
                   do ir=1,norb
                      Dij(ir,ic) = work(ir+(ic-1)*norb)
                   end do
                end do
                deallocate( work )
! ---
                j=count( Dij /= 0.0d0 )
                k=0 !------> check single or multi reference
                do l=1,norb
                   if ( Dij(l,l) /= 0.0d0 ) k=k+1
                end do
                if ( j == k ) then
                   do l=1,norb
                      anorm(l) = Dij(l,l)
                   end do
                   Dij(:,:)=0.0d0
                end if !------> check single or multi reference (end)
! ---
                end if

             else if ( ckey(1:9) == "</PP_DIJ>" ) then
                !write(*,*) ckey,"exit"
                exit
             end if

          end do ! i

       end if ! </PP_NONLOCAL>

       if ( ckey(1:12) == "<PP_RHOATOM " ) then

          !write(*,*) ckey(1:12)

          read(g,*) cdd(1:nrr)
          !exit

       end if

       if ( ckey(1:13) == "<PP_SPIN_ORB>" ) then

          flag_spin_orb = .true.

          do j=1,norb
             read(g,*) cbuf
             write(*,*) j, lo(j),cbuf
          end do

       end if

    end do ! loop
10 continue

    pp%lref(ik) = ubound(pp%vpp,2)

    pp%zion = Zps
    pp%zps(ik) = Zps
    pp%mlps(ik) = maxval( lo(1:norb) )

    pp%nproj(:,ik)=0
    do i=1,norb
       l=lo(i)
       pp%nproj(l,ik) = pp%nproj(l,ik) + 1
    end do

    if( flag_spin_orb )then
       do ll=1,pp%mlps(ik)
          pp%nproj(ll,ik) = pp%nproj(ll,ik)/2
       end do
    end if

    if( flag_spin_orb )then
       allocate( i2l(norb) ); i2l=0
       i=0
       l0=0
       do ll=0,pp%mlps(ik)
       do l=l0,l0+pp%nproj(ll,ik)-1
          i=i+1
          i2l(i)=l
          if( ll > 0 )then
             i=i+1
             i2l(i)=l
          end if
       end do
       l0=l
       end do
    end if

    if( flag_spin_orb )then
       pp%anorm(:,ik)=0.0d0
       do i=1,norb
          l=i2l(i)
          if( pp%anorm(l,ik) == 0.0d0 )then
             pp%anorm(l,ik)=anorm(i)
          else
             pp%anorm_so(l,ik)=anorm(i)
          end if
       end do
       l0=0
       do ll=0,pp%mlps(ik)
       do l=l0,l0+pp%nproj(ll,ik)-1
          write(*,*) l,ll,pp%anorm(l,ik), pp%anorm_so(l,ik)
       end do
       l0=l
       end do
    else
       do i=1,norb
          pp%anorm(i-1,ik) = anorm(i)
       end do
    end if

    if ( rr(1) == 0.0d0 ) then
       pp%mr(ik) = nrr
       do i=1,nrr
          pp%rad(i,ik)=rr(i)
       end do
       if( flag_spin_orb )then
          pp%vpp(:,:)=0.0d0
          do i=1,norb
             l=i2l(i)
             if( all(pp%vpp(:,l)==0.0d0) )then
                do j=1,nrr
                   pp%vpp(j-1,l) = sqrt(0.5d0)*viod(j,i)
                end do
             else
                do j=1,nrr
                   pp%vpp_so(j-1,l) = sqrt(0.5d0)*viod(j,i)
                end do
             end if
          end do
       else
          do j=1,norb
             do i=1,nrr
                pp%vpp(i-1,j-1) = sqrt(0.5d0)*viod(i,j)
             end do
          end do
       end if
       do i=1,nrr
          pp%vpp(i-1,pp%lref(ik)) = 0.5d0*vql(i)
       end do
       do i=1,nrr
          rhor_nlcc(i-1,0) = cdc(i)
       end do
       do i=1,nrr
          pp%rho_pp_tbl(i,ik) = cdd(i)
       end do
    else
       pp%mr(ik) = nrr + 1
       do i=2,nrr+1
          pp%rad(i,ik) = rr(i-1)
       end do
       pp%rad(1,ik)=0.0d0
       if( flag_spin_orb )then
          pp%vpp(:,:)=0.0d0
          do i=1,norb
             l=i2l(i)
             if( all(pp%vpp(:,l)==0.0d0) )then
                do j=1,nrr
                   pp%vpp(j,l) = sqrt(0.5d0)*viod(j,i)
                end do
                pp%vpp(0,l) = pp%vpp(1,l)
             else
                do j=1,nrr
                   pp%vpp_so(j,l) = sqrt(0.5d0)*viod(j,i)
                end do
                pp%vpp_so(0,l) = pp%vpp_so(1,l)
             end if
          end do
       else
          do j=1,norb
             do i=1,nrr
                pp%vpp(i,j-1) = sqrt(0.5d0)*viod(i,j)
             end do
             pp%vpp(0,j-1)=pp%vpp(1,j-1)
          end do
       end if
       do i=1,nrr
          pp%vpp(i,pp%lref(ik)) = 0.5d0*vql(i)
       end do
       pp%vpp(0,pp%lref(ik))=pp%vpp(1,pp%lref(ik))
       do i=1,nrr
          rhor_nlcc(i,0) = cdc(i)
       end do
       rhor_nlcc(0,0) = rhor_nlcc(1,0)
       do i=1,nrr
          pp%rho_pp_tbl(i+1,ik) = cdd(i)
       end do
       pp%rho_pp_tbl(1,ik) = cdd(1)
    end if

    if ( all(rx(2:nrr)==rx(1)) ) then
       dx = rr(2) - rr(1)
       do i=pp%mr(ik)+1,pp%nrmax
          pp%rad(i,ik) = (i-1)*dx
       end do
    else
       dx = log(rr(2)) - log(rr(1))
       x1 = log(rr(1))
       do i=pp%mr(ik)+1,pp%nrmax
          x = x1 + (i-2)*dx
          pp%rad(i,ik) = exp(x)
       end do
    end if

    !i=0
    !l0=0
    !do ll=0,pp%mlps(ik)
    !   do l=l0,l0+pp%nproj(ll,ik)-1
    !      i=i+1
    !      write(*,'(1x,3i4,3f20.15)') i,l,ll,pp%anorm(l,ik),sum(pp%vpp(:,l)**2)*2.0d0,sum(viod(:,i)**2) 
    !      if( ll > 0 )then
    !         i=i+1
    !         write(*,'(1x,3i4,3f20.15)') i,l,ll,pp%anorm_so(l,ik),sum(pp%vpp_so(:,l)**2)*2.0d0,sum(viod(:,i)**2) 
    !      end if
    !   end do
    !   l0=l
    !end do

    if( flag_spin_orb )then
       lo=0
       i=0
       l0=0
       do ll=0,pp%mlps(ik)
       do l=l0,l0+pp%nproj(ll,ik)-1
          lo(l+1)=ll
       end do
       l0=l
       end do
       rrc=0.0d0
       l0=0
       do ll=0,pp%mlps(ik)
       do l=l0,l0+pp%nproj(ll,ik)-1
          do i=pp%mr(ik),1,-1
             if ( abs(pp%vpp(i-1,l)) > 1.0d-6 ) then
                rrc(ll) = max( rrc(ll), pp%rad(i,ik) )
                exit
             end if
          end do
          do i=pp%mr(ik),1,-1
             if ( abs(pp%vpp_so(i-1,l)) > 1.0d-6 ) then
                rrc(ll) = max( rrc(ll), pp%rad(i,ik) )
                exit
             end if
          end do
       end do
       l0=l
       end do
    else
       rrc=0.0d0
       do j=1,norb
          l=lo(j)
          do i=pp%mr(ik),1,-1
             if ( abs(pp%vpp(i-1,j-1)) > 1.0d-6 ) then
                rrc(l) = max( rrc(l), pp%rad(i,ik) )
                exit
             end if
          end do
       end do
    end if

    if ( allocated(Dij)   ) deallocate( Dij )
    if ( allocated(anorm) ) deallocate( anorm )
    deallocate( viod, no, lo )
    deallocate( cdd, vql, cdc, rx, rr )

    write(*,*) "*** Unified Pseudopotenetial Format (UPF) ***"
    write(*,*) "Znuc=",pp%zps(ik)
    write(*,*) "max angular momentum =",pp%mlps(ik)
    write(*,*) "# of orbitals =",norb
    write(*,*) "               ",(pp%nproj(l,ik),l=0,pp%mlps(ik))
    write(*,*) "cut off radius =",(rrc(l),l=0,pp%mlps(ik))
    write(*,*) "# of radial mesh points =",pp%mr(ik)
    write(*,*) "uVu integral (anorm) ="
    write(*,'(1x,8f10.5)') ( pp%anorm(i,ik),i=0,norb-1 )
    if( flag_spin_orb )then
       write(*,*) "uVu integral (anorm_so) ="
       write(*,'(1x,8f10.5)') ( pp%anorm_so(i,ik),i=0,norb-1 )
    end if
    r=0.0d0
    do i=1,pp%mr(ik)
       r=r+rhor_nlcc(i,0)*pp%rad(i+1,ik)**2*(pp%rad(i+1,ik)-pp%rad(i,ik))
    end do
    write(*,*) "sum(rhor_nlcc)=",r, r*4.0d0*acos(-1.0d0)

    r=0.0d0
    do i=1,pp%mr(ik)
       r=r+pp%rho_pp_tbl(i,ik)*(pp%rad(i+1,ik)-pp%rad(i,ik))
    end do
    write(*,*) "sum(rho_pp_tbl)=",r
    !rewind 100
    !do i=1,pp%mr(ik)
    !  if ( pp%rad(i,ik) == 0.0d0 ) cycle
    !  write(100,*) pp%rad(i,ik), pp%rho_pp_tbl(i,ik)/(4.0d0*acos(-1.0d0)*pp%rad(i,ik)**2)
    !end do

    return
  end subroutine read_ps_upf_ver201


  subroutine get_num_from_string( buf, str, n )
    implicit none
    character(*),intent(IN) :: buf, str
    integer,intent(OUT) :: n
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


end module read_ps_upf_module
