!
!  Copyright 2018-2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module prep_pp_sub
  implicit none

contains

subroutine init_ps(lg,mg,system,info,fg,poisson,pp,ppg,Vpsl)
  use structures
  use hamiltonian, only: update_kvector_nonlocalpt
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root
  use salmon_global, only: iperiodic,natom,quiet,yn_spinorbit
  use prep_pp_so_sub, only: calc_uv_so
  use prep_pp_plusU_sub, only: calc_uv_plusU, PLUS_U_ON
  use timer
  implicit none
  type(s_rgrid)           ,intent(in) :: lg,mg
  type(s_dft_system)      ,intent(in) :: system
  type(s_parallel_info)   ,intent(in) :: info
  type(s_reciprocal_grid) ,intent(in) :: fg
  type(s_poisson)                     :: poisson
  type(s_pp_info)         ,intent(in) :: pp
  type(s_pp_grid)                     :: ppg
  type(s_scalar)                      :: Vpsl
  !
  character(17) :: property
  real(8) :: matrix_a(3,3),al(3,3),rshift(3),hvol,hgs(3)
  integer :: ia_s,ia_e
  logical :: flag_cuboid
  real(8) :: Rion_min(3), Rion_max(3), rps_max
  integer :: nc(3), ixyz, n

  call timer_begin(LOG_INIT_PS_TOTAL)

  if(allocated(ppg%save_udVtbl_a)) then
    property='update'
  else
    property='initial'
    if(comm_is_root(nproc_id_global))then
      if (.not. quiet) then
      write(*,*) ''
      write(*,*) '============init_ps=============='
      end if
    endif
    allocate(ppg%mps(natom))
    if (.not. allocated(ppg%jxyz_max)) then
      allocate(ppg%jxyz_max(1:3,natom))
      allocate(ppg%jxyz_min(1:3,natom))
      allocate(ppg%jxyz_changed(natom))
      ppg%jxyz_max = 0
      ppg%jxyz_min = ppg%nps
      ppg%jxyz_changed = .false.
    end if
    n=maxval(pp%nproj)*(pp%lmax+1)**2
    allocate(ppg%lma_tbl(n,natom))
    allocate(ppg%ia_tbl(n*natom))
    allocate(ppg%rinv_uvu(n*natom))
  endif
  
  ia_s = info%ia_s
  ia_e = info%ia_e
  matrix_a = system%rmatrix_A
  al = system%primitive_a
  hvol = system%hvol
  hgs = system%Hgs

  flag_cuboid = .true.
  if( abs(al(1,2)).ge.1d-10 .or. abs(al(1,3)).ge.1d-10.or. &
         abs(al(2,3)).ge.1d-10 )  flag_cuboid=.false.

  if( flag_cuboid ) then
     rps_max   = maxval( pp%rps(:)) + 1.5d0*maxval(hgs) + 1d-2
  endif

  if(iperiodic==0)then
    nc(:)=0
  else if(iperiodic==3)then
    if( flag_cuboid ) then
       do ixyz=1,3
          Rion_min(ixyz) = minval(system%Rion(ixyz,:))
          Rion_max(ixyz) = maxval(system%Rion(ixyz,:))
          if( Rion_min(ixyz) + 2d0*al(ixyz,ixyz) .gt. al(ixyz,ixyz) + rps_max .and. &
              Rion_max(ixyz) - 2d0*al(ixyz,ixyz) .lt.               - rps_max ) then
             nc(ixyz) = 1
          else
             nc(ixyz) = 2
          endif
       enddo
    else
       nc(:)=2
    endif
  end if

  if(iperiodic==0)then
    if(mod(lg%num(1),2)==1)then
      rshift(1)=0.d0
    else
      rshift(1)=-0.5d0*hgs(1)
    end if
    if(mod(lg%num(2),2)==1)then
      rshift(2)=0.d0
    else
      rshift(2)=-0.5d0*hgs(2)
    end if
    if(mod(lg%num(3),2)==1)then
      rshift(3)=0.d0
    else
      rshift(3)=-0.5d0*hgs(3)
    end if
  else if(iperiodic==3)then
    rshift(1)=-hgs(1)
    rshift(2)=-hgs(2)
    rshift(3)=-hgs(3)
  end if

  if (.not. allocated(ppg%Rion_old)) then
    call cache_jxyz(ppg,system%Rion)
  end if

  call timer_begin(LOG_INIT_PS_CALC_NPS)
  call calc_nps
  call timer_end(LOG_INIT_PS_CALC_NPS)

  call timer_begin(LOG_INIT_PS_CALC_JXYZ)
  call calc_jxyz
  call timer_end(LOG_INIT_PS_CALC_JXYZ)

  call cache_jxyz(ppg,system%Rion)

  call timer_begin(LOG_INIT_PS_LMA_UV)
  call set_lma
  call calc_uv
  if ( yn_spinorbit=='y' ) then
    call calc_uv_so(pp,ppg,lg%num,hgs,hvol,property)
  end if
  if ( PLUS_U_ON ) then
    call calc_uv_plusU( pp, ppg, property )
  end if
  call timer_end(LOG_INIT_PS_LMA_UV)

  call timer_begin(LOG_INIT_PS_CALC_VPSL)
  select case(iperiodic)
  case(0)
    call calc_Vpsl_isolated(lg,mg,system,info,pp,fg,Vpsl,ppg,property)
  case(3)
    call calc_Vpsl_periodic(lg,mg,system,info,pp,fg,poisson,Vpsl,ppg,property)
  end select
  call timer_end(LOG_INIT_PS_CALC_VPSL)

  call timer_begin(LOG_INIT_PS_UVPSI)
  call init_uvpsi_summation(ppg,info%icomm_r)
  call init_uvpsi_table(ppg)
  call timer_end(LOG_INIT_PS_UVPSI)

  if(iperiodic==3) then
    call update_kvector_nonlocalpt(info%ik_s,info%ik_e,system,ppg)
  end if
  
  if(comm_is_root(nproc_id_global) .and. property=='initial' .and. (.not. quiet)) write(*,*)'end init_ps'

  call timer_end(LOG_INIT_PS_TOTAL)
  return
  
contains

!-----------------------------------------------------------------------------------------------------------------------------------

  subroutine calc_nps
    use salmon_global,only : kion
    use communication, only: comm_get_max,comm_get_groupinfo,comm_logical_or
    implicit none
    !
    integer :: ia,ik,i1,i2,i3,j1,j2,j3,j,ixyz
    integer :: mps_tmp
    real(8) :: tmpx,tmpy,tmpz
    real(8) :: x,y,z,r,u,v,w
    real(8) :: xyz(3)

    mps_tmp = 0
    ppg%jxyz_changed(:) = .false.
!$omp parallel do default(none) &
!$omp private(ia,ik,j,i1,i2,i3,j1,j2,j3,tmpx,tmpy,tmpz,x,y,z,r,u,v,w,xyz) &
!$omp shared(ia_s,ia_e,kion,nc,flag_cuboid,system,al,hgs,rps_max,mg,rshift,matrix_a,pp,ppg) &
!$omp reduction(max:mps_tmp)
    do ia=ia_s,ia_e
      ik=kion(ia)
      j=0
      do j1=-nc(1),nc(1)
        if( flag_cuboid ) then
          xyz(1) = system%Rion(1,ia) + j1*al(1,1)
          if( xyz(1) .le. mg%is(1)* hgs(1) - rps_max  .or. &
              xyz(1) .ge. mg%ie(1)* hgs(1) + rps_max ) cycle
        endif
      do j2=-nc(2),nc(2)
        if( flag_cuboid ) then
          xyz(2) = system%Rion(2,ia) + j2*al(2,2)
          if( xyz(2) .le. mg%is(2)* hgs(2) - rps_max  .or. &
              xyz(2) .ge. mg%ie(2)* hgs(2) + rps_max ) cycle
        endif
      do j3=-nc(3),nc(3)
        if( flag_cuboid ) then
          xyz(3) = system%Rion(3,ia) + j3*al(3,3)
          if( xyz(3) .le. mg%is(3)* hgs(3) - rps_max  .or. &
              xyz(3) .ge. mg%ie(3)* hgs(3) + rps_max ) cycle
        endif

        tmpx = system%Rion(1,ia) + j1*al(1,1) + j2*al(1,2) + j3*al(1,3)
        tmpy = system%Rion(2,ia) + j1*al(2,1) + j2*al(2,2) + j3*al(2,3)
        tmpz = system%Rion(3,ia) + j1*al(3,1) + j2*al(3,2) + j3*al(3,3)
        do i3=mg%is(3),mg%ie(3)
        do i2=mg%is(2),mg%ie(2)
        do i1=mg%is(1),mg%ie(1)
          u = i1*hgs(1) + rshift(1)
          v = i2*hgs(2) + rshift(2)
          w = i3*hgs(3) + rshift(3)
          x = u*matrix_a(1,1) + v*matrix_a(1,2) + w*matrix_a(1,3) - tmpx
          y = u*matrix_a(2,1) + v*matrix_a(2,2) + w*matrix_a(2,3) - tmpy
          z = u*matrix_a(3,1) + v*matrix_a(3,2) + w*matrix_a(3,3) - tmpz
          r = sqrt(x*x+y*y+z*z)
          if (r<pp%rps(ik)+1.d-12) then
            j=j+1
            if (ppg%mps_old(ia) < j) then
              ppg%jxyz_changed(ia) = .true.
            else
              ppg%jxyz_changed(ia) = ppg%jxyz_changed(ia)          .or. &
                                     i1 /= ppg%jxyz_old(1,j,ia) .or. &
                                     i2 /= ppg%jxyz_old(2,j,ia) .or. &
                                     i3 /= ppg%jxyz_old(3,j,ia)
            end if
          end if
        end do
        end do
        end do
      end do
      end do
      end do
      ppg%jxyz_changed(ia) = ppg%jxyz_changed(ia) .or. (ppg%mps_old(ia) /= j)
      mps_tmp = max(mps_tmp,j)
    end do
!$omp end parallel do

    ppg%nps=mps_tmp
    if (allocated(ppg%jxyz_old)) then
      ppg%nps = max(ppg%nps, size(ppg%jxyz_old,2))
    end if
    call comm_get_max(ppg%nps,info%icomm_ko)
    call comm_logical_or(ppg%jxyz_changed,info%icomm_ko)

  end subroutine calc_nps

!-----------------------------------------------------------------------------------------------------------------------------------

  subroutine calc_jxyz
    use salmon_global,only : kion
    use communication,only: comm_get_groupinfo,comm_summation
    implicit none
    !
    integer :: ia,i,ik,i1,i2,i3,j1,j2,j3,j
    integer :: ixyz
    real(8) :: tmpx,tmpy,tmpz
    real(8) :: r,x,y,z,u,v,w
    real(8) :: xyz(3)
    
    allocate(ppg%jxyz(3,ppg%nps,natom))
    allocate(ppg%rxyz(3,ppg%nps,natom))

    ppg%jxyz = 0
    ppg%rxyz = 0d0
    ppg%mps  = 0

!$omp parallel do default(none) &
!$omp private(ia,ik,j,i,i1,i2,i3,j1,j2,j3,tmpx,tmpy,tmpz,x,y,z,r,u,v,w,xyz) &
!$omp shared(ia_s,ia_e,natom,kion,nc,al,system,hgs,rshift,matrix_a,pp,ppg,mg,flag_cuboid,rps_max)
    do ia=ia_s,ia_e
      if (ppg%jxyz_changed(ia)) then
        ik=kion(ia)
        j=0
        do j1=-nc(1),nc(1)
          if( flag_cuboid ) then
            xyz(1) = system%Rion(1,ia) + j1*al(1,1)
            if( xyz(1) .le. mg%is(1)* hgs(1) - rps_max  .or. &
                xyz(1) .ge. mg%ie(1)* hgs(1) + rps_max ) cycle
          endif
        do j2=-nc(2),nc(2)
          if( flag_cuboid ) then
            xyz(2) = system%Rion(2,ia) + j2*al(2,2)
            if( xyz(2) .le. mg%is(2)* hgs(2) - rps_max  .or. &
                xyz(2) .ge. mg%ie(2)* hgs(2) + rps_max ) cycle
          endif
        do j3=-nc(3),nc(3)
          if( flag_cuboid ) then
            xyz(3) = system%Rion(3,ia) + j3*al(3,3)
            if( xyz(3) .le. mg%is(3)* hgs(3) - rps_max  .or. &
                xyz(3) .ge. mg%ie(3)* hgs(3) + rps_max ) cycle
          endif

          tmpx = system%Rion(1,ia) + j1*al(1,1) + j2*al(1,2) + j3*al(1,3)
          tmpy = system%Rion(2,ia) + j1*al(2,1) + j2*al(2,2) + j3*al(2,3)
          tmpz = system%Rion(3,ia) + j1*al(3,1) + j2*al(3,2) + j3*al(3,3)
          do i3=mg%is(3),mg%ie(3)
          do i2=mg%is(2),mg%ie(2)
          do i1=mg%is(1),mg%ie(1)
            u = i1*hgs(1) + rshift(1)
            v = i2*hgs(2) + rshift(2)
            w = i3*hgs(3) + rshift(3)
            x = u*matrix_a(1,1) + v*matrix_a(1,2) + w*matrix_a(1,3) - tmpx
            y = u*matrix_a(2,1) + v*matrix_a(2,2) + w*matrix_a(2,3) - tmpy
            z = u*matrix_a(3,1) + v*matrix_a(3,2) + w*matrix_a(3,3) - tmpz
            r = sqrt(x*x+y*y+z*z)
            if (r<pp%rps(ik)+1.d-12) then
              j = j + 1
              if (j<=ppg%nps) then
                ppg%jxyz(1,j,ia) = i1
                ppg%jxyz(2,j,ia) = i2
                ppg%jxyz(3,j,ia) = i3
                ppg%rxyz(1,j,ia) = x
                ppg%rxyz(2,j,ia) = y
                ppg%rxyz(3,j,ia) = z
              end if
            end if
          end do
          end do
          end do
        end do
        end do
        end do
        ppg%mps(ia) = j
      else
        i = ppg%mps_old(ia)
        ppg%mps(ia) = i
        ppg%jxyz(1:3,1:i,ia) = ppg%jxyz_old(1:3,1:i,ia)
        do j=1,i
          ppg%rxyz(1:3,j,ia) = ppg%rxyz_old(1:3,j,ia) - (system%Rion(1:3,ia) - ppg%rion_old(1:3,ia))
        end do
      end if
    end do
!$omp end parallel do

    call comm_summation(ppg%jxyz,info%icomm_ko)
    call comm_summation(ppg%rxyz,info%icomm_ko)
    call comm_summation(ppg%mps, info%icomm_ko)

  end subroutine calc_jxyz

!-----------------------------------------------------------------------------------------------------------------------------------

  subroutine set_lma
    use salmon_global,only : kion
    implicit none
    integer :: lma,lm,ia,ik,m,l,ll,l0

    lma=0
    do ia=1,natom
      ik=kion(ia)
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        if(pp%inorm(l,ik)==0) cycle
        do m=-ll,ll
          lma=lma+1
        enddo
      enddo
      l0=l
      enddo
    enddo
    ppg%nlma=lma
    
    ppg%lma_tbl=0
    ppg%ia_tbl=0
    ppg%rinv_uvu=0.0d0
    
    allocate(ppg%uv(ppg%nps,ppg%nlma),ppg%duv(ppg%nps,ppg%nlma,3))

    lma=0
    do ia=1,natom
      ik=kion(ia)
      lm=0
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        if(pp%inorm(l,ik)==0) cycle
        do m=-ll,ll
          lm=lm+1
          lma=lma+1
          ppg%lma_tbl(lm,ia)=lma
          ppg%ia_tbl(lma)=ia
        enddo
      enddo
      l0=l
      enddo
    enddo
    return
  end subroutine set_lma
  
!-----------------------------------------------------------------------------------------------------------------------------------

  subroutine calc_uv
    use salmon_global,only : kion,nelem
    use math_constants,only : pi
    use salmon_math,only : ylm,dylm,spline
    implicit none
    integer :: ia,ik,j,l,lm,m,ll,l0
    integer :: ilma,intr,ir,lma
    real(8),allocatable :: xn(:),yn(:),an(:),bn(:),cn(:),dn(:)
    real(8) :: uvr(0:2*pp%lmax+1), r,x,y,z, xx

    if( property == 'initial' ) then
      allocate( ppg%save_udVtbl_a(pp%nrmax,0:2*pp%lmax+1,nelem) )
      allocate( ppg%save_udVtbl_b(pp%nrmax,0:2*pp%lmax+1,nelem) )
      allocate( ppg%save_udVtbl_c(pp%nrmax,0:2*pp%lmax+1,nelem) )
      allocate( ppg%save_udVtbl_d(pp%nrmax,0:2*pp%lmax+1,nelem) )

      do ik=1,nelem
        allocate(xn(0:pp%nrps(ik)-1),yn(0:pp%nrps(ik)-1),an(0:pp%nrps(ik)-2) &
                ,bn(0:pp%nrps(ik)-2),cn(0:pp%nrps(ik)-2),dn(0:pp%nrps(ik)-2))
    
        xn(0:pp%nrps(ik)-1) = pp%radnl(1:pp%nrps(ik),ik)
        l0=0
        do ll=0,pp%mlps(ik)
        do l=l0,l0+pp%nproj(ll,ik)-1
          yn(0:pp%nrps(ik)-1) = pp%udvtbl(1:pp%nrps(ik),l,ik)
          call spline(pp%nrps(ik),xn,yn,an,bn,cn,dn)
          ppg%save_udvtbl_a(1:pp%nrps(ik)-1,l,ik) = an(0:pp%nrps(ik)-2)
          ppg%save_udvtbl_b(1:pp%nrps(ik)-1,l,ik) = bn(0:pp%nrps(ik)-2)
          ppg%save_udvtbl_c(1:pp%nrps(ik)-1,l,ik) = cn(0:pp%nrps(ik)-2)
          ppg%save_udvtbl_d(1:pp%nrps(ik)-1,l,ik) = dn(0:pp%nrps(ik)-2)
        end do
        l0=l
        end do
        deallocate(xn,yn,an,bn,cn,dn)
      enddo
    end if
    
    do ia=1,natom
       ik=kion(ia)

    !!$omp parallel
    !!$omp do private(j,x,y,z,r,ir,intr,xx,l,lm,m,uvr,ilma,l0,ll)
       do j=1,ppg%mps(ia)
         x=ppg%rxyz(1,j,ia)
         y=ppg%rxyz(2,j,ia)
         z=ppg%rxyz(3,j,ia)
         r=sqrt(x*x+y*y+z*z)+1d-50
         do ir=1,pp%nrps(ik)
           if(pp%radnl(ir,ik).gt.r) exit
         enddo
         intr=ir-1
         if(intr.lt.0.or.intr.ge.pp%nrps(ik)) stop 'bad intr at prep_ps'
         xx = r - pp%radnl(intr,ik)

         l0=0
         do ll=0,pp%mlps(ik)
         do l=l0,l0+pp%nproj(ll,ik)-1
            uvr(l)=   ppg%save_udvtbl_a(intr,l,ik)*xx**3 &
                    + ppg%save_udvtbl_b(intr,l,ik)*xx**2 &
                    + ppg%save_udvtbl_c(intr,l,ik)*xx    &
                    + ppg%save_udvtbl_d(intr,l,ik)
         enddo
         l0=l
         enddo
   
         lm=0
         l0=0
         do ll=0,pp%mlps(ik)
         do l=l0,l0+pp%nproj(ll,ik)-1
           if(pp%inorm(l,ik)==0) cycle
           do m=-ll,ll
             lm=lm+1
             ilma=ppg%lma_tbl(lm,ia)
             ppg%uv(j,ilma)   = uvr(l)* ylm(x,y,z,ll,m)
           enddo
         enddo
         l0=l
         enddo
   
       enddo
   !!$omp end do
   !!$omp end parallel
   
    enddo

    lma=0
    do ia=1,natom
      ik=kion(ia)
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        if(pp%inorm(l,ik)==0) cycle
        do m=-ll,ll
          lma=lma+1
          ppg%rinv_uvu(lma)=dble(pp%inorm(l,ik))*hvol
        enddo
      enddo
      l0=l
      enddo
    enddo

  end subroutine calc_uv

end subroutine init_ps

!===================================================================================================================================

SUBROUTINE dealloc_init_ps(ppg)
  use structures, only: s_pp_grid
  implicit none
  type(s_pp_grid) :: ppg

  deallocate(ppg%jxyz, ppg%rxyz, ppg%uv, ppg%duv)
  if(allocated(ppg%zekr_uV)) deallocate(ppg%zekr_uV)

  if (allocated(ppg%irange_atom))    deallocate(ppg%irange_atom)
  if (allocated(ppg%ireferred_atom)) deallocate(ppg%ireferred_atom)
  if (allocated(ppg%ilocal_nlma2ilma)) deallocate(ppg%ilocal_nlma2ilma)
  if (allocated(ppg%ilocal_nlma2ia))   deallocate(ppg%ilocal_nlma2ia)
  
  return
END SUBROUTINE dealloc_init_ps

!===================================================================================================================================

SUBROUTINE calc_Vpsl_isolated(lg,mg,system,info,pp,fg,vpsl,ppg,property)
  use structures
  use salmon_global,only : natom, kion, quiet, method_poisson, nelem, yn_ffte
#ifdef USE_FFTW
  use salmon_global,only : yn_fftw
#endif
  use math_constants,only : pi,zi
  use parallelization, only: nproc_id_global
  use communication, only: comm_summation
  implicit none
  type(s_rgrid)          ,intent(in) :: lg,mg
  type(s_dft_system)     ,intent(in) :: system
  type(s_parallel_info)  ,intent(in) :: info
  type(s_pp_info)        ,intent(in) :: pp
  type(s_reciprocal_grid),intent(in) :: fg
  type(s_scalar)                     :: vpsl
  type(s_pp_grid)                    :: ppg
  character(17)          ,intent(in) :: property
  !
  integer :: ix,iy,iz,ak,ik
  integer :: ia,i,j,a,intr
  real(8) :: ratio1,ratio2,r
  integer :: ifgx_s,ifgx_e
  integer :: ifgy_s,ifgy_e
  integer :: ifgz_s,ifgz_e
  real(8) :: g(3),gd,s,g2sq,r1,dr,vloc_av
  complex(8) :: tmp_exp
  complex(8),allocatable :: vtmp1(:,:,:,:),vtmp2(:,:,:,:)

  if(.not.allocated(ppg%Vpsl_ion)) then
    allocate(ppg%Vpsl_ion(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:natom))
  end if
  
  Vpsl%f=0.d0

  do a=1,natom
    ak=Kion(a)
    do j=1,3
      if(abs(system%Rion(j,a))<lg%num(j)*system%Hgs(j))then
        continue
      else
        if (.not. quiet) &
        write(*,*) "Rion error",nproc_id_global,a,j,system%Rion(j,a)
      end if
    end do
    do ix=mg%is(1),mg%ie(1)
    do iy=mg%is(2),mg%ie(2)
    do iz=mg%is(3),mg%ie(3)
      r=sqrt( (lg%coordinate(ix,1)-system%Rion(1,a))**2      &
             +(lg%coordinate(iy,2)-system%Rion(2,a))**2      &
             +(lg%coordinate(iz,3)-system%Rion(3,a))**2 )+1.d-50
      call bisection(r,intr,ak,pp%nrmax,pp%rad)
      ratio1=(r-pp%rad(intr,ak))/(pp%rad(intr+1,ak)-pp%rad(intr,ak)) ; ratio2=1.d0-ratio1
      if(intr>0.and.intr<=pp%nrmax)then
        continue
      else
        write(*,*) "intr error",nproc_id_global,intr,r
      end if

      Vpsl%f(ix,iy,iz)=Vpsl%f(ix,iy,iz)      &
                  +ratio1*pp%vpp_f(intr,pp%Lref(ak),ak)      &
                  +ratio2*pp%vpp_f(intr-1,pp%Lref(ak),ak)  !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate
      ppg%Vpsl_ion(ix,iy,iz,a) = ratio1*pp%vpp_f(intr,pp%Lref(ak),ak) + ratio2*pp%vpp_f(intr-1,pp%Lref(ak),ak)
    end do
    end do
    end do
  end do

  allocate(ppg%zekr_uV(ppg%nps,ppg%nlma,1))
  ppg%zekr_uV(:,:,1) = dcmplx(ppg%uV)

  if(method_poisson=='ft')then
#ifdef USE_FFTW
    if(yn_fftw=='n')then
#endif
      if(yn_ffte=='n')then
        ifgx_s = (mg%is(1)-lg%is(1))*2+1
        ifgx_e = (mg%is(1)-lg%is(1))*2+mg%num(1)*2
        ifgy_s = (mg%is(2)-lg%is(2))*2+1
        ifgy_e = (mg%is(2)-lg%is(2))*2+mg%num(2)*2
        ifgz_s = (mg%is(3)-lg%is(3))*2+1
        ifgz_e = (mg%is(3)-lg%is(3))*2+mg%num(3)*2
      else
        if(mod(info%nporbital,4)==0)then
          ! start and end point of reciprocal grids for x, y, z
          ifgx_s = 1
          ifgx_e = 2*lg%num(1)
          if(info%id_y_isolated_ffte >= info%isize_y_isolated_ffte/2) then
            ifgy_s = mg%is(2)-lg%is(2)+1+lg%num(2)
          else
            ifgy_s = mg%is(2)-lg%is(2)+1
          end if
          ifgy_e = ifgy_s+mg%num(2)-1
          if(info%id_z_isolated_ffte >= info%isize_z_isolated_ffte/2) then
            ifgz_s = mg%is(3)-lg%is(3)+1+lg%num(3)
          else
            ifgz_s = mg%is(3)-lg%is(3)+1
          end if
          ifgz_e = ifgz_s+mg%num(3)-1
        else
          ! start and end point of reciprocal grids for x, y, z
          ifgx_s = 1
          ifgx_e = 2*lg%num(1)
          ifgy_s = 1
          ifgy_e = 2*lg%num(2)
          ifgz_s = 1
          ifgz_e = 2*lg%num(3)
        end if
      end if
#ifdef USE_FFTW
    else if(yn_fftw=='y')then
      if(mod(info%nporbital,2)==0)then
        ! start and end point of reciprocal grids for x, y, z
        ifgx_s = 1
        ifgx_e = 2*lg%num(1)
        ifgy_s = 1
        ifgy_e = 2*lg%num(2)
        if(info%iaddress_isolated_fftw(4)==1) then
          ifgz_s = mg%is(3)-lg%is(3)+1+lg%num(3)
        else
          ifgz_s = mg%is(3)-lg%is(3)+1
        end if
        ifgz_e = ifgz_s+mg%num(3)-1
      else
        ! start and end point of reciprocal grids for x, y, z
        ifgx_s = 1
        ifgx_e = 2*lg%num(1)
        ifgy_s = 1
        ifgy_e = 2*lg%num(2)
        ifgz_s = 1
        ifgz_e = 2*lg%num(3)
      end if
    endif
#endif

    if( property == 'initial' ) then
      allocate(ppg%zrhoG_ion(ifgx_s:ifgx_e,ifgy_s:ifgy_e,ifgz_s:ifgz_e)  & ! rho_ion(G)
            & ,ppg%zVG_ion  (ifgx_s:ifgx_e,ifgy_s:ifgy_e,ifgz_s:ifgz_e,nelem)) ! V_ion(G)

      ppg%zVG_ion = 0d0


  !$omp parallel
  !$omp do private(ik,ix,iy,iz,g,g2sq,s,r1,dr,i,vloc_av) collapse(3)
      do ik=1,nelem
        do iz=ifgz_s,ifgz_e
        do iy=ifgy_s,ifgy_e
        do ix=ifgx_s,ifgx_e
          g(1) = fg%vec_G(1,ix,iy,iz)
          g(2) = fg%vec_G(2,ix,iy,iz)
          g(3) = fg%vec_G(3,ix,iy,iz)
          g2sq = sqrt(g(1)**2+g(2)**2+g(3)**2)
          s=0.d0
          if (fg%if_Gzero(ix,iy,iz)) then
            do i=2,pp%nrloc(ik)
              r1=0.5d0*(pp%rad(i,ik)+pp%rad(i-1,ik))
              dr=pp%rad(i,ik)-pp%rad(i-1,ik)
              vloc_av = 0.5d0*(pp%vloctbl(i,ik)+pp%vloctbl(i-1,ik))
              s=s+4d0*pi*(r1**2*vloc_av+r1*pp%zps(ik))*dr
            end do
          else
            do i=2,pp%nrloc(ik)
              r1=0.5d0*(pp%rad(i,ik)+pp%rad(i-1,ik))
              dr=pp%rad(i,ik)-pp%rad(i-1,ik)
              vloc_av = 0.5d0*(pp%vloctbl(i,ik)+pp%vloctbl(i-1,ik))
              s=s+4d0*pi*sin(g2sq*r1)/g2sq*(r1*vloc_av+pp%zps(ik))*dr !Vloc - coulomb
            end do
          end if
          ppg%zVG_ion(ix,iy,iz,ik) = s
        end do
        end do
        end do
      end do
  !$omp end do
  !$omp end parallel

    end if

#ifdef USE_FFTW
    if(yn_fftw=='n')then
#endif
      if(yn_ffte=='y')then
        allocate(vtmp1(ifgx_s:ifgx_e,ifgy_s:ifgy_e,ifgz_s:ifgz_e,1:2))
        allocate(vtmp2(ifgx_s:ifgx_e,ifgy_s:ifgy_e,ifgz_s:ifgz_e,1:2))
      else
        allocate(vtmp1((mg%is(1)-lg%is(1))*2+1:(mg%is(1)-lg%is(1))*2+mg%num(1)*2, &
                       (mg%is(2)-lg%is(2))*2+1:(mg%is(2)-lg%is(2))*2+mg%num(2)*2, &
                       (mg%is(3)-lg%is(3))*2+1:(mg%is(3)-lg%is(3))*2+mg%num(3)*2 ,1:2))
        allocate(vtmp2((mg%is(1)-lg%is(1))*2+1:(mg%is(1)-lg%is(1))*2+mg%num(1)*2, &
                       (mg%is(2)-lg%is(2))*2+1:(mg%is(2)-lg%is(2))*2+mg%num(2)*2, &
                       (mg%is(3)-lg%is(3))*2+1:(mg%is(3)-lg%is(3))*2+mg%num(3)*2 ,1:2))
      end if
#ifdef USE_FFTW
    else if(yn_fftw=='y')then
      allocate(vtmp1(ifgx_s:ifgx_e,ifgy_s:ifgy_e,ifgz_s:ifgz_e,1:2))
      allocate(vtmp2(ifgx_s:ifgx_e,ifgy_s:ifgy_e,ifgz_s:ifgz_e,1:2))
    end if
#endif

! vtmp(:,:,:,1)=V_ion(G): local part of the pseudopotential in the G space
    vtmp1 = 0d0
  !$omp parallel do collapse(2) private(ix,iy,iz,g,ia,ik,gd,tmp_exp)
    do iz=ifgz_s,ifgz_e
    do iy=ifgy_s,ifgy_e
    do ix=ifgx_s,ifgx_e
      g(1) = fg%vec_G(1,ix,iy,iz)
      g(2) = fg%vec_G(2,ix,iy,iz)
      g(3) = fg%vec_G(3,ix,iy,iz)
      do ia=info%ia_s,info%ia_e
        ik=kion(ia)
        gd = g(1)*system%Rion(1,ia) + g(2)*system%Rion(2,ia) + g(3)*system%Rion(3,ia)
        tmp_exp = exp(-zi*gd)/system%det_A
        vtmp1(ix,iy,iz,1) = vtmp1(ix,iy,iz,1) + ( ppg%zVG_ion(ix,iy,iz,ik) - fg%coef(ix,iy,iz)*pp%zps(ik) ) *tmp_exp ! V_ion(G)
        vtmp1(ix,iy,iz,2) = vtmp1(ix,iy,iz,2) + pp%zps(ik)*tmp_exp ! rho_ion(G)
      end do
      end do
      end do
    end do
  !$omp end parallel do
 
#ifdef USE_FFTW
    if(yn_fftw=='n')then
#endif
      if(yn_ffte=='y')then
        ppg%zrhoG_ion = vtmp1(:,:,:,2)
      else
        call comm_summation(vtmp1,vtmp2,(ifgx_e-ifgx_s+1)*(ifgy_e-ifgy_s+1)*(ifgz_e-ifgz_s+1)*2,info%icomm_ko)
        ppg%zrhoG_ion = vtmp2(:,:,:,2)
      end if
#ifdef USE_FFTW
    else if(yn_fftw=='y')then
      ppg%zrhoG_ion = vtmp1(:,:,:,2)
    end if
#endif

    deallocate(vtmp1,vtmp2)

  end if

  return
END SUBROUTINE calc_Vpsl_isolated

!===================================================================================================================================

subroutine calc_vpsl_periodic(lg,mg,system,info,pp,fg,poisson,Vpsl,ppg,property)
  use salmon_global,only : nelem, kion, yn_ffte
  use communication, only: comm_summation
  use math_constants,only : pi,zi
  use structures
  implicit none
  type(s_rgrid)          ,intent(in) :: lg,mg
  type(s_dft_system)     ,intent(in) :: system
  type(s_parallel_info)  ,intent(in) :: info
  type(s_pp_info)        ,intent(in) :: pp
  type(s_reciprocal_grid),intent(in) :: fg
  type(s_poisson)                    :: poisson
  type(s_scalar)                     :: Vpsl
  type(s_pp_grid)                    :: ppg
  character(17)          ,intent(in) :: property
  !
  integer :: ia,i,ik,ix,iy,iz,kx,ky,kz,iiy,iiz
  real(8) :: g(3),gd,s,g2sq,r1,dr,vloc_av
  complex(8) :: tmp_exp
  complex(8) :: vtmp1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:2)
  complex(8) :: vtmp2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:2)

  if( property == 'initial' ) then
  
    allocate(ppg%zrhoG_ion(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)) & ! rho_ion(G)
          & ,ppg%zVG_ion  (mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),nelem)) ! V_ion(G)

    ppg%zVG_ion = 0d0
  !$omp parallel
  !$omp do private(ik,ix,iy,iz,g,g2sq,s,r1,dr,i,vloc_av) collapse(3)
    do ik=1,nelem
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        g(1) = fg%vec_G(1,ix,iy,iz)
        g(2) = fg%vec_G(2,ix,iy,iz)
        g(3) = fg%vec_G(3,ix,iy,iz)
        g2sq = sqrt(g(1)**2+g(2)**2+g(3)**2)
        s=0.d0
        if (fg%if_Gzero(ix,iy,iz)) then
          do i=2,pp%nrloc(ik)
            r1=0.5d0*(pp%rad(i,ik)+pp%rad(i-1,ik))
            dr=pp%rad(i,ik)-pp%rad(i-1,ik)
            vloc_av = 0.5d0*(pp%vloctbl(i,ik)+pp%vloctbl(i-1,ik))
            s=s+4d0*pi*(r1**2*vloc_av+r1*pp%zps(ik))*dr
          end do
        else
          do i=2,pp%nrloc(ik)
            r1=0.5d0*(pp%rad(i,ik)+pp%rad(i-1,ik))
            dr=pp%rad(i,ik)-pp%rad(i-1,ik)
            vloc_av = 0.5d0*(pp%vloctbl(i,ik)+pp%vloctbl(i-1,ik))
            s=s+4d0*pi*sin(g2sq*r1)/g2sq*(r1*vloc_av+pp%zps(ik))*dr !Vloc - coulomb
          end do
        end if
        ppg%zVG_ion(ix,iy,iz,ik) = s
      end do
      end do
      end do
    end do
  !$omp end do
  !$omp end parallel

  endif

! vtmp(:,:,:,1)=V_ion(G): local part of the pseudopotential in the G space
  vtmp1 = 0d0
  !$omp parallel do collapse(2) private(ix,iy,iz,g,ia,ik,gd,tmp_exp)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    g(1) = fg%vec_G(1,ix,iy,iz)
    g(2) = fg%vec_G(2,ix,iy,iz)
    g(3) = fg%vec_G(3,ix,iy,iz)
    do ia=info%ia_s,info%ia_e
      ik=kion(ia)
      gd = g(1)*system%Rion(1,ia) + g(2)*system%Rion(2,ia) + g(3)*system%Rion(3,ia)
      tmp_exp = exp(-zi*gd)/system%det_A
      vtmp1(ix,iy,iz,1) = vtmp1(ix,iy,iz,1) + ( ppg%zVG_ion(ix,iy,iz,ik) - fg%coef(ix,iy,iz)*pp%zps(ik) ) *tmp_exp ! V_ion(G)
      vtmp1(ix,iy,iz,2) = vtmp1(ix,iy,iz,2) + pp%zps(ik)*tmp_exp ! rho_ion(G)
    end do
    end do
    end do
  end do
  !$omp end parallel do
  
  call comm_summation(vtmp1,vtmp2,mg%num(1)*mg%num(2)*mg%num(3)*2,info%icomm_ko)
  ppg%zrhoG_ion = vtmp2(:,:,:,2)

! Vpsl=V_ion(r): local part of the pseudopotential in the r space

  if(yn_ffte=='n') then
  ! cf. poisson_periodic.f90

  !$omp workshare
    poisson%ff1x = 0d0
  !$omp end workshare
  
  !$omp workshare
    poisson%ff1y = 0d0
  !$omp end workshare
  
  !$omp workshare
    poisson%ff1z = 0d0
  !$omp end workshare

  !$OMP parallel do private(kz,ky,kx)
    do kz = mg%is(3),mg%ie(3)
    do ky = mg%is(2),mg%ie(2)
    do kx = mg%is(1),mg%ie(1)
      poisson%ff1z(kx,ky,kz) = vtmp2(kx,ky,kz,1) ! V_ion(G)
    end do
    end do
    end do
    call comm_summation(poisson%ff1z,poisson%ff2z,mg%num(1)*mg%num(2)*lg%num(3),info%icomm_z)

  !$OMP parallel do private(iz,ky,kx)
    do iz = mg%is(3),mg%ie(3)
    do ky = mg%is(2),mg%ie(2)
    do kx = mg%is(1),mg%ie(1)
      poisson%ff1y(kx,ky,iz)=sum(fg%egz(:,iz)*poisson%ff2z(kx,ky,:))
    end do
    end do
    end do
    call comm_summation(poisson%ff1y,poisson%ff2y,mg%num(1)*lg%num(2)*mg%num(3),info%icomm_y)

  !$OMP parallel do private(iz,iy,kx)
    do iz = mg%is(3),mg%ie(3)
    do iy = mg%is(2),mg%ie(2)
    do kx = mg%is(1),mg%ie(1)
      poisson%ff1x(kx,iy,iz)=sum(fg%egy(:,iy)*poisson%ff2y(kx,:,iz))
    end do
    end do
    end do
    call comm_summation(poisson%ff1x,poisson%ff2x,lg%num(1)*mg%num(2)*mg%num(3),info%icomm_x)

  !$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz = mg%is(3),mg%ie(3)
    do iy = mg%is(2),mg%ie(2)
    do ix = mg%is(1),mg%ie(1)
      Vpsl%f(ix,iy,iz) = sum(fg%egx(:,ix)*poisson%ff2x(:,iy,iz))
    end do
    end do
    end do
  
  else
  ! cf. poisson_ffte.f90
  
    poisson%b_ffte=0.d0
  !$OMP parallel do private(iz,iy,ix,iiz,iiy) collapse(2)
    do iz = 1,mg%num(3)
    do iy = 1,mg%num(2)
    do ix = mg%is(1),mg%ie(1)
      iiz=iz+mg%is(3)-1
      iiy=iy+mg%is(2)-1
      poisson%b_ffte(ix,iy,iz) = vtmp2(ix,iiy,iiz,1) ! V_ion(G)
    end do
    end do
    end do
    call comm_summation(poisson%b_ffte,poisson%a_ffte,size(poisson%a_ffte),info%icomm_x)

    CALL PZFFT3DV_MOD(poisson%a_ffte,poisson%b_ffte,lg%num(1),lg%num(2),lg%num(3),   &
                      info%isize_y,info%isize_z,1, &
                      info%icomm_y,info%icomm_z)

  !$OMP parallel do private(iz,iy,iiz,iiy) collapse(2)
    do iz=1,mg%num(3)
    do iy=1,mg%num(2)
      iiz=iz+mg%is(3)-1
      iiy=iy+mg%is(2)-1
      Vpsl%f(mg%is(1):mg%ie(1),iiy,iiz) = poisson%b_ffte(mg%is(1):mg%ie(1),iy,iz)*system%ngrid
    end do
    end do
    
  end if

  return
end subroutine calc_vpsl_periodic

!===================================================================================================================================

subroutine cache_jxyz(ppg,sysRion)
  use salmon_global, only: natom
  use structures
  implicit none
  type(s_pp_grid)    :: ppg
  real(8),intent(in) :: sysRion(:,:)

  if (allocated(ppg%mps_old))  deallocate(ppg%mps_old)
  if (allocated(ppg%rion_old)) deallocate(ppg%rion_old)
  if (allocated(ppg%jxyz_old)) deallocate(ppg%jxyz_old)
  if (allocated(ppg%rxyz_old)) deallocate(ppg%rxyz_old)

  allocate(ppg%mps_old(natom))
  allocate(ppg%rion_old(size(sysRion,1),size(sysRion,2)))
  ppg%rion_old = sysRion

  if (allocated(ppg%jxyz)) then
    ppg%mps_old  = ppg%mps

    allocate(ppg%jxyz_old(size(ppg%jxyz,1),size(ppg%jxyz,2),size(ppg%jxyz,3)))
    allocate(ppg%rxyz_old(size(ppg%rxyz,1),size(ppg%rxyz,2),size(ppg%rxyz,3)))
    ppg%jxyz_old = ppg%jxyz
    ppg%rxyz_old = ppg%rxyz
  else
    ppg%mps_old  = -1

    allocate(ppg%jxyz_old(3,1,natom))
    allocate(ppg%rxyz_old(3,1,natom))
    ppg%jxyz_old = -1
    ppg%rxyz_old = 0d0
  end if
end subroutine

!===================================================================================================================================

subroutine bisection(xx,inode,iak,nr,rad_psl)
  use salmon_global,only : nelem
  implicit none
  integer,intent(out) :: inode
  integer,intent(in)  :: iak
  integer,intent(in)  :: nr
  real(8),intent(in)  :: rad_psl(nr,nelem)
  real(8),intent(in)  :: xx
  integer :: imin,imax
  
  imin=1
  imax=nr
  do while (imax-imin>1)
    inode=(imin+imax)/2
    if(xx>rad_psl(inode,iak))then
      imin=inode
    else
      imax=inode
    end if
  end do
  inode=imin

end subroutine bisection

!===================================================================================================================================

subroutine init_uvpsi_summation(ppg,icomm_r)
  use structures,    only: s_pp_grid
  use salmon_global, only: natom
  use communication, only: comm_get_groupinfo &
                          ,comm_allgather &
                          ,comm_create_group_byid &
                          ,comm_group_null &
                          ,comm_free_group
  implicit none
  type(s_pp_grid),intent(inout) :: ppg
  integer,intent(in) :: icomm_r

  integer :: ilma,ia
  integer :: irank_r,isize_r,n,i
  integer :: nlma
  logical :: t,u
  logical,allocatable :: ireferred_atom_comm_r(:,:),iupdated(:)
  integer,allocatable :: iranklist(:)

  call comm_get_groupinfo(icomm_r, irank_r, isize_r)

  nlma = ppg%Nlma

  allocate(iranklist(isize_r))
  allocate(ireferred_atom_comm_r(natom,isize_r))
  allocate(iupdated(natom))
  allocate(ppg%irange_atom(2,natom))
  allocate(ppg%ireferred_atom(natom))
  if (.not. allocated(ppg%ireferred_atom_comm_r)) then
    allocate(ppg%ireferred_atom_comm_r(natom,isize_r))
    ppg%ireferred_atom_comm_r = .false.
  end if
  if (.not. allocated(ppg%icomm_atom)) then
    allocate(ppg%icomm_atom(natom))
    ppg%icomm_atom(:) = comm_group_null
  end if

  ppg%irange_atom(1,:) = 1
  ppg%irange_atom(2,:) = 0
#ifdef USE_OPENACC
!$acc parallel loop private(ia,ilma)
#else
!$omp parallel do private(ia,ilma)
#endif
  do ia=1,natom
    ppg%ireferred_atom(ia) = (ppg%mps(ia) > 0)

    ! forward search
    do ilma=1,nlma
      if (ppg%ia_tbl(ilma) == ia) then
        ppg%irange_atom(1,ia) = ilma
        exit
      end if
    end do

    ! backward search
    do ilma=nlma,1,-1
      if (ppg%ia_tbl(ilma) == ia) then
        ppg%irange_atom(2,ia) = ilma
        exit
      end if
    end do
  end do
#ifdef USE_OPENACC
!$acc end parallel
#else
!$omp end parallel do
#endif

  call comm_allgather(ppg%ireferred_atom, ireferred_atom_comm_r, icomm_r)

  iupdated = .false.
#ifdef USE_OPENACC
!$acc parallel loop private(ia,i,t,u)
#else
!$omp parallel do private(ia,i,t,u)
#endif
  do ia=1,natom
    do i=1,isize_r
      t = ppg%ireferred_atom_comm_r(ia,i)
      u = ireferred_atom_comm_r(ia,i)
      if (t .neqv. u) then
        iupdated(ia) = .true.
        exit
      end if
    end do
  end do
#ifdef USE_OPENACC
!$acc end parallel
#else
!$omp end parallel do
#endif

  ppg%ireferred_atom_comm_r = ireferred_atom_comm_r

  do ia=1,natom
    if (iupdated(ia)) then
      call comm_free_group(ppg%icomm_atom(ia))
      n = 0
      do i=1,isize_r
        if (ppg%ireferred_atom_comm_r(ia,i)) then
          n = n + 1
          iranklist(n) = i - 1
        end if
      end do
      ppg%icomm_atom(ia) = comm_create_group_byid(icomm_r, iranklist(1:n))
    end if
  end do
end subroutine init_uvpsi_summation

!===================================================================================================================================

subroutine init_uvpsi_table(ppg)
  use structures,    only: s_pp_grid
  implicit none
  type(s_pp_grid),intent(inout) :: ppg
  integer :: ilma,ia,ilocal,ilocal_nlma

  ilocal_nlma = 0
#ifdef USE_OPENACC
!$acc parallel loop private(ilma,ia) reduction(+:ilocal_nlma)
#else
!$omp parallel do private(ilma,ia) reduction(+:ilocal_nlma)
#endif
  do ilma=1,ppg%nlma
    ia = ppg%ia_tbl(ilma)
    if (ppg%ireferred_atom(ia)) ilocal_nlma = ilocal_nlma + 1
  end do
#ifdef USE_OPENACC
!$acc end parallel
#else
!$omp end parallel do
#endif
  ppg%ilocal_nlma = ilocal_nlma

  allocate(ppg%ilocal_nlma2ilma(ppg%ilocal_nlma))
  allocate(ppg%ilocal_nlma2ia  (ppg%ilocal_nlma))
  ilocal = 0
  do ilma=1,ppg%nlma
    ia = ppg%ia_tbl(ilma)
    if (ppg%ireferred_atom(ia)) then
      ilocal = ilocal + 1
      ppg%ilocal_nlma2ilma(ilocal) = ilma
      ppg%ilocal_nlma2ia  (ilocal) = ia
    end if
  end do
end subroutine init_uvpsi_table

end module prep_pp_sub
