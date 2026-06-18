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
  use nvtx_wrapper
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

  call nvtxStartRange('init_ps', __LINE__)
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

  if (property == 'initial') then
    call write_vloctbl_derivative_diagnostics(pp)
    call write_local_sr_u_diagnostics(pp)
  end if

  call timer_begin(LOG_INIT_PS_UVPSI)
  call init_uvpsi_summation(ppg,info%icomm_r)
  call init_uvpsi_table(ppg)
  call timer_end(LOG_INIT_PS_UVPSI)

#ifdef USE_OPENACC
  if(info%if_divide_rspace) then
  else
    allocate(ppg%uVpsibox(ppg%nlma, &
                          system%nspin, &
                          info%io_s:info%io_e, &
                          info%ik_s:info%ik_e, &
                          info%im_s:info%im_e))
    call init_uvpsi_blocking(ppg,mg)
  end if
#endif

  if(iperiodic==3) then
    call update_kvector_nonlocalpt(info%ik_s,info%ik_e,system,ppg)
  end if
  
  if(comm_is_root(nproc_id_global) .and. property=='initial' .and. (.not. quiet)) write(*,*)'end init_ps'

  call timer_end(LOG_INIT_PS_TOTAL)
  call nvtxEndRange
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
    
    allocate(ppg%uv(ppg%nps,ppg%nlma))

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
         if( intr <= 0 .or. pp%nrps(ik) <= intr ) stop &
         & 'Invalid r-grid sampling: The atomic positions should be shifted slightly by parallel translation.'
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

subroutine build_local_sr_shared_u(pp, ik)
  use salmon_global, only: method_loc_sr_origin, probe_loc_sr_fit_du12_scale, probe_loc_sr_fit_u1_scale
  use structures, only: s_pp_info
  implicit none
  type(s_pp_info), intent(inout) :: pp
  integer, intent(in) :: ik
  integer :: i, i1, i2, nr, switch_idx
  real(8) :: a1, a3, denom, dr, r1, r2, u1_fit, u1_raw, u2_fit, u2_raw, v1, v2

  nr = pp%nrps(ik)
  if (nr < 3) stop 'build_local_sr_shared_u: need at least 3 radial points'

  pp%u_sr_tbl(:,ik) = 0d0
  pp%du_sr_seg(:,ik) = 0d0
  pp%u_sr_origin_coef(:,:,ik) = 0d0
  pp%r_sr_origin_switch(ik) = 0d0

  do i = 1, nr
    pp%u_sr_tbl(i,ik) = pp%rad(i,ik) * pp%vloctbl(i,ik) + dble(pp%zps(ik))
  end do

  do i = 1, nr - 1
    dr = pp%rad(i+1,ik) - pp%rad(i,ik)
    if (dr > 0d0) then
      pp%du_sr_seg(i,ik) = (pp%u_sr_tbl(i+1,ik) - pp%u_sr_tbl(i,ik)) / dr
    else
      pp%du_sr_seg(i,ik) = 0d0
    end if
  end do
  if (nr > 1) pp%du_sr_seg(nr,ik) = pp%du_sr_seg(nr-1,ik)

  i1 = 0
  i2 = 0
  do i = 1, nr
    if (pp%rad(i,ik) <= 0d0) cycle
    if (i1 == 0) then
      i1 = i
    else
      i2 = i
      exit
    end if
  end do
  if (i1 == 0 .or. i2 == 0) stop 'build_local_sr_shared_u: unable to fit origin polynomial'

  r1 = pp%rad(i1,ik)
  r2 = pp%rad(i2,ik)
  u1_raw = pp%u_sr_tbl(i1,ik)
  u2_raw = pp%u_sr_tbl(i2,ik)
  u1_fit = probe_loc_sr_fit_u1_scale * u1_raw
  u2_fit = u1_raw + probe_loc_sr_fit_du12_scale * (u2_raw - u1_raw)
  select case(trim(method_loc_sr_origin))
  case('poly3')
    switch_idx = 3
    if (nr < switch_idx) stop 'build_local_sr_shared_u: invalid origin switch'
    pp%r_sr_origin_switch(ik) = pp%rad(switch_idx,ik)
    denom = r1 * r2 * (r2*r2 - r1*r1)
    if (abs(denom) <= tiny(1d0)) stop 'build_local_sr_shared_u: degenerate origin fit'
    a3 = (u2_fit * r1 - u1_fit * r2) / denom
    a1 = (u1_fit - a3 * r1**3) / r1
  case('vsr_linear')
    if (abs(r2 - r1) <= tiny(1d0)) stop 'build_local_sr_shared_u: degenerate vsr extrapolation'
    v1 = u1_fit / r1
    v2 = u2_fit / r2
    a1 = v1 - (v2 - v1) * r1 / (r2 - r1)
    a3 = 0d0
    pp%r_sr_origin_switch(ik) = r1
  case default
    stop 'build_local_sr_shared_u: unsupported method_loc_sr_origin'
  end select
  pp%u_sr_origin_coef(1,1,ik) = a1
  pp%u_sr_origin_coef(2,1,ik) = a3
end subroutine build_local_sr_shared_u

!-----------------------------------------------------------------------------------------------------------------------------------

subroutine build_local_sr_shared_u_stress_spline(pp, ik)
  use salmon_math, only: spline
  use structures, only: s_pp_info
  implicit none
  type(s_pp_info), intent(inout) :: pp
  integer, intent(in) :: ik
  integer :: i, nk
  real(8) :: xn(0:pp%nrmax), yn(0:pp%nrmax)
  real(8) :: an(0:pp%nrmax-1), bn(0:pp%nrmax-1), cn(0:pp%nrmax-1), dn(0:pp%nrmax-1)

  pp%nr_u_sr_stress(ik) = 0
  pp%rad_u_sr_stress(:,ik) = 0d0
  pp%u_sr_stress_tbl(:,ik) = 0d0
  pp%u_sr_stress_a(:,ik) = 0d0
  pp%u_sr_stress_b(:,ik) = 0d0
  pp%u_sr_stress_c(:,ik) = 0d0
  pp%u_sr_stress_d(:,ik) = 0d0

  xn = 0d0
  yn = 0d0
  an = 0d0
  bn = 0d0
  cn = 0d0
  dn = 0d0
  xn(0) = 0d0
  yn(0) = 0d0
  nk = 1

  do i = 1, pp%nrloc(ik)
    if (pp%rad(i,ik) <= 0d0) cycle
    if (pp%rad(i,ik) <= xn(nk-1)) stop 'build_local_sr_shared_u_stress_spline: radial grid must be strictly increasing'
    xn(nk) = pp%rad(i,ik)
    yn(nk) = pp%rad(i,ik) * pp%vloctbl(i,ik) + dble(pp%zps(ik))
    nk = nk + 1
  end do
  if (nk < 3) stop 'build_local_sr_shared_u_stress_spline: need origin plus 2 positive local radial points'

  pp%nr_u_sr_stress(ik) = nk
  pp%rad_u_sr_stress(0:nk-1,ik) = xn(0:nk-1)
  pp%u_sr_stress_tbl(0:nk-1,ik) = yn(0:nk-1)
  call spline(nk, xn(0:nk-1), yn(0:nk-1), an(0:nk-2), bn(0:nk-2), cn(0:nk-2), dn(0:nk-2))
  pp%u_sr_stress_a(0:nk-2,ik) = an(0:nk-2)
  pp%u_sr_stress_b(0:nk-2,ik) = bn(0:nk-2)
  pp%u_sr_stress_c(0:nk-2,ik) = cn(0:nk-2)
  pp%u_sr_stress_d(0:nk-2,ik) = dn(0:nk-2)
end subroutine build_local_sr_shared_u_stress_spline

!-----------------------------------------------------------------------------------------------------------------------------------

pure subroutine eval_local_sr_shared_u(pp, ik, r, u_r, du_r, intr)
  use structures, only: s_pp_info
  implicit none
  type(s_pp_info), intent(in) :: pp
  integer, intent(in) :: ik
  real(8), intent(in) :: r
  real(8), intent(out) :: u_r, du_r
  integer, intent(in), optional :: intr
  integer :: i, seg
  real(8) :: a1, a3, rsw

  a1 = pp%u_sr_origin_coef(1,1,ik)
  a3 = pp%u_sr_origin_coef(2,1,ik)
  rsw = pp%r_sr_origin_switch(ik)

  if (r < rsw) then
    u_r = a1*r + a3*r*r*r
    du_r = a1 + 3d0*a3*r*r
    return
  end if

  if (present(intr)) then
    seg = intr
  else
    seg = 1
    do i = 1, pp%nrps(ik) - 1
      if (r < pp%rad(i+1,ik)) exit
      seg = i
    end do
  end if
  if (pp%nrps(ik) > 1) then
    seg = max(1, min(seg, pp%nrps(ik)-1))
  else
    seg = 1
  end if

  u_r = pp%u_sr_tbl(seg,ik) + pp%du_sr_seg(seg,ik) * (r - pp%rad(seg,ik))
  du_r = pp%du_sr_seg(seg,ik)
end subroutine eval_local_sr_shared_u

!-----------------------------------------------------------------------------------------------------------------------------------

pure subroutine eval_local_sr_shared_u_stress(pp, ik, r, u_r, du_r, intr)
  use structures, only: s_pp_info
  implicit none
  type(s_pp_info), intent(in) :: pp
  integer, intent(in) :: ik
  real(8), intent(in) :: r
  real(8), intent(out) :: u_r, du_r
  integer, intent(in), optional :: intr
  integer :: i, nk, seg
  real(8) :: dx

  nk = pp%nr_u_sr_stress(ik)
  if (nk < 2) then
    u_r = 0d0
    du_r = 0d0
    return
  end if

  if (r <= 0d0) then
    u_r = 0d0
    du_r = pp%u_sr_stress_c(0,ik)
    return
  end if

  if (present(intr)) then
    seg = max(1, min(intr, nk-1))
    if (r < pp%rad_u_sr_stress(seg-1,ik) .or. (seg < nk-1 .and. r >= pp%rad_u_sr_stress(seg,ik))) then
      seg = nk - 1
      do i = 1, nk - 1
        if (r < pp%rad_u_sr_stress(i,ik)) then
          seg = i
          exit
        end if
      end do
    end if
  else
    seg = nk - 1
    do i = 1, nk - 1
      if (r < pp%rad_u_sr_stress(i,ik)) then
        seg = i
        exit
      end if
    end do
  end if

  dx = r - pp%rad_u_sr_stress(seg-1,ik)
  u_r = ((pp%u_sr_stress_a(seg-1,ik) * dx + pp%u_sr_stress_b(seg-1,ik)) * dx + pp%u_sr_stress_c(seg-1,ik)) * dx &
    & + pp%u_sr_stress_d(seg-1,ik)
  du_r = (3d0 * pp%u_sr_stress_a(seg-1,ik) * dx + 2d0 * pp%u_sr_stress_b(seg-1,ik)) * dx &
    & + pp%u_sr_stress_c(seg-1,ik)
end subroutine eval_local_sr_shared_u_stress

!-----------------------------------------------------------------------------------------------------------------------------------

real(8) function local_sr_aux_target_dr(gmag, npt_loc_sr_aux_2pi) result(dr_target)
  use math_constants, only: pi
  implicit none
  real(8), intent(in) :: gmag
  integer, intent(in) :: npt_loc_sr_aux_2pi

  if (npt_loc_sr_aux_2pi <= 0 .or. gmag <= 0d0) then
    dr_target = huge(1d0)
  else
    dr_target = 2d0*pi / (dble(npt_loc_sr_aux_2pi) * gmag)
  end if
end function local_sr_aux_target_dr

!-----------------------------------------------------------------------------------------------------------------------------------

real(8) function integrate_local_sr_vg_shell(pp, ik, gmag, npt_loc_sr_aux_2pi) result(vg)
  use structures, only: s_pp_info
  use math_constants, only: pi
  implicit none
  type(s_pp_info), intent(in) :: pp
  integer, intent(in) :: ik
  real(8), intent(in) :: gmag
  integer, intent(in) :: npt_loc_sr_aux_2pi
  integer :: i, isub, nsub
  real(8) :: r_lo, r_hi, r_mid, dr_seg, dr_sub, dr_target, u_r, du_r

  vg = 0d0
  do i = 2, pp%nrloc(ik)
    r_lo = pp%rad(i-1,ik)
    r_hi = pp%rad(i,ik)
    dr_seg = r_hi - r_lo
    if (dr_seg <= 0d0) cycle
    nsub = 1
    if (npt_loc_sr_aux_2pi > 0 .and. gmag > 0d0) then
      dr_target = local_sr_aux_target_dr(gmag, npt_loc_sr_aux_2pi)
      nsub = max(1, ceiling(dr_seg / dr_target))
    end if
    dr_sub = dr_seg / dble(nsub)
    do isub = 1, nsub
      r_mid = r_lo + (dble(isub) - 0.5d0) * dr_sub
      call eval_local_sr_shared_u(pp,ik,r_mid,u_r,du_r,i-1)
      if (gmag == 0d0) then
        vg = vg + 4d0*pi*r_mid*u_r*dr_sub
      else
        vg = vg + 4d0*pi*sin(gmag*r_mid)/gmag*u_r*dr_sub
      end if
    end do
  end do
end function integrate_local_sr_vg_shell

!-----------------------------------------------------------------------------------------------------------------------------------

real(8) function integrate_local_sr_dvg_dg2_shell(pp, ik, gmag, npt_loc_sr_aux_2pi) result(dvg)
  use structures, only: s_pp_info
  use math_constants, only: pi
  implicit none
  type(s_pp_info), intent(in) :: pp
  integer, intent(in) :: ik
  real(8), intent(in) :: gmag
  integer, intent(in) :: npt_loc_sr_aux_2pi
  integer :: i, isub, nsub
  real(8) :: r_lo, r_hi, r_mid, dr_seg, dr_sub, dr_target, u_r, du_r, x

  dvg = 0d0
  if (gmag <= 0d0) return

  do i = 2, pp%nrloc(ik)
    r_lo = pp%rad(i-1,ik)
    r_hi = pp%rad(i,ik)
    dr_seg = r_hi - r_lo
    if (dr_seg <= 0d0) cycle
    nsub = 1
    if (npt_loc_sr_aux_2pi > 0) then
      dr_target = local_sr_aux_target_dr(gmag, npt_loc_sr_aux_2pi)
      nsub = max(1, ceiling(dr_seg / dr_target))
    end if
    dr_sub = dr_seg / dble(nsub)
    do isub = 1, nsub
      r_mid = r_lo + (dble(isub) - 0.5d0) * dr_sub
      call eval_local_sr_shared_u(pp,ik,r_mid,u_r,du_r,i-1)
      x = gmag * r_mid
      if (x < 1d-2) then
        dvg = dvg + 4d0*pi*u_r * r_mid**3 * (-1d0/6d0 + x**2/60d0 - x**4/1680d0) * dr_sub
      else
        dvg = dvg + 4d0*pi*u_r * (x*cos(x) - sin(x)) / (2d0*gmag**3) * dr_sub
      end if
    end do
  end do
end function integrate_local_sr_dvg_dg2_shell

!===================================================================================================================================

real(8) function integrate_local_sr_vg_shell_stress(pp, ik, gmag, npt_loc_sr_aux_2pi) result(vg)
  use structures, only: s_pp_info
  use math_constants, only: pi
  implicit none
  type(s_pp_info), intent(in) :: pp
  integer, intent(in) :: ik
  real(8), intent(in) :: gmag
  integer, intent(in) :: npt_loc_sr_aux_2pi
  integer :: i, isub, nk, nsub
  real(8) :: r_lo, r_hi, r_mid, dr_seg, dr_sub, dr_target, u_r, du_r

  vg = 0d0
  nk = pp%nr_u_sr_stress(ik)
  if (nk < 2) return

  do i = 1, nk - 1
    r_lo = pp%rad_u_sr_stress(i-1,ik)
    r_hi = pp%rad_u_sr_stress(i,ik)
    dr_seg = r_hi - r_lo
    if (dr_seg <= 0d0) cycle
    nsub = 1
    if (npt_loc_sr_aux_2pi > 0 .and. gmag > 0d0) then
      dr_target = local_sr_aux_target_dr(gmag, npt_loc_sr_aux_2pi)
      nsub = max(1, ceiling(dr_seg / dr_target))
    end if
    dr_sub = dr_seg / dble(nsub)
    do isub = 1, nsub
      r_mid = r_lo + (dble(isub) - 0.5d0) * dr_sub
      call eval_local_sr_shared_u_stress(pp,ik,r_mid,u_r,du_r,i)
      if (gmag == 0d0) then
        vg = vg + 4d0*pi*r_mid*u_r*dr_sub
      else
        vg = vg + 4d0*pi*sin(gmag*r_mid)/gmag*u_r*dr_sub
      end if
    end do
  end do
end function integrate_local_sr_vg_shell_stress

!-----------------------------------------------------------------------------------------------------------------------------------

real(8) function integrate_local_sr_dvg_dg2_shell_stress(pp, ik, gmag, npt_loc_sr_aux_2pi) result(dvg)
  use structures, only: s_pp_info
  use math_constants, only: pi
  implicit none
  type(s_pp_info), intent(in) :: pp
  integer, intent(in) :: ik
  real(8), intent(in) :: gmag
  integer, intent(in) :: npt_loc_sr_aux_2pi
  integer :: i, isub, nk, nsub
  real(8) :: r_lo, r_hi, r_mid, dr_seg, dr_sub, dr_target, u_r, du_r, x

  dvg = 0d0
  if (gmag <= 0d0) return
  nk = pp%nr_u_sr_stress(ik)
  if (nk < 2) return

  do i = 1, nk - 1
    r_lo = pp%rad_u_sr_stress(i-1,ik)
    r_hi = pp%rad_u_sr_stress(i,ik)
    dr_seg = r_hi - r_lo
    if (dr_seg <= 0d0) cycle
    nsub = 1
    if (npt_loc_sr_aux_2pi > 0) then
      dr_target = local_sr_aux_target_dr(gmag, npt_loc_sr_aux_2pi)
      nsub = max(1, ceiling(dr_seg / dr_target))
    end if
    dr_sub = dr_seg / dble(nsub)
    do isub = 1, nsub
      r_mid = r_lo + (dble(isub) - 0.5d0) * dr_sub
      call eval_local_sr_shared_u_stress(pp,ik,r_mid,u_r,du_r,i)
      x = gmag * r_mid
      if (x < 1d-2) then
        dvg = dvg + 4d0*pi*u_r * r_mid**3 * (-1d0/6d0 + x**2/60d0 - x**4/1680d0) * dr_sub
      else
        dvg = dvg + 4d0*pi*u_r * (x*cos(x) - sin(x)) / (2d0*gmag**3) * dr_sub
      end if
    end do
  end do
end function integrate_local_sr_dvg_dg2_shell_stress

!===================================================================================================================================

subroutine write_local_sr_u_diagnostics(pp)
  use communication, only: comm_is_root
  use filesystem, only: open_filehandle
  use parallelization, only: nproc_id_global
  use salmon_global, only: base_directory, ps_format, sysname
  use structures, only: s_pp_info
  implicit none
  type(s_pp_info), intent(in) :: pp
  integer :: fh, i, i1, i2, i3, ik, nlcc_flag, nr, seg
  real(8) :: a1, a3, rsw, r1, r2, r3, r_at_max_du_jump, dr_jump
  real(8) :: u_r1, u_r2, u_r3, u_switch_poly, du_switch_poly
  real(8) :: u_switch_seg, du_switch_seg, delta_u_switch, delta_du_switch
  real(8) :: u_at_rps, du_tail, max_abs_du_jump
  logical :: has_nlcc
  character(256) :: filename

  if (.not. comm_is_root(nproc_id_global)) return

  filename = trim(base_directory)//trim(sysname)//'_pp_local_sr_u_summary.data'
  fh = open_filehandle(filename, status="replace")
  write(fh, '("#",1X,A)') 'Local-SR shared-u property summary'
  write(fh, '("#",1X,A)') 'u(r) = r * V_sr(r) = r * V_loc(r) + Z'
  write(fh, '("#",99(1X,I0,":",A,"[",A,"]"))') &
    & 1, "ik", "none", &
    & 2, "symbol", "none", &
    & 3, "ps_format", "none", &
    & 4, "nlcc", "none", &
    & 5, "nrps", "none", &
    & 6, "nrloc", "none", &
    & 7, "rps", "a.u.", &
    & 8, "r1", "a.u.", &
    & 9, "r2", "a.u.", &
    & 10, "r3", "a.u.", &
    & 11, "r_switch", "a.u.", &
    & 12, "a1", "Ha", &
    & 13, "a3", "Ha/a.u.^2", &
    & 14, "u_r1", "Ha*a.u.", &
    & 15, "u_r2", "Ha*a.u.", &
    & 16, "u_r3", "Ha*a.u.", &
    & 17, "delta_u_switch", "Ha*a.u.", &
    & 18, "delta_du_switch", "Ha", &
    & 19, "u_at_rps", "Ha*a.u.", &
    & 20, "du_tail", "Ha", &
    & 21, "max_abs_du_jump", "Ha", &
    & 22, "r_at_max_du_jump", "a.u."

  do ik = 1, size(pp%atom_symbol)
    nr = pp%nrps(ik)
    if (nr < 3) cycle

    i1 = 0
    i2 = 0
    i3 = 0
    do i = 1, nr
      if (pp%rad(i,ik) <= 0d0) cycle
      if (i1 == 0) then
        i1 = i
      else if (i2 == 0) then
        i2 = i
      else
        i3 = i
        exit
      end if
    end do
    if (i1 == 0 .or. i2 == 0 .or. i3 == 0) cycle

    r1 = pp%rad(i1,ik)
    r2 = pp%rad(i2,ik)
    r3 = pp%rad(i3,ik)
    u_r1 = pp%u_sr_tbl(i1,ik)
    u_r2 = pp%u_sr_tbl(i2,ik)
    u_r3 = pp%u_sr_tbl(i3,ik)

    a1 = pp%u_sr_origin_coef(1,1,ik)
    a3 = pp%u_sr_origin_coef(2,1,ik)
    rsw = pp%r_sr_origin_switch(ik)
    u_switch_poly = a1*rsw + a3*rsw**3
    du_switch_poly = a1 + 3d0*a3*rsw**2
    call eval_local_sr_shared_u(pp, ik, rsw, u_switch_seg, du_switch_seg)
    delta_u_switch = u_switch_seg - u_switch_poly
    delta_du_switch = du_switch_seg - du_switch_poly

    seg = max(1, min(pp%nrloc(ik)-1, nr-1))
    u_at_rps = pp%u_sr_tbl(pp%nrloc(ik),ik)
    du_tail = pp%du_sr_seg(seg,ik)

    max_abs_du_jump = 0d0
    r_at_max_du_jump = pp%rad(i2,ik)
    do i = 2, nr - 1
      dr_jump = abs(pp%du_sr_seg(i,ik) - pp%du_sr_seg(i-1,ik))
      if (dr_jump > max_abs_du_jump) then
        max_abs_du_jump = dr_jump
        r_at_max_du_jump = pp%rad(i,ik)
      end if
    end do

    has_nlcc = maxval(abs(pp%rho_nlcc_tbl(:,ik))) > 0d0
    nlcc_flag = 0
    if (has_nlcc) nlcc_flag = 1

    write(fh,'(1X,I4,1X,A8,1X,A10,1X,I1,1X,I6,1X,I6,16(1X,ES23.15E3))') &
      & ik, trim(pp%atom_symbol(ik)), trim(ps_format(ik)), nlcc_flag, pp%nrps(ik), pp%nrloc(ik), &
      & pp%rps(ik), r1, r2, r3, rsw, a1, a3, u_r1, u_r2, u_r3, delta_u_switch, delta_du_switch, &
      & u_at_rps, du_tail, max_abs_du_jump, r_at_max_du_jump
  end do

  close(fh)
end subroutine write_local_sr_u_diagnostics

!===================================================================================================================================

subroutine write_vloctbl_derivative_diagnostics(pp)
  use communication, only: comm_is_root
  use filesystem, only: open_filehandle
  use parallelization, only: nproc_id_global
  use salmon_global, only: base_directory, ps_format, quiet, sysname
  use structures, only: s_pp_info
  implicit none
  type(s_pp_info), intent(in) :: pp
  integer :: fh, i, ik, nlcc_flag, nr
  real(8), allocatable :: d_num(:)
  real(8) :: d_mid, d_seg, dr, dr_max, dr_min, max_abs_node_diff
  real(8) :: max_abs_seg_end_diff, max_abs_seg_mid_diff, mean_abs_seg_mid_diff
  real(8) :: r1, r2, r3, r4, r_at_max_node, r_at_max_seg_end, r_at_max_seg_mid
  real(8) :: r_mid, slope, sum_abs_seg_mid_diff, v_lo, v_hi
  logical :: has_nlcc
  character(256) :: filename
  character(7) :: nlcc_label

  if (.not. comm_is_root(nproc_id_global)) return

  filename = trim(base_directory)//trim(sysname)//'_pp_local_derivative_check.data'
  fh = open_filehandle(filename, status="replace")
  write(fh, '("#",1X,A)') 'Local-potential derivative diagnostics'
  write(fh, '("#",1X,A)') 'Compare stored dvloctbl against fresh node-wise finite differences and interval slopes'
  write(fh, '("#",99(1X,I0,":",A,"[",A,"]"))') &
    & 1, "ik", "none", &
    & 2, "symbol", "none", &
    & 3, "ps_format", "none", &
    & 4, "nlcc", "none", &
    & 5, "nrloc", "none", &
    & 6, "rps", "a.u.", &
    & 7, "dr_min", "a.u.", &
    & 8, "dr_max", "a.u.", &
    & 9, "max_abs_node_diff", "Ha/a.u.", &
    & 10, "r_at_max_node", "a.u.", &
    & 11, "max_abs_seg_end_diff", "Ha/a.u.", &
    & 12, "r_at_max_seg_end", "a.u.", &
    & 13, "max_abs_seg_mid_diff", "Ha/a.u.", &
    & 14, "r_at_max_seg_mid", "a.u.", &
    & 15, "mean_abs_seg_mid_diff", "Ha/a.u."

  do ik = 1, size(pp%atom_symbol)
    nr = pp%nrloc(ik)
    if (nr < 3) cycle

    allocate(d_num(nr))
    d_num = 0d0

    do i = 2, nr - 1
      r1 = pp%rad(i+1,ik) - pp%rad(i,ik)
      r2 = pp%rad(i+1,ik) - pp%rad(i+2,ik)
      r3 = pp%rad(i+2,ik) - pp%rad(i,ik)
      r4 = r1 / r2
      d_num(i) = (r4 + 1d0) * (pp%vloctbl(i,ik) - pp%vloctbl(i-1,ik)) / r1 &
               - (pp%vloctbl(i+1,ik) - pp%vloctbl(i-1,ik)) / r3 * r4
    end do
    d_num(1) = d_num(2) - (d_num(3) - d_num(2)) / (pp%rad(3,ik) - pp%rad(2,ik)) &
             * (pp%rad(2,ik) - pp%rad(1,ik))
    if (nr < ubound(pp%rad,1)) then
      d_num(nr) = d_num(nr-1) + (d_num(nr-1) - d_num(nr-2)) / (pp%rad(nr,ik) - pp%rad(nr-1,ik)) &
                * (pp%rad(nr+1,ik) - pp%rad(nr,ik))
    else
      d_num(nr) = d_num(nr-1)
    end if

    dr_min = huge(1d0)
    dr_max = 0d0
    max_abs_node_diff = 0d0
    max_abs_seg_end_diff = 0d0
    max_abs_seg_mid_diff = 0d0
    r_at_max_node = pp%rad(1,ik)
    r_at_max_seg_end = 0.5d0 * (pp%rad(1,ik) + pp%rad(2,ik))
    r_at_max_seg_mid = r_at_max_seg_end
    sum_abs_seg_mid_diff = 0d0

    do i = 1, nr
      if (abs(pp%dvloctbl(i,ik) - d_num(i)) > max_abs_node_diff) then
        max_abs_node_diff = abs(pp%dvloctbl(i,ik) - d_num(i))
        r_at_max_node = pp%rad(i,ik)
      end if
    end do

    do i = 1, nr - 1
      dr = pp%rad(i+1,ik) - pp%rad(i,ik)
      if (dr <= 0d0) cycle
      dr_min = min(dr_min, dr)
      dr_max = max(dr_max, dr)
      v_lo = pp%vloctbl(i,ik)
      v_hi = pp%vloctbl(i+1,ik)
      slope = (v_hi - v_lo) / dr
      d_seg = max(abs(pp%dvloctbl(i,ik) - slope), abs(pp%dvloctbl(i+1,ik) - slope))
      d_mid = abs(0.5d0 * (pp%dvloctbl(i,ik) + pp%dvloctbl(i+1,ik)) - slope)
      r_mid = 0.5d0 * (pp%rad(i,ik) + pp%rad(i+1,ik))
      if (d_seg > max_abs_seg_end_diff) then
        max_abs_seg_end_diff = d_seg
        r_at_max_seg_end = r_mid
      end if
      if (d_mid > max_abs_seg_mid_diff) then
        max_abs_seg_mid_diff = d_mid
        r_at_max_seg_mid = r_mid
      end if
      sum_abs_seg_mid_diff = sum_abs_seg_mid_diff + d_mid
    end do

    if (dr_min == huge(1d0)) dr_min = 0d0
    mean_abs_seg_mid_diff = sum_abs_seg_mid_diff / max(1, nr - 1)
    has_nlcc = maxval(abs(pp%rho_nlcc_tbl(:,ik))) > 0d0
    nlcc_flag = 0
    if (has_nlcc) nlcc_flag = 1
    if (has_nlcc) then
      nlcc_label = 'nlcc'
    else
      nlcc_label = 'no-nlcc'
    end if

    write(fh, '(I4,1X,A2,1X,A16,1X,I1,1X,I8,1X,10(ES23.15E3,1X))') &
      & ik, trim(pp%atom_symbol(ik)), trim(ps_format(ik)), nlcc_flag, nr, &
      & pp%rps(ik), dr_min, dr_max, max_abs_node_diff, r_at_max_node, &
      & max_abs_seg_end_diff, r_at_max_seg_end, max_abs_seg_mid_diff, &
      & r_at_max_seg_mid, mean_abs_seg_mid_diff

    if (.not. quiet) then
      write(*, '(A,1X,I0,1X,A,1X,A,1X,A,1X,3(A,ES12.4E3))') &
        & 'pp-deriv-check', ik, trim(pp%atom_symbol(ik)), trim(ps_format(ik)), nlcc_label, &
        & 'node=', max_abs_node_diff, ' seg_mid=', max_abs_seg_mid_diff, &
        & ' seg_end=', max_abs_seg_end_diff
    end if

    deallocate(d_num)
  end do

  close(fh)
  if (.not. quiet) write(*, '(A,1X,A)') 'wrote local-potential derivative diagnostics:', trim(filename)

end subroutine write_vloctbl_derivative_diagnostics

!===================================================================================================================================

SUBROUTINE dealloc_init_ps(ppg)
  use structures, only: s_pp_grid
  implicit none
  type(s_pp_grid) :: ppg

  deallocate(ppg%jxyz, ppg%rxyz, ppg%uv)
  if(allocated(ppg%zekr_uV)) deallocate(ppg%zekr_uV)

  if (allocated(ppg%irange_atom))    deallocate(ppg%irange_atom)
  if (allocated(ppg%ireferred_atom)) deallocate(ppg%ireferred_atom)
  if (allocated(ppg%ilocal_nlma2ilma)) deallocate(ppg%ilocal_nlma2ilma)
  if (allocated(ppg%ilocal_nlma2ia))   deallocate(ppg%ilocal_nlma2ia)
  if (allocated(ppg%uVpsibox))       deallocate(ppg%uVpsibox)
  
  return
END SUBROUTINE dealloc_init_ps

!===================================================================================================================================

SUBROUTINE calc_Vpsl_isolated(lg,mg,system,info,pp,fg,vpsl,ppg,property)
  use structures
  use salmon_global,only : natom, kion, quiet, method_poisson, nelem, yn_ffte, yn_out_stress, yn_spinorbit, npt_loc_sr_aux_2pi
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
  real(8) :: g(3),gd,g2sq
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
  if(yn_spinorbit == 'y') then
    allocate(ppg%zekr_uV_so(ppg%nps,size(ppg%ia_tbl_so),1,2,1))
    ppg%zekr_uV_so(:,:,1,1:2,1) = dcmplx(ppg%uv_so(:,:,1:2,1))
  end if

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
      if(yn_out_stress == 'y') then
        allocate(ppg%zVG_ion_stress(ifgx_s:ifgx_e,ifgy_s:ifgy_e,ifgz_s:ifgz_e,nelem))
        ppg%zVG_ion_stress = 0d0
        allocate(ppg%dVG_ion_dG2(ifgx_s:ifgx_e,ifgy_s:ifgy_e,ifgz_s:ifgz_e,nelem))
        ppg%dVG_ion_dG2 = 0d0
      end if

      ppg%zVG_ion = 0d0


  !$omp parallel
  !$omp do private(ik,ix,iy,iz,g,g2sq) collapse(3)
      do ik=1,nelem
        do iz=ifgz_s,ifgz_e
        do iy=ifgy_s,ifgy_e
        do ix=ifgx_s,ifgx_e
          g(1) = fg%vec_G(1,ix,iy,iz)
          g(2) = fg%vec_G(2,ix,iy,iz)
          g(3) = fg%vec_G(3,ix,iy,iz)
          g2sq = sqrt(g(1)**2+g(2)**2+g(3)**2)
          ppg%zVG_ion(ix,iy,iz,ik) = integrate_local_sr_vg_shell(pp, ik, g2sq, npt_loc_sr_aux_2pi)
        end do
        end do
        end do
      end do
  !$omp end do
  !$omp end parallel

      if(yn_out_stress == 'y') then
  !$omp parallel do private(ik,ix,iy,iz,g,g2sq) collapse(3)
        do ik=1,nelem
          do iz=ifgz_s,ifgz_e
          do iy=ifgy_s,ifgy_e
          do ix=ifgx_s,ifgx_e
            g(1) = fg%vec_G(1,ix,iy,iz)
            g(2) = fg%vec_G(2,ix,iy,iz)
            g(3) = fg%vec_G(3,ix,iy,iz)
            g2sq = sqrt(g(1)**2+g(2)**2+g(3)**2)
            ppg%zVG_ion_stress(ix,iy,iz,ik) = integrate_local_sr_vg_shell_stress(pp, ik, g2sq, npt_loc_sr_aux_2pi)
          end do
          end do
          end do
        end do
  !$omp end parallel do

  !$omp parallel do private(ik,ix,iy,iz,g,g2sq) collapse(3)
        do ik=1,nelem
          do iz=ifgz_s,ifgz_e
          do iy=ifgy_s,ifgy_e
          do ix=ifgx_s,ifgx_e
            if(fg%if_Gzero(ix,iy,iz)) cycle
            g(1) = fg%vec_G(1,ix,iy,iz)
            g(2) = fg%vec_G(2,ix,iy,iz)
            g(3) = fg%vec_G(3,ix,iy,iz)
            g2sq = sqrt(g(1)**2+g(2)**2+g(3)**2)
            ppg%dVG_ion_dG2(ix,iy,iz,ik) = integrate_local_sr_dvg_dg2_shell_stress(pp, ik, g2sq, npt_loc_sr_aux_2pi)
          end do
          end do
          end do
        end do
  !$omp end parallel do
      end if

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
  use salmon_global,only : nelem, kion, yn_ffte, yn_out_stress, npt_loc_sr_aux_2pi
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
  real(8) :: g(3),gd,g2sq
  complex(8) :: tmp_exp
  complex(8) :: vtmp1(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:2)
  complex(8) :: vtmp2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:2)

  if( property == 'initial' ) then
  
    allocate(ppg%zrhoG_ion(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)) & ! rho_ion(G)
          & ,ppg%zVG_ion  (mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),nelem)) ! V_ion(G)
    if(yn_out_stress == 'y') then
      allocate(ppg%zVG_ion_stress(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),nelem))
      ppg%zVG_ion_stress = 0d0
      allocate(ppg%dVG_ion_dG2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),nelem))
      ppg%dVG_ion_dG2 = 0d0
    end if

    ppg%zVG_ion = 0d0
  !$omp parallel
  !$omp do private(ik,ix,iy,iz,g,g2sq) collapse(3)
    do ik=1,nelem
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        g(1) = fg%vec_G(1,ix,iy,iz)
        g(2) = fg%vec_G(2,ix,iy,iz)
        g(3) = fg%vec_G(3,ix,iy,iz)
        g2sq = sqrt(g(1)**2+g(2)**2+g(3)**2)
        ppg%zVG_ion(ix,iy,iz,ik) = integrate_local_sr_vg_shell(pp, ik, g2sq, npt_loc_sr_aux_2pi)
      end do
      end do
      end do
    end do
  !$omp end do
  !$omp end parallel

    if(yn_out_stress == 'y') then
  !$omp parallel do private(ik,ix,iy,iz,g,g2sq) collapse(3)
      do ik=1,nelem
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          g(1) = fg%vec_G(1,ix,iy,iz)
          g(2) = fg%vec_G(2,ix,iy,iz)
          g(3) = fg%vec_G(3,ix,iy,iz)
          g2sq = sqrt(g(1)**2+g(2)**2+g(3)**2)
          ppg%zVG_ion_stress(ix,iy,iz,ik) = integrate_local_sr_vg_shell_stress(pp, ik, g2sq, npt_loc_sr_aux_2pi)
        end do
        end do
        end do
      end do
  !$omp end parallel do

  !$omp parallel do private(ik,ix,iy,iz,g,g2sq) collapse(3)
      do ik=1,nelem
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          if(fg%if_Gzero(ix,iy,iz)) cycle
          g(1) = fg%vec_G(1,ix,iy,iz)
          g(2) = fg%vec_G(2,ix,iy,iz)
          g(3) = fg%vec_G(3,ix,iy,iz)
          g2sq = sqrt(g(1)**2+g(2)**2+g(3)**2)
          ppg%dVG_ion_dG2(ix,iy,iz,ik) = integrate_local_sr_dvg_dg2_shell_stress(pp, ik, g2sq, npt_loc_sr_aux_2pi)
        end do
        end do
        end do
      end do
  !$omp end parallel do
    end if

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
!$acc parallel loop private(ia,ilma) copyin(ppg)
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
!$omp parallel do private(ilma,ia) reduction(+:ilocal_nlma)
  do ilma=1,ppg%nlma
    ia = ppg%ia_tbl(ilma)
    if (ppg%ireferred_atom(ia)) ilocal_nlma = ilocal_nlma + 1
  end do
!$omp end parallel do
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

!===================================================================================================================================

subroutine init_uvpsi_blocking(ppg,lg)
  use structures, only:s_pp_grid,s_rgrid
  implicit none
  type(s_pp_grid), intent(inout) :: ppg
  type(s_rgrid),   intent(in)    :: lg

  integer, allocatable :: v2nlma(:), i2vi(:)
  integer :: nl,ilma,i,j,ix,iy,iz,max_nlma,vi,n,ia

  nl = product(lg%num)
  allocate(v2nlma(0:nl-1))

  v2nlma = 0
  do ilma=1,ppg%nlma
    ia = ppg%ia_tbl(ilma)

    do j=1,ppg%mps(ia)
      ix = ppg%jxyz(1,j,ia) - lg%is(1)
      iy = ppg%jxyz(2,j,ia) - lg%is(2)
      iz = ppg%jxyz(3,j,ia) - lg%is(3)
      i  = ix + iy*lg%num(1) + iz*lg%num(1)*lg%num(2)
      v2nlma(i) = v2nlma(i) + 1
    end do
  end do

  allocate(i2vi(0:nl-1))
  allocate(ppg%v2j(3, 0:nl-1))

  max_nlma = 0
  vi = 0
  do iz=0,lg%num(3)-1
  do iy=0,lg%num(2)-1
  do ix=0,lg%num(1)-1
    i = ix + iy*lg%num(1) + iz*lg%num(1)*lg%num(2)

    max_nlma = max(max_nlma, v2nlma(i))

    i2vi(i) = -1
    if (v2nlma(i) > 0) then
      ppg%v2j(1,vi) = ix + lg%is(1)
      ppg%v2j(2,vi) = iy + lg%is(2)
      ppg%v2j(3,vi) = iz + lg%is(3)
      i2vi(i)       = vi
      vi            = vi + 1
    end if
  end do
  end do
  end do
  ppg%max_vi = vi

  allocate(ppg%v2nlma(0:ppg%max_vi-1))
  allocate(ppg%k2ilma(0:ppg%max_vi-1, max_nlma))
  allocate(ppg%k2j(0:ppg%max_vi-1, max_nlma))
  
  ppg%v2nlma = 0

  do ilma=1,ppg%nlma
    ia = ppg%ia_tbl(ilma)

    do j=1,ppg%mps(ia)
      ix = ppg%jxyz(1,j,ia) - lg%is(1)
      iy = ppg%jxyz(2,j,ia) - lg%is(2)
      iz = ppg%jxyz(3,j,ia) - lg%is(3)
      i  = ix + iy*lg%num(1) + iz*lg%num(1)*lg%num(2)

      vi = i2vi(i)

      ppg%v2nlma(vi) = ppg%v2nlma(vi) + 1
      n = ppg%v2nlma(vi)

      ppg%k2ilma(vi,n) = ilma
      ppg%k2j   (vi,n) = j
    end do
  end do

end subroutine init_uvpsi_blocking

end module prep_pp_sub
