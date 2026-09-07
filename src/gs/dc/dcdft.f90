!
!  Copyright 2019-2024 SALMON developers
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
module dcdft
  use dc_fragment_geometry, only: get_fragment_domain, optimize_fragment_geometry
  implicit none
contains

  subroutine init_dcdft(dc,pp,mixing,ewald)
    use structures
    use salmon_global, only: nproc_k, nproc_ob, nproc_rgrid, nproc_rgrid_tot &
    & , nstate, nelec, yn_dc, nstate_frag
    implicit none
    type(s_dcdft)        ,intent(inout) :: dc
    type(s_pp_info)      ,intent(inout) :: pp
    type(s_mixing)       ,intent(inout) :: mixing
    type(s_ewald_ion_ion),intent(inout) :: ewald
    !
    integer :: nproc_ob_tmp, nproc_rgrid_tmp(3)
    
    if(nproc_k/=1) stop "DC method (yn_dc=y): nproc_k must be 1 for both the total system and fragments."
    nproc_ob_tmp = nproc_ob
    nproc_rgrid_tmp = nproc_rgrid
    dc%nstate_tot = nstate
    dc%nstate_frag = nstate_frag
    
  ! total system
    nproc_ob = 1 ! override
    nproc_rgrid = nproc_rgrid_tot ! override
    yn_dc = 't' ! override !!!!!! future work: remove
    call init_total
    
  ! fragment
    nproc_ob = nproc_ob_tmp ! override
    nproc_rgrid = nproc_rgrid_tmp ! override
    nstate = dc%nstate_frag ! override
    yn_dc = 'y' ! override !!!!!! future work: remove
    call init_comm_frag
    call init_fragment
    
  contains
  
    subroutine init_total
      use parallelization, only: nproc_group_global, nproc_id_global, nproc_size_global
      use initialization_sub, only: init_dft, init_nion_div
      use sendrecv_grid, only: dealloc_cache
      use mixing_sub, only: init_mixing
      use salmon_pp, only: read_pslfile
      use prep_pp_sub, only: init_ps
      use salmon_global, only: num_fragment, nelec, base_directory, method_init_density
      use filesystem, only: atomic_create_directory
      use read_gs, only: read_dns_cube
      use Total_Energy, only: init_ewald
      implicit none
      integer :: i
      type(s_stencil) :: stencil_dummy
      type(s_ofile) :: ofile_dummy
      
    ! MPI for the total system
      dc%icomm_tot = nproc_group_global
      dc%id_tot = nproc_id_global
      dc%isize_tot = nproc_size_global
      
    ! base_directory for the total system
      if(base_directory /= './') stop "DC method (yn_dc=y): base_directory must be default."
      dc%base_directory = './data_dcdft/total/'
      call atomic_create_directory(dc%base_directory,dc%icomm_tot,dc%id_tot)
      base_directory = dc%base_directory ! override
    
      dc%n_frag = num_fragment(1)*num_fragment(2)*num_fragment(3) ! # of the fragments
      dc%elec_num_tot = dble(nelec) ! # of total electrons
     
    ! initialization for the total system
      call init_dft(dc%icomm_tot,dc%info_tot,dc%lg_tot,dc%mg_tot,dc%system_tot, &
      & stencil_dummy,dc%fg_tot,dc%poisson_tot,dc%srg_tot,dc%srg_scalar_tot,ofile_dummy)
      deallocate(dc%system_tot%rocc)
      
      call allocate_scalar(dc%mg_tot,dc%rho_tot)
      call allocate_scalar(dc%mg_tot,dc%vh_tot)
      call allocate_scalar(dc%mg_tot,dc%vpsl_tot)
      allocate(dc%rho_tot_s(dc%system_tot%nspin),dc%vloc_tot(dc%system_tot%nspin),dc%vxc_tot(dc%system_tot%nspin))
      do i=1,dc%system_tot%nspin
        call allocate_scalar(dc%mg_tot,dc%rho_tot_s(i))
        call allocate_scalar(dc%mg_tot,dc%vloc_tot(i))
        call allocate_scalar(dc%mg_tot,dc%vxc_tot(i))
      end do
      
    ! mixing
      mixing%num_rho_stock = 21
      call init_mixing(dc%system_tot%nspin,dc%mg_tot,mixing)
      
    ! Vpsl
      call read_pslfile(dc%system_tot,pp)
      call init_ps(dc%lg_tot,dc%mg_tot,dc%system_tot,dc%info_tot,dc%fg_tot,dc%poisson_tot, &
      & pp,dc%ppg_tot,dc%vpsl_tot)
      
      if(method_init_density=='read_dns_cube') then
      ! read the initial density for the total system
        call read_dns_cube(dc%lg_tot,dc%mg_tot,dc%system_tot,dc%info_tot,dc%rho_tot,dc%rho_tot_s)
      end if
      
    ! Ewald
      call init_nion_div(dc%system_tot,dc%lg_tot,dc%mg_tot,dc%info_tot)
      call init_ewald(dc%system_tot,dc%info_tot,ewald)
    
    end subroutine init_total
  
    subroutine init_comm_frag
      use parallelization, only: nproc_group_global, nproc_id_global, nproc_size_global
      use communication, only: comm_create_group,comm_get_groupinfo
      use filesystem, only: atomic_create_directory
      use salmon_global, only: base_directory
      implicit none
      integer :: icomm_frag,isize_frag,id_frag
      integer :: npg,i,j,k,m
      
    ! set dc%i_frag (fragment index)
      npg = dc%isize_tot / dc%n_frag
      m = mod(dc%isize_tot,dc%n_frag) ! nproc = npg*dc%n_frag + m
      k=0
      do j=0,dc%n_frag-1
      do i=0,npg-1
        if(j*npg+i==dc%id_tot) then
          dc%i_frag=j
          k=1
          exit
          exit
        end if
      end do
      end do
      if(k==0) dc%i_frag = dc%id_tot-npg*dc%n_frag
      dc%i_frag = dc%i_frag + 1 ! = 1:dc%n_frag
      
    ! split communicator
      icomm_frag = comm_create_group(dc%icomm_tot,dc%i_frag,dc%id_tot) ! dc%i_frag : color, dc%id_tot : key
      call comm_get_groupinfo(icomm_frag, id_frag, isize_frag)
      
    ! MPI for the fragment
      dc%icomm_frag = icomm_frag
      dc%id_frag = id_frag
      dc%isize_frag = isize_frag
      
    ! Override global variables
      nproc_group_global = icomm_frag
      nproc_id_global = id_frag
      nproc_size_global = isize_frag
      write(base_directory, '(a, i6.6, a)') './data_dcdft/fragments/', dc%i_frag, '/'
      
    ! base_directory for the fragment
      call atomic_create_directory(base_directory,icomm_frag,id_frag)
      
!write(*,'(a,5i10)') " i_frag,id_F,isize_F,id,isize",dc%i_frag,id_frag,isize_frag,dc%id_tot,dc%isize_tot
      
    end subroutine init_comm_frag
    
    subroutine init_fragment
      use salmon_global, only: num_rgrid_buffer, kion, rion, natom, num_rgrid, al, num_fragment, &
      & yn_dc_fragment_optimization
      implicit none
      integer :: i_frag,n,i,j,k,ii,jj,kk
      integer :: iatom,iatom_frag,max_natom_frag
      integer, allocatable :: kion_frag(:,:),natom_frag(:)
      integer :: nxyz_domain(3), nxyz_domain_frag(3)
      real(8) :: dr
      real(8) :: r1(3),r2(3),r(3)
      real(8) :: ldomain(3),lbuffer(3),ldomain_frag(3)
      real(8), allocatable :: rion_frag(:,:,:)
    
    ! length of domain
      ldomain(1:3) = al(1:3) / dble(num_fragment(1:3))
      dc%optimized_fragment_geometry = .false.
      max_natom_frag = 27*natom
      allocate(kion_frag(max_natom_frag,dc%n_frag),natom_frag(dc%n_frag))
      allocate(rion_frag(3,max_natom_frag,dc%n_frag))
      kion_frag = 0
      natom_frag = 0
      rion_frag = 0d0
      
      do n=1,3 ! x,y,z
      ! rion --> rion = [0:al] (total system)
        do i=1,natom
          rion(n,i) = r_periodic(rion(n,i),al(n))
          if(rion(n,i) < 0d0 .or. rion(n,i) > al(n)) stop "DC method (yn_dc=y): rion"
        end do
      ! dc%nxyz_domain: # of grid points for each domain
      ! dc%nxyz_buffer: # of grid points for the buffer region
        if(mod(num_rgrid(n),num_fragment(n))==0) then
          dc%nxyz_domain(n) = num_rgrid(n) / num_fragment(n)
          dc%nxyz_buffer(n) = num_rgrid_buffer(n)
          dr = al(n)/dble(num_rgrid(n))
          lbuffer(n) = dr * dc%nxyz_buffer(n) ! length of the buffer region
        else
          stop "DC method (yn_dc=y): mod(num_rgrid,num_fragment) /= 0"
        end if
      end do ! n=x,y,z

      if (yn_dc_fragment_optimization == 'y') then
        call build_optimized_fragment_geometry(num_fragment, num_rgrid, al, natom, rion)
      else
        call build_uniform_fragment_geometry(num_fragment, ldomain)
      end if
      
    ! variables for each fragment
      i_frag = 1
      do i=1,num_fragment(1)
      do j=1,num_fragment(2)
      do k=1,num_fragment(3)
      ! boundaries of the fragment i_frag
        call get_fragment_domain(dc, i_frag, nxyz_domain_frag)
        ldomain_frag(1:3) = al(1:3) * dble(nxyz_domain_frag(1:3)) / dble(num_rgrid(1:3))
        r1 = dc%rxyz_frag(:,i_frag) - lbuffer
        r2 = dc%rxyz_frag(:,i_frag) + ldomain_frag + lbuffer
      ! atom count
        iatom_frag = 0
        do iatom=1,natom
          do ii=-1,1
          do jj=-1,1
          do kk=-1,1
            r(1:3) = rion(1:3,iatom) ! r = [0:al]
            r(1) = r(1) + dble(ii)*al(1)
            r(2) = r(2) + dble(jj)*al(2)
            r(3) = r(3) + dble(kk)*al(3)
            if( r1(1) <= r(1) .and. r(1) < r2(1)  .and. &
            &   r1(2) <= r(2) .and. r(2) < r2(2)  .and. &
            &   r1(3) <= r(3) .and. r(3) < r2(3)  ) then
              iatom_frag = iatom_frag + 1
              if(iatom_frag > max_natom_frag) stop "DC method (yn_dc=y): fragment atom count exceeds periodic image buffer"
              rion_frag(1:3,iatom_frag,i_frag) = r(1:3) - dc%rxyz_frag(1:3,i_frag)
              kion_frag(iatom_frag,i_frag) = kion(iatom)
            end if
          end do
          end do
          end do
        end do
        natom_frag(i_frag) = iatom_frag
        i_frag = i_frag + 1
      end do
      end do
      end do
    
    ! set variables for own fragment
      call get_fragment_domain(dc, dc%i_frag, nxyz_domain)
    
    ! nelec (total system) --> nelec (fragment)
      nelec = nelec * natom_frag(dc%i_frag) / natom ! initial guess
    
    ! al, num_rgrid (total system) --> al, num_rgrid (fragment)
      ldomain_frag(1:3) = al(1:3) * dble(nxyz_domain(1:3)) / dble(num_rgrid(1:3))
      al = ldomain_frag + 2d0*lbuffer
      num_rgrid = nxyz_domain + 2*dc%nxyz_buffer
      
    ! natom, rion, kion (total system) --> natom, rion, kion (fragment)
      natom = natom_frag(dc%i_frag)
      deallocate(rion,kion)
      allocate(rion(3,natom),kion(natom))
      rion(1:3,1:natom) = rion_frag(1:3,1:natom,dc%i_frag)
      kion(1:natom) = kion_frag(1:natom,dc%i_frag)
      do i=1,natom
        do n=1,3 ! x,y,z
          rion(n,i) = r_periodic(rion(n,i),al(n))
        end do
      end do
      
    ! dc%jxyz_tot: r-grid (fragment) --> r-grid (total)
      allocate(dc%jxyz_tot(maxval(num_rgrid),3))
      do n=1,3 ! x,y,z
        do i=1,num_rgrid(n) ! r-grid (fragment)
          if(i <= nxyz_domain(n) + dc%nxyz_buffer(n)) then
            j = dc%ixyz_frag(n,dc%i_frag) + i
          else
            j = dc%ixyz_frag(n,dc%i_frag) + ( i - num_rgrid(n) ) ! minus region
          end if
          j = mod(j+dc%lg_tot%num(n)-1,dc%lg_tot%num(n))+1
          dc%jxyz_tot(i,n) = j ! r-grid (total)
        end do
      end do
      
      if(dc%id_frag==0) then
        write(*,'(a,6i10)') "fragment, natom, nelec, ixyz_frag: ",dc%i_frag, natom, nelec, dc%ixyz_frag(1:3,dc%i_frag)
      end if
    
    end subroutine init_fragment

    subroutine build_uniform_fragment_geometry(num_fragment, ldomain)
      implicit none
      integer, intent(in) :: num_fragment(3)
      real(8), intent(in) :: ldomain(3)
      integer :: i_frag, i, j, k

      if (allocated(dc%ixyz_frag)) deallocate(dc%ixyz_frag)
      if (allocated(dc%rxyz_frag)) deallocate(dc%rxyz_frag)
      if (allocated(dc%nxyz_domain_frag)) deallocate(dc%nxyz_domain_frag)
      allocate(dc%ixyz_frag(3,dc%n_frag),dc%rxyz_frag(3,dc%n_frag),dc%nxyz_domain_frag(3,dc%n_frag))

      i_frag = 1
      do i=1,num_fragment(1)
      do j=1,num_fragment(2)
      do k=1,num_fragment(3)
        dc%ixyz_frag(1,i_frag) = (i-1)*dc%nxyz_domain(1)
        dc%ixyz_frag(2,i_frag) = (j-1)*dc%nxyz_domain(2)
        dc%ixyz_frag(3,i_frag) = (k-1)*dc%nxyz_domain(3)
        dc%rxyz_frag(1,i_frag) = dble(i-1)*ldomain(1)
        dc%rxyz_frag(2,i_frag) = dble(j-1)*ldomain(2)
        dc%rxyz_frag(3,i_frag) = dble(k-1)*ldomain(3)
        dc%nxyz_domain_frag(1:3,i_frag) = dc%nxyz_domain(1:3)
        i_frag = i_frag + 1
      end do
      end do
      end do
    end subroutine build_uniform_fragment_geometry

    subroutine build_optimized_fragment_geometry(num_fragment, num_rgrid, al, natom, rion)
      implicit none
      integer, intent(in) :: num_fragment(3), num_rgrid(3), natom
      real(8), intent(in) :: al(3), rion(3, natom)
      integer :: axis, ifrag, ix_frag, iy_frag, iz_frag
      integer :: axis_offset(3, maxval(num_fragment))
      integer :: widths(3)
      integer :: i

      call optimize_fragment_geometry(dc, num_fragment, num_rgrid, al, natom, rion)
      if (allocated(dc%ixyz_frag)) deallocate(dc%ixyz_frag)
      if (allocated(dc%rxyz_frag)) deallocate(dc%rxyz_frag)
      allocate(dc%ixyz_frag(3,dc%n_frag),dc%rxyz_frag(3,dc%n_frag))

      axis_offset(:, :) = 0
      do axis = 1, 3
        do i = 2, num_fragment(axis)
          axis_offset(axis, i) = axis_offset(axis, i - 1) + dc%nxyz_domain_frag(axis, fragment_id_for_axis(axis, i - 1, num_fragment))
        end do
      end do

      ifrag = 0
      do ix_frag = 1, num_fragment(1)
      do iy_frag = 1, num_fragment(2)
      do iz_frag = 1, num_fragment(3)
        ifrag = ifrag + 1
        widths(1:3) = dc%nxyz_domain_frag(1:3, ifrag)
        dc%ixyz_frag(1, ifrag) = axis_offset(1, ix_frag)
        dc%ixyz_frag(2, ifrag) = axis_offset(2, iy_frag)
        dc%ixyz_frag(3, ifrag) = axis_offset(3, iz_frag)
        dc%rxyz_frag(1, ifrag) = al(1) * dble(dc%ixyz_frag(1, ifrag)) / dble(num_rgrid(1))
        dc%rxyz_frag(2, ifrag) = al(2) * dble(dc%ixyz_frag(2, ifrag)) / dble(num_rgrid(2))
        dc%rxyz_frag(3, ifrag) = al(3) * dble(dc%ixyz_frag(3, ifrag)) / dble(num_rgrid(3))
      end do
      end do
      end do
    end subroutine build_optimized_fragment_geometry

    integer function fragment_id_for_axis(axis, iseg, num_fragment) result(ifrag_axis)
      integer, intent(in) :: axis, iseg, num_fragment(3)

      select case(axis)
      case(1)
        ifrag_axis = ((iseg - 1) * num_fragment(2)) * num_fragment(3) + 1
      case(2)
        ifrag_axis = (iseg - 1) * num_fragment(3) + 1
      case default
        ifrag_axis = iseg
      end select
    end function fragment_id_for_axis
    
    function r_periodic(r,a) ! r --> r_periodic in [0,a]
      implicit none
      real(8) :: r_periodic
      real(8),intent(in) :: r,a
      r_periodic = r
      do while (r_periodic < 0d0)
        r_periodic = r_periodic + a
      end do
      do while (r_periodic > a)
        r_periodic = r_periodic - a
      end do
    end function r_periodic

  end subroutine init_dcdft
  
!===================================================================================================================================

  subroutine finalize_dcdft(dc)
    use structures
    implicit none
    type(s_dcdft), intent(inout) :: dc
    integer :: i

    ! DC allocates a second full-system workspace inside main_dft.  Release it
    ! explicitly so compiler/runtime automatic cleanup at subroutine return has
    ! no nested DG/DC buffers left to guess about.
    call deallocate_dft_system(dc%system_tot)
    call finalize_parallel_info(dc%info_tot)
    call deallocate_rgrid(dc%lg_tot)
    call deallocate_rgrid(dc%mg_tot)
    call deallocate_pp_grid(dc%ppg_tot)
    call finalize_reciprocal_grid(dc%fg_tot)
    call finalize_poisson(dc%poisson_tot)
    call finalize_sendrecv_grid_storage(dc%srg_tot)
    call finalize_sendrecv_grid_storage(dc%srg_scalar_tot)

    call deallocate_scalar(dc%vpsl_tot)
    call deallocate_scalar(dc%vh_tot)
    call deallocate_scalar(dc%rho_tot)

    if (allocated(dc%rho_tot_s)) then
      do i = 1, size(dc%rho_tot_s)
        call deallocate_scalar(dc%rho_tot_s(i))
      end do
      deallocate(dc%rho_tot_s)
    end if

    if (allocated(dc%vloc_tot)) then
      do i = 1, size(dc%vloc_tot)
        call deallocate_scalar(dc%vloc_tot(i))
      end do
      deallocate(dc%vloc_tot)
    end if

    if (allocated(dc%vxc_tot)) then
      do i = 1, size(dc%vxc_tot)
        call deallocate_scalar(dc%vxc_tot(i))
      end do
      deallocate(dc%vxc_tot)
    end if

    if (allocated(dc%nxyz_domain_frag)) deallocate(dc%nxyz_domain_frag)
    if (allocated(dc%ixyz_frag)) deallocate(dc%ixyz_frag)
    if (allocated(dc%rxyz_frag)) deallocate(dc%rxyz_frag)
    if (allocated(dc%jxyz_tot)) deallocate(dc%jxyz_tot)

  contains

    subroutine finalize_parallel_info(info)
      type(s_parallel_info), intent(inout) :: info
      if (allocated(info%imap)) deallocate(info%imap)
      if (allocated(info%imap_isolated_ffte)) deallocate(info%imap_isolated_ffte)
      if (allocated(info%ia_mg)) deallocate(info%ia_mg)
      if (allocated(info%irank_io)) deallocate(info%irank_io)
      if (allocated(info%io_s_all)) deallocate(info%io_s_all)
      if (allocated(info%io_e_all)) deallocate(info%io_e_all)
      if (allocated(info%numo_all)) deallocate(info%numo_all)
#ifdef USE_SCALAPACK
      if (allocated(info%ndiv)) deallocate(info%ndiv)
      if (allocated(info%i_tbl)) deallocate(info%i_tbl)
      if (allocated(info%j_tbl)) deallocate(info%j_tbl)
      if (allocated(info%iloc_tbl)) deallocate(info%iloc_tbl)
      if (allocated(info%jloc_tbl)) deallocate(info%jloc_tbl)
#endif
#ifdef USE_FFTW
      if (allocated(info%imap_isolated_fftw)) deallocate(info%imap_isolated_fftw)
#endif
    end subroutine finalize_parallel_info

    subroutine finalize_reciprocal_grid(fg)
      type(s_reciprocal_grid), intent(inout) :: fg
      if (allocated(fg%if_Gzero)) deallocate(fg%if_Gzero)
      if (allocated(fg%vec_G)) deallocate(fg%vec_G)
      if (allocated(fg%coef)) deallocate(fg%coef)
      if (allocated(fg%exp_ewald)) deallocate(fg%exp_ewald)
      if (allocated(fg%egx)) deallocate(fg%egx)
      if (allocated(fg%egxc)) deallocate(fg%egxc)
      if (allocated(fg%egy)) deallocate(fg%egy)
      if (allocated(fg%egyc)) deallocate(fg%egyc)
      if (allocated(fg%egz)) deallocate(fg%egz)
      if (allocated(fg%egzc)) deallocate(fg%egzc)
      if (allocated(fg%coef_nabla)) deallocate(fg%coef_nabla)
      if (allocated(fg%coef_gxgy0)) deallocate(fg%coef_gxgy0)
      if (allocated(fg%cos_cGdt)) deallocate(fg%cos_cGdt)
      if (allocated(fg%sin_cGdt)) deallocate(fg%sin_cGdt)
    end subroutine finalize_reciprocal_grid

    subroutine finalize_poisson(poisson)
      type(s_poisson), intent(inout) :: poisson
      if (allocated(poisson%ipole_tbl)) deallocate(poisson%ipole_tbl)
      if (allocated(poisson%ig_num)) deallocate(poisson%ig_num)
      if (allocated(poisson%ig)) deallocate(poisson%ig)
      if (allocated(poisson%ig_bound)) deallocate(poisson%ig_bound)
      if (allocated(poisson%wkbound)) deallocate(poisson%wkbound)
      if (allocated(poisson%wkbound2)) deallocate(poisson%wkbound2)
      if (allocated(poisson%zrhoG_ele)) deallocate(poisson%zrhoG_ele)
      if (allocated(poisson%ff1x)) deallocate(poisson%ff1x)
      if (allocated(poisson%ff1y)) deallocate(poisson%ff1y)
      if (allocated(poisson%ff1z)) deallocate(poisson%ff1z)
      if (allocated(poisson%ff2x)) deallocate(poisson%ff2x)
      if (allocated(poisson%ff2y)) deallocate(poisson%ff2y)
      if (allocated(poisson%ff2z)) deallocate(poisson%ff2z)
      if (allocated(poisson%ff1)) deallocate(poisson%ff1)
      if (allocated(poisson%ff2)) deallocate(poisson%ff2)
      if (allocated(poisson%ff3x)) deallocate(poisson%ff3x)
      if (allocated(poisson%ff3y)) deallocate(poisson%ff3y)
      if (allocated(poisson%ff3z)) deallocate(poisson%ff3z)
      if (allocated(poisson%ff4x)) deallocate(poisson%ff4x)
      if (allocated(poisson%ff4y)) deallocate(poisson%ff4y)
      if (allocated(poisson%ff4z)) deallocate(poisson%ff4z)
      if (allocated(poisson%dgf)) deallocate(poisson%dgf)
      if (allocated(poisson%a_ffte)) deallocate(poisson%a_ffte)
      if (allocated(poisson%b_ffte)) deallocate(poisson%b_ffte)
#ifdef USE_FFTW
      if (allocated(poisson%fftw1)) deallocate(poisson%fftw1)
      if (allocated(poisson%fftw2)) deallocate(poisson%fftw2)
#endif
    end subroutine finalize_poisson

    subroutine finalize_sendrecv_grid_storage(srg)
      type(s_sendrecv_grid), intent(inout) :: srg
      integer :: idir, iside, itype

      do idir = 1, 3
        do iside = 1, 2
          do itype = 1, 2
            if (allocated(srg%cache(itype, iside, idir)%dbuf)) deallocate(srg%cache(itype, iside, idir)%dbuf)
            if (allocated(srg%cache(itype, iside, idir)%zbuf)) deallocate(srg%cache(itype, iside, idir)%zbuf)
          end do
        end do
      end do
      srg%if_pcomm_real8_initialized = .false.
      srg%if_pcomm_complex8_initialized = .false.
    end subroutine finalize_sendrecv_grid_storage

  end subroutine finalize_dcdft

!===================================================================================================================================
  
  ! rho_s (fragment) --> dc%rho_tot_s (total system)
  subroutine calc_rho_total_dcdft(nspin,lg,mg,info,rho_s,dc)
    use structures
    use communication, only: comm_summation
    implicit none
    integer,              intent(in) :: nspin
    type(s_rgrid),        intent(in) :: lg,mg
    type(s_parallel_info),intent(in) :: info
    type(s_scalar),       intent(in) :: rho_s(nspin)
    type(s_dcdft)                    :: dc
    !
    integer :: ix,iy,iz,ispin,ix_tot,iy_tot,iz_tot
    real(8),dimension(lg%num(1),lg%num(2),lg%num(3),nspin) :: frg_tmp,frg
    real(8),dimension(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),nspin) :: tot_tmp,tot
    
  ! rho_s (fragment)
    frg_tmp = 0d0
    do ispin=1,nspin
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      frg_tmp(ix,iy,iz,ispin) = rho_s(ispin)%f(ix,iy,iz)
    end do
    end do
    end do
    end do
    call comm_summation(frg_tmp,frg,lg%num(1)*lg%num(2)*lg%num(3)*nspin,info%icomm_r)
    
  ! rho_s (total)
    tot_tmp = 0d0
    if(info%id_rko==0) then ! info%id_rko == 0 : representative process of each fragment
      do ispin=1,nspin
      do iz=1,dc%nxyz_domain(3); iz_tot = dc%jxyz_tot(iz,3)
      do iy=1,dc%nxyz_domain(2); iy_tot = dc%jxyz_tot(iy,2)
      do ix=1,dc%nxyz_domain(1); ix_tot = dc%jxyz_tot(ix,1)
        tot_tmp(ix_tot,iy_tot,iz_tot,ispin) = frg(ix,iy,iz,ispin)
      end do
      end do
      end do
      end do
    end if
    call comm_summation(tot_tmp,tot,dc%lg_tot%num(1)*dc%lg_tot%num(2)*dc%lg_tot%num(3)*nspin,dc%icomm_tot)
    do ispin=1,nspin
    do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)
    do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
    do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
      dc%rho_tot_s(ispin)%f(ix,iy,iz) = tot(ix,iy,iz,ispin)
    end do
    end do
    end do
    end do

    if(dc%id_tot==0) then
      write(*,*) "integral(rho_tot)=",sum(tot)*dc%system_tot%hvol," Ne=",dc%elec_num_tot
    end if
    
  end subroutine calc_rho_total_dcdft
  
!===================================================================================================================================
  
  ! dc%vloc_tot (total system) --> v_local (fragment)
  subroutine calc_vlocal_fragment_dcdft(nspin,mg,vloc,dc)
    use structures
    use communication, only: comm_summation
    implicit none
    integer,      intent(in) :: nspin
    type(s_rgrid),intent(in) :: mg
    type(s_scalar)           :: vloc(nspin)
    type(s_dcdft)            :: dc
    !
    integer :: ix,iy,iz,ispin,ix_tot,iy_tot,iz_tot
    real(8),dimension(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3),nspin) :: tot_tmp,tot
    
  ! vloc (total)
    tot_tmp = 0d0
    do ispin=1,nspin
    do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)
    do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
    do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
      tot_tmp(ix,iy,iz,ispin) = dc%vloc_tot(ispin)%f(ix,iy,iz)
    end do
    end do
    end do
    end do
    call comm_summation(tot_tmp,tot,dc%lg_tot%num(1)*dc%lg_tot%num(2)*dc%lg_tot%num(3)*nspin,dc%icomm_tot)
    
  ! vloc (fragment)
    do ispin=1,nspin
    do iz=mg%is(3),mg%ie(3) ; iz_tot = dc%jxyz_tot(iz,3)
    do iy=mg%is(2),mg%ie(2) ; iy_tot = dc%jxyz_tot(iy,2)
    do ix=mg%is(1),mg%ie(1) ; ix_tot = dc%jxyz_tot(ix,1)
      vloc(ispin)%f(ix,iy,iz) = tot(ix_tot,iy_tot,iz_tot,ispin)
    end do
    end do
    end do
    end do
    
  end subroutine calc_vlocal_fragment_dcdft

!===================================================================================================================================

  subroutine capture_dg_dc_seed_payload_dcdft(system,energy,spsi,dc,residual,&
      iteration,payload,ok,message)
    use structures,only:s_dft_system,s_dft_energy,s_orbital,s_dcdft
    use dg_dc_seed_checkpoint,only:s_dg_dc_seed_payload
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    implicit none
    type(s_dft_system),intent(in)::system
    type(s_dft_energy),intent(in)::energy
    type(s_orbital),intent(in)::spsi
    type(s_dcdft),intent(in)::dc
    real(8),intent(in)::residual
    integer,intent(in)::iteration
    type(s_dg_dc_seed_payload),intent(out)::payload
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::allocation_status

    ok=.false.;message=''
    if(system%nspin/=1.or..not.system%if_real_orbital.or..not.allocated(spsi%rwf).or.&
       .not.allocated(dc%rho_tot_s).or..not.allocated(dc%vloc_tot).or.&
       .not.allocated(dc%rho_tot_s(1)%f).or..not.allocated(dc%vloc_tot(1)%f).or.&
       .not.allocated(energy%esp).or..not.allocated(system%rocc))then
      message='DG DC seed capture requires a complete real single-spin state';return
    endif
    if(.not.ieee_is_finite(system%mu).or..not.ieee_is_finite(residual).or.&
       residual<0d0.or.iteration<0)then
      message='DG DC seed capture received invalid scalar provenance';return
    endif
    allocation_status=0
    allocate(payload%rwf(lbound(spsi%rwf,1):ubound(spsi%rwf,1),&
      lbound(spsi%rwf,2):ubound(spsi%rwf,2),lbound(spsi%rwf,3):ubound(spsi%rwf,3),&
      lbound(spsi%rwf,4):ubound(spsi%rwf,4),lbound(spsi%rwf,5):ubound(spsi%rwf,5),&
      lbound(spsi%rwf,6):ubound(spsi%rwf,6),lbound(spsi%rwf,7):ubound(spsi%rwf,7)),&
      payload%rho_tot(lbound(dc%rho_tot_s(1)%f,1):ubound(dc%rho_tot_s(1)%f,1),&
      lbound(dc%rho_tot_s(1)%f,2):ubound(dc%rho_tot_s(1)%f,2),&
      lbound(dc%rho_tot_s(1)%f,3):ubound(dc%rho_tot_s(1)%f,3)),&
      payload%vloc_tot(lbound(dc%vloc_tot(1)%f,1):ubound(dc%vloc_tot(1)%f,1),&
      lbound(dc%vloc_tot(1)%f,2):ubound(dc%vloc_tot(1)%f,2),&
      lbound(dc%vloc_tot(1)%f,3):ubound(dc%vloc_tot(1)%f,3)),&
      payload%esp(lbound(energy%esp,1):ubound(energy%esp,1),&
      lbound(energy%esp,2):ubound(energy%esp,2),lbound(energy%esp,3):ubound(energy%esp,3)),&
      payload%rocc(lbound(system%rocc,1):ubound(system%rocc,1),&
      lbound(system%rocc,2):ubound(system%rocc,2),lbound(system%rocc,3):ubound(system%rocc,3)),&
      stat=allocation_status)
    if(allocation_status/=0)then
      message='cannot allocate DG DC seed capture payload';return
    endif
    payload%rwf=spsi%rwf;payload%rho_tot=dc%rho_tot_s(1)%f
    payload%vloc_tot=dc%vloc_tot(1)%f;payload%esp=energy%esp
    payload%rocc=system%rocc;payload%mu=system%mu
    payload%residual=residual;payload%iteration=iteration
    ok=all(ieee_is_finite(payload%rwf)).and.all(ieee_is_finite(payload%rho_tot)).and.&
      all(ieee_is_finite(payload%vloc_tot)).and.all(ieee_is_finite(payload%esp)).and.&
      all(ieee_is_finite(payload%rocc))
    if(.not.ok)message='DG DC seed capture contains non-finite state'
  end subroutine capture_dg_dc_seed_payload_dcdft

  subroutine rebuild_dg_dc_seed_derived_state_dcdft(mg,info,system,spsi,rho,&
      rho_s,vlocal,dc,ok,message)
    use structures,only:s_rgrid,s_parallel_info,s_dft_system,s_orbital,s_scalar,s_dcdft
    use density_matrix,only:calc_density
    use communication,only:comm_summation
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    implicit none
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_dft_system),intent(in)::system
    type(s_orbital),intent(in)::spsi
    type(s_scalar),intent(inout)::rho,rho_s(system%nspin),vlocal(system%nspin)
    type(s_dcdft),intent(inout)::dc
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::ispin,local_bad,global_bad

    ok=.false.;message='';local_bad=0
    if(system%nspin/=1.or..not.system%if_real_orbital.or..not.allocated(spsi%rwf).or.&
       .not.allocated(dc%rho_tot_s).or..not.allocated(dc%vloc_tot).or.&
       .not.allocated(dc%rho_tot_s(1)%f).or..not.allocated(dc%vloc_tot(1)%f).or.&
       .not.allocated(dc%rho_tot%f).or..not.allocated(rho%f).or.&
       .not.allocated(rho_s(1)%f).or..not.allocated(vlocal(1)%f))local_bad=1
    call comm_summation(local_bad,global_bad,dc%icomm_tot)
    if(global_bad/=0)then
      message='incomplete arrays for restored DG DC seed state';return
    endif

    dc%rho_tot%f=0d0
    do ispin=1,system%nspin
      dc%rho_tot%f=dc%rho_tot%f+dc%rho_tot_s(ispin)%f
    enddo
    call calc_density(system,rho_s,spsi,info,mg)
    rho%f=0d0
    do ispin=1,system%nspin
      rho%f=rho%f+rho_s(ispin)%f
    enddo
    call calc_vlocal_fragment_dcdft(system%nspin,mg,vlocal,dc)
    local_bad=merge(0,1,all(ieee_is_finite(dc%rho_tot%f)).and.&
      all(ieee_is_finite(rho%f)).and.all(ieee_is_finite(rho_s(1)%f)).and.&
      all(ieee_is_finite(vlocal(1)%f)))
    call comm_summation(local_bad,global_bad,dc%icomm_tot)
    if(global_bad/=0)then
      message='non-finite restored DG DC derived state';return
    endif
    ok=.true.
  end subroutine rebuild_dg_dc_seed_derived_state_dcdft
  
!===================================================================================================================================

! cf. src/gs/occupation.f90
  SUBROUTINE ne2mu_dcdft(mg,info,energy,spsi,dc,system)
    use structures
    use communication, only: comm_summation
    use occupation_kernel, only: solve_weighted_state_occupations
    use salmon_global, only: temperature,yn_spinorbit
    implicit none
    type(s_rgrid),        intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_dft_energy),   intent(in) :: energy
    type(s_orbital),      intent(in) :: spsi
    type(s_dcdft),        intent(in) :: dc
    type(s_dft_system)               :: system
    !
    integer :: state_count
    real(8) :: wspin,electron_count
    real(8) :: ne_each(system%no,system%nspin)
    real(8),dimension(system%no,system%nspin,dc%n_frag) :: rocc,esp,ne_frag_orb,wrk1,wrk2
    real(8),allocatable :: state_weights(:),solved_occupations(:)
    logical :: occupation_ok
    character(256) :: occupation_message

    if(system%nspin==1) then
      wspin = 2d0
    else if(system%nspin==2) then
      wspin = 1d0
    end if
    
    call calc_ne_each
    wrk1 = 0d0
    wrk2 = 0d0
    if(info%id_rko==0) then ! info%id_rko == 0 : representative process of each fragment
      wrk1(1:system%no,1:system%nspin,dc%i_frag) = energy%esp(1:system%no,1,1:system%nspin)
      wrk2(1:system%no,1:system%nspin,dc%i_frag) = ne_each(1:system%no,1:system%nspin)
    end if
    call comm_summation(wrk1,esp,        system%no*system%nspin*dc%n_frag,dc%icomm_tot)
    call comm_summation(wrk2,ne_frag_orb,system%no*system%nspin*dc%n_frag,dc%icomm_tot)

    state_count=system%no*system%nspin*dc%n_frag
    allocate(state_weights(state_count))
    state_weights=reshape(ne_frag_orb,[state_count])
    if(yn_spinorbit=='y')state_weights=0.5d0*state_weights
    call solve_weighted_state_occupations(reshape(esp,[state_count]),state_weights,&
      dc%elec_num_tot,max(0d0,temperature),wspin,solved_occupations,system%mu,electron_count,&
      occupation_ok,occupation_message)
    if(.not.occupation_ok)then
      if(dc%id_tot==0)write(0,'(2a)')'DC occupation solve failed: ',trim(occupation_message)
      error stop 'authoritative DC occupation kernel failed'
    endif
    rocc=reshape(solved_occupations,[system%no,system%nspin,dc%n_frag])
    system%rocc(1:system%no,1,1:system%nspin)=rocc(:,:,dc%i_frag)

    return
    
  contains
  
    subroutine calc_ne_each
      implicit none
      integer :: io,ispin,ix,iy,iz
      real(8) :: wrk(system%no,system%nspin)
      
      wrk = 0d0
      if(allocated(spsi%rwf)) then
        do ispin=1,system%nspin
        do io=info%io_s,info%io_e
          do iz=mg%is(3),min(mg%ie(3),dc%nxyz_domain(3)) ! core region only
          do iy=mg%is(2),min(mg%ie(2),dc%nxyz_domain(2)) ! core region only
          do ix=mg%is(1),min(mg%ie(1),dc%nxyz_domain(1)) ! core region only
            wrk(io,ispin) = wrk(io,ispin) + (abs(spsi%rwf(ix,iy,iz,ispin,io,1,1))**2) * system%hvol
          end do
          end do
          end do
        end do
        end do
      else if(allocated(spsi%zwf)) then
        do ispin=1,system%nspin
        do io=info%io_s,info%io_e
          do iz=mg%is(3),min(mg%ie(3),dc%nxyz_domain(3)) ! core region only
          do iy=mg%is(2),min(mg%ie(2),dc%nxyz_domain(2)) ! core region only
          do ix=mg%is(1),min(mg%ie(1),dc%nxyz_domain(1)) ! core region only
            wrk(io,ispin) = wrk(io,ispin) + (abs(spsi%zwf(ix,iy,iz,ispin,io,1,1))**2) * system%hvol
          end do
          end do
          end do
        end do
        end do
      else
        stop "ne2mu_dcdft: neither rwf nor zwf is allocated."
      end if
      call comm_summation(wrk,ne_each,system%no*system%nspin,info%icomm_rko)
      
    end subroutine calc_ne_each

  END SUBROUTINE ne2mu_dcdft
  
!===================================================================================================================================

  subroutine calc_total_energy_dcdft(mg,system,info,v_local,spsi,shpsi,sttpsi,ewald,pp,rion_update,dc,energy)
    use structures
    use communication, only: comm_summation
    use Total_Energy, only: calc_Total_Energy_periodic
    use salmon_global, only: kion !!!!!! future work: remove (kion --> system%kion)
    implicit none
    type(s_rgrid),        intent(in) :: mg
    type(s_dft_system),   intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_scalar)       ,intent(in) :: v_local(system%nspin)
    type(s_orbital),      intent(in) :: spsi,shpsi,sttpsi
    type(s_ewald_ion_ion),intent(in) :: ewald
    type(s_pp_info)      ,intent(in) :: pp
    logical              ,intent(in) :: rion_update
    type(s_dcdft),        intent(in) :: dc
    type(s_dft_energy)               :: energy
    !
    integer :: ispin,io
    integer,dimension(3) :: is,ie
    real(8) :: E_tmp,E_local(2),E_sum(2)
    complex(8) :: ztmp
    
    is(1:3) = mg%is(1:3)
    ie(1:3) = min(mg%ie(1:3),dc%nxyz_domain(1:3)) ! core region only
    
  ! kinetic energy (E_kin)
    E_tmp = 0d0
    if(allocated(spsi%rwf)) then
      do ispin=1,system%Nspin
      do io=info%io_s,info%io_e
        E_tmp = E_tmp + system%rocc(io,1,ispin) &
                    * sum(  spsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1) &
                        * sttpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1) ) * system%Hvol
      end do
      end do
    else if(allocated(spsi%zwf)) then
      do ispin=1,system%Nspin
      do io=info%io_s,info%io_e
        ztmp = sum( conjg(spsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1)) &
            &     * sttpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1) )
        E_tmp = E_tmp + system%rocc(io,1,ispin) * dble(ztmp) * system%Hvol
      end do
      end do
    else
      stop "calc_total_energy_dcdft: neither rwf nor zwf is allocated (E_kin)."
    end if
    E_local(1) = E_tmp

  ! nonlocal part (E_ion_nloc)
    E_tmp = 0d0
    if(allocated(spsi%rwf)) then
      do ispin=1,system%Nspin
      do io=info%io_s,info%io_e
        E_tmp = E_tmp + system%rocc(io,1,ispin) * system%hvol &
          * sum( spsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1) &
             * (shpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1) &
            - (sttpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1) &
       + V_local(ispin)%f(is(1):ie(1),is(2):ie(2),is(3):ie(3)) &
               * spsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1) &
              ) &
            ) &
          )
      end do
      end do
    else if(allocated(spsi%zwf)) then
      do ispin=1,system%Nspin
      do io=info%io_s,info%io_e
        ztmp = sum( conjg(spsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1)) &
            &     * (shpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1) &
            &     - (sttpsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1) &
            &     + V_local(ispin)%f(is(1):ie(1),is(2):ie(2),is(3):ie(3)) &
            &       * spsi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1)) ) )
        E_tmp = E_tmp + system%rocc(io,1,ispin) * dble(ztmp) * system%hvol
      end do
      end do
    else
      stop "calc_total_energy_dcdft: neither rwf nor zwf is allocated (E_ion_nloc)."
    end if
    E_local(2) = E_tmp
    
  ! summation in each fragment
    call comm_summation(E_local,E_sum,2,info%icomm_rko)
    
  ! summation over the total system
    E_local = 0d0
    if(info%id_rko == 0) E_local = E_sum ! info%id_rko == 0 : representative process of each fragment
    call comm_summation(E_local,E_sum,2,dc%icomm_tot)
      
    energy%E_kin = E_sum(1)
    energy%E_ion_nloc = E_sum(2)
    
  ! override (fragment --> total)
    deallocate(kion)
    allocate(kion(dc%system_tot%nion))
    kion = dc%system_tot%kion
    
    call calc_Total_Energy_periodic(dc%mg_tot,ewald,dc%system_tot,dc%info_tot,pp,dc%ppg_tot &
      & ,dc%fg_tot,dc%poisson_tot,rion_update,energy)
      
  ! override (total --> fragment)
    deallocate(kion)
    allocate(kion(system%nion))
    kion = system%kion
    
  end subroutine calc_total_energy_dcdft

!===================================================================================================================================
  
  subroutine write_total_dcdft(system,dc)
    use structures
    use communication, only: comm_summation
    use salmon_global, only: natom, kion, rion, base_directory, yn_out_dns
    use writefield, only: write_dns
    implicit none
    type(s_dcdft) :: dc
    type(s_dft_system),intent(in) :: system
    !
    character(256) :: dir_tmp
    
  ! override (fragment --> total)
    natom = dc%system_tot%nion
    deallocate(kion,rion)
    allocate(kion(natom),rion(3,natom))
    kion = dc%system_tot%kion
    rion = dc%system_tot%rion
    dir_tmp = base_directory
    base_directory = dc%base_directory
    
    if(yn_out_dns =='y') call write_dns(dc%lg_tot,dc%mg_tot,dc%system_tot,dc%info_tot,dc%rho_tot_s)
    
  ! override (total --> fragment)
    natom = system%nion
    deallocate(kion,rion)
    allocate(kion(natom),rion(3,natom))
    kion = system%kion
    rion = system%rion
    base_directory = dir_tmp
    
  end subroutine write_total_dcdft

  subroutine prepare_dg_hybrid_divided_dc_controls(dc,convergence_mode,density_threshold,&
      initial_total_density,ok,message)
    use structures,only:s_dcdft
    use salmon_global,only:convergence,threshold
    use mpi,only:MPI_Allreduce,MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_SUCCESS
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    implicit none
    type(s_dcdft),intent(in)::dc
    character(16),intent(out)::convergence_mode
    real(8),intent(out)::density_threshold
    real(8),allocatable,intent(out)::initial_total_density(:,:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::allocation_status,ix,iy,iz,ierr

    ok=.false.;message='';convergence_mode=convergence;density_threshold=threshold
    if(trim(convergence_mode)/='rho_dne'.and.trim(convergence_mode)/='norm_rho'.and.&
        trim(convergence_mode)/='norm_rho_dng')then
      message='divided Hybrid SCF requires a density convergence quantity';return
    endif
    if(.not.ieee_is_finite(density_threshold).or.density_threshold<=0d0)then
      message='divided Hybrid SCF requires the existing positive DC threshold';return
    endif
    allocate(initial_total_density(dc%lg_tot%num(1),dc%lg_tot%num(2),dc%lg_tot%num(3)),&
      stat=allocation_status)
    if(allocation_status/=0)then
      message='divided Hybrid SCF could not snapshot dc%rho_tot';return
    endif
    initial_total_density=0d0
    do iz=dc%mg_tot%is(3),dc%mg_tot%ie(3)
    do iy=dc%mg_tot%is(2),dc%mg_tot%ie(2)
    do ix=dc%mg_tot%is(1),dc%mg_tot%ie(1)
      initial_total_density(ix,iy,iz)=dc%rho_tot%f(ix,iy,iz)
    enddo
    enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,initial_total_density,size(initial_total_density),&
      MPI_DOUBLE_PRECISION,MPI_SUM,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)then
      message='divided Hybrid SCF could not gather dc%rho_tot';return
    endif
    if(.not.all(ieee_is_finite(initial_total_density)))then
      message='divided Hybrid SCF received a non-finite dc%rho_tot';return
    endif
    ok=.true.
  end subroutine prepare_dg_hybrid_divided_dc_controls

  subroutine load_dg_hybrid_distributed_dc_density(dc,point_ids,values,ok,message)
    use structures,only:s_dcdft
    use mpi
    use,intrinsic::iso_fortran_env,only:int64
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    implicit none
    type(s_dcdft),intent(inout)::dc
    integer(int64),intent(in)::point_ids(:)
    real(8),intent(in)::values(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::rank,nproc,ierr,p,r,owner,local_bad,global_bad,gx,gy,gz
    integer,allocatable::bounds(:,:),send_counts(:),recv_counts(:),send_displs(:),recv_displs(:),cursor(:),&
      point_multiplicity(:,:,:)
    integer(int64),allocatable::send_ids(:),recv_ids(:)
    real(8),allocatable::send_values(:),recv_values(:)

    ok=.false.;message='';local_bad=0
    call MPI_Comm_rank(dc%icomm_tot,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(dc%icomm_tot,nproc,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(size(point_ids)/=size(values).or.any(point_ids<1_int64).or.&
        any(point_ids>product(int(dc%lg_tot%num,int64))).or..not.all(ieee_is_finite(values)))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid distributed DC density payload';return;endif
    allocate(bounds(6,nproc),send_counts(nproc),recv_counts(nproc),send_displs(nproc),recv_displs(nproc),cursor(nproc))
    call MPI_Allgather([dc%mg_tot%is,dc%mg_tot%ie],6,MPI_INTEGER,bounds,6,MPI_INTEGER,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS)then;message='distributed DC density layout exchange failed';return;endif
    send_counts=0
    do p=1,size(point_ids)
      call decode_density_point(point_ids(p),dc%lg_tot%num,gx,gy,gz);owner=-1
      do r=1,nproc
        if(gx>=bounds(1,r).and.gx<=bounds(4,r).and.gy>=bounds(2,r).and.gy<=bounds(5,r).and.&
            gz>=bounds(3,r).and.gz<=bounds(6,r))then
          if(owner/=-1)local_bad=1
          owner=r-1
        endif
      enddo
      if(owner<0)then;local_bad=1;else;send_counts(owner+1)=send_counts(owner+1)+1;endif
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='distributed DC density point has no unique owner';return;endif
    call MPI_Alltoall(send_counts,1,MPI_INTEGER,recv_counts,1,MPI_INTEGER,dc%icomm_tot,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='distributed DC density count exchange failed';return
    endif
    send_displs(1)=0;recv_displs(1)=0
    do r=2,nproc
      send_displs(r)=send_displs(r-1)+send_counts(r-1);recv_displs(r)=recv_displs(r-1)+recv_counts(r-1)
    enddo
    allocate(send_ids(size(point_ids)),send_values(size(values)),recv_ids(sum(recv_counts)),recv_values(sum(recv_counts)))
    cursor=send_displs
    do p=1,size(point_ids)
      call decode_density_point(point_ids(p),dc%lg_tot%num,gx,gy,gz);owner=0
      do r=1,nproc
        if(gx>=bounds(1,r).and.gx<=bounds(4,r).and.gy>=bounds(2,r).and.gy<=bounds(5,r).and.&
            gz>=bounds(3,r).and.gz<=bounds(6,r))owner=r-1
      enddo
      cursor(owner+1)=cursor(owner+1)+1
      send_ids(cursor(owner+1))=point_ids(p);send_values(cursor(owner+1))=values(p)
    enddo
    call MPI_Alltoallv(send_ids,send_counts,send_displs,MPI_INTEGER8,recv_ids,recv_counts,recv_displs,&
      MPI_INTEGER8,dc%icomm_tot,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='distributed DC density ID exchange failed';return;endif
    call MPI_Alltoallv(send_values,send_counts,send_displs,MPI_DOUBLE_PRECISION,&
      recv_values,recv_counts,recv_displs,MPI_DOUBLE_PRECISION,dc%icomm_tot,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='distributed DC density value exchange failed';return;endif
    allocate(point_multiplicity(dc%mg_tot%is(1):dc%mg_tot%ie(1),dc%mg_tot%is(2):dc%mg_tot%ie(2),&
      dc%mg_tot%is(3):dc%mg_tot%ie(3)));point_multiplicity=0
    do p=1,size(recv_ids)
      call decode_density_point(recv_ids(p),dc%lg_tot%num,gx,gy,gz)
      point_multiplicity(gx,gy,gz)=point_multiplicity(gx,gy,gz)+1
    enddo
    local_bad=merge(0,1,.not.any(point_multiplicity/=1))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,dc%icomm_tot,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='distributed DC density catalog is not exactly once';return
    endif
    dc%rho_tot_s(1)%f=0d0
    do p=1,size(recv_ids)
      call decode_density_point(recv_ids(p),dc%lg_tot%num,gx,gy,gz)
      dc%rho_tot_s(1)%f(gx,gy,gz)=recv_values(p)
    enddo
    dc%rho_tot%f=dc%rho_tot_s(1)%f
    ok=.true.
  contains
    pure subroutine decode_density_point(point_id,grid_size,x,y,z)
      integer(int64),intent(in)::point_id
      integer,intent(in)::grid_size(3)
      integer,intent(out)::x,y,z
      x=int(modulo(point_id-1_int64,int(grid_size(1),int64)))+1
      y=int(modulo((point_id-1_int64)/int(grid_size(1),int64),int(grid_size(2),int64)))+1
      z=int((point_id-1_int64)/int(grid_size(1)*grid_size(2),int64))+1
    end subroutine decode_density_point
  end subroutine load_dg_hybrid_distributed_dc_density

end module dcdft
