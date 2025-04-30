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
      type(s_sendrecv_grid) :: srg_dummy
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
      & stencil_dummy,dc%fg_tot,dc%poisson_tot,srg_dummy,dc%srg_scalar_tot,ofile_dummy)
      deallocate(dc%system_tot%rocc)
      call dealloc_cache(srg_dummy)
      
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
      use salmon_global, only: num_rgrid_buffer, kion, rion, natom, num_rgrid, al, num_fragment
      implicit none
      integer :: i_frag,n,i,j,k,ii,jj,kk
      integer :: iatom,iatom_frag
      integer :: kion_frag(natom,dc%n_frag),natom_frag(dc%n_frag)
      real(8) :: dr
      real(8) :: r1(3),r2(3),r(3)
      real(8) :: ldomain(3),lbuffer(3)
      real(8) :: rion_frag(3,natom,dc%n_frag)
    
    ! length of domain
      ldomain(1:3) = al(1:3) / dble(num_fragment(1:3))
      
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

    ! position of each fragment
      allocate(dc%ixyz_frag(3,dc%n_frag),dc%rxyz_frag(3,dc%n_frag))
      i_frag = 1
      do i=1,num_fragment(1)
      do j=1,num_fragment(2)
      do k=1,num_fragment(3)
      ! ix_total = ix_fragment + dc%ixyz_frag(1,dc%i_frag), ix_fragment=[1:dc%nxyz_domain(1)], etc.
        dc%ixyz_frag(1,i_frag) = (i-1)*dc%nxyz_domain(1)
        dc%ixyz_frag(2,i_frag) = (j-1)*dc%nxyz_domain(2)
        dc%ixyz_frag(3,i_frag) = (k-1)*dc%nxyz_domain(3)
        dc%rxyz_frag(1,i_frag) = dble(i-1)*ldomain(1)
        dc%rxyz_frag(2,i_frag) = dble(j-1)*ldomain(2)
        dc%rxyz_frag(3,i_frag) = dble(k-1)*ldomain(3)
        i_frag = i_frag + 1
      end do
      end do
      end do
      
    ! variables for each fragment
      i_frag = 1
      do i=1,num_fragment(1)
      do j=1,num_fragment(2)
      do k=1,num_fragment(3)
      ! boundaries of the fragment i_frag
        r1 = dc%rxyz_frag(:,i_frag) - lbuffer
        r2 = dc%rxyz_frag(:,i_frag) + ldomain + lbuffer
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
    
    ! nelec (total system) --> nelec (fragment)
      nelec = nelec * natom_frag(dc%i_frag) / natom ! initial guess
    
    ! al, num_rgrid (total system) --> al, num_rgrid (fragment)
      al = ldomain + 2d0*lbuffer
      num_rgrid = dc%nxyz_domain + 2*dc%nxyz_buffer
      
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
          if(i <= dc%nxyz_domain(n) + dc%nxyz_buffer(n)) then
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

! cf. src/gs/occupation.f90
  SUBROUTINE ne2mu_dcdft(mg,info,energy,spsi,dc,system)
    use structures
    use communication, only: comm_summation
    implicit none
    type(s_rgrid),        intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_dft_energy),   intent(in) :: energy
    type(s_orbital),      intent(in) :: spsi
    type(s_dcdft),        intent(in) :: dc
    type(s_dft_system)               :: system
    !
    real(8) :: wspin,emax,emin
    real(8) :: ne_each(system%no,system%nspin)
    real(8),dimension(system%no,system%nspin,dc%n_frag) :: rocc,esp,ne_frag_orb,wrk1,wrk2

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

    emin = minval(esp)
    emax = maxval(esp)
    call ne2mu_core(dc%elec_num_tot,emax,emin,system%mu)

    system%rocc(1:system%no,1,1:system%nspin) = rocc(1:system%no,1:system%nspin,dc%i_frag)

    return
    
  contains
  
    subroutine calc_ne_each
      implicit none
      integer :: io,ispin,ix,iy,iz
      real(8) :: wrk(system%no,system%nspin)
      
      wrk = 0d0
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
      call comm_summation(wrk,ne_each,system%no*system%nspin,info%icomm_rko)
      
    end subroutine calc_ne_each

    SUBROUTINE mu2ne(muin,neout)
      use salmon_global, only: temperature, yn_spinorbit
      implicit none
      real(8), intent(in)  :: muin
      real(8), intent(out) :: neout
      !
      integer :: i_frag,ispin,io
      real(8) :: fact
      
      neout=0d0

      if(temperature==0d0) then
        do i_frag=1,dc%n_frag
        do ispin=1,system%nspin
        do io=1,system%no
           fact = esp(io,ispin,i_frag) - muin
           if(fact > 0d0) then
              rocc(io,ispin,i_frag) = 0d0
           else
              rocc(io,ispin,i_frag) = wspin
           endif
           neout = neout + rocc(io,ispin,i_frag) * ne_frag_orb(io,ispin,i_frag)
        end do
        end do
        end do
      else
        do i_frag=1,dc%n_frag
        do ispin=1,system%nspin
        do io=1,system%no
           fact = (esp(io,ispin,i_frag)-muin)/temperature
           if(fact.ge.40.d0) then
              rocc(io,ispin,i_frag) = 0d0
           else
              rocc(io,ispin,i_frag) = wspin/( 1d0 + exp( fact ) )
           endif
           neout = neout + rocc(io,ispin,i_frag) * ne_frag_orb(io,ispin,i_frag)
        end do
        end do
        end do
      end if
      
      if(yn_spinorbit=='y') then
        neout = neout*0.5d0 !!! For the SO mode, ispin=2 components are duplicate copy of jspin=1.
      end if

    END SUBROUTINE mu2ne
    
    subroutine ne2mu_core(nein,emax,emin,muout)
      implicit none
      real(8),intent(in)  :: nein,emax,emin
      real(8),intent(out) :: muout
      !
      integer :: iter,ii,p5,p1,p2,nc
      real(8) :: mu1,mu2,mu3,ne1,ne3,ne3o,diff_ne,diff_mu,muo,diff_ne2,wspin
      
      mu1 = emin
      mu2 = emax
      nc=0

      ITERATION: do iter=1,100000
        diff_ne = 100.d0
        diff_ne2 = 100.d0

        ne1  = 0d0
        ne3  = 0d0
        ne3o = 0d0

        diff_mu =100.d0
        muo=0.0d0

        call mu2ne(mu1,ne1)

        do ii=1,1000
          if ( ii .eq. 1000 ) then
            if ( nc .le. 50)  then
               nc= nc + 1
               mu1 = emin - 0.2d0*dble(nc)
               mu2 = emax + 0.2d0*dble(nc)
               cycle ITERATION
            else
               !print *,'=================================='
               !print *,'Const Ne does not converged!!!!!!!'
               !print *,'=================================='
               exit ITERATION
            endif
          endif
          if ( diff_ne < 1d-10  .and. diff_mu < 1d-9  &
             .and. diff_ne2 < 1.d-9 ) exit ITERATION
          mu3 = mu1 + (mu2-mu1)/2.d0
          ne3=0
          call mu2ne(mu3,ne3)
          diff_ne = abs((ne3-ne3o)/ne3)
          diff_ne2 = abs(nein-ne3)
          if ( (ne1-nein)*(ne3-nein) > 0 ) then
            mu1=mu3
            ne1=ne3
          else
            mu2=mu3
          end if

          ne3o = ne3
          diff_mu = mu3 - muo
          muo  = mu3
        end do
      end do ITERATION

      muout = mu3
      call mu2ne(muout,ne3)
    
    end subroutine ne2mu_core
    
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
    
    is(1:3) = mg%is(1:3)
    ie(1:3) = min(mg%ie(1:3),dc%nxyz_domain(1:3)) ! core region only
    
  ! kinetic energy (E_kin)
    E_tmp = 0d0
    do ispin=1,system%Nspin
    do io=info%io_s,info%io_e
      E_tmp = E_tmp + system%rocc(io,1,ispin) &
                  * sum(  spsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1) &
                      * sttpsi%rwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),ispin,io,1,1) ) * system%Hvol
    end do
    end do
    E_local(1) = E_tmp

  ! nonlocal part (E_ion_nloc)
    E_tmp = 0d0
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
  
end module dcdft
