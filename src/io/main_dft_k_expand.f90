!
!  Copyright 2020 SALMON developers
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
!=======================================================================

#include "config.h"

subroutine main_dft_k_expand
  use salmon_global, only: method_wf_distributor
  implicit none

  if(method_wf_distributor=='single') then
     call main_dft_k_expand_single
  else if(method_wf_distributor=='slice') then
     call main_dft_k_expand_slice
  else
     stop "error in method_wf_distributor option"
  endif
end subroutine main_dft_k_expand

subroutine main_dft_k_expand_single
  use math_constants, only: pi, zi
  use structures
  use inputoutput
  use salmon_global
  use parallelization, only: nproc_id_global,nproc_group_global
  use communication, only: comm_is_root, comm_summation, comm_bcast, comm_sync_all
  use salmon_xc
  use timer
  use scf_iteration_sub
  use density_matrix, only: calc_density
  use writefield
  use salmon_pp, only: calc_nlcc
  use hartree_sub, only: hartree
  use write_sub
  use read_gs
  use code_optimization
  use initialization_sub
  use occupation
  use prep_pp_sub
  use mixing_sub
  use checkpoint_restart_sub
  use filesystem
  use hamiltonian
  use total_energy
  use initialization_dft
  implicit none
  type(s_rgrid) :: lg
  type(s_rgrid) :: mg
  type(s_parallel_info) :: info
  type(s_sendrecv_grid) :: srg, srg_scalar
  type(s_orbital) :: spsi,shpsi,sttpsi
  type(s_dft_system) :: system
  type(s_poisson) :: poisson
  type(s_stencil) :: stencil
  type(s_scalar) :: rho,Vh,Vpsl
  type(s_scalar),allocatable :: V_local(:),rho_s(:),Vxc(:)
  type(s_reciprocal_grid) :: fg
  type(s_pp_info) :: pp
  type(s_pp_grid) :: ppg
  type(s_pp_nlcc) :: ppn
  type(s_dft_energy) :: energy
  type(s_mixing) :: mixing
  type(s_ofile)  :: ofl
  type(s_k_expand) :: kex
  
  integer :: i,ix,iy,iz,ik
  integer :: Miter,nspin
  character(256) :: gdir, wdir0
  character(256),allocatable :: wdir(:)

  call timer_begin(LOG_TOTAL)
  call timer_begin(LOG_INIT_GS)
  
  !change theory input keyword : behave like dft(gs) calc
  theory = 'dft'
  calc_mode='GS'
  
  !check condition
  if(yn_restart /= 'y') stop "error: yn_restart must be y"
  if(process_allocation /= 'orbital_sequential') stop "error: process_allocation must be orbital_sequential"

  call init_dft(nproc_group_global,info,lg,mg,system,stencil,fg,poisson,srg,srg_scalar,ofl)
  allocate( rho_s(system%nspin),V_local(system%nspin),Vxc(system%nspin) )


  call initialization1_dft( system, energy, stencil, fg, poisson,  &
                            lg, mg,  &
                            info,  &
                            srg, srg_scalar,  &
                            rho, rho_s, Vh, V_local, Vpsl, Vxc,  &
                            spsi, shpsi, sttpsi,  &
                            pp, ppg, ppn,  &
                            ofl )

  nspin = system%nspin
  spsi%update_zwf_overlap = .false.
  mixing%num_rho_stock = 21
  call init_mixing(nspin,mg,mixing)  !maybe not necessary


  if(system%nspin /= 1) stop "error: nspin must be 1"
  if(nproc_k /= system%nk) stop "error: nproc_k must be # of k-points"
  if(mod(nstate,nproc_ob)/=0.or.mod(nelec/2,(nstate/nproc_ob))/=0) &
    stop "error: must be mod(nstate,nproc_ob)==0.and.mod(nelec/2,(nstate/nproc_ob))==0"

  ! read restart data
  call restart_gs(lg,mg,system,info,spsi,Miter,mixing=mixing)

  ! initialization for k-expand
  call init_k_expand(system%nk,kex)
  call get_print_rank_numbers(kex,info)

  !(prepare directory)
  call init_dir_out_restart(ofl)
  allocate(wdir(kex%nk))
  call generate_restart_directory_name_k_expand(kex,ofl%dir_out_restart,gdir,wdir)
  do ik=1,kex%nk
    !write(*,*) trim(wdir(ik))
     call create_directory(wdir(ik))
  enddo

  do ik=1,kex%nk
     if(kex%myrank(ik)==0) then
        wdir0 = wdir(ik)
        call write_info_bin(     wdir0,kex,system,Miter)
        call write_atomic_coor(  wdir0,kex,system)
        call write_atomic_vel(   wdir0,kex,system)
        call write_occ_bin(      wdir0,kex,system)
     endif
  enddo
  call write_rho_inout_bin(wdir0,kex,system,lg,mg,info,mixing)
  call write_wfn_bin(wdir,kex,system,spsi,lg,mg,info,mixing)


  !(log)
  if(comm_is_root(nproc_id_global)) then
     write(*,*)
     write(*,'(a)')         "  New parameters after expanding:"
     write(*,'(a, i6)')     "    natom =", kex%natom
     write(*,'(a, i6)')     "    nelec =", kex%nelec
     write(*,'(a, i6)')     "    nstate=", kex%nstate
     write(*,'(a,3i6)')     "    num_rgrid=", kex%num_rgrid(1:3)
     write(*,'(a,3f22.12)') "    al[A]    =", kex%al(1:3)*au_length_aa

     write(*,'(a)')         "    process_allocation = orbital_sequential"
     write(*,'(a, i6)')     "    nproc_k  =" , kex%nk_new
     write(*,'(a, i6)')     "    nproc_ob =" , nproc_ob * kex%nk
     write(*,'(a,3i6)')     "    nproc_rgrid =", nproc_rgrid(1)*kex%nkx, nproc_rgrid(1)*kex%nky, nproc_rgrid(2)*kex%nkz
  endif


  call timer_end(LOG_TOTAL)

contains

subroutine init_k_expand(mk,kex)
  use parallelization, only: nproc_id_global, nproc_group_global
  use communication, only: comm_is_root, comm_bcast
  use salmon_global, only: al,natom,nelec,nstate,num_rgrid
  use inputoutput, only: au_length_aa
  implicit none
  type(s_k_expand) :: kex
  integer :: nkxyz(3), ix,iy,iz,ixyz, iu1,iu2, mk, i,ik
  real(8) :: primitive_b(3)
  character(256) :: ifile_k_data

  kex%nk = mk
  allocate( kex%isupercell(3,kex%nk), kex%k_vec(3,kex%nk) )

  if(comm_is_root(nproc_id_global)) then
     write(*,*)"  Convert wavefunction from GS/k-points into supercell/gamma-point"
     write(*,*)"  assuming :"
     write(*,*)"  (1) # of k-points in GS is over 2 and that in RT is 1 (gamma only)"
     write(*,*)"  (2) k-points in GS is Monkhorst-pac (using user-defined k-points)"
     write(*,*)"  (3) spacial grid size does not change"
     write(*,*)"  (4) # of k-points in GS and # of blocks in super-cell in RT is the same"
     write(*,*)"  (5) cuboid cell"
     write(*,*)"  (6) process_allocation = orbital_sequential"
     write(*,*)"  (7) nproc_rgrid(1,2,3) --> x Nkx, x Nky, x Nkz"
     write(*,*)"  (8) nproc_k = nk in the input data file"
     write(*,*)"  (9) mod(nstate,nproc_ob)==0 and mod(nelec/2,(nstate/nproc_ob))==0"
     write(*,*)"  (10)nspin=1, periodic system, ...."
     write(*,*)

     !(read information from "info_k_expand.dat")
     iu1=91
     open(iu1,file="info_k_expand.dat",status="old")
     read(iu1,'(a)') ifile_k_data  !file name of k-point information printed by SALMON (xxx_k.data)
     read(iu1,*) kex%nkx, kex%nky, kex%nkz   ! super cell size ; must correspond to k-point

     nkxyz(1) = kex%nkx
     nkxyz(2) = kex%nky
     nkxyz(3) = kex%nkz
     if(kex%nk .ne. kex%nkx * kex%nky * kex%nkz) stop "Error the super cell size or k-points in gs"

     ixyz=0
     do ix=1,nkxyz(1)
     do iy=1,nkxyz(2)
     do iz=1,nkxyz(3)
        ixyz=ixyz+1
        kex%isupercell(1,ixyz)=ix
        kex%isupercell(2,ixyz)=iy
        kex%isupercell(3,ixyz)=iz
     enddo
     enddo
     enddo
     close(iu1)

     !(read information of k-point in GS)
     iu2=92
     open(iu2,file=trim(ifile_k_data),status="old")
     do i=1,5
        read(iu2,*)
     enddo
     do ik=1,kex%nk
        read(iu2,*) i, kex%k_vec(1:3,ik)  ! k=[-0.5:0.5], weight is assumed to 1
     enddo
     read(iu2,*)
     read(iu2,*) primitive_b(1:3)   !=2*pi/a

     write(*,*) "  k_vec [a.u.]"
     do ik=1,kex%nk
        kex%k_vec(:,ik) = kex%k_vec(:,ik) * primitive_b(:)
        write(*,*) ik, real(kex%k_vec(:,ik))
     enddo

     close(iu2)

  endif

  call comm_bcast(kex%nkx,       nproc_group_global)
  call comm_bcast(kex%nky,       nproc_group_global)
  call comm_bcast(kex%nkz,       nproc_group_global)
  call comm_bcast(kex%isupercell,nproc_group_global)
  call comm_bcast(kex%k_vec,     nproc_group_global)

  kex%natom    = natom * kex%nk
  kex%nelec    = nelec * kex%nk
  kex%nstate   = nstate* kex%nk
  kex%num_rgrid(:) = num_rgrid(:)* nkxyz(:)
  kex%al(:)    = al(:) * nkxyz(:)

end subroutine init_k_expand

subroutine get_print_rank_numbers(kex,info)
  use communication, only: comm_is_root, comm_bcast
  implicit none
  type(s_k_expand) :: kex
  type(s_parallel_info),intent(in) :: info
  integer :: icnt,nl,i1,i2,i3,i4,i5,k1,k2,k3,nkx,nky,nkz
  integer :: nproc_k, nproc_o, nproc_r(3)

   !!(before)
   !! (in subroutine init_communicator_dft in init_communicator.f90)
   ! nl=-1
   ! do i3=0,nproc_d_o(3)-1
   ! do i2=0,nproc_d_o(2)-1
   ! do i1=0,nproc_d_o(1)-1
   ! do i5=0,nproc_k-1
   ! do i4=0,nproc_ob-1
   !   nl = nl + 1
   !   info%imap(i1,i2,i3,i4,i5) = nl
   !   if (nl == myrank) then
   !     info%iaddress = [i1,i2,i3,i4,i5]
   !   end if
   ! end do
   ! end do
   ! end do
   ! end do
   ! end do

    nkx = kex%nkx
    nky = kex%nky
    nkz = kex%nkz

    nproc_k  = info%npk
    nproc_o  = info%nporbital
    nproc_r  = info%nprgrid

    !(after)
    allocate( kex%myrank(kex%nk) )
    allocate( kex%iaddress(8,kex%nk) )

    icnt=0
    nl=-1
    do k3= 1,nkz
    do i3= 0,nproc_r(3)-1
    do k2= 1,nky
    do i2= 0,nproc_r(2)-1
    do k1= 1,nkx
    do i1= 0,nproc_r(1)-1
    do i4= 0,nproc_o-1
    do i5= 0,nproc_k-1
       nl= nl + 1
       if ( info%iaddress(1)==i1 .and. &
            info%iaddress(2)==i2 .and. &
            info%iaddress(3)==i3 .and. &
            info%iaddress(4)==i4 .and. &
            info%iaddress(5)==i5  ) then
            icnt = icnt + 1
            kex%myrank(icnt) = nl
            kex%iaddress(1,icnt) = i1
            kex%iaddress(2,icnt) = k1
            kex%iaddress(3,icnt) = i2
            kex%iaddress(4,icnt) = k2
            kex%iaddress(5,icnt) = i3
            kex%iaddress(6,icnt) = k3
            kex%iaddress(7,icnt) = i4
            kex%iaddress(8,icnt) = i5
       endif
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do

   !write(*,'(a,100i5)') "myrank  =", kex%myrank(:)

  end subroutine get_print_rank_numbers

  subroutine generate_restart_directory_name_k_expand(kex,basedir,gdir,pdir)
    implicit none
    type(s_k_expand) :: kex
    character(*),  intent(in)  :: basedir
    character(256),intent(out) :: gdir
    character(256),intent(out) :: pdir(kex%nk)
    integer :: ik, nproc_id_kex

    ! global directory
    write(gdir,'(A,I6.6,A)')   trim(basedir)

    ! process private directory for k-expand
    do ik=1,kex%nk
       nproc_id_kex = kex%myrank(ik)
       write(pdir(ik),'(A,A,I6.6,A)') trim(gdir),'rank_',nproc_id_kex,'/'
    enddo

  end subroutine generate_restart_directory_name_k_expand


  subroutine write_info_bin(odir,kex,system,iter)
    use parallelization, only: nproc_size_global
    implicit none
    type(s_k_expand) :: kex
    type(s_dft_system) :: system
    character(*)    :: odir
    character(1024) :: dir_file_out
    logical :: if_real_orbital_tmp
    integer :: iu3, iter, nprocs

    iu3 = 93
    dir_file_out = trim(odir)//"/info.bin"
    open(iu3,file=dir_file_out,form='unformatted')

    kex%nk_new = 1
    kex%no_new = system%no * kex%nk
    ! should comm  (kex%nk_new, kex%no_new)
    !iter       = 0
    nprocs     = nproc_size_global * kex%nk
    if_real_orbital_tmp = .false.
  
    write(iu3) kex%nk_new
    write(iu3) kex%no_new
    write(iu3) iter
    write(iu3) nprocs
    write(iu3) if_real_orbital_tmp

    write(*,*) "info_bin:", kex%nk_new, kex%no_new, iter, nprocs
    
    close(iu3)

  end subroutine write_info_bin

  subroutine write_atomic_coor(odir,kex,system)
    use salmon_global, only: al,natom,atom_name,kion,unit_length
    use inputoutput, only: au_length_aa
    implicit none
    type(s_k_expand) :: kex
    type(s_dft_system) :: system
    character(*)    :: odir
    character(1024) :: dir_file_out
    integer :: iu1, ia, ixyz
    real(8) :: uconv, Rion_new(3)

    if(unit_length=='AA')then ; uconv = au_length_aa
    else                      ; uconv = 1d0   !au
    endif

    iu1 = 91
    dir_file_out = trim(odir)//"/atomic_coor.txt"
    open(iu1,file=dir_file_out,status="unknown")

    do ia = 1,natom
       atom_name(ia) = adjustl(atom_name(ia))
       atom_name(ia) = "'"//trim(atom_name(ia))//"'"
    enddo
    do ixyz= 1,kex%nk
       do ia = 1,natom
          Rion_new(:) = system%Rion(:,ia) + al(:) * (kex%isupercell(:,ixyz)-1)
          write(iu1,7000) trim(atom_name(ia)), Rion_new(1:3)*uconv, kion(ia)
       enddo
    enddo
7000 format(" ",a6,"  ",3f18.10,i4)
    close(iu1)

  end subroutine write_atomic_coor

  subroutine write_atomic_vel(odir,kex,system)
    use salmon_global, only: natom
    implicit none
    type(s_k_expand) :: kex
    type(s_dft_system) :: system
    character(*)    :: odir
    character(1024) :: dir_file_out
    integer :: iu1, ia, ixyz

    iu1 = 91
    dir_file_out = trim(odir)//"/atomic_vel.txt"
    open(iu1,file=dir_file_out,status="unknown")

    do ixyz= 1,kex%nk
       do ia = 1,natom
          write(iu1,7100) system%Velocity(1:3,ia)
       enddo
    enddo
7100 format(3f18.10)
    close(iu1)

  end subroutine write_atomic_vel

  subroutine write_occ_bin(odir,kex,system)
    implicit none
    type(s_k_expand) :: kex
    type(s_dft_system) :: system
    character(*)    :: odir
    character(1024) :: dir_file_out
    integer :: iu1, io,ik, io_new
    real(8),allocatable :: rocc_new(:,:,:)

    iu1 = 91
    dir_file_out = trim(odir)//"/occupation.bin"
    open(iu1,file=dir_file_out,form='unformatted')

    allocate(rocc_new(1:kex%no_new, 1:kex%nk_new, 1:system%nspin))

    rocc_new(:,:,:) = 0d0
    io_new=0
    do io=1,system%no
    do ik=1,kex%nk
       io_new = io_new + 1
       rocc_new(io_new,1,nspin) = system%rocc(io,ik,nspin)
      !write(*,*) "rocc", io_new, rocc_new(io_new,1,nspin)
    enddo
    enddo
    write(iu1) rocc_new(1:kex%no_new,1:kex%nk_new,1:system%nspin)
    close(iu1)

  deallocate(rocc_new)

end subroutine write_occ_bin

subroutine write_rho_inout_bin(odir,kex,system,lg,mg,info,mixing)
  implicit none
  type(s_k_expand) :: kex
  type(s_dft_system) :: system
  type(s_rgrid) :: lg,mg
  type(s_parallel_info) :: info
  type(s_mixing) :: mixing
  character(*)    :: odir
  character(1024) :: dir_file_out  !,ofile_cube
  integer :: iu1,ik, num_rho_stock
  integer :: lg_is_new(3), lg_ie_new(3), ixyz, nval
  integer :: ix_add, iy_add, iz_add,iix,iiy,iiz
  real(8) :: norm
  real(8),allocatable :: rho_in_new(:,:,:,:), rho_out_new(:,:,:,:)
  real(8),allocatable :: tmp1_rho_in_new(:,:,:), tmp2_rho_in_new(:,:,:)
  real(8),allocatable :: tmp1_rho_out_new(:,:,:),tmp2_rho_out_new(:,:,:)

  lg_is_new(:) = lg%is(:)
  lg_ie_new(1) = lg%ie(1) * kex%nkx
  lg_ie_new(2) = lg%ie(2) * kex%nky
  lg_ie_new(3) = lg%ie(3) * kex%nkz

  nval = lg_ie_new(1)*lg_ie_new(2)*lg_ie_new(3)
  num_rho_stock = mixing%num_rho_stock

  allocate( rho_in_new( lg_ie_new(1),lg_ie_new(2),lg_ie_new(3),num_rho_stock+1))
  allocate( rho_out_new(lg_ie_new(1),lg_ie_new(2),lg_ie_new(3),num_rho_stock)  )

  !$omp parallel do collapse(3) private(i,ix,iy,iz)
  do i=1,num_rho_stock
  do ix=lg_is_new(1),lg_ie_new(1)
  do iy=lg_is_new(2),lg_ie_new(2)
  do iz=lg_is_new(3),lg_ie_new(3)
     rho_in_new( ix,iy,iz,i) = 0d0
     rho_out_new(ix,iy,iz,i) = 0d0
  enddo
  enddo
  enddo
  enddo

  allocate( tmp1_rho_in_new(lg_ie_new(1),lg_ie_new(2),lg_ie_new(3)) )
  allocate( tmp2_rho_in_new(lg_ie_new(1),lg_ie_new(2),lg_ie_new(3)) )
  allocate( tmp1_rho_out_new(lg_ie_new(1),lg_ie_new(2),lg_ie_new(3)) )
  allocate( tmp2_rho_out_new(lg_ie_new(1),lg_ie_new(2),lg_ie_new(3)) )

  tmp1_rho_in_new  = 0d0
  tmp2_rho_in_new  = 0d0
  tmp1_rho_out_new = 0d0
  tmp2_rho_out_new = 0d0


  do i=1,num_rho_stock+1
     do ixyz=1,kex%nk

        ix_add = lg%ie(1) * (kex%isupercell(1,ixyz)-1)
        iy_add = lg%ie(2) * (kex%isupercell(2,ixyz)-1)
        iz_add = lg%ie(3) * (kex%isupercell(3,ixyz)-1)

        !$omp parallel do collapse(2) private(ix,iy,iz,iix,iiy,iiz)
        do iz= mg%is(3),mg%ie(3)
        do iy= mg%is(2),mg%ie(2)
        do ix= mg%is(1),mg%ie(1)

           iix = ix + ix_add
           iiy = iy + iy_add
           iiz = iz + iz_add
           tmp1_rho_in_new(iix,iiy,iiz)  = mixing%rho_in(i)%f(ix,iy,iz)

           if(i.le.num_rho_stock) &
           tmp1_rho_out_new(iix,iiy,iiz) = mixing%rho_out(i)%f(ix,iy,iz)

        enddo
        enddo
        enddo

     enddo

     call comm_summation(tmp1_rho_in_new, tmp2_rho_in_new,  nval,info%icomm_r)

     norm = sum( tmp2_rho_in_new(:,:,:) )*system%hvol
     write(*,*) "norm of rho in density", i,real(norm)
     !tmp2_rho_in_new(:,:,:) = tmp2_rho_in_new(:,:,:)/norm

     rho_in_new(1:lg_ie_new(1),1:lg_ie_new(2),1:lg_ie_new(3),i) = &
            tmp2_rho_in_new(1:lg_ie_new(1),1:lg_ie_new(2),1:lg_ie_new(3))


     if(i.le.num_rho_stock) then
        call comm_summation(tmp1_rho_out_new,tmp2_rho_out_new, nval,info%icomm_r)

        norm = sum( tmp2_rho_out_new(:,:,:) )*system%hvol
        write(*,*) "norm of rho out density", i, real(norm)

        rho_out_new(1:lg_ie_new(1),1:lg_ie_new(2),1:lg_ie_new(3),i) = &
             tmp2_rho_out_new(1:lg_ie_new(1),1:lg_ie_new(2),1:lg_ie_new(3))
     endif


     !for check
     !ofile_cube="rho.cube"
     !call  write_cube(tmp2_rho_in_new,lg_ie_new(1),lg_ie_new(2),lg_ie_new(3),system%hgs,ofile_cube)

  enddo ! i


  !only print process only from here
  do ik=1,kex%nk
     if(kex%myrank(ik)==0) then

        iu1 = 91
        dir_file_out = trim(odir)//"/rho_inout.bin"
        open(iu1,file=dir_file_out,form='unformatted')
        do i=1,num_rho_stock+1
           write(iu1) rho_in_new(1:lg_ie_new(1),1:lg_ie_new(2),1:lg_ie_new(3),i)
        end do
        do i=1,num_rho_stock
           write(iu1) rho_out_new(1:lg_ie_new(1),1:lg_ie_new(2),1:lg_ie_new(3),i)
        end do

        close(iu1)
     endif
  enddo

  deallocate(rho_in_new,rho_out_new, tmp1_rho_in_new,tmp2_rho_in_new)

end subroutine write_rho_inout_bin

subroutine write_wfn_bin(odir,kex,system,spsi,lg,mg,info,mixing)
  implicit none
  type(s_k_expand) :: kex
  type(s_dft_system) :: system
  type(s_orbital) :: spsi
  type(s_rgrid) :: lg,mg
  type(s_parallel_info) :: info
  type(s_mixing)  :: mixing
  character(*)    :: odir(kex%nk)
  character(1024) :: dir_file_out    !,ofile_cube
 !real(8) :: tmp_norm(system%nk,system%no), norm_all(system%nk,system%no)
  real(8) :: scale, rshift(3), r(3)  !,norm
  integer :: ix_add, iy_add, iz_add  !,mx,my,mz
  integer :: iu,id,ip_x,ip_y,ip_z, inkx,inky,inkz, ip_o,ip_k
  integer :: im,ik,is,ix,iy,iz,iix,iiy,iiz !, io
 !real(8),allocatable :: orb_real(:,:,:)
  complex(8) :: ai,ekr
  complex(8),allocatable :: zwf_ekr(:,:,:,:)

  !!check norm -> later normalize to improbe accuracy
  !im = 1
  !is = 1 
  !tmp_norm(system%nk,system%no) = 0d0
  !
  !do ik= info%ik_s,info%ik_e
  !do io= info%io_s,info%io_e
  !    tmp_norm(ik,io) = sum(abs(spsi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),is,io,ik,im))**2)
  !enddo
  !enddo
  !
  !call comm_summation(tmp_norm,norm_all,system%nk*system%no,info%icomm_r)
  !
  !do ik= info%ik_s,info%ik_e
  !do io= info%io_s,info%io_e
  !   norm_all(ik,io) = norm_all(ik,io) * system%hvol
  !   scale = 1d0/sqrt(norm_all(ik,io))
  !   write(*,'(a,2i6,f20.12)') "norm of orbital   ", ik,io, norm_all(ik,io)
  !   spsi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),is,io,ik,im)= &
  !   spsi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),is,io,ik,im)*scale
  !enddo
  !enddo
!--------------------

  ai        = ( 0d0, 1d0 )
  rshift(:) = -system%hgs(:)    !probably only if(yn_domain_parallel=='y')
  scale = 1d0/sqrt(dble(kex%nk))

  allocate( zwf_ekr(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),info%io_s:info%io_e) )

  do id =1,kex%nk

     ![i1,k1,i2,k2,i3,k3,i4,i5]
     ip_x = kex%iaddress(1,id)
     inkx = kex%iaddress(2,id)
     ip_y = kex%iaddress(3,id)
     inky = kex%iaddress(4,id)
     ip_z = kex%iaddress(5,id)
     inkz = kex%iaddress(6,id)
     ip_o = kex%iaddress(7,id)
     ip_k = kex%iaddress(8,id)

     !now assuming ik_s=ik_e, im_s=im_e, nspin=1
     if(info%im_s /= info%im_e) stop "error im_s /= im_e"
     if(info%ik_s /= info%ik_e) stop "error ik_s /= ik_e"

     im = info%im_s
     ik = info%ik_s
     is = 1 

     iu = 1000 + kex%myrank(id)
     dir_file_out = trim(odir(id))//"/wfn.bin"
     open(iu,file=dir_file_out,form='unformatted',access='stream')

     ix_add = (lg%ie(1)-lg%is(1)+1) * (inkx-1)
     iy_add = (lg%ie(2)-lg%is(2)+1) * (inky-1)
     iz_add = (lg%ie(3)-lg%is(3)+1) * (inkz-1)

     do iz = mg%is(3), mg%ie(3)
        iiz = iz + iz_add
        r(3) = iiz*system%hgs(3) + rshift(3)

        do iy = mg%is(2), mg%ie(2)
           iiy = iy + iy_add
           r(2) = iiy*system%hgs(2) + rshift(2)

           do ix = mg%is(1), mg%ie(1)
              iix = ix + ix_add
              r(1) = iix*system%hgs(1) + rshift(1)

              ekr = exp(ai * sum(kex%k_vec(:,ik)*r(:)) ) * scale
              zwf_ekr(ix,iy,iz,info%io_s:info%io_e) = &
                     spsi%zwf(ix,iy,iz,is,info%io_s:info%io_e,ik,im)* ekr

           enddo
        enddo
     enddo

     !write
     write(iu) zwf_ekr(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),info%io_s:info%io_e)

    ! !for check
    ! if(ip_o==20) then
    !    mx=mg%ie(1)
    !    my=mg%ie(2)
    !    mz=mg%ie(3)
    !    allocate(orb_real(mx,my,mz))
    !    orb_real(1:mx,1:my,1:mz) = real(zwf_ekr(1:mx,1:my,1:mz,ip_o))
    !    ofile_cube="orb_20.cube"
    !    call  write_cube(orb_real,mx,my,mz,system%hgs,ofile_cube)
    ! endif

     close(iu)
  enddo
  if(allocated(zwf_ekr)) deallocate( zwf_ekr )

end subroutine write_wfn_bin

subroutine  write_cube(orb,mx,my,mz,hgs,ofl)
  implicit none
  integer :: mx,my,mz,fp,i,j,n,ix,iy,iz
  real(8) :: orb(mx,my,mz), hgs(3), crd(3,2)
  character(*) :: ofl

  n=2
  crd(:,1)=0d0
  crd(1,2)=mx*hgs(1)
  crd(2,2)=my*hgs(2)
  crd(3,2)=mz*hgs(3)

  fp=1
  open(fp,file=trim(ofl),status="unknown")
  write(fp,*) "distribution"
  write(fp,*) "All values here are in a.u."
  write(fp,'(i5,3f12.6)') n, 0d0, 0d0, 0d0
  write(fp,'(i5,3f12.6)') mx,hgs(1),0.d0,0.d0
  write(fp,'(i5,3f12.6)') my,0.d0,hgs(2),0.d0
  write(fp,'(i5,3f12.6)') mz,0.d0,0.d0,hgs(3)
  do i=1,n
     !ik=Kion(iatom)
     !write(fp,'(i5,4f12.6)') izatom(ik),dble(izatom(ik)),(Rion(j,iatom),j=1,3)
     write(fp,'(i5,4f12.6)') 8,dble(8),(crd(j,i),j=1,3)
  end do
  
  do ix=1,mx
  do iy=1,my
     write(fp,'(6(1X,E23.15E3))', advance="yes") (orb(ix,iy,iz),iz=1,mz)
  end do
  end do
  close(fp)

  return
end subroutine write_cube

end subroutine main_dft_k_expand_single
