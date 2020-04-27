!
!  Copyright 2020-2020 SALMON developers
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

subroutine main_dft_k_expand_slice
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
  use input_pp_sub
  use prep_pp_sub
  use mixing_sub
  use checkpoint_restart_sub
  use filesystem
  use hamiltonian
  use salmon_total_energy
  use initialization_dft
  implicit none
  type(s_rgrid) :: lg
  type(s_rgrid) :: mg
  type(s_parallel_info) :: info
  type(s_sendrecv_grid) :: srg, srg_ng
  type(s_orbital) :: spsi,shpsi,sttpsi
  type(s_dft_system) :: system
  type(s_poisson) :: poisson
  type(s_stencil) :: stencil
  type(s_scalar) :: srho,sVh,sVpsl
  type(s_scalar),allocatable :: V_local(:),srho_s(:),sVxc(:)
  type(s_reciprocal_grid) :: fg
  type(s_pp_info) :: pp
  type(s_pp_grid) :: ppg
  type(s_pp_nlcc) :: ppn
  type(s_dft_energy) :: energy
  type(s_mixing) :: mixing
  type(s_ofile)  :: ofl
  type(s_k_expand) :: kex
  
  integer :: i,ix,iy,iz,ndir_r, ndir_w
  integer :: nspin, nblock_orbital
  character(256) :: rgdir,wgdir, rdir0,wdir0
  character(256),allocatable :: wdir(:), rdir(:)

  call timer_begin(LOG_TOTAL)
  call timer_begin(LOG_INIT_GS)
  
  !change theory input keyword : behave like dft(gs) calc
  theory = 'dft'
  calc_mode='GS'
  
  !check condition
  if(yn_restart /= 'y') stop "error: yn_restart must be y"
  if(method_wf_distributor /= 'slice') stop "error: method_wf_distributor must be slice"

  call init_dft(nproc_group_global,info,lg,mg,system,stencil,fg,poisson,srg,srg_ng,ofl)
  allocate( srho_s(system%nspin),V_local(system%nspin),sVxc(system%nspin) )


  call initialization1_dft( system, energy, stencil, fg, poisson,  &
                            lg, mg,  &
                            info,  &
                            srg, srg_ng,  &
                            srho, srho_s, sVh, V_local, sVpsl, sVxc,  &
                            spsi, shpsi, sttpsi,  &
                            pp, ppg, ppn,  &
                            ofl )

  nspin = system%nspin
  spsi%update_zwf_overlap = .false.
  mixing%num_rho_stock = 21
  call init_mixing(nspin,mg,mixing)  !maybe not necessary


  if(system%nspin /= 1) stop "error: nspin must be 1"
  if(nproc_k /= system%nk) stop "error: nproc_k must be # of k-points"


  ! initialization for k-expand
  call init_k_expand(system%nk,kex)
  call assign_rank_mumber_to_read_write_files(kex,system,info)

  nblock_orbital = min(system%no, nblock_wf_distribute)
  if(mod(system%no,nblock_orbital)==0) then
     ndir_r =  int(system%no/nblock_orbital) * system%nk
  else
     ndir_r = (int(system%no/nblock_orbital)+1) * system%nk
  endif
  allocate(rdir(ndir_r))
  call get_restart_reading_directory_name(system,nblock_orbital,rgdir,rdir)


  !(prepare directory)
  call init_dir_out_restart(ofl)
  if(mod(kex%no_new,nblock_orbital)==0) then
     ndir_w = int(kex%no_new/nblock_orbital)
  else
     ndir_w = int(kex%no_new/nblock_orbital)+1
  endif
  allocate(wdir(ndir_w))
  call gen_restart_writing_directory_name_k_expand(kex,nblock_orbital,ofl%dir_out_restart,wgdir,wdir)

  if(comm_is_root(nproc_id_global)) then
     rdir0 = rgdir
     wdir0 = wgdir
     call read_write_info_bin(     rdir0,wdir0,kex,system)
     call read_write_atomic_coor(  rdir0,wdir0,kex,system)
     call read_write_atomic_vel(   rdir0,wdir0,kex,system)
     call read_write_occ_bin(      rdir0,wdir0,kex,system)
     call read_write_rho_inout_bin(rdir0,wdir0,kex,system,lg,mixing)
  endif
  call read_write_wfn_bin(rdir,wdir,kex,nblock_orbital,system,lg)


  !(log)
  if(comm_is_root(nproc_id_global)) then
     write(*,*)
     write(*,'(a)')         "  New parameters after expanding:"
     write(*,'(a, i6)')     "    natom =", kex%natom
     write(*,'(a, i6)')     "    nelec =", kex%nelec
     write(*,'(a, i6)')     "    nstate=", kex%nstate
     write(*,'(a,3i6)')     "    num_rgrid=", kex%num_rgrid(1:3)
     write(*,'(a,3f22.12)') "    al[A]    =", kex%al(1:3)*au_length_aa
    !write(*,'(a)')         "    process_allocation = orbital_sequential"
    !write(*,'(a, i6)')     "    nproc_k  =" , kex%nk_new
    !write(*,'(a, i6)')     "    nproc_ob =" , nproc_ob * kex%nk
    !write(*,'(a,3i6)')     "    nproc_rgrid =", nproc_rgrid(1)*kex%nkx, nproc_rgrid(1)*kex%nky, nproc_rgrid(2)*kex%nkz
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
     write(*,*)"  (0) method_wf_distributor = slice"
     write(*,*)"  (1) # of k-points in GS is over 2 and that in RT is 1 (gamma only)"
     write(*,*)"  (2) k-points in GS is Monkhorst-pac (using user-defined k-points)"
     write(*,*)"  (3) spacial grid size does not change"
     write(*,*)"  (4) # of k-points in GS and # of blocks in super-cell in RT is the same"
     write(*,*)"  (5) cuboid cell"
     write(*,*)"  (6) process_allocation = orbital_sequential"
     write(*,*)"  (7) nproc_rgrid(1,2,3) --> x Nkx, x Nky, x Nkz"
     write(*,*)"  (8) nproc_k = nk in the input data file"
     write(*,*)"  (9)nspin=1, periodic system, ...."
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

subroutine assign_rank_mumber_to_read_write_files(kex,system,info)
  use parallelization, only: nproc_size_global
  use communication, only: comm_is_root, comm_bcast
  implicit none
  type(s_k_expand) :: kex
  type(s_dft_system) :: system
  type(s_parallel_info),intent(in) :: info
  integer :: nfile,io_new, myrank, icnt,nl,i4,i5

  myrank   = info%id_rko

  kex%nk_new = 1
  kex%no_new = system%no * kex%nk

  nfile = system%nk * system%no

  if( nproc_size_global >= nfile ) then
     kex%nmax=1
  else
     kex%nmax = nfile / nproc_size_global + 1
  endif
  
  allocate( kex%myrank(kex%nmax) )
  allocate( kex%iaddress(2,kex%nmax), kex%iaddress_new(2,kex%nmax) )

  kex%myrank(:) = -1
  kex%iaddress(:,:) = -1
  kex%iaddress_new(:,:) = -1

  !!(for reading)
  icnt=0
  nl=-1
  do i5=1,system%nk
  do i4=1,system%no
     nl = nl + 1
     if(mod(nl,nproc_size_global) == myrank) then
        icnt = icnt + 1
        kex%myrank(icnt) = nl
        kex%iaddress(1,icnt) = i4  !o
        kex%iaddress(2,icnt) = i5  !k
     end if
  end do
  end do

  io_new = 0
  do i4=1,system%no
  do i5=1,system%nk
     io_new = io_new + 1
     do icnt=1,kex%nmax
        if(kex%iaddress(1,icnt) == i4 .and. &
           kex%iaddress(2,icnt) == i5) then
           kex%iaddress_new(1,icnt) = io_new !orbital
           kex%iaddress_new(2,icnt) = 1      !k
        end if
     enddo
  end do
  end do

 !write(*,'(a,100i5)') "myrank  =", kex%myrank(:)

  end subroutine

  subroutine get_restart_reading_directory_name(system,nblock,gdir,pdir)
    use salmon_global, only: directory_read_data
    implicit none
    type(s_dft_system) :: system
    character(256),intent(out) :: gdir
    character(256),intent(out) :: pdir(:)
    integer :: ik,io,nblock, i

    write(gdir,'(A,I6.6,A)')   trim(directory_read_data)

    i=0
    do ik=1,system%nk 
    do io=1,system%no, nblock
       i=i+1
       write (pdir(i),'(A,I3.3,A,I6.6)') trim(gdir)//'k_',ik,'_ob_',io
    enddo
    enddo

  end subroutine get_restart_reading_directory_name

  subroutine gen_restart_writing_directory_name_k_expand(kex,nblock,odir,gdir,pdir)
    implicit none
    type(s_k_expand) :: kex
    character(*),  intent(in)  :: odir
    character(256),intent(out) :: gdir
    character(256),intent(out) :: pdir(:)
    integer :: icnt, i,ik_new,io_new, nblock

    write(gdir,'(A,I6.6,A)')   trim(odir)

    ik_new=1
    i=0
    do io_new= 1,kex%no_new, nblock
       i=i+1
       write(pdir(i),'(A,I3.3,A,I6.6)') trim(gdir)//'k_',ik_new,'_ob_',io_new
       do icnt=1,kex%nmax
          if(kex%iaddress_new(1,icnt)==io_new) then
             call create_directory(pdir(i))
          endif
       enddo
    enddo

  end subroutine gen_restart_writing_directory_name_k_expand


  subroutine read_write_info_bin(idir,odir,kex,system)
    implicit none
    type(s_k_expand) :: kex
    type(s_dft_system) :: system
    character(*)    :: idir,odir
    character(1024) :: dir_file_in, dir_file_out
    logical :: if_real_orbital_tmp
    integer :: iu1,iu3, iter, mk,mo,nprocs

    iu1 = 91
    iu3 = 93
    dir_file_in  = trim(idir)//"/info.bin"
    dir_file_out = trim(odir)//"/info.bin"
    open(iu1,file=dir_file_in, form='unformatted')
    open(iu3,file=dir_file_out,form='unformatted')


    read(iu1) mk
    read(iu1) mo
    read(iu1) iter
    read(iu1) nprocs
    read(iu1) if_real_orbital_tmp
    close(iu1)

    if(system%nk /= mk) stop "error in reading info.bin: mk"
    if(system%no /= mo) stop "error in reading info.bin: mo"
    if(if_real_orbital_tmp) stop "error in reading info.bin:if_real_orbital"
   !iter   = 0
    nprocs = nprocs * kex%nk
  
    write(iu3) kex%nk_new
    write(iu3) kex%no_new
    write(iu3) iter
    write(iu3) nprocs
    write(iu3) if_real_orbital_tmp
    close(iu3)

   !write(*,*) "info_bin:", kex%nk_new, kex%no_new, iter, nprocs
    
  end subroutine read_write_info_bin

  subroutine read_write_atomic_coor(idir,odir,kex,system)
    use salmon_global, only: al,natom,atom_name,kion,unit_length
    use inputoutput, only: au_length_aa
    implicit none
    type(s_k_expand) :: kex
    type(s_dft_system) :: system
    character(*)    :: idir,odir
    character(1024) :: dir_file_out  !,dir_file_in
    integer :: iu3, ia, ixyz  !,iu1
    real(8) :: uconv, Rion_new(3)

    if(unit_length=='AA')then ; uconv = au_length_aa
    else                      ; uconv = 1d0   !au
    endif

    !atomic coordinate is read from input file, not from restart data
   !iu1 = 91
    iu3 = 93
   !dir_file_in  = trim(idir)//"/atomic_coor.txt"
    dir_file_out = trim(odir)//"/atomic_coor.txt"
   !open(iu1,file=dir_file_in, status="old")
    open(iu3,file=dir_file_out,status="unknown")

    do ia = 1,natom
       atom_name(ia) = adjustl(atom_name(ia))
       atom_name(ia) = "'"//trim(atom_name(ia))//"'"
    enddo
    do ixyz= 1,kex%nk
       do ia = 1,natom
          Rion_new(:) = system%Rion(:,ia) + al(:) * (kex%isupercell(:,ixyz)-1)
          write(iu3,7000) trim(atom_name(ia)), Rion_new(1:3)*uconv, kion(ia)
       enddo
    enddo
7000 format(" ",a6,"  ",3f18.10,i4)
    close(iu3)

  end subroutine read_write_atomic_coor

  subroutine read_write_atomic_vel(idir,odir,kex,system)
    use salmon_global, only: natom
    implicit none
    type(s_k_expand) :: kex
    type(s_dft_system) :: system
    character(*)    :: idir,odir
    character(1024) :: dir_file_out  !,dir_file_in
    integer :: iu3, ia, ixyz  !,iu1

    !atomic velocity is not read from restart data
   !iu1 = 91
    iu3 = 93
   !dir_file_in  = trim(idir)//"/atomic_vel.txt"
    dir_file_out = trim(odir)//"/atomic_vel.txt"
   !open(iu1,file=dir_file_in, status="old")
    open(iu3,file=dir_file_out,status="unknown")

    do ixyz= 1,kex%nk
       do ia = 1,natom
          write(iu3,7100) system%Velocity(1:3,ia)
       enddo
    enddo
7100 format(3f18.10)
    close(iu3)

  end subroutine read_write_atomic_vel

  subroutine read_write_occ_bin(idir,odir,kex,system)
    implicit none
    type(s_k_expand) :: kex
    type(s_dft_system) :: system
    character(*)    :: idir,odir
    character(1024) :: dir_file_in,dir_file_out
    integer :: iu1,iu3, io,ik, io_new
    real(8),allocatable :: rocc_new(:,:,:)

    iu1 = 91
    iu3 = 93
    dir_file_in  = trim(idir)//"/occupation.bin"
    dir_file_out = trim(odir)//"/occupation.bin"
    open(iu1,file=dir_file_in, form='unformatted')
    open(iu3,file=dir_file_out,form='unformatted')


    read(iu1) system%rocc(1:system%no,1:system%nk,1:system%nspin)
    close(iu1)


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
    write(iu3) rocc_new(1:kex%no_new,1:kex%nk_new,1:system%nspin)
    close(iu3)

  deallocate(rocc_new)

end subroutine read_write_occ_bin

subroutine read_write_rho_inout_bin(idir,odir,kex,system,lg,mixing)
  implicit none
  type(s_k_expand) :: kex
  type(s_dft_system) :: system
  type(s_rgrid) :: lg
  type(s_mixing) :: mixing
  character(*)    :: idir,odir
  character(1024) :: dir_file_in,dir_file_out  !,ofile_cube
  integer :: iu1,iu3,num_rho_stock, num_stock_inout
  integer :: lg_is_new(3), lg_ie_new(3), ixyz
  integer :: ix_add, iy_add, iz_add,iix,iiy,iiz
 !real(8) :: norm
  real(8),allocatable :: rho(:,:,:), rho_new(:,:,:)

  iu1 = 91
  iu3 = 93
  dir_file_in  = trim(idir)//"/rho_inout.bin"
  dir_file_out = trim(odir)//"/rho_inout.bin"
  open(iu1,file=dir_file_in, form='unformatted')  !read
  open(iu3,file=dir_file_out,form='unformatted')  !write

  lg_is_new(:) = lg%is(:)
  lg_ie_new(1) = lg%ie(1) * kex%nkx
  lg_ie_new(2) = lg%ie(2) * kex%nky
  lg_ie_new(3) = lg%ie(3) * kex%nkz


  allocate( rho( lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3) ) )
  allocate( rho_new( lg_ie_new(1),lg_ie_new(2),lg_ie_new(3) ) )

  num_rho_stock = mixing%num_rho_stock
  num_stock_inout = num_rho_stock + num_rho_stock+1

  do i=1,num_stock_inout

     read(iu1) rho( lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3) )

     do ixyz=1,kex%nk

        ix_add = lg%ie(1) * (kex%isupercell(1,ixyz)-1)
        iy_add = lg%ie(2) * (kex%isupercell(2,ixyz)-1)
        iz_add = lg%ie(3) * (kex%isupercell(3,ixyz)-1)

        !$omp parallel do collapse(2) private(ix,iy,iz,iix,iiy,iiz)
        do iz= lg%is(3),lg%ie(3)
        do iy= lg%is(2),lg%ie(2)
        do ix= lg%is(1),lg%ie(1)

           iix = ix + ix_add
           iiy = iy + iy_add
           iiz = iz + iz_add
           rho_new(iix,iiy,iiz)  = rho(ix,iy,iz)

        enddo
        enddo
        enddo

     enddo

     write(iu3) rho_new( 1:lg_ie_new(1),1:lg_ie_new(2),1:lg_ie_new(3) )

     !!(check)--norm
     !norm = sum( rho_new(:,:,:) )*system%hvol
     !write(*,*) "norm of rho in density", i,norm
     !rho_new(:,:,:) = rho_new(:,:,:)/norm


     !!(check)--print cube file
     !if(i==22) then
     !   ofile_cube="rho.cube"
     !   call write_cube(rho_new,lg_ie_new(1),lg_ie_new(2),lg_ie_new(3),system%hgs,ofile_cube)
     !endif

  enddo ! i
  deallocate(rho, rho_new)

end subroutine read_write_rho_inout_bin

subroutine read_write_wfn_bin(idir,odir,kex,nblock,system,lg)
  implicit none
  type(s_k_expand) :: kex
  type(s_dft_system) :: system
  type(s_rgrid) :: lg
  character(*)    :: idir(:), odir(:)
  character(1024) :: dir_file_in,dir_file_out  !,ofile_cube
  integer :: i,icnt,io,ik,io_new, io_r,io_w, ik_r,ik_w  
  integer :: myrank, i_rdir,i_wdir, nblock
  integer :: lg_is_new(3), lg_ie_new(3)
  integer :: ix_add, iy_add, iz_add
  integer :: iu,ixyz, ix,iy,iz,iix,iiy,iiz   !,mx,my,mz
  real(8) :: scale, rshift(3), r(3)  !,norm
 !real(8),allocatable :: orb_real(:,:,:)
  complex(8) :: ai,ekr
  complex(8),allocatable :: zwf(:,:,:), zwf_ekr(:,:,:)

  ai        = ( 0d0, 1d0 )
  rshift(:) = -system%hgs(:)    !probably only if(yn_domain_parallel=='y')
  scale = 1d0/sqrt(dble(kex%nk))

  lg_is_new(:) = lg%is(:)
  lg_ie_new(1) = lg%ie(1) * kex%nkx
  lg_ie_new(2) = lg%ie(2) * kex%nky
  lg_ie_new(3) = lg%ie(3) * kex%nkz

  allocate( zwf(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)) )
  allocate( zwf_ekr(1:lg_ie_new(1),1:lg_ie_new(2),1:lg_ie_new(3)) )


  do icnt=1,kex%nmax

     myrank = kex%myrank(icnt)
     if(myrank .le. -1) cycle
     io_r = kex%iaddress(1,icnt)
     ik_r = kex%iaddress(2,icnt)
     io_w = kex%iaddress_new(1,icnt)
     ik_w = kex%iaddress_new(2,icnt)

     i=0
     do ik=1,system%nk 
     do io=1,system%no, nblock
        i=i+1
        if(ik_r==ik .and. io_r.ge.io .and. io_r.lt.io+nblock) then
           i_rdir = i
           exit
        endif
     enddo
     enddo

     i=0
     do io_new= 1,kex%no_new, nblock
        i=i+1
        if(io_w.ge.io_new .and. io_w.lt.io_new+nblock) then
           i_wdir = i
           exit
        endif
     enddo

     iu = 1000 + myrank

     write(dir_file_in, '(A,I6.6,A)') trim(idir(i_rdir))//'/wfn_ob_',io_r,'.dat'
     write(dir_file_out,'(A,I6.6,A)') trim(odir(i_wdir))//'/wfn_ob_',io_w,'.dat'
    !write(*,'(3a)') trim(dir_file_in),"  ", trim(dir_file_out)

     open(iu,file=dir_file_in,form='unformatted',access='stream',status='old')
     read(iu) zwf(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
     close(iu)

     !!(check)-norm
     !norm = sum(abs(zwf(:,:,:))**2) * system%hvol
     !write(*,*) "norm=", norm

     do ixyz=1,kex%nk

        ix_add = lg%ie(1) * (kex%isupercell(1,ixyz)-1)
        iy_add = lg%ie(2) * (kex%isupercell(2,ixyz)-1)
        iz_add = lg%ie(3) * (kex%isupercell(3,ixyz)-1)

        do iz = lg%is(3), lg%ie(3)
           iiz = iz + iz_add
           r(3) = iiz*system%hgs(3) + rshift(3)

           do iy = lg%is(2), lg%ie(2)
              iiy = iy + iy_add
              r(2) = iiy*system%hgs(2) + rshift(2)

              do ix = lg%is(1), lg%ie(1)
                 iix = ix + ix_add
                 r(1) = iix*system%hgs(1) + rshift(1)

                 ekr = exp(ai * sum(kex%k_vec(:,ik_r)*r(:)) ) * scale
                 zwf_ekr(iix,iiy,iiz) = zwf(ix,iy,iz)* ekr
                 
              enddo
           enddo
        enddo
     enddo

     !!(check)--norm
     !norm = sum(abs(zwf_ekr(:,:,:))**2) * system%hvol
     !write(*,*) "norm(before)=", norm
     !zwf_ekr(:,:,:) = zwf_ekr(:,:,:) / sqrt(norm)
     !norm = sum(abs(zwf_ekr(:,:,:))**2) * system%hvol
     !write(*,*) "norm(after)=", norm

     !write
     open(iu,file=dir_file_out,form='unformatted',access='stream')
     write(iu) zwf_ekr(1:lg_ie_new(1),1:lg_ie_new(2),1:lg_ie_new(3))
     close(iu)

     !!(check)--print cube file
     !if(io_w==20) then
     !   mx=lg_ie_new(1)
     !   my=lg_ie_new(2)
     !   mz=lg_ie_new(3)
     !   allocate(orb_real(mx,my,mz))
     !   orb_real(1:mx,1:my,1:mz) = real(zwf_ekr(1:mx,1:my,1:mz))
     !   ofile_cube="orb_20.cube"
     !   call  write_cube(orb_real,mx,my,mz,system%hgs,ofile_cube)
     !endif


  enddo
  deallocate( zwf, zwf_ekr )

end subroutine read_write_wfn_bin

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
     !write(fp,'(i5,4f12.6)') izatom(ik),dble(izatom(ik)),(rion(j,iatom),j=1,3)
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

end subroutine main_dft_k_expand_slice
