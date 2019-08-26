!
!  Copyright 2019 SALMON developers
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

MODULE global_variables_scf

use inputoutput
use scf_data
use allocate_mat_sub
use deallocate_mat_sub
use new_world_sub
use init_sendrecv_sub
use change_order_sub
use read_pslfile_sub
use allocate_psl_sub
use persistent_comm
use structure_opt_sub
use calc_allob_sub
use salmon_total_energy
use hpsi_sub
implicit none

END MODULE global_variables_scf

!=======================================================================

subroutine Real_Space_DFT
use structures
use salmon_parallel, only: nproc_id_global, nproc_size_global, nproc_group_global, &
                           nproc_group_h, nproc_id_kgrid
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use salmon_xc, only: init_xc, finalize_xc
use timer
use calc_iobnum_sub
use check_mg_sub
use check_ng_sub
use scf_iteration_sub
use rmmdiis_sub
use density_matrix, only: calc_density
use writefield
use global_variables_scf
use sendrecv_grid, only: s_sendrecv_grid, init_sendrecv_grid
use salmon_pp, only: calc_nlcc
use force_sub
use calc_iroot_sub
use gram_schmidt_orth, only: gram_schmidt 
use print_sub
use read_gs
use code_optimization
use salmon_initialization
implicit none
integer :: ix,iy,iz,ik,ikoa,is,i,j
integer :: iter,iatom,iob,p1,p2,p5,ii,jj,iflag,jspin
real(8) :: sum0,sum1
character(100) :: file_atoms_coo, comment_line
complex(8),allocatable :: zpsi_tmp(:,:,:,:,:)
real(8) :: rNebox1,rNebox2
integer :: itmg,nspin,n,nn
integer :: neig(1:3, 1:2)
integer :: neig_ng(1:3, 1:2)

type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_rgrid) :: ng
type(s_orbital_parallel) :: info
type(s_field_parallel) :: info_field
type(s_sendrecv_grid) :: srg, srg_ng, srg_ob_1
type(s_orbital) :: spsi,shpsi,sttpsi
type(s_dft_system) :: system
type(s_poisson_cg) :: poisson_cg
type(s_stencil) :: stencil
type(s_scalar) :: srho
type(s_scalar) :: sVh,sVpsl
type(s_scalar),allocatable :: V_local(:),srho_s(:),sVxc(:)
type(s_reciprocal_grid) :: fg
type(s_pp_nlcc) :: ppn
type(s_dft_energy) :: energy
type(s_cg)  :: cg
type(s_mixing) :: mixing

logical :: rion_update

call init_xc(xc_func, ispin, cval, xcname=xc, xname=xname, cname=cname)

iSCFRT=1
ihpsieff=0
iflag_comm_rho=1

iblacsinit=0

call timer_begin(LOG_TOTAL)


call timer_begin(LOG_INIT_GS)
inumcpu_check=0

call setbN(bnmat)
call setcN(cnmat)

call check_dos_pdos

call convert_input_scf(info,info_field,file_atoms_coo,mixing,poisson_cg)

call init_dft(lg,system,stencil)
if(stencil%if_orthogonal) then
  if(comm_is_root(nproc_id_global)) write(*,*) "orthogonal cell: using al"
else
  if(comm_is_root(nproc_id_global)) write(*,*) "non-orthogonal cell: using al_vec[1,2,3]"
end if
allocate(system%mass(1:nelem))

call set_filename

call setk(k_sta, k_end, k_num, num_kpoints_rd, nproc_k, info%id_k)

call calc_iobnum(itotMST,nproc_id_kgrid,iobnum,nproc_ob)

if(iflag_opt==1)then
   call structure_opt_ini(MI)
   flag_opt_conv=.false.
   write(comment_line,10) 0
   call write_xyz(comment_line,"new","r  ",system)
10 format("#opt iteration step=",i5)
end if
call timer_end(LOG_INIT_GS)


Structure_Optimization_Iteration : do iopt=1,iter_opt
Multigrid_Iteration : do img=1,ntmg

if(iopt==1)then
  select case( IC )
  case default ! New calculation

    call timer_begin(LOG_INIT_GS)

    Hvol = system%Hvol
    Hgs = system%Hgs
    Miter = 0        ! Miter: Iteration counter set to zero
    itmg=img
    call set_imesh_oddeven(itmg)
    call init_mesh(lg,mg)
    call set_gridcoo
    call init_mesh_s(ng)
    call check_mg(mg)
    call check_ng(ng)

  case(1,3) ! Continue the previous calculation

    call IN_data(lg,mg,ng,info,info_field,system,stencil,mixing)

  end select

  call init_updown(info)
  call init_itype
  call init_sendrecv_matrix
  select case(iperiodic)
  case(0)
    if(layout_multipole==2.or.layout_multipole==3) call make_corr_pole(ng,poisson_cg)
  end select
  call set_ig_bound(ng,poisson_cg)

  call allocate_mat(ng)
  call set_icoo1d
  call allocate_sendrecv
  call init_persistent_requests(info)
  call init_code_optimization

  ! Initialization of s_sendrecv_grid structure (experimental implementation)
  neig(1, 1) = iup_array(1)
  neig(1, 2) = idw_array(1)
  neig(2, 1) = jup_array(1)
  neig(2, 2) = jdw_array(1)
  neig(3, 1) = kup_array(1)
  neig(3, 2) = kdw_array(1)
  call init_sendrecv_grid(srg, mg, iobnum * k_num, info%icomm_r, neig)
  call init_sendrecv_grid(srg_ob_1, mg, 1, info%icomm_r, neig)

  neig_ng(1, 1) = iup_array(2)
  neig_ng(1, 2) = idw_array(2)
  neig_ng(2, 1) = jup_array(2)
  neig_ng(2, 2) = jdw_array(2)
  neig_ng(3, 1) = kup_array(2)
  neig_ng(3, 2) = kdw_array(2)
  call init_sendrecv_grid(srg_ng, ng, 1, &
    & nproc_group_global, neig_ng)

  if(ispin==0)then
    nspin=1
  else
    nspin=2
  end if

  system%rocc(1:itotMST,1:system%nk,1) = rocc(1:itotMST,1:system%nk)

  allocate(energy%esp(system%no,system%nk,system%nspin))

  info%im_s = 1
  info%im_e = 1
  info%numm = 1
  info%ik_s = k_sta
  info%ik_e = k_end
  info%numk = k_num
  info%io_s = 1
  info%io_e = iobnum/nspin
  info%numo = iobnum/nspin

  info%if_divide_rspace = nproc_d_o_mul.ne.1
  info%if_divide_orbit  = nproc_ob.ne.1
  info%icomm_rko  = nproc_group_global
  allocate(info%occ(info%io_s:info%io_e, info%ik_s:info%ik_e, 1:system%nspin,1) &
            ,info%io_tbl(info%io_s:info%io_e), info%jo_tbl(1:system%no) &
            ,info%irank_jo(1:system%no))

  info%jo_tbl(:) = 0 !(initial value)
  do iob=info%io_s,info%io_e
    call calc_allob(iob,jj,itotmst,mst,iobnum)
    info%io_tbl(iob) = jj
    info%jo_tbl(jj)  = iob
  end do

  do jspin=1,system%nspin
    do ik=info%ik_s,info%ik_e
      do iob=info%io_s,info%io_e
        jj = info%io_tbl(iob)
        info%occ(iob,ik,jspin,1) = system%rocc(jj,ik,jspin)*system%wtk(ik)
      end do
    end do
  end do

  do jj=1, system%no
    call calc_iroot(jj,info%irank_jo(jj),ilsda,nproc_ob,itotmst,mst)
  end do

  allocate(srho_s(system%nspin),V_local(system%nspin),sVxc(system%nspin))

  call allocate_scalar(mg,srho)
  call allocate_scalar(mg,sVh)
  call allocate_scalar(mg,sVpsl)
  do jspin=1,system%nspin
    call allocate_scalar(mg,srho_s(jspin))
    call allocate_scalar(mg,V_local(jspin))
    call allocate_scalar(mg,sVxc(jspin))
  end do

  select case(iperiodic)
  case(0)
    call allocate_orbital_real(system%nspin,mg,info,spsi)
    call allocate_orbital_real(system%nspin,mg,info,shpsi)
  case(3)
    call allocate_orbital_complex(system%nspin,mg,info,spsi)
    call allocate_orbital_complex(system%nspin,mg,info,shpsi)
    call allocate_orbital_complex(system%nspin,mg,info,sttpsi)
  end select

  if(iperiodic==3)then
    allocate (zpsi_tmp(mg%is_overlap(1):mg%ie_overlap(1) &
    &                 ,mg%is_overlap(2):mg%ie_overlap(2) &
    &                 ,mg%is_overlap(3):mg%ie_overlap(3) &
    &                 ,1:iobnum,k_sta:k_end))
  end if

  if(iperiodic==3 .and. iflag_hartree==4)then
    call prep_poisson_fft(ng)
  end if

  if(.not. allocated(Vpsl)) allocate( Vpsl(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
  if(.not. allocated(Vpsl_atom)) allocate( Vpsl_atom(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),MI) )
  if(iflag_ps.eq.0)then
    Vpsl=0d0
  else
    call read_pslfile(system)
    call allocate_psl
    call init_ps(ng,system%primitive_a,system%primitive_b,stencil%rmatrix_A,info%icomm_r)
  end if
  sVpsl%f = Vpsl

  if(iperiodic==3) then
    allocate(stencil%vec_kAc(3,info%ik_s:info%ik_e))
    stencil%vec_kAc(:,info%ik_s:info%ik_e) = system%vec_k(:,info%ik_s:info%ik_e)
    call update_kvector_nonlocalpt(ppg,stencil%vec_kAc,info%ik_s,info%ik_e)
  end if

  if(iperiodic==3) call get_fourier_grid_G(fg)

  select case( IC )
  case default ! New calculation

    if(iobnum >= 1)then
      select case(iperiodic)
      case(0)
        allocate( psi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum,k_sta:k_end) )
      case(3)
        allocate( ttpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
        allocate( zpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum,k_sta:k_end) )
      end select
    end if
    if(iswitch_orbital_mesh==1.or.iflag_subspace_diag==1)then
      select case(iperiodic)
      case(0)
        allocate( psi_mesh(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3),1:itotMST,1) )
      case(3)
        allocate( zpsi_mesh(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3),1:itotMST,num_kpoints_rd) )
      end select
    end if

    if(read_gs_wfn_k=='n') then
      call init_wf_ns(1)
      ! Store to psi/zpsi
      select case(iperiodic)
      case(0)
        do ik=k_sta,k_end
        do iob=1,info%numo
          do is=1,nspin
            spsi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),is,iob,ik,1) = &
            & psi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),iob+(is-1)*info%numo,ik)
          end do
        end do
        end do
      case(3)
        do ik=k_sta,k_end
        do iob=1,info%numo
          do is=1,nspin
            spsi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),is,iob,ik,1) = &
            & zpsi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),iob+(is-1)*info%numo,ik)
          end do
        end do
      end do
      end select
    else
      if(iperiodic==0) stop "error: read_gs_wfn_k='y' & iperiodic=0"
      call read_wfn(lg,mg,spsi,info,system)
    end if

    call gram_schmidt(system, mg, info, spsi)

    allocate( rho(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
    if(ilsda == 1)then
      allocate( rho_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2) )
    end if

    allocate(mixing%srho_in(1:mixing%num_rho_stock+1))
    allocate(mixing%srho_out(1:mixing%num_rho_stock+1))
    do i=1,mixing%num_rho_stock+1
      allocate(mixing%srho_in(i)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
      allocate(mixing%srho_out(i)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
      mixing%srho_in(i)%f(:,:,:)=0.d0
      mixing%srho_out(i)%f(:,:,:)=0.d0
    end do

    if(ilsda==1)then
      allocate(mixing%srho_s_in(1:mixing%num_rho_stock+1,2))
      allocate(mixing%srho_s_out(1:mixing%num_rho_stock+1,2))
      do j=1,2
        do i=1,mixing%num_rho_stock+1
          allocate(mixing%srho_s_in(i,j)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
          allocate(mixing%srho_s_out(i,j)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
          mixing%srho_s_in(i,j)%f(:,:,:)=0.d0
          mixing%srho_s_out(i,j)%f(:,:,:)=0.d0
        end do
      end do
    end if

    if(read_gs_dns_cube == 'n') then
      call calc_density(srho_s,spsi,info,mg,nspin)
    else
      if(ispin/=0) stop "read_gs_dns_cube=='n' & ispin/=0"
      call read_dns(lg,mg,srho_s(1)%f) ! cube file only
    end if

    srho%f = 0d0
    do jspin=1,nspin
      srho%f = srho%f + srho_s(jspin)%f
    end do
    rho = srho%f

    allocate (Vlocal(mg_sta(1):mg_end(1),  &
                mg_sta(2):mg_end(2),  &
                mg_sta(3):mg_end(3),nspin))

    allocate( Vh(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
    Vh=0.d0

    call Hartree_ns(lg,mg,ng,info_field,system,poisson_cg,srg_ng,stencil,srho,sVh,fg)
    Vh = sVh%f

    if(ilsda == 0) then
      allocate( Vxc(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
    else if(ilsda == 1) then
      allocate( Vxc_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2) )
    end if
    allocate( esp(itotMST,num_kpoints_rd) )

    call exc_cor_ns(ng, srg_ng, system%nspin, srho_s, ppn, sVxc, energy%E_xc)

    call allgatherv_vlocal(ng,info,system%nspin,sVh,sVpsl,sVxc,V_local)
    do jspin=1,system%nspin
      Vlocal(:,:,:,jspin) = V_local(jspin)%f
    end do

    call calc_eigen_energy(energy,spsi,shpsi,sttpsi,system,info,mg,V_local,stencil,srg,ppg)
    select case(iperiodic)
    case(0)
      call calc_Total_Energy_isolated(energy,system,info,ng,pp,srho_s,sVh,sVxc)
    case(3)
      rion_update = .true. ! it's first calculation
      call calc_Total_Energy_periodic(energy,system,pp,fg,rion_update)
    end select
    esp = energy%esp(:,:,1) !++++++++

  case(1,3) ! Continue the previous calculation

    select case(iperiodic)
    case(0)
      do ik=k_sta,k_end
      do iob=1,info%numo
        do is=1,nspin
          spsi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),is,iob,ik,1) = &
          & psi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),iob+(is-1)*info%numo,ik)
        end do
      end do
      end do
    case(3)
      do ik=k_sta,k_end
      do iob=1,info%numo
        do is=1,nspin
          spsi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),is,iob,ik,1) = &
          & zpsi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),iob+(is-1)*info%numo,ik)
        end do
      end do
    end do
    end select
    srho%f = rho
    if(ilsda == 1)then
      srho_s(1)%f = rho_s(:,:,:,1)
      srho_s(2)%f = rho_s(:,:,:,2)
    end if
    do jspin=1,nspin
      V_local(jspin)%f = Vlocal(:,:,:,jspin)
    end do

  end select

  call timer_end(LOG_INIT_GS)

else if(iopt>=2)then
  call timer_begin(LOG_INIT_GS)
  Miter = 0        ! Miter: Iteration counter set to zero
  if(iflag_ps/=0) then
    call dealloc_init_ps(ppg,ppg_all,ppn)
    call init_ps(ng,system%primitive_a,system%primitive_b,stencil%rmatrix_A,info%icomm_r)
    if(iperiodic==3) call get_fourier_grid_G(fg)

  end if
  call timer_end(LOG_INIT_GS)
end if


if(comm_is_root(nproc_id_global)) then
  write(*,*) '-----------------------------------------------'
  select case(iperiodic)
  case(0)
    write(*,'(1x,"iter =",i6,5x,"Total Energy =",f19.8,5x,"Vh iteration =",i4)') Miter,energy%E_tot*2d0*Ry,iterVh
  case(3)
    write(*,'(1x,"iter =",i6,5x,"Total Energy =",f19.8)') Miter,energy%E_tot*2d0*Ry
  end select
  do ik=1,num_kpoints_rd
    if(ik<=3)then
      if(iperiodic==3) write(*,*) "k=",ik
      do p5=1,(itotMST+3)/4
        p1=4*(p5-1)+1
        p2=4*p5 ; if ( p2 > itotMST ) p2=itotMST
        write(*,'(1x,4(i5,f15.4,2x))') (iob,energy%esp(iob,ik,1)*2d0*Ry,iob=p1,p2)
      end do
      if(iperiodic==3) write(*,*)
    end if
  end do
end if


!---------------------------------------- Iteration

call timer_begin(LOG_INIT_GS_ITERATION)
iflag=1
iterVh=1000
sum1=1.0d9

iflag_diisjump=0

allocate(idiis_sd(itotMST))
idiis_sd=0

if(img==1.and.iopt==1) allocate(norm_diff_psi_stock(itotMST,1))
norm_diff_psi_stock=1.0d9

if(img>=2.or.iopt>=2) deallocate(rho_stock,Vlocal_stock)
allocate(rho_stock(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3),1))
if(ilsda==0)then
  allocate(Vlocal_stock(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3),1))
else
  allocate(Vlocal_stock(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3),1:2))
end if

if(ilsda==0)then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    rho_stock(ix,iy,iz,1)=rho(ix,iy,iz)
    Vlocal_stock(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)
  end do
  end do
  end do
else
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    rho_stock(ix,iy,iz,1)=rho(ix,iy,iz)
    Vlocal_stock(ix,iy,iz,1:2)=Vlocal(ix,iy,iz,1:2)
  end do
  end do
  end do
end if

! Setup NLCC term from pseudopotential
call calc_nlcc(pp, system, mg, ppn)

if (comm_is_root(nproc_id_global)) then
  write(*, '(1x, a, es23.15e3)') "Maximal rho_NLCC=", maxval(ppn%rho_nlcc)
  write(*, '(1x, a, es23.15e3)') "Maximal tau_NLCC=", maxval(ppn%tau_nlcc)
end if    
call timer_end(LOG_INIT_GS_ITERATION)

call timer_begin(LOG_GS_ITERATION)
DFT_Iteration : do iter=1,iDiter(img)


  if(sum1<threshold) cycle DFT_Iteration

  Miter=Miter+1

  ! for calc_total_energy_periodic
  rion_update = check_rion_update() .or. (iter == 1)

  if(temperature_k>=0.d0.and.Miter>iditer_notemperature) then
    if(iperiodic.eq.3) then
      call ne2mu_p
    else
      call ne2mu
    endif
  else
    call calc_occupation
  endif

  system%rocc(1:itotMST,1:system%nk,1) = rocc(1:itotMST,1:system%nk)

  do jspin=1,system%nspin
    do ik=info%ik_s,info%ik_e
      do iob=info%io_s,info%io_e
        jj = info%io_tbl(iob)
        info%occ(iob,ik,jspin,1) = system%rocc(jj,ik,jspin)*system%wtk(ik)
      end do
    end do
  end do

  call copy_density(system%nspin,ng,srho_s,mixing)

  if(iscf_order==1)then

    call scf_iteration(mg,ng,system,info,stencil,srg,spsi,shpsi,srho,srho_s,itotmst,mst, &
                       cg,ppg,V_local,  &
                       iflag_diisjump,energy, &
                       norm_diff_psi_stock, &
                       Miter,iDiterYBCG,   &
                       iflag_subspace_diag,iditer_nosubspace_diag,ifmst,mixing,iter)

    call timer_begin(LOG_CALC_HARTREE)
    if(imesh_s_all==1.or.(imesh_s_all==0.and.nproc_id_global<nproc_d_o_mul*nproc_d_g_mul_dm))then
      call Hartree_ns(lg,mg,ng,info_field,system,poisson_cg,srg_ng,stencil,srho,sVh,fg)
    end if
    call timer_end(LOG_CALC_HARTREE)

    call timer_begin(LOG_CALC_EXC_COR)
    if(imesh_s_all==1.or.(imesh_s_all==0.and.nproc_id_global<nproc_d_o_mul*nproc_d_g_mul_dm))then
      call exc_cor_ns(ng, srg_ng, system%nspin, srho_s, ppn, sVxc, energy%E_xc)
    end if
    call timer_end(LOG_CALC_EXC_COR)

    call allgatherv_vlocal(ng,info,system%nspin,sVh,sVpsl,sVxc,V_local)

    call timer_begin(LOG_CALC_TOTAL_ENERGY)
    call calc_eigen_energy(energy,spsi,shpsi,sttpsi,system,info,mg,V_local,stencil,srg,ppg)
    select case(iperiodic)
    case(0)
      call calc_Total_Energy_isolated(energy,system,info,ng,pp,srho_s,sVh,sVxc)
    case(3)
      call calc_Total_Energy_periodic(energy,system,pp,fg,rion_update)
    end select
    esp = energy%esp(:,:,1) !++++++++
    call timer_end(LOG_CALC_TOTAL_ENERGY)


    call timer_begin(LOG_CALC_CHANGE_ORDER)
    if(iperiodic==0)then  
      call change_order(psi,info)
    end if
    call timer_end(LOG_CALC_CHANGE_ORDER)

  end if

  call timer_begin(LOG_WRITE_GS_RESULTS)

  select case(convergence)
    case('rho_dne')
      sum0=0.d0
!$OMP parallel do reduction(+:sum0) private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3) 
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        sum0=sum0+abs(srho%f(ix,iy,iz)-rho_stock(ix,iy,iz,1))
      end do
      end do
      end do
      call comm_summation(sum0,sum1,nproc_group_h)
      if(ispin==0)then
        sum1=sum1*Hvol/(dble(ifMST(1))*2.d0)
      else if(ispin==1)then
        sum1=sum1*Hvol/dble(ifMST(1)+ifMST(2))
      end if
    case('norm_rho','norm_rho_dng')
      sum0=0.d0
!$OMP parallel do reduction(+:sum0) private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3) 
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        sum0=sum0+(srho%f(ix,iy,iz)-rho_stock(ix,iy,iz,1))**2
      end do
      end do
      end do
      call comm_summation(sum0,sum1,nproc_group_h)
      if(convergence=='norm_rho_dng')then
        sum1=sum1/dble(lg_num(1)*lg_num(2)*lg_num(3))
      end if
    case('norm_pot','norm_pot_dng')
      sum0=0.d0
!$OMP parallel do reduction(+:sum0) private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3) 
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        sum0=sum0+(V_local(1)%f(ix,iy,iz)-Vlocal_stock(ix,iy,iz,1))**2
      end do
      end do
      end do
      call comm_summation(sum0,sum1,nproc_group_h)
      if(convergence=='norm_pot_dng')then
        sum1=sum1/dble(lg_num(1)*lg_num(2)*lg_num(3))
      end if
  end select 

  if(comm_is_root(nproc_id_global)) then
    write(*,*) '-----------------------------------------------'
    select case(iperiodic)
    case(0)
      if(iflag_diisjump == 1) then
        write(*,'("Diisjump occured. Steepest descent was used.")')
      end if
      write(*,'(1x,"iter =",i6,5x,"Total Energy =",f19.8,5x,"Vh iteration =",i4)') Miter,energy%E_tot*2d0*Ry,iterVh
    case(3)
      write(*,'(1x,"iter =",i6,5x,"Total Energy =",f19.8,5x)') Miter,energy%E_tot*2d0*Ry
    end select
    do ik=1,num_kpoints_rd
      if(ik<=3)then
        if(iperiodic==3) write(*,*) "k=",ik
        do p5=1,(itotMST+3)/4
          p1=4*(p5-1)+1
          p2=4*p5 ; if ( p2 > itotMST ) p2=itotMST
          write(*,'(1x,4(i5,f15.4,2x))') (iob,energy%esp(iob,ik,1)*2d0*Ry,iob=p1,p2)
        end do
        if(iperiodic==3) write(*,*) 
      end if
    end do

    select case(convergence)
      case('rho_dne')
        write(*,'("iter and int_x|rho_i(x)-rho_i-1(x)|dx/nelec     = ",i6,e15.8)') Miter,sum1
      case('norm_rho')
        write(*,'("iter and ||rho_i(ix)-rho_i-1(ix)||**2              = ",i6,e15.8)') Miter,sum1/a_B**6
      case('norm_rho_dng')
        write(*,'("iter and ||rho_i(ix)-rho_i-1(ix)||**2/(# of grids) = ",i6,e15.8)') Miter,sum1/a_B**6
      case('norm_pot')
        write(*,'("iter and ||Vlocal_i(ix)-Vlocal_i-1(ix)||**2              = ",i6,e15.8)') Miter,     &
                                                                         sum1*(2.d0*Ry)**2/a_B**6
      case('norm_pot_dng')
        write(*,'("iter and ||Vlocal_i(ix)-Vlocal_i-1(ix)||**2/(# of grids) = ",i6,e15.8)') Miter,     &
                                                                         sum1*(2.d0*Ry)**2/a_B**6
    end select
  end if 
  rNebox1=0.d0 
!$OMP parallel do reduction(+:rNebox1) private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    rNebox1=rNebox1+srho%f(ix,iy,iz)
  end do
  end do
  end do
  call comm_summation(rNebox1,rNebox2,nproc_group_global)
  if(comm_is_root(nproc_id_global))then
    write(*,*) "Ne=",rNebox2*Hvol
  end if
  call timer_end(LOG_WRITE_GS_RESULTS)


if(ilsda==0)then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    rho_stock(ix,iy,iz,1)=srho%f(ix,iy,iz)
    Vlocal_stock(ix,iy,iz,1)=V_local(1)%f(ix,iy,iz)
  end do
  end do
  end do
else if(ilsda==1)then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    rho_stock(ix,iy,iz,1)=srho%f(ix,iy,iz)
    Vlocal_stock(ix,iy,iz,1)=V_local(1)%f(ix,iy,iz)
    Vlocal_stock(ix,iy,iz,2)=V_local(2)%f(ix,iy,iz)
  end do
  end do
  end do
end if

end do DFT_Iteration

! for OUT_data
Vh = sVh%f
rho = srho%f
if(ilsda == 1) then
  do jspin=1,system%nspin
    Vxc_s(:,:,:,jspin) = sVxc(jspin)%f
  end do
else
  Vxc = sVxc(1)%f
end if
Exc = energy%E_xc

! Store to psi/zpsi
select case(iperiodic)
case(0)
  do ik=k_sta,k_end
  do iob=1,info%numo
    do is=1,nspin
      !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        psi(ix,iy,iz,iob+(is-1)*info%numo,ik)=spsi%rwf(ix,iy,iz,is,iob,ik,1)
      end do
      end do
      end do
    end do
  end do
  end do
case(3)
  do ik=k_sta,k_end
  do iob=1,info%numo
    do is=1,nspin
      !$OMP parallel do private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        zpsi(ix,iy,iz,iob+(is-1)*info%numo,ik)=spsi%zwf(ix,iy,iz,is,iob,ik,1)
      end do
      end do
      end do
    end do
  end do
end do
end select

! output the wavefunctions for next GS calculations
if(write_gs_wfn_k == 'y') then
  if(iperiodic==3) then
    call write_wfn(lg,mg,spsi,info,system)
    
    ! Experimental Implementation of Inner-Product Outputs:
    call write_prod_dk_data(lg, mg, system, info, spsi) 
  else
    write(*,*) "error: write_gs_wfn_k='y' & iperiodic=0"
  end if
end if

! output transition moment
if(yn_out_tm  == 'y') then
  if(iperiodic==3) then
    call write_k_data(system,stencil)
    call write_tm_data(spsi,system,info,mg,stencil,srg,ppg)
  else
    write(*,*) "error: yn_out_tm='y' & iperiodic=0"
  end if
end if

! force
!if(iflag_opt==1) then
if (iperiodic == 3 .and. iflag_hartree == 4) then
  ! NOTE: calc_force_salmon hangs under this configuration due to ppg%vpsl_atom
  ! does not allocate.
else
   call calc_force_salmon(system,pp,fg,info,mg,stencil,srg,ppg,spsi)
   if(comm_is_root(nproc_id_global))then
      write(*,*) "===== force ====="
      do iatom=1,MI
         select case(unit_system)
         case('au','a.u.')
            write(*,'(i6,3e16.8)') iatom,(system%Force(ix,iatom),ix=1,3)
         case('A_eV_fs')
            write(*,'(i6,3e16.8)') iatom,(system%Force(ix,iatom)*2.d0*Ry/a_B,ix=1,3)
         end select
      end do
   end if
end if
!end if

if(iperiodic==3) deallocate(stencil%vec_kAc,ppg%zekr_uV)

deallocate(idiis_sd)
call timer_end(LOG_GS_ITERATION)

call timer_begin(LOG_DEINIT_GS_ITERATION)
if(iflag_opt==1) then
  call structure_opt_check(MI,iopt,flag_opt_conv,system%Force)
  if(.not.flag_opt_conv) call structure_opt(MI,iopt,system)
  !! Rion is old variables to be removed 
  !! but currently it is used in many subroutines.
  Rion(:,:) = system%Rion(:,:) 

  write(comment_line,10) iopt
  call write_xyz(comment_line,"add","r  ",system)

  if(comm_is_root(nproc_id_global))then
    write(*,*) "atomic coordinate"
    do iatom=1,MI
       write(*,20) "'"//trim(AtomName(Kion(iatom)))//"'",  &
                   (system%Rion(jj,iatom)*ulength_from_au,jj=1,3), &
                   Kion(iatom), flag_opt_atom(iatom)
    end do
20  format(a5,3f16.8,i3,a3)
  end if

  if(flag_opt_conv) then
    call structure_opt_fin
    exit Multigrid_Iteration
  end if

else
   select case(iperiodic)
   case(0) ; deallocate(spsi%rwf)
   case(3) ; deallocate(spsi%zwf)
   end select
end if
call timer_end(LOG_DEINIT_GS_ITERATION)


end do Multigrid_Iteration
if(flag_opt_conv)then
  exit Structure_Optimization_Iteration
end if
end do Structure_Optimization_Iteration


!---------------------------------------- Output
call timer_begin(LOG_WRITE_GS_RESULTS)

call band_information

call write_eigen

if(yn_out_psi=='y') then
  call writepsi(lg)
end if

if(yn_out_dns=='y') then
  call writedns(lg,mg,ng,rho,matbox_m,matbox_m2,icoo1d,hgs,igc_is,igc_ie,gridcoo,iscfrt)
end if

if(yn_out_dos=='y') then
  call calc_dos(info)
end if

if(yn_out_pdos=='y') then
  call calc_pdos(info)
end if

if(OC==2)then
  call prep_ini
end if

if(yn_out_elf=='y')then
  allocate(elf(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),      &
               lg_sta(3):lg_end(3)))
  call calcELF(mg,ng,info,srho,0,srg_ob_1)
  call writeelf(lg,elf,icoo1d,hgs,igc_is,igc_ie,gridcoo,iscfrt)
  deallocate(elf)
end if
call timer_end(LOG_WRITE_GS_RESULTS)


call timer_begin(LOG_WRITE_LDA_DATA)
! LDA data
! subroutines in scf_data.f90
if ( OC==1.or.OC==2.or.OC==3 ) then
  call OUT_data(ng,mixing)
end if
call timer_end(LOG_WRITE_LDA_DATA)


! LDA information
call timer_begin(LOG_WRITE_LDA_INFOS)
if(comm_is_root(nproc_id_global)) then
  open(1,file=LDA_info)

  write(1,*) "Total number of iteration = ", Miter
  write(1,*)
  select case (ilsda)
  case(0)
    write(1,*) "Number of states = ", nstate
    write(1,*) "Number of electrons = ", ifMST(1)*2
  case(1)
    write(1,*) "Number of states = ", (nstate_spin(is),is=1,2)
    write(1,*) "Number of electrons = ", (nelec_spin(is),is=1,2)
  end select
  write(1,*)
  write(1,*) "Total energy (eV) = ", energy%E_tot*2d0*Ry
  write(1,*) "1-particle energies (eV)"
  select case (ilsda)
  case(0)
    do p5=1,(nstate+3)/4
      p1=4*(p5-1)+1
      p2=4*p5 ; if ( p2 > nstate ) p2=nstate
      write(1,'(1x,4(i5,f15.4,2x))') (iob,energy%esp(iob,1,1)*2d0*Ry,iob=p1,p2)
    end do
  case(1)
    do is=1,2
      select case(is)
      case(1)
        write(1,*) "for up-spin"
        do p5=1,(nstate_spin(is)+3)/4
          p1=4*(p5-1)+1
          p2=4*p5 ; if ( p2 > nstate_spin(1) ) p2=nstate_spin(1)
          write(1,'(1x,4(i5,f15.4,2x))') (iob,energy%esp(iob,1,1)*2d0*Ry,iob=p1,p2)
        end do
      case(2)
        write(1,*) "for down-spin"
        do p5=1,(nstate_spin(is)+3)/4
          p1=4*(p5-1)+1+nstate_spin(1)
          p2=4*p5+nstate_spin(1) ; if ( p2 > nstate_spin(1)+nstate_spin(2) ) p2=nstate_spin(1)+nstate_spin(2)
          write(1,'(1x,4(i5,f15.4,2x))') (iob-nstate_spin(1),energy%esp(iob,1,1)*2d0*Ry,iob=p1,p2)
        end do
      end select
    end do
  end select
  write(1,*)

  do ii=1,ntmg
    write(1,'(1x,a,3f14.8)') "Size of the box (A) = ", rLsize(:,ii)*a_B
  end do

  write(1,'(1x,a,3f14.8)')   "Grid spacing (A)    = ", (Hgs(jj)*a_B,jj=1,3)
  write(1,*)
  write(1,'(1x,"Number of atoms = ",i8)') MI
  do ik=1,MKI
    write(1,'(1x,"iZatom(",i3,")     = ",i8)') ik, iZatom(ik)
  end do
  write(1,*)
  write(1,*) "Ref. and max angular momentum",      &
             " and pseudo-core radius of PP (A)"
  do ikoa=1,MKI
     write(1,'(1x,"(",i3,")  "," Ref, Max, Rps =",2i4,f8.3)')      &
                              ikoa,Lref(ikoa),Mlps(ikoa),Rps(ikoa)*a_B
  end do

  close(1)

end if

call timer_end(LOG_WRITE_LDA_INFOS)

deallocate(Vlocal)
call finalize_xc(xc_func)

call timer_end(LOG_TOTAL)

contains

subroutine band_information
  implicit none
  integer :: ik
  real(8),dimension(num_kpoints_rd) :: esp_vb_min,esp_vb_max,esp_cb_min,esp_cb_max
  if(comm_is_root(nproc_id_global) .and. itotfMST<itotMST) then
    do ik=1,num_kpoints_rd
      esp_vb_min(ik)=minval(energy%esp(1:itotfMST,ik,1))
      esp_vb_max(ik)=maxval(energy%esp(1:itotfMST,ik,1))
      esp_cb_min(ik)=minval(energy%esp(itotfMST+1:itotMST,ik,1))
      esp_cb_max(ik)=maxval(energy%esp(itotfMST+1:itotMST,ik,1))
    end do
    write(*,*) 'band information-----------------------------------------'
    write(*,*) 'Bottom of VB',minval(esp_vb_min(:))
    write(*,*) 'Top of VB',maxval(esp_vb_max(:))
    write(*,*) 'Bottom of CB',minval(esp_cb_min(:))
    write(*,*) 'Top of CB',maxval(esp_cb_max(:))
    write(*,*) 'Fundamental gap',minval(esp_cb_min(:))-maxval(esp_vb_max(:))
    write(*,*) 'BG between same k-point',minval(esp_cb_min(:)-esp_vb_max(:))
    write(*,*) 'Physicaly upper bound of CB for DOS',minval(esp_cb_max(:))
    write(*,*) 'Physicaly upper bound of CB for eps(omega)',minval(esp_cb_max(:)-esp_vb_min(:))
    write(*,*) '---------------------------------------------------------'
    write(*,*) 'Bottom of VB[eV]',minval(esp_vb_min(:))*2.0*Ry
    write(*,*) 'Top of VB[eV]',maxval(esp_vb_max(:))*2.0*Ry
    write(*,*) 'Bottom of CB[eV]',minval(esp_cb_min(:))*2.0*Ry
    write(*,*) 'Top of CB[eV]',maxval(esp_cb_max(:))*2.0*Ry
    write(*,*) 'Fundamental gap[eV]',(minval(esp_cb_min(:))-maxval(esp_vb_max(:)))*2.0*Ry
    write(*,*) 'BG between same k-point[eV]',(minval(esp_cb_min(:)-esp_vb_max(:)))*2.0*Ry
    write(*,*) '---------------------------------------------------------'
  end if
  return
end subroutine band_information

subroutine get_fourier_grid_G(fg)
  use structures, only: s_reciprocal_grid
  implicit none
  type(s_reciprocal_grid) :: fg

  if(allocated(fg%Gx))       deallocate(fg%Gx,fg%Gy,fg%Gz)
  if(allocated(fg%zrhoG_ion)) deallocate(fg%zrhoG_ion,fg%zrhoG_ele,fg%zdVG_ion)

  jj = system%ngrid/nproc_size_global
  fg%ig_s = nproc_id_global*jj+1
  fg%ig_e = (nproc_id_global+1)*jj
  if(nproc_id_global==nproc_size_global-1) fg%ig_e = system%ngrid
  fg%icomm_G = nproc_group_global
  fg%ng = system%ngrid
  allocate(fg%Gx(fg%ng),fg%Gy(fg%ng),fg%Gz(fg%ng))
  allocate(fg%zrhoG_ion(fg%ng),fg%zrhoG_ele(fg%ng),fg%zdVG_ion(fg%ng,nelem))
  if(iflag_hartree==2)then
     fg%iGzero = nGzero
     fg%Gx = Gx
     fg%Gy = Gy
     fg%Gz = Gz
     fg%zrhoG_ion = rhoion_G
     fg%zdVG_ion = dVloc_G
  else if(iflag_hartree==4)then
     fg%iGzero = 1
     fg%Gx = 0.d0
     fg%Gy = 0.d0
     fg%Gz = 0.d0
     fg%zrhoG_ion = 0.d0
     fg%zdVG_ion = 0.d0
     do iz=1,lg_num(3)/NPUZ
     do iy=1,lg_num(2)/NPUY
     do ix=ng%is(1)-lg%is(1)+1,ng%ie(1)-lg%is(1)+1
        n=(iz-1)*lg_num(2)/NPUY*lg_num(1)+(iy-1)*lg_num(1)+ix
        nn=ix-(ng%is(1)-lg%is(1)+1)+1+(iy-1)*ng%num(1)+(iz-1)*lg%num(2)/NPUY*ng%num(1)+fg%ig_s-1
        fg%Gx(nn) = Gx(n)
        fg%Gy(nn) = Gy(n)
        fg%Gz(nn) = Gz(n)
        fg%zrhoG_ion(nn) = rhoion_G(n)
        fg%zdVG_ion(nn,:) = dVloc_G(n,:)
     enddo
     enddo
     enddo
  end if

end subroutine get_fourier_grid_G

subroutine init_code_optimization
  implicit none
  integer :: ignum(3)

  call switch_stencil_optimization(mg%num)
  call switch_openmp_parallelization(mg%num)

  if(iperiodic==3 .and. nproc_d_o(1)*nproc_d_o(2)*nproc_d_o(3)==1) then
    ignum = mg%num
  else
    ignum = mg%num + (nd*2)
  end if
  call set_modulo_tables(ignum)

  if (comm_is_root(nproc_id_global)) then
    call optimization_log(nproc_k, nproc_ob, nproc_d_o, nproc_d_g)
  end if
end subroutine

END subroutine Real_Space_DFT

!=======================================================================
!========================================= Grid generation and labelling

SUBROUTINE init_mesh(lg,mg)
use structures, only: s_rgrid
use salmon_parallel, only: nproc_id_global, nproc_size_global
use salmon_communication, only: comm_is_root
use inputoutput, only: iperiodic
use global_variables_scf
implicit none
type(s_rgrid) :: lg
type(s_rgrid),intent(out) :: mg

if(comm_is_root(nproc_id_global))      &
    print *,"----------------------------------- init_mesh"

lg_sta(1:3) = lg%is(1:3)
lg_end(1:3) = lg%ie(1:3)
lg_num(1:3) = lg%num(1:3)
!call setlg(lg,lg_sta,lg_end,lg_num,ista_Mx_ori,iend_Mx_ori,inum_Mx_ori,    &
!           Hgs,Nd,rLsize1,imesh_oddeven,iperiodic)
call check_fourier

allocate(ista_Mxin(3,0:nproc_size_global-1),iend_Mxin(3,0:nproc_size_global-1))
allocate(inum_Mxin(3,0:nproc_size_global-1))

call setmg(mg,mg_sta,mg_end,mg_num,ista_Mxin,iend_Mxin,inum_Mxin,  &
           lg_sta,lg_num,nproc_size_global,nproc_id_global,nproc_d_o,nproc_k,nproc_ob,iscfrt)

if(comm_is_root(nproc_id_global)) write(*,*) "Mx     =", iend_Mx_ori

if(iperiodic==3 .and. nproc_d_o(1)*nproc_d_o(2)*nproc_d_o(3)==1) then
  if(comm_is_root(nproc_id_global)) write(*,*) "r-space parallelization: off"
  mg%is(1:3)=lg%is(1:3)
  mg%ie(1:3)=lg%ie(1:3)
  mg%num(1:3)=lg%num(1:3)
  mg%is_overlap(1:3)=lg%is_overlap(1:3)
  mg%ie_overlap(1:3)=lg%ie_overlap(1:3)
  mg%is_array(1:3)=lg%is_array(1:3)
  mg%ie_array(1:3)=lg%ie_array(1:3)
  if(allocated(mg%idx)) deallocate(mg%idx)
  if(allocated(mg%idy)) deallocate(mg%idy)
  if(allocated(mg%idz)) deallocate(mg%idz)
  allocate(mg%idx(mg%is_overlap(1):mg%ie_overlap(1)) &
          ,mg%idy(mg%is_overlap(2):mg%ie_overlap(2)) &
          ,mg%idz(mg%is_overlap(3):mg%ie_overlap(3)))
  mg%idx = lg%idx
  mg%idy = lg%idy
  mg%idz = lg%idz
end if

return

END SUBROUTINE init_mesh

