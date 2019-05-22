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
use copy_psi_mesh_sub
use calc_density_sub
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
use structures!, only: s_rgrid, s_wf_info, s_wavefunction, s_system, s_stencil
use salmon_parallel, only: nproc_id_global, nproc_size_global, nproc_group_global, &
                           nproc_group_h, nproc_id_kgrid, nproc_id_orbitalgrid, &
                           nproc_group_korbital, nproc_id_korbital, nproc_group_rho, &
                           nproc_group_kgrid, nproc_group_k
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use salmon_xc, only: init_xc, finalize_xc
use timer
use write_performance_results, only: write_gs_performance
use iso_fortran_env, only: output_unit
use calc_iobnum_sub
use check_mg_sub
use check_ng_sub
use dtcg_sub
use gscg_sub
use dtcg_periodic_sub
use gscg_periodic_sub
use rmmdiis_sub
use subspace_diag_sub
use subspace_diag_periodic_sub
use writefield
use global_variables_scf
use lattice
use sendrecv_grid, only: s_sendrecv_grid, init_sendrecv_grid
use salmon_pp, only: calc_nlcc
use force_sub
use calc_iroot_sub
use gram_schmidt_orth, only: gram_schmidt 
use md_ground_state, only: write_xyz   !temporary
implicit none

integer :: ix,iy,iz,ik,ikoa, is
integer :: iter,iatom,iob,p1,p2,p5,ii,jj,iflag,jspin
real(8) :: sum0,sum1
character(100) :: file_atoms_coo, comment_line
complex(8),allocatable :: zpsi_tmp(:,:,:,:,:)
real(8) :: rNebox1,rNebox2
integer :: nspin,n,nn,itmg
integer :: neig(1:3, 1:2)
integer :: neig_ng(1:3, 1:2)

type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_rgrid) :: ng
type(s_wf_info) :: info_ob
type(s_wf_info) :: info
type(s_sendrecv_grid) :: srg, srg_ob_1, srg_ob, srg_ng
type(s_wavefunction) :: spsi,shpsi,sttpsi
type(s_system) :: system
type(s_stencil) :: stencil
type(s_scalar) :: sVh
type(s_scalar),allocatable :: V_local(:),srho(:),sVxc(:)
type(s_fourier_grid) :: fg
type(s_pp_nlcc) :: ppn
type(s_energy) :: energy
type(s_force)  :: force


call init_xc(xc_func, ispin, cval, xcname=xc, xname=xname, cname=cname)

iSCFRT=1
ihpsieff=0
iflag_comm_rho=1

iblacsinit=0

call timer_begin(LOG_TOTAL)


call timer_begin(LOG_INIT_GS)
inumcpu_check=0

call setbN
call setcN

call check_dos_pdos

call convert_input_scf(file_atoms_coo)

if(al_vec1(2)==0d0 .and. al_vec1(3)==0d0 .and. al_vec2(1)==0d0 .and. &
   al_vec2(3)==0d0 .and. al_vec3(1)==0d0 .and. al_vec3(2)==0d0) then
  if(comm_is_root(nproc_id_global)) write(*,*) "orthogonal cell: using al"
  stencil%if_orthogonal = .true.
  system%al = 0d0
  system%al(1,1) = al(1)
  system%al(2,2) = al(2)
  system%al(3,3) = al(3)
else
  if(comm_is_root(nproc_id_global)) write(*,*) "non-orthogonal cell: using al_vec[1,2,3]"
  stencil%if_orthogonal = .false.
  system%al(1:3,1) = al_vec1
  system%al(1:3,2) = al_vec2
  system%al(1:3,3) = al_vec3
  rLsize(1,1) = sqrt(sum(al_vec1**2))
  rLsize(2,1) = sqrt(sum(al_vec2**2))
  rLsize(3,1) = sqrt(sum(al_vec3**2))
  if(sum(abs(num_rgrid)) /= 0 .and. sum(abs(dl)) == 0d0) then
    Harray(1:3,1) = rLsize(1:3,1) / dble(num_rgrid(1:3))
  end if
end if

call set_filename

call setk(k_sta, k_end, k_num, num_kpoints_rd, nproc_k, nproc_id_orbitalgrid)

call calc_iobnum(itotMST,nproc_ob,nproc_id_kgrid,iobnum,nproc_ob,iparaway_ob)

if(iflag_opt==1)then
  call structure_opt_ini(MI)
  iopt_tranc=0
end if
call timer_end(LOG_INIT_GS)


Structure_Optimization_Iteration : do iopt=1,iter_opt
Multigrid_Iteration : do img=1,ntmg

if(iopt==1)then
  select case( IC )

!------------------------------ New calculation

  case default

    call timer_begin(LOG_INIT_GS)

    Hgs(1:3)=Harray(1:3,1)
    Hvol=Hgs(1)*Hgs(2)*Hgs(3)
    Miter = 0        ! Miter: Iteration counter set to zero
    itmg=img
    call set_imesh_oddeven(itmg)
    call init_mesh(lg,mg,itmg)
    call set_gridcoo
    call init_mesh_s(ng)
    call check_mg(mg)
    call check_ng(ng)
    lg%ndir = 3
    mg%ndir = 3
    ng%ndir = 3
    system%ngrid = lg_num(1)*lg_num(2)*lg_num(3)

    call init_lattice(system,stencil,lg)
    Hvol = system%Hvol
    Hgs = system%Hgs
    
    call init_updown
    call init_itype
    call init_sendrecv_matrix
    select case(iperiodic)
    case(0)
      if(MEO==2.or.MEO==3) call make_corr_pole
    end select
    call make_icoobox_bound
        
    call allocate_mat
    call set_icoo1d
    call allocate_sendrecv
    call init_persistent_requests

    ! Initialization of s_sendrecv_grid structure (experimental implementation)
    neig(1, 1) = iup_array(1)
    neig(1, 2) = idw_array(1)
    neig(2, 1) = jup_array(1)
    neig(2, 2) = jdw_array(1)
    neig(3, 1) = kup_array(1)
    neig(3, 2) = kdw_array(1)
    call init_sendrecv_grid(srg, mg, iobnum * k_num, &
      & nproc_group_korbital, nproc_id_korbital, neig)
    call init_sendrecv_grid(srg_ob, mg, nspin, &
      & nproc_group_korbital, nproc_id_korbital, neig)
    call init_sendrecv_grid(srg_ob_1, mg, 1, &
      & nproc_group_korbital, nproc_id_korbital, neig)

    neig_ng(1, 1) = iup_array(2)
    neig_ng(1, 2) = idw_array(2)
    neig_ng(2, 1) = jup_array(2)
    neig_ng(2, 2) = jdw_array(2)
    neig_ng(3, 1) = kup_array(2)
    neig_ng(3, 2) = kdw_array(2)
    call init_sendrecv_grid(srg_ng, ng, 1, &
      & nproc_group_global, nproc_id_global, neig_ng)

    if(ispin==0)then
      nspin=1
    else
      nspin=2
    end if
    
    system%iperiodic = iperiodic
    system%ngrid = lg_num(1)*lg_num(2)*lg_num(3)
    system%nspin = nspin
    system%no = itotMST
    system%nk = num_kpoints_rd
    system%nion = MI
    
    allocate(system%Rion(3,system%nion) &
            ,system%wtk(system%nk) &
            ,system%rocc(system%no,system%nk,system%nspin))
    system%wtk = wtk
    system%rion = rion
    
    allocate(energy%esp(system%no,system%nk,system%nspin))

    info%im_s = 1
    info%im_e = 1
    info%numm = 1
    info%ik_s=k_sta
    info%ik_e=k_end
    info%numk=k_num
    info%io_s=1
    info%io_e=iobnum/nspin
    info%numo=iobnum/nspin
    
    info%if_divide_rspace = nproc_mxin_mul.ne.1
    info%irank_r(1) = iup_array(1)
    info%irank_r(2) = idw_array(1)
    info%irank_r(3) = jup_array(1)
    info%irank_r(4) = jdw_array(1)
    info%irank_r(5) = kup_array(1)
    info%irank_r(6) = kdw_array(1)
    info%icomm_r = nproc_group_korbital
    info%icomm_o = nproc_group_kgrid
    info%icomm_ko = nproc_group_rho
    info%icomm_ro = nproc_group_k
    info%icomm_rko = nproc_group_global
    allocate(info%occ(info%io_s:info%io_e, info%ik_s:info%ik_e, 1:system%nspin) &
              ,info%io_tbl(info%io_s:info%io_e), info%jo_tbl(1:system%no) &
              ,info%irank_jo(1:system%no))
 
    info%jo_tbl(:) = 0 !(initial value)
    do iob=info%io_s,info%io_e
      call calc_allob(iob,jj,iparaway_ob,itotmst,mst,iobnum)
      info%io_tbl(iob) = jj
      info%jo_tbl(jj) = iob
    end do
    
    do jj=1, system%no
      call calc_iroot(jj,info%irank_jo(jj),ilsda,nproc_ob,iparaway_ob,itotmst,mst)
    end do
    
    if(stencil%if_orthogonal) then
      stencil%lap0 = -0.5d0*cNmat(0,Nd)*(1.d0/Hgs(1)**2+1.d0/Hgs(2)**2+1.d0/Hgs(3)**2)
    else
      if(info%if_divide_rspace) stop "error: nonorthogonal lattice and r-space parallelization"
      stencil%lap0 = -0.5d0*cNmat(0,Nd)*( stencil%coef_F(1)/Hgs(1)**2 + stencil%coef_F(2)/Hgs(2)**2 + stencil%coef_F(3)/Hgs(3)**2 )
    end if
    do jj=1,3
      do ii=1,4
        stencil%lapt(ii,jj) = cnmat(ii,4)/hgs(jj)**2
        stencil%nabt(ii,jj) = bnmat(ii,4)/hgs(jj)
      end do
    end do
    
    info_ob%im_s = 1
    info_ob%im_e = 1
    info_ob%numm = 1
    info_ob%ik_s = 1
    info_ob%ik_e = 1
    info_ob%numk = 1
    info_ob%io_s = 1
    info_ob%io_e = 1
    info_ob%numo = 1
    info_ob%if_divide_rspace = nproc_mxin_mul.ne.1
    info_ob%irank_r(1) = iup_array(1)
    info_ob%irank_r(2) = idw_array(1)
    info_ob%irank_r(3) = jup_array(1)
    info_ob%irank_r(4) = jdw_array(1)
    info_ob%irank_r(5) = kup_array(1)
    info_ob%irank_r(6) = kdw_array(1)
    info_ob%icomm_r = nproc_group_korbital
    
    allocate(V_local(system%nspin),srho(system%nspin),sVxc(system%nspin))
    do jspin=1,system%nspin
      allocate(V_local(jspin)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
      allocate(srho(jspin)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
      allocate(sVxc(jspin)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    end do
    allocate(sVh%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    
    select case(iperiodic)
    case(0)
      allocate(spsi%rwf(mg%is_array(1):mg%ie_array(1),  &
                        mg%is_array(2):mg%ie_array(2),  &
                        mg%is_array(3):mg%ie_array(3),  &
                        nspin,  &
                        info%io_s:info%io_e,  &
                        info%ik_s:info%ik_e,  &
                        1))
      allocate(shpsi%rwf(mg%is_array(1):mg%ie_array(1),  &
                         mg%is_array(2):mg%ie_array(2),  &
                         mg%is_array(3):mg%ie_array(3),  &
                         nspin,  &
                         info%io_s:info%io_e,  &
                         info%ik_s:info%ik_e,  &
                         1))
    !$OMP parallel do private(ik,iob,iz,iy,ix) collapse(4)
      do ik=info%ik_s,info%ik_e
      do iob=info%io_s,info%io_e
        do is=1,nspin
          do iz=mg%is_array(3),mg%ie_array(3)
          do iy=mg%is_array(2),mg%ie_array(2)
          do ix=mg%is_array(1),mg%ie_array(1)
            spsi%rwf(ix,iy,iz,is,iob,ik,1)=0.d0
          end do
          end do
          end do
        end do
      end do
      end do
    case(3)
      allocate(spsi%zwf(mg%is_array(1):mg%ie_array(1),  &
                        mg%is_array(2):mg%ie_array(2),  &
                        mg%is_array(3):mg%ie_array(3),  &
                        nspin,  &
                        info%io_s:info%io_e,  &
                        info%ik_s:info%ik_e,  &
                        1))
      allocate(shpsi%zwf(mg%is_array(1):mg%ie_array(1),  &
                         mg%is_array(2):mg%ie_array(2),  &
                         mg%is_array(3):mg%ie_array(3),  &
                         nspin,  &
                         info%io_s:info%io_e,  &
                         info%ik_s:info%ik_e,  &
                         1))
      allocate(sttpsi%zwf(mg%is_array(1):mg%ie_array(1),  &
                         mg%is_array(2):mg%ie_array(2),  &
                         mg%is_array(3):mg%ie_array(3),  &
                         nspin,  &
                         info%io_s:info%io_e,  &
                         info%ik_s:info%ik_e,  &
                         1))
    !$OMP parallel do private(ik,iob,iz,iy,ix) collapse(4)
      do ik=info%ik_s,info%ik_e
      do iob=info%io_s,info%io_e
        do is=1,nspin
          do iz=mg%is_array(3),mg%ie_array(3)
          do iy=mg%is_array(2),mg%ie_array(2)
          do ix=mg%is_array(1),mg%ie_array(1)
            spsi%zwf(ix,iy,iz,is,iob,ik,1)=0.d0
          end do
          end do
          end do
        end do
      end do
      end do
    end select

    if(iperiodic==3)then
      allocate (zpsi_tmp(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
                 1:iobnum,k_sta:k_end))
      allocate(k_rd(3,num_kpoints_rd),ksquare(num_kpoints_rd))
      call init_k_rd(k_rd,ksquare,1,system%brl)
    end if

    allocate( Vpsl(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )
    if(icalcforce==1)then
      allocate( Vpsl_atom(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),MI) )
    end if
    
    if(iperiodic==3 .and. iflag_hartree==4)then
      call prep_poisson_fft
    end if

    if(iflag_ps.eq.0)then
      Vpsl=0d0
    else
      call read_pslfile
      call allocate_psl
      call init_ps(system%al,system%brl,stencil%matrix_A)
    end if

    if(iperiodic==3)then
      jj = system%ngrid/nproc_size_global
      fg%ig_s = nproc_id_global*jj+1
      fg%ig_e = (nproc_id_global+1)*jj
      if(nproc_id_global==nproc_size_global-1) fg%ig_e = system%ngrid
      fg%icomm_fourier = nproc_group_global
      fg%ng = system%ngrid
      allocate(fg%Gx(fg%ng),fg%Gy(fg%ng),fg%Gz(fg%ng))
      allocate(fg%rhoG_ion(fg%ng),fg%rhoG_elec(fg%ng),fg%dVG_ion(fg%ng,nelem))
      if(iflag_hartree==2)then
        fg%iGzero = nGzero
        fg%Gx = Gx
        fg%Gy = Gy
        fg%Gz = Gz
        fg%rhoG_ion = rhoion_G
        fg%dVG_ion = dVloc_G
      else if(iflag_hartree==4)then
        fg%iGzero = 1
        fg%Gx = 0.d0
        fg%Gy = 0.d0
        fg%Gz = 0.d0
        fg%rhoG_ion = 0.d0
        fg%dVG_ion = 0.d0
        do iz=1,lg_num(3)/NPUZ
        do iy=1,lg_num(2)/NPUY
        do ix=ng%is(1)-lg%is(1)+1,ng%ie(1)-lg%is(1)+1
          n=(iz-1)*lg_num(2)/NPUY*lg_num(1)+(iy-1)*lg_num(1)+ix
          nn=ix-(ng%is(1)-lg%is(1)+1)+1+(iy-1)*ng%num(1)+(iz-1)*lg%num(2)/NPUY*ng%num(1)+fg%ig_s-1
          fg%Gx(nn) = Gx(n)
          fg%Gy(nn) = Gy(n)
          fg%Gz(nn) = Gz(n)
          fg%rhoG_ion(nn) = rhoion_G(n)
          fg%dVG_ion(nn,:) = dVloc_G(n,:)
        enddo
        enddo
        enddo
      end if
    end if

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
        allocate( psi_mesh(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:itotMST,1) ) 
      case(3)
        allocate( zpsi_mesh(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:itotMST,num_kpoints_rd) )
      end select
    end if

    call init_wf_ns(1)
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
            spsi%rwf(ix,iy,iz,is,iob,ik,1)=psi(ix,iy,iz,iob+(is-1)*info%numo,ik)
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
            spsi%zwf(ix,iy,iz,is,iob,ik,1)=zpsi(ix,iy,iz,iob+(is-1)*info%numo,ik)
          end do
          end do
          end do
        end do
      end do
    end do
    end select 

    
    call timer_begin(LOG_CALC_GRAM_SCHMIDT)
    call gram_schmidt(system, mg, info, spsi)
    call timer_end(LOG_CALC_GRAM_SCHMIDT)


    allocate( rho(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )  
    allocate( rho_in(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:num_rho_stock+1) )  
    allocate( rho_out(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:num_rho_stock+1) ) 
    rho_in=0.d0
    rho_out=0.d0
                                
    if(ilsda == 1)then
      allocate( rho_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2) )  
      allocate( rho_s_in(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:num_rho_stock+1,2) )  
      allocate( rho_s_out(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:num_rho_stock+1,2) )  
      rho_s_in=0.d0
      rho_s_out=0.d0
    end if
    rho=0.d0 


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

    select case(iperiodic)
    case(0)
      call calc_density(psi)
    case(3)
      call calc_density(zpsi)
    end select

    if(ilsda==0)then
      allocate (Vlocal(mg_sta(1):mg_end(1),  &
              mg_sta(2):mg_end(2),  &
              mg_sta(3):mg_end(3),1))
    else if(ilsda==1)then
      allocate (Vlocal(mg_sta(1):mg_end(1),  &
              mg_sta(2):mg_end(2),  &
              mg_sta(3):mg_end(3),2))
    end if

    allocate( Vh(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )  
    Vh=0.d0

    call Hartree_ns(lg,mg,ng,system%Brl,srg_ng,stencil)
    
    if(ilsda == 0) then
      allocate( Vxc(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)) )  
    else if(ilsda == 1) then
      allocate( Vxc_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),2) )  
    end if
    allocate( esp(itotMST,num_kpoints_rd) )

    call exc_cor_ns(ppn)

    call allgatherv_vlocal

    if(iperiodic==3) then
      allocate(stencil%kAc(k_sta:k_end,3))
      do jj=1,3
        stencil%kAc(k_sta:k_end,jj) = k_rd(jj,k_sta:k_end)
      end do
      call update_kvector_nonlocalpt(ppg,stencil%kAc,k_sta,k_end)
    end if
    do jspin=1,system%nspin
      V_local(jspin)%f = Vlocal(:,:,:,jspin)
    end do
    if(ilsda == 1) then
      do jspin=1,system%nspin
        srho(jspin)%f = rho_s(:,:,:,jspin)
        sVxc(jspin)%f = Vxc_s(:,:,:,jspin)
      end do
    else
      srho(1)%f = rho
      sVxc(1)%f = Vxc
    end if
    sVh%f = Vh
    energy%E_xc = Exc
    call calc_eigen_energy(energy,spsi,shpsi,sttpsi,system,info,mg,V_local,stencil,srg,ppg)
    select case(iperiodic)
    case(0)
      call calc_Total_Energy_isolated(energy,system,info,ng,pp,srho,sVh,sVxc)
    case(3)
      if(iflag_hartree==2)then
        fg%rhoG_elec = rhoe_G
      else if(iflag_hartree==4)then
        do iz=1,lg_num(3)/NPUZ
        do iy=1,lg_num(2)/NPUY
        do ix=ng%is(1)-lg%is(1)+1,ng%ie(1)-lg%is(1)+1
          n=(iz-1)*lg_num(2)/NPUY*lg_num(1)+(iy-1)*lg_num(1)+ix
          nn=ix-(ng%is(1)-lg%is(1)+1)+1+(iy-1)*ng%num(1)+(iz-1)*lg%num(2)/NPUY*ng%num(1)+fg%ig_s-1
          fg%rhoG_elec(nn) = rhoe_G(n)
        enddo
        enddo
        enddo
      end if
      call calc_Total_Energy_periodic(energy,system,pp,fg)
    end select
    esp = energy%esp(:,:,1) !++++++++
    if(iperiodic==3) deallocate(stencil%kAc,ppg%ekr_uV)

    call timer_end(LOG_INIT_GS)
        
!------------------------------ Continue the previous calculation

  case(1,3)
    call timer_begin(LOG_INIT_GS_RESTART)

    call IN_data(lg,mg,ng,system,info,stencil)

    call allocate_mat
    call set_icoo1d
    call allocate_sendrecv
    call init_persistent_requests

    ! Initialization of s_sendrecv_grid structure (experimental implementation)
    neig(1, 1) = iup_array(1)
    neig(1, 2) = idw_array(1)
    neig(2, 1) = jup_array(1)
    neig(2, 2) = jdw_array(1)
    neig(3, 1) = kup_array(1)
    neig(3, 2) = kdw_array(1)
    call init_sendrecv_grid(srg, mg, iobnum * k_num, &
      & nproc_group_korbital, nproc_id_korbital, neig)
    call init_sendrecv_grid(srg_ob, mg, nspin, &
      & nproc_group_korbital, nproc_id_korbital, neig)
    call init_sendrecv_grid(srg_ob_1, mg, 1, &
      & nproc_group_korbital, nproc_id_korbital, neig)

    neig_ng(1, 1) = iup_array(2)
    neig_ng(1, 2) = idw_array(2)
    neig_ng(2, 1) = jup_array(2)
    neig_ng(2, 2) = jdw_array(2)
    neig_ng(3, 1) = kup_array(2)
    neig_ng(3, 2) = kdw_array(2)
    call init_sendrecv_grid(srg_ng, ng, 1, &
      & nproc_group_global, nproc_id_global, neig_ng)

    if(ispin==0)then
      nspin=1
    else
      nspin=2
    end if
    
    system%iperiodic = iperiodic
    system%ngrid = lg_num(1)*lg_num(2)*lg_num(3)
    system%nspin = nspin
    system%no = itotMST
    system%nk = num_kpoints_rd
    system%nion = MI
    
    allocate(system%Rion(3,system%nion) &
            ,system%wtk(system%nk) &
            ,system%rocc(system%no,system%nk,system%nspin))
    system%wtk = wtk
    system%rion = rion
    
    allocate(energy%esp(system%no,system%nk,system%nspin))

    info%im_s = 1
    info%im_e = 1
    info%numm = 1
    info%ik_s=k_sta
    info%ik_e=k_end
    info%numk=k_num
    info%io_s=1
    info%io_e=iobnum/nspin
    info%numo=iobnum/nspin
    
    info%if_divide_rspace = nproc_mxin_mul.ne.1
    info%irank_r(1) = iup_array(1)
    info%irank_r(2) = idw_array(1)
    info%irank_r(3) = jup_array(1)
    info%irank_r(4) = jdw_array(1)
    info%irank_r(5) = kup_array(1)
    info%irank_r(6) = kdw_array(1)
    info%icomm_r = nproc_group_korbital
    info%icomm_o = nproc_group_kgrid
    info%icomm_ko = nproc_group_rho
    info%icomm_ro = nproc_group_k
    info%icomm_rko = nproc_group_global
    allocate(info%occ(info%io_s:info%io_e, info%ik_s:info%ik_e, 1:system%nspin) &
              ,info%io_tbl(info%io_s:info%io_e), info%jo_tbl(1:system%no) &
              ,info%irank_jo(1:system%no))
 
    info%jo_tbl(:) = 0 !(initial value)
    do iob=info%io_s,info%io_e
      call calc_allob(iob,jj,iparaway_ob,itotmst,mst,iobnum)
      info%io_tbl(iob) = jj
      info%jo_tbl(jj) = iob
    end do
    
    do jj=1, system%no
      call calc_iroot(jj,info%irank_jo(jj),ilsda,nproc_ob,iparaway_ob,itotmst,mst)
    end do
    
    if(stencil%if_orthogonal) then
      stencil%lap0 = -0.5d0*cNmat(0,Nd)*(1.d0/Hgs(1)**2+1.d0/Hgs(2)**2+1.d0/Hgs(3)**2)
    else
      if(info%if_divide_rspace) stop "error: nonorthogonal lattice and r-space parallelization"
      stencil%lap0 = -0.5d0*cNmat(0,Nd)*( stencil%coef_F(1)/Hgs(1)**2 + stencil%coef_F(2)/Hgs(2)**2 + stencil%coef_F(3)/Hgs(3)**2 )
    end if
    do jj=1,3
      do ii=1,4
        stencil%lapt(ii,jj) = cnmat(ii,4)/hgs(jj)**2
        stencil%nabt(ii,jj) = bnmat(ii,4)/hgs(jj)
      end do
    end do
    
    info_ob%im_s = 1
    info_ob%im_e = 1
    info_ob%numm = 1
    info_ob%ik_s = 1
    info_ob%ik_e = 1
    info_ob%numk = 1
    info_ob%io_s = 1
    info_ob%io_e = 1
    info_ob%numo = 1
    info_ob%if_divide_rspace = nproc_mxin_mul.ne.1
    info_ob%irank_r(1) = iup_array(1)
    info_ob%irank_r(2) = idw_array(1)
    info_ob%irank_r(3) = jup_array(1)
    info_ob%irank_r(4) = jdw_array(1)
    info_ob%irank_r(5) = kup_array(1)
    info_ob%irank_r(6) = kdw_array(1)
    info_ob%icomm_r = nproc_group_korbital
    
    allocate(V_local(system%nspin),srho(system%nspin),sVxc(system%nspin))
    do jspin=1,system%nspin
      allocate(V_local(jspin)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
      allocate(srho(jspin)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
      allocate(sVxc(jspin)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    end do
    allocate(sVh%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))

    select case(iperiodic)
    case(0)
      allocate(spsi%rwf(mg%is_array(1):mg%ie_array(1),  &
                        mg%is_array(2):mg%ie_array(2),  &
                        mg%is_array(3):mg%ie_array(3),  &
                        nspin,  &
                        info%io_s:info%io_e,  &
                        info%ik_s:info%ik_e,  &
                        1))
      allocate(shpsi%rwf(mg%is_array(1):mg%ie_array(1),  &
                         mg%is_array(2):mg%ie_array(2),  &
                         mg%is_array(3):mg%ie_array(3),  &
                         nspin,  &
                         info%io_s:info%io_e,  &
                         info%ik_s:info%ik_e,  &
                         1))
    !$OMP parallel do private(ik,iob,iz,iy,ix) collapse(4)
      do ik=info%ik_s,info%ik_e
      do iob=info%io_s,info%io_e
        do is=1,nspin
          do iz=mg%is_array(3),mg%ie_array(3)
          do iy=mg%is_array(2),mg%ie_array(2)
          do ix=mg%is_array(1),mg%ie_array(1)
            spsi%rwf(ix,iy,iz,is,iob,ik,1)=0.d0
          end do
          end do
          end do
        end do
      end do
      end do
    case(3)
      allocate(spsi%zwf(mg%is_array(1):mg%ie_array(1),  &
                        mg%is_array(2):mg%ie_array(2),  &
                        mg%is_array(3):mg%ie_array(3),  &
                        nspin,  &
                        info%io_s:info%io_e,  &
                        info%ik_s:info%ik_e,  &
                        1))
      allocate(shpsi%zwf(mg%is_array(1):mg%ie_array(1),  &
                         mg%is_array(2):mg%ie_array(2),  &
                         mg%is_array(3):mg%ie_array(3),  &
                         nspin,  &
                         info%io_s:info%io_e,  &
                         info%ik_s:info%ik_e,  &
                         1))
      allocate(sttpsi%zwf(mg%is_array(1):mg%ie_array(1),  &
                         mg%is_array(2):mg%ie_array(2),  &
                         mg%is_array(3):mg%ie_array(3),  &
                         nspin,  &
                         info%io_s:info%io_e,  &
                         info%ik_s:info%ik_e,  &
                         1))
    !$OMP parallel do private(ik,iob,iz,iy,ix) collapse(4)
      do ik=info%ik_s,info%ik_e
      do iob=info%io_s,info%io_e
        do is=1,nspin
          do iz=mg%is_array(3),mg%ie_array(3)
          do iy=mg%is_array(2),mg%ie_array(2)
          do ix=mg%is_array(1),mg%ie_array(1)
            spsi%zwf(ix,iy,iz,is,iob,ik,1)=0.d0
          end do
          end do
          end do
        end do
      end do
      end do
    end select

    if(iperiodic==3)then
      allocate (zpsi_tmp(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
                 1:iobnum,k_sta:k_end))
      allocate(k_rd(3,num_kpoints_rd),ksquare(num_kpoints_rd))
      call init_k_rd(k_rd,ksquare,1)
    end if

    if(iperiodic==3)then
      call prep_poisson_fft
    end if

    if(iflag_ps/=0) then
      call read_pslfile
      call allocate_psl
      call init_ps(system%al,system%brl,stencil%matrix_A)
    end if

    if(iperiodic==3)then
      jj = system%ngrid/nproc_size_global
      fg%ig_s = nproc_id_global*jj+1
      fg%ig_e = (nproc_id_global+1)*jj
      if(nproc_id_global==nproc_size_global-1) fg%ig_e = system%ngrid
      fg%icomm_fourier = nproc_group_global
      fg%ng = system%ngrid
      allocate(fg%Gx(fg%ng),fg%Gy(fg%ng),fg%Gz(fg%ng))
      allocate(fg%rhoG_ion(fg%ng),fg%rhoG_elec(fg%ng),fg%dVG_ion(fg%ng,nelem))
      if(iflag_hartree==2)then
        fg%iGzero = nGzero
        fg%Gx = Gx
        fg%Gy = Gy
        fg%Gz = Gz
        fg%rhoG_ion = rhoion_G
        fg%dVG_ion = dVloc_G
      else if(iflag_hartree==4)then
        fg%iGzero = 1
        fg%Gx = 0.d0
        fg%Gy = 0.d0
        fg%Gz = 0.d0
        fg%rhoG_ion = 0.d0
        fg%dVG_ion = 0.d0
        do iz=1,lg_num(3)/NPUZ
        do iy=1,lg_num(2)/NPUY
        do ix=ng%is(1)-lg%is(1)+1,ng%ie(1)-lg%is(1)+1
          n=(iz-1)*lg_num(2)/NPUY*lg_num(1)+(iy-1)*lg_num(1)+ix
          nn=ix-(ng%is(1)-lg%is(1)+1)+1+(iy-1)*ng%num(1)+(iz-1)*lg%num(2)/NPUY*ng%num(1)+fg%ig_s-1
          fg%Gx(nn) = Gx(n)
          fg%Gy(nn) = Gy(n)
          fg%Gz(nn) = Gz(n)
          fg%rhoG_ion(nn) = rhoion_G(n)
          fg%dVG_ion(nn,:) = dVloc_G(n,:)
        enddo
        enddo
        enddo
      end if
    end if

    call init_updown
    call init_itype
    call init_sendrecv_matrix
    if(MEO==2.or.MEO==3) call make_corr_pole
    call make_icoobox_bound
  end select

  call timer_end(LOG_INIT_GS_RESTART)

else if(iopt>=2)then
  call timer_begin(LOG_INIT_GS)
  Miter = 0        ! Miter: Iteration counter set to zero
  if(iflag_ps/=0) then
    call dealloc_init_ps(ppg,ppg_all,ppn)
    call init_ps(system%al,system%brl,stencil%matrix_A)
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
allocate(rho_stock(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1))
if(ilsda==0)then
  allocate(Vlocal_stock(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1))
else
  allocate(Vlocal_stock(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),1:2))
end if

if(ilsda==0)then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    rho_stock(ix,iy,iz,1)=rho(ix,iy,iz)
    Vlocal_stock(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)
  end do
  end do
  end do
else
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
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
  select case(convergence)
    case('rho_dne')
      if(sum1<threshold) cycle DFT_Iteration
    case('norm_rho','norm_rho_dng')
      if(sum1<threshold_norm_rho) cycle DFT_Iteration
    case('norm_pot','norm_pot_dng')
      if(sum1<threshold_norm_pot) cycle DFT_Iteration
  end select 

  Miter=Miter+1

  if(temperature_k>=0.d0.and.Miter>iditer_notemperature) then
    if(iperiodic.eq.3) then
      call ne2mu_p
    else
      call ne2mu
    endif
  else
    call calc_occupation
  endif

  system%rocc(:,:,1) = rocc

  do jspin=1,system%nspin
    do ik=info%ik_s,info%ik_e
      do iob=info%io_s,info%io_e
        jj = info%io_tbl(iob)+(jspin-1)*mst(1)
        info%occ(iob,ik,jspin) = system%rocc(jj,ik,jspin)*system%wtk(ik)
      end do
    end do
  end do

  call copy_density



  if(iscf_order==1)then
    call timer_begin(LOG_CALC_MINIMIZATION)
    if( amin_routine == 'cg' .or.       &
      ( amin_routine == 'cg-diis' .and. Miter <= iDiterYBCG) ) then
      select case(iperiodic)
      case(0)
        select case(gscg)
        case('y')
          call sgscg(mg,nspin,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,hvol,ilsda,nproc_ob,iparaway_ob, &
                     rxk_ob,rhxk_ob,rgk_ob,rpk_ob,   &
                     info_ob,bnmat,cnmat,hgs,ppg,vlocal)
        case('n')
          call dtcg(mg,nspin,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,hvol,ilsda,nproc_ob,iparaway_ob,   &
                    info_ob,bnmat,cnmat,hgs,ppg,vlocal)
        end select
      case(3)
        select case(gscg)
        case('y')
          call gscg_periodic(mg,nspin,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,hvol,ilsda,nproc_ob,iparaway_ob,   &
                             zxk_ob,zhxk_ob,zgk_ob,zpk_ob,zpko_ob,zhtpsi_ob,  &
                             info_ob,bnmat,cnmat,hgs,ppg,vlocal,num_kpoints_rd,k_rd)
        case('n')
          call dtcg_periodic(mg,nspin,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,hvol,ilsda,nproc_ob,iparaway_ob,   &
                             info_ob,bnmat,cnmat,hgs,ppg,vlocal,num_kpoints_rd,k_rd)
        end select
      end select
    else if( amin_routine  == 'diis' .or. amin_routine == 'cg-diis' ) then
      select case(iperiodic)
      case(0)
        call rmmdiis(mg,nspin,info,stencil,srg_ob_1,spsi,itotmst,mst,num_kpoints_rd,hvol,iflag_diisjump,energy%esp, &
                     norm_diff_psi_stock,info_ob,bnmat,cnmat,hgs,ppg,vlocal,iparaway_ob)
      case(3)
        stop "rmmdiis method is not implemented for periodic systems."
      end select
    end if
    call timer_end(LOG_CALC_MINIMIZATION)

    call timer_begin(LOG_CALC_GRAM_SCHMIDT)
    call gram_schmidt(system, mg, info, spsi)
    call timer_end(LOG_CALC_GRAM_SCHMIDT)




    !call debug_var_dump(system, mg, info, spsi, iter)  !uemoto!
    if(iflag_subspace_diag==1)then
      if(Miter>iDiter_nosubspace_diag)then
        select case(iperiodic)
        case(0)      
          call subspace_diag(mg,info,stencil,srg_ob_1,spsi,ilsda,nproc_ob,iparaway_ob,iobnum,itotmst,k_sta,k_end,mst,ifmst,hvol,  &
                info_ob,bnmat,cnmat,hgs,ppg,vlocal)

        case(3)

          call subspace_diag_periodic(mg,info,stencil,srg_ob_1,spsi,ilsda,nproc_ob,iparaway_ob,  &
                                      iobnum,itotmst,k_sta,k_end,mst,ifmst,hvol,   &
                                      info_ob,bnmat,cnmat,hgs,ppg,vlocal,num_kpoints_rd,k_rd)
        end select
      end if
    end if
  
    ! Store to psi/zpsi for calc_density
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

    call timer_begin(LOG_CALC_RHO)
    select case(iperiodic)
    case(0)
      call calc_density(psi)
    case(3)
      call calc_density(zpsi)
    end select

    select case(amixing)
      case ('simple')
        call simple_mixing(1.d0-rmixrate,rmixrate)
      case ('broyden')
        call buffer_broyden_ns(iter)
    end select
    call timer_end(LOG_CALC_RHO)
    
  
    call timer_begin(LOG_CALC_HARTREE)
    if(imesh_s_all==1.or.(imesh_s_all==0.and.nproc_id_global<nproc_Mxin_mul*nproc_Mxin_mul_s_dm))then
      call Hartree_ns(lg,mg,ng,system%Brl,srg_ng,stencil)
    end if
    call timer_end(LOG_CALC_HARTREE)
  
  
    call timer_begin(LOG_CALC_EXC_COR)
    if(imesh_s_all==1.or.(imesh_s_all==0.and.nproc_id_global<nproc_Mxin_mul*nproc_Mxin_mul_s_dm))then
      call exc_cor_ns(ppn)
    end if
    call timer_end(LOG_CALC_EXC_COR)
   

    call timer_begin(LOG_ALLGATHERV_TOTAL)
    call allgatherv_vlocal
    call timer_end(LOG_ALLGATHERV_TOTAL)
    
  
    call timer_begin(LOG_CALC_TOTAL_ENERGY)
    if(iperiodic==3) then
      allocate(stencil%kAc(k_sta:k_end,3))
      do jj=1,3
        stencil%kAc(k_sta:k_end,jj) = k_rd(jj,k_sta:k_end)
      end do
      call update_kvector_nonlocalpt(ppg,stencil%kAc,k_sta,k_end)
    end if
    do jspin=1,system%nspin
      V_local(jspin)%f = Vlocal(:,:,:,jspin)
    end do
    if(ilsda == 1) then
      do jspin=1,system%nspin
        srho(jspin)%f = rho_s(:,:,:,jspin)
        sVxc(jspin)%f = Vxc_s(:,:,:,jspin)
      end do
    else
      srho(1)%f = rho
      sVxc(1)%f = Vxc
    end if
    sVh%f = Vh
    energy%E_xc = Exc
    call calc_eigen_energy(energy,spsi,shpsi,sttpsi,system,info,mg,V_local,stencil,srg,ppg)
    select case(iperiodic)
    case(0)
      call calc_Total_Energy_isolated(energy,system,info,ng,pp,srho,sVh,sVxc)
    case(3)
      if(iflag_hartree==2)then
        fg%rhoG_elec = rhoe_G
      else if(iflag_hartree==4)then
        do iz=1,lg_num(3)/NPUZ
        do iy=1,lg_num(2)/NPUY
        do ix=ng%is(1)-lg%is(1)+1,ng%ie(1)-lg%is(1)+1
          n=(iz-1)*lg_num(2)/NPUY*lg_num(1)+(iy-1)*lg_num(1)+ix
          nn=ix-(ng%is(1)-lg%is(1)+1)+1+(iy-1)*ng%num(1)+(iz-1)*lg%num(2)/NPUY*ng%num(1)+fg%ig_s-1
          fg%rhoG_elec(nn) = rhoe_G(n)
        enddo
        enddo
        enddo
      end if
      call calc_Total_Energy_periodic(energy,system,pp,fg)
    end select
    esp = energy%esp(:,:,1) !++++++++
    if(iperiodic==3) deallocate(stencil%kAc,ppg%ekr_uV)
    call timer_end(LOG_CALC_TOTAL_ENERGY)


    call timer_begin(LOG_CALC_CHANGE_ORDER)
    if(iperiodic==0)then  
      call change_order(psi)
    end if
    call timer_end(LOG_CALC_CHANGE_ORDER)
  
  else if(iscf_order==2)then

    call timer_begin(LOG_CALC_GRAM_SCHMIDT)
    call gram_schmidt(system, mg, info, spsi)
    call timer_end(LOG_CALC_GRAM_SCHMIDT)


    call timer_begin(LOG_CALC_SUBSPACE_DIAG)
    if(Miter>iDiter_nosubspace_diag)then
      select case(iperiodic)
      case(0)

        call subspace_diag(mg,info,stencil,srg_ob_1,spsi,ilsda,nproc_ob,iparaway_ob,iobnum,itotmst,k_sta,k_end,mst,ifmst,hvol,  &
                info_ob,bnmat,cnmat,hgs,ppg,vlocal)

      case(3)

        call subspace_diag_periodic(mg,info,stencil,srg_ob_1,spsi,ilsda,nproc_ob,iparaway_ob,  &
                                    iobnum,itotmst,k_sta,k_end,mst,ifmst,hvol,   &
                                    info_ob,bnmat,cnmat,hgs,ppg,vlocal,num_kpoints_rd,k_rd)

      end select
    end if
    call timer_end(LOG_CALC_SUBSPACE_DIAG)

    call timer_begin(LOG_CALC_GRAM_SCHMIDT)
    call gram_schmidt(system, mg, info, spsi)    
    call timer_end(LOG_CALC_GRAM_SCHMIDT)

    call timer_begin(LOG_CALC_MINIMIZATION)
    if( amin_routine == 'cg' .or. (amin_routine == 'cg-diis' .and. Miter <= iDiterYBCG) ) then
      select case(iperiodic)
      case(0)
        select case(gscg)
        case('y')
          call sgscg(mg,nspin,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,hvol,ilsda,nproc_ob,iparaway_ob, &
                     rxk_ob,rhxk_ob,rgk_ob,rpk_ob,   &
                     info_ob,bnmat,cnmat,hgs,ppg,vlocal)
        case('n')
          call dtcg(mg,nspin,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,hvol,ilsda,nproc_ob,iparaway_ob,  &
                    info_ob,bnmat,cnmat,hgs,ppg,vlocal)
        end select
      case(3)
        select case(gscg)
        case('y')
          call gscg_periodic(mg,nspin,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,hvol,ilsda,nproc_ob,iparaway_ob,   &
                             zxk_ob,zhxk_ob,zgk_ob,zpk_ob,zpko_ob,zhtpsi_ob,   &
                             info_ob,bnmat,cnmat,hgs,ppg,vlocal,num_kpoints_rd,k_rd)
        case('n')
          call dtcg_periodic(mg,nspin,info,stencil,srg_ob_1,spsi,iflag,itotmst,mst,hvol,ilsda,nproc_ob,iparaway_ob,   &
                             info_ob,bnmat,cnmat,hgs,ppg,vlocal,num_kpoints_rd,k_rd)
        end select
      end select
    else if( amin_routine == 'diis' .or. amin_routine == 'cg-diis' ) then
      select case(iperiodic)
      case(0)
        call rmmdiis(mg,nspin,info,stencil,srg_ob_1,spsi,itotmst,mst,num_kpoints_rd,hvol,iflag_diisjump,energy%esp, &
                     norm_diff_psi_stock,info_ob,bnmat,cnmat,hgs,ppg,vlocal,iparaway_ob)
      case(3)
        stop "rmmdiis method is not implemented for periodic systems."
      end select
    end if
    call timer_end(LOG_CALC_MINIMIZATION)


    call timer_begin(LOG_CALC_GRAM_SCHMIDT)
    call gram_schmidt(system, mg, info, spsi)    
    call timer_end(LOG_CALC_GRAM_SCHMIDT)





    ! Store to psi/zpsi for calc_density
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

    call timer_begin(LOG_CALC_RHO)
    select case(iperiodic)
    case(0)
      call calc_density(psi)
    case(3)
      call calc_density(zpsi)
    end select
    call timer_end(LOG_CALC_RHO)


    select case(amixing)
      case ('simple')
        call simple_mixing(1.d0-rmixrate,rmixrate)
      case ('broyden')
        call buffer_broyden_ns(iter)
    end select
 

    call timer_begin(LOG_CALC_HARTREE)
    if(imesh_s_all==1.or.(imesh_s_all==0.and.nproc_id_global<nproc_Mxin_mul*nproc_Mxin_mul_s_dm))then
      call Hartree_ns(lg,mg,ng,system%brl,srg_ng,stencil)
    end if
    call timer_end(LOG_CALC_HARTREE)

  
    call timer_begin(LOG_CALC_EXC_COR)
    if(imesh_s_all==1.or.(imesh_s_all==0.and.nproc_id_global<nproc_Mxin_mul*nproc_Mxin_mul_s_dm))then
      call exc_cor_ns(ppn)
    end if
    call timer_end(LOG_CALC_EXC_COR)
   

    call timer_begin(LOG_ALLGATHERV_TOTAL)
    call allgatherv_vlocal
    call timer_end(LOG_ALLGATHERV_TOTAL)

    
    call timer_begin(LOG_CALC_RHO)
    select case(iperiodic)
    case(0)
      call calc_density(psi)
    case(3)
      call calc_density(zpsi)
    end select
    call timer_end(LOG_CALC_RHO)


    call timer_begin(LOG_CALC_TOTAL_ENERGY)
    if(iperiodic==3) then
      allocate(stencil%kAc(k_sta:k_end,3))
      do jj=1,3
        stencil%kAc(k_sta:k_end,jj) = k_rd(jj,k_sta:k_end)
      end do
      call update_kvector_nonlocalpt(ppg,stencil%kAc,k_sta,k_end)
    end if
    do jspin=1,system%nspin
      V_local(jspin)%f = Vlocal(:,:,:,jspin)
    end do
    if(ilsda == 1) then
      do jspin=1,system%nspin
        srho(jspin)%f = rho_s(:,:,:,jspin)
        sVxc(jspin)%f = Vxc_s(:,:,:,jspin)
      end do
    else
      srho(1)%f = rho
      sVxc(1)%f = Vxc
    end if
    sVh%f = Vh
    energy%E_xc = Exc
    call calc_eigen_energy(energy,spsi,shpsi,sttpsi,system,info,mg,V_local,stencil,srg,ppg)
    select case(iperiodic)
    case(0)
      call calc_Total_Energy_isolated(energy,system,info,ng,pp,srho,sVh,sVxc)
    case(3)
      if(iflag_hartree==2)then
        fg%rhoG_elec = rhoe_G
      else if(iflag_hartree==4)then
        do iz=1,lg_num(3)/NPUZ
        do iy=1,lg_num(2)/NPUY
        do ix=ng%is(1)-lg%is(1)+1,ng%ie(1)-lg%is(1)+1
          n=(iz-1)*lg_num(2)/NPUY*lg_num(1)+(iy-1)*lg_num(1)+ix
          nn=ix-(ng%is(1)-lg%is(1)+1)+1+(iy-1)*ng%num(1)+(iz-1)*lg%num(2)/NPUY*ng%num(1)+fg%ig_s-1
          fg%rhoG_elec(nn) = rhoe_G(n)
        enddo
        enddo
        enddo
      end if
      call calc_Total_Energy_periodic(energy,system,pp,fg)
    end select
    esp = energy%esp(:,:,1) !++++++++
    if(iperiodic==3) deallocate(stencil%kAc,ppg%ekr_uV)
    call timer_end(LOG_CALC_TOTAL_ENERGY)
  end if

  call timer_begin(LOG_WRITE_RESULTS)

  select case(convergence)
    case('rho_dne')
      sum0=0.d0
!$OMP parallel do reduction(+:sum0) private(iz,iy,ix)
      do iz=ng_sta(3),ng_end(3) 
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        sum0=sum0+abs(rho(ix,iy,iz)-rho_stock(ix,iy,iz,1))
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
      do iz=ng_sta(3),ng_end(3) 
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        sum0=sum0+(rho(ix,iy,iz)-rho_stock(ix,iy,iz,1))**2
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
      do iz=ng_sta(3),ng_end(3) 
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        sum0=sum0+(Vlocal(ix,iy,iz,1)-Vlocal_stock(ix,iy,iz,1))**2
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
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    rNebox1=rNebox1+rho(ix,iy,iz)
  end do
  end do
  end do
  call comm_summation(rNebox1,rNebox2,nproc_group_global)
  if(comm_is_root(nproc_id_global))then
    write(*,*) "Ne=",rNebox2*Hvol
  end if
  call timer_end(LOG_WRITE_RESULTS)


if(ilsda==0)then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    rho_stock(ix,iy,iz,1)=rho(ix,iy,iz)
    Vlocal_stock(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)
  end do
  end do
  end do
else if(ilsda==1)then
!$OMP parallel do private(iz,iy,ix)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    rho_stock(ix,iy,iz,1)=rho(ix,iy,iz)
    Vlocal_stock(ix,iy,iz,1:2)=Vlocal(ix,iy,iz,1:2)
  end do
  end do
  end do
end if

end do DFT_Iteration


deallocate(idiis_sd)
call timer_end(LOG_GS_ITERATION)


call timer_begin(LOG_DEINIT_GS_ITERATION)
if(icalcforce==1) then
  call calc_force_salmon(force,system,pp,fg,info,mg,stencil,srg,ppg,spsi)
  if(comm_is_root(nproc_id_global))then
    write(*,*) "===== force ====="
    do iatom=1,MI
      select case(unit_system)
      case('au','a.u.')
        write(*,'(i6,3e16.8)') iatom,(force%f(ix,iatom),ix=1,3)
      case('A_eV_fs')
        write(*,'(i6,3e16.8)') iatom,(force%f(ix,iatom)*2.d0*Ry/a_B,ix=1,3)
      end select
    end do
  end if
end if

if(iflag_opt==1) then
  call structure_opt_check(MI,iopt,iopt_tranc,force)
  if(iopt_tranc/=1) call structure_opt(MI,iopt,force,system%Rion)
  !! Rion is old variables to be removed 
  !! but currently it is used in many subroutines.
  Rion(:,:) = system%Rion(:,:) 
  if(comm_is_root(nproc_id_global))then
    write(*,*) "atomic coordinate"
    do iatom=1,MI
       write(*,10) "'"//trim(AtomName(Kion(iatom)))//"'",  &
                   (system%Rion(jj,iatom)*ulength_from_au,jj=1,3), &
                   Kion(iatom), flag_geo_opt_atom(iatom)
    end do
10  format(a5,3f16.8,i3,a3)
  end if
  if(iopt_tranc==1) then
    call structure_opt_fin
    exit Multigrid_Iteration
  end if

else
   select case(iperiodic)
   case(0)
      deallocate(spsi%rwf)
   case(3)
      deallocate(spsi%zwf)
   end select
end if
call timer_end(LOG_DEINIT_GS_ITERATION)


end do Multigrid_Iteration
if(iopt_tranc==1)then
  exit Structure_Optimization_Iteration
end if
end do Structure_Optimization_Iteration


!---------------------------------------- Output
call timer_begin(LOG_WRITE_RESULTS)

call band_information

call write_eigen

if(out_psi=='y') then
  call writepsi(lg)
end if

if(out_dns=='y') then
  call writedns(lg,mg,ng,rho,matbox_m,matbox_m2,icoo1d,hgs,igc_is,igc_ie,gridcoo,iscfrt)
end if

if(out_dos=='y') then
  call calc_dos
end if

if(out_pdos=='y') then
  call calc_pdos
end if

if(OC==2)then
  call prep_ini
end if

if(out_elf=='y')then
  allocate(elf(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),      &
               lg_sta(3):lg_end(3)))
  call calcELF
  call writeelf(lg,elf,icoo1d,hgs,igc_is,igc_ie,gridcoo,iscfrt)
  deallocate(elf)
end if
call timer_end(LOG_WRITE_RESULTS)


call timer_begin(LOG_WRITE_LDA_DATA)
! LDA data
! subroutines in scf_data.f90
if ( OC==1.or.OC==2.or.OC==3 ) then
  call OUT_data
end if
call timer_end(LOG_WRITE_LDA_DATA)


! LDA information
call timer_begin(LOG_WRITE_INFOS)
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
call timer_end(LOG_WRITE_INFOS)

deallocate(Vlocal)
call finalize_xc(xc_func)

call timer_end(LOG_TOTAL)


if(comm_is_root(nproc_id_global))then
  call write_gs_performance(output_unit)
end if

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

END subroutine Real_Space_DFT

!=======================================================================
!========================================= Grid generation and labelling

SUBROUTINE init_mesh(lg,mg,itmg)
use structures, only: s_rgrid
use salmon_parallel, only: nproc_id_global, nproc_size_global
use salmon_communication, only: comm_is_root
use inputoutput, only: iperiodic
use global_variables_scf
implicit none
type(s_rgrid),intent(out) :: lg
type(s_rgrid),intent(out) :: mg
integer,intent(in) :: itmg
integer :: j
real(8) :: rLsize1(3)

if(comm_is_root(nproc_id_global))      &
    print *,"----------------------------------- init_mesh"

rLsize1(:)=rLsize(:,itmg)
call setlg(lg,lg_sta,lg_end,lg_num,ista_Mx_ori,iend_Mx_ori,inum_Mx_ori,    &
           Hgs,Nd,rLsize1,imesh_oddeven,iperiodic,iscfrt)
call check_fourier

allocate(ista_Mxin(3,0:nproc_size_global-1),iend_Mxin(3,0:nproc_size_global-1))
allocate(inum_Mxin(3,0:nproc_size_global-1))

call setmg(mg,mg_sta,mg_end,mg_num,ista_Mxin,iend_Mxin,inum_Mxin,  &
           lg_sta,lg_num,nproc_size_global,nproc_id_global,nproc_Mxin,nproc_k,nproc_ob,isequential,iscfrt)

if(comm_is_root(nproc_id_global)) write(*,*) "Mx     =", iend_Mx_ori

if(iperiodic==3 .and. nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin(3)==1) then
  if(comm_is_root(nproc_id_global)) write(*,*) "r-space parallelization: off"
  lg%is(1:3)=lg_sta(1:3)
  lg%ie(1:3)=lg_end(1:3)
  lg%num(1:3)=lg_num(1:3)
  lg%is_overlap(1:3)=lg_sta(1:3)-nd
  lg%ie_overlap(1:3)=lg_end(1:3)+nd
  lg%is_array(1:3)=lg_sta(1:3)
  lg%ie_array(1:3)=lg_end(1:3)

  if(allocated(lg%idx)) deallocate(lg%idx)
  if(allocated(lg%idy)) deallocate(lg%idy)
  if(allocated(lg%idz)) deallocate(lg%idz)
  allocate(lg%idx(lg%is_overlap(1):lg%ie_overlap(1)) &
          ,lg%idy(lg%is_overlap(2):lg%ie_overlap(2)) &
          ,lg%idz(lg%is_overlap(3):lg%ie_overlap(3)))
  do j=lg%is_overlap(1),lg%ie_overlap(1)
    lg%idx(j) = mod(j+lg%num(1)-1,lg%num(1))+1
  end do
  do j=lg%is_overlap(2),lg%ie_overlap(2)
    lg%idy(j) = mod(j+lg%num(2)-1,lg%num(2))+1
  end do
  do j=lg%is_overlap(3),lg%ie_overlap(3)
    lg%idz(j) = mod(j+lg%num(3)-1,lg%num(3))+1
  end do

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

