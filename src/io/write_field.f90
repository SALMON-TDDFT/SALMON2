!
!  Copyright 2019-2020 SALMON developers
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
module writefield
  implicit none

contains

!===================================================================================================================================

subroutine write_dns(lg,mg,system,info,rho_s,rho0_s,itt)
  use inputoutput, only: au_length_aa
  use salmon_global, only: format_voxel_data,theory
  use structures
  implicit none
  type(s_rgrid)     ,intent(in)          :: lg
  type(s_rgrid)     ,intent(in)          :: mg
  type(s_dft_system),intent(in)          :: system
  type(s_parallel_info),intent(in)       :: info
  type(s_scalar)    ,intent(in)          :: rho_s (system%nspin)
  type(s_scalar)    ,intent(in),optional :: rho0_s(system%nspin)
  integer           ,intent(in),optional :: itt
  !
  integer :: ispin,ix,iy,iz
  character(60) :: suffix,suffix_u,suffix_d
  character(30) :: phys_quantity
  character(10) :: filenum
  real(8),dimension(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:system%nspin) :: wrk1,wrk2,dns
  
  call calc_wrk2(rho_s)
  !$omp workshare
  dns = wrk2
  !$omp end workshare

  select case(theory)
  case('dft','dft_band','dft_md') 
    suffix = "dns"
    suffix_u = "dns_u"
    suffix_d = "dns_d"
  case('tddft_response','tddft_pulse','single_scale_maxwell_tddft','multi_scale_maxwell_tddft')
    write(filenum, '(i6.6)') itt
    suffix = "dns_"//adjustl(filenum)
    suffix_u = "dns_u_"//adjustl(filenum)
    suffix_d = "dns_d_"//adjustl(filenum)
  case default
    stop 'invalid theory @ write_dns'
  end select
  phys_quantity = "dns"
  
  if(system%nspin==1) then
    call write_dns_core(suffix,dns(:,:,:,1))
  else if(system%nspin==2) then
    call write_dns_core(suffix_u,dns(:,:,:,1))
    call write_dns_core(suffix_d,dns(:,:,:,2))
    !$omp workshare
    wrk1(:,:,:,1) = dns(:,:,:,1) + dns(:,:,:,2)
    !$omp end workshare
    call write_dns_core(suffix,wrk1(:,:,:,1))
  else
    stop 'invalid nspin'
  end if
  
  if(.not.present(rho0_s)) then
    return
  end if

  select case(theory)
  case('tddft_response','tddft_pulse','single_scale_maxwell_tddft','multi_scale_maxwell_tddft')

    call calc_wrk2(rho0_s)
    !$omp workshare
    dns = dns - wrk2
    !$omp end workshare

    write(filenum, '(i6.6)') itt
    suffix = "dnsdiff_"//adjustl(filenum)
    suffix_u = "dnsdiff_u_"//adjustl(filenum)
    suffix_d = "dnsdiff_d_"//adjustl(filenum)
    phys_quantity = "dnsdiff"
    
    if(system%nspin==1) then
      call write_dns_core(suffix,dns(:,:,:,1))
    else if(system%nspin==2) then
      call write_dns_core(suffix_u,dns(:,:,:,1))
      call write_dns_core(suffix_d,dns(:,:,:,2))
      !$omp workshare
      wrk1(:,:,:,1) = dns(:,:,:,1) + dns(:,:,:,2)
      !$omp end workshare
      call write_dns_core(suffix,wrk1(:,:,:,1))
    else
      stop 'invalid nspin'
    end if
    
  case default
    stop 'invalid theory @ write_dns'
  end select
  
contains

  subroutine calc_wrk2(trho)
    use communication, only: comm_summation
    implicit none
    type(s_scalar),intent(in) :: trho(system%nspin)
    !$omp workshare
    wrk1 = 0.d0
    !$omp end workshare
    do ispin=1,system%nspin
      !$OMP parallel do collapse(2) private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        wrk1(ix,iy,iz,ispin) = trho(ispin)%f(ix,iy,iz)
      end do
      end do
      end do
    end do
    if(format_voxel_data=='avs')then
      !$omp workshare
      wrk1 = wrk1 /(au_length_aa**3)
      !$omp end workshare
    end if
    call comm_summation(wrk1,wrk2,lg%num(1)*lg%num(2)*lg%num(3)*system%nspin,info%icomm_r)
  end subroutine calc_wrk2
  
  subroutine write_dns_core(tsuffix,tdns)
    use write_file3d, only: write_avs, write_cube, write_vtk
    implicit none
    real(8)      ,intent(in) :: tdns(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
    character(60),intent(in) :: tsuffix
    character(20)            :: header_unit
    if(format_voxel_data=='avs')then
      header_unit='A**(-3)'
      call write_avs(lg,103,tsuffix,header_unit,tdns)
    else if(format_voxel_data=='cube')then
      call write_cube(lg,103,tsuffix,phys_quantity,tdns,system)
    else if(format_voxel_data=='vtk')then
      call write_vtk(lg,103,tsuffix,tdns,system%hgs)
    end if
  end subroutine write_dns_core

end subroutine write_dns

!===================================================================================================================================
!subroutine write_dns_ac_je(info,mg,system,rho,j_e,itt,action)
subroutine write_dns_ac_je(info,mg,system,rho,fw,itt,action)
  use structures, only: s_dft_system,s_parallel_info,s_rgrid,s_scalar,allocate_scalar,deallocate_scalar,s_singlescale  !,s_vector
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root
  use salmon_global, only: kion,izatom
  use filesystem, only: create_directory
  implicit none
  type(s_parallel_info),intent(in) :: info
  type(s_dft_system),intent(in) :: system
  type(s_rgrid),intent(in)  :: mg
 !type(s_vector),intent(in) :: j_e
  type(s_singlescale) :: fw
  real(8),intent(in) :: rho(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  integer,intent(in) :: itt
  character(3),intent(in) :: action
  integer :: fp
  real(8) :: dummy(1:3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  character(10) :: filenum1, filenum2
  character(256):: wdir1, wdir2, ofile

  wdir1 = "output_ion/"
  wdir2 = "output_dns_ac_je"
  if(comm_is_root(info%id_ko)) then
     write(filenum2,'(i6.6)') info%id_r
     wdir2 = "output_dns_ac_je_idr"//trim(filenum2)//"/"
  endif


  if(action=='new') then

     if(comm_is_root(nproc_id_global)) then
        call create_directory(wdir1)
     endif

     if(comm_is_root(info%id_ko)) then
        call create_directory(wdir2)  
     endif

  else

     if(comm_is_root(nproc_id_global)) then

        fp = 299
        write(filenum1,'(i6.6)') itt
        filenum1=adjustl(filenum1)
        ofile =trim(wdir1)//"it"//trim(filenum1)//"_ion.bin"
        open(fp,file=trim(ofile),form='unformatted',access='stream')
        write(fp) system%nion
        write(fp) system%Rion(1:3,1:system%nion)
        write(fp) izatom(kion(1:system%nion))
        close(fp)

     endif

     if(comm_is_root(info%id_ko)) then

        fp = 300 + info%id_r
        write(filenum1,'(i6.6)') itt
        write(filenum2,'(i6.6)') info%id_r
        filenum1=adjustl(filenum1)
        filenum2=adjustl(filenum2)
       !ofile= trim(wdir2)//"it"//trim(filenum1)//"_idr"//trim(filenum2)//".bin"
        ofile= trim(wdir2)//"it"//trim(filenum1)//".bin"
        open(fp,file=trim(ofile),form='unformatted',access='stream')
        write(fp) mg%is(1), mg%ie(1),mg%is(2), mg%ie(2),mg%is(3), mg%ie(3)
        write(fp) system%hgs(1:3)
        write(fp) rho(      mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
        if(itt==0) then
        dummy(1:3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))=0d0
        write(fp) dummy(1:3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
        write(fp) dummy(1:3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
        else
       !write(fp) j_e%v(1:3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
        write(fp) fw%curr(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:3)
        write(fp) system%Ac_micro%v( &
                        1:3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
        endif
        close(fp)
        
     endif

  endif  !action

end subroutine write_dns_ac_je

!===================================================================================================================================

subroutine write_elf(itt,lg,mg,system,info,stencil,rho,srg,srg_scalar,tpsi)
  use salmon_global       ,only: format_voxel_data,theory
  use structures
  use math_constants      ,only: pi
  use communication       ,only: comm_summation
  use misc_routines       ,only: get_wtime
  use sendrecv_grid       ,only: update_overlap_complex8,update_overlap_real8
  use stencil_sub         ,only: calc_gradient_field
  use write_file3d
  implicit none
  integer                 ,intent(in) :: itt
  type(s_rgrid)           ,intent(in) :: lg,mg
  type(s_dft_system)      ,intent(in) :: system
  type(s_parallel_info)   ,intent(in) :: info
  type(s_stencil)         ,intent(in) :: stencil
  type(s_scalar)          ,intent(in) :: rho
  type(s_sendrecv_grid)               :: srg,srg_scalar
  type(s_orbital)                     :: tpsi
  !
  real(8):: elf(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  character(60) :: suffix
  character(30) :: phys_quantity
  character(10) :: filenum
  character(20) :: header_unit
  integer :: nspin,no,nk,ik_s,ik_e,io_s,io_e,is(3),ie(3)
  integer :: io,ik,ispin,ix,iy,iz
  integer,parameter :: Nd=4
  !
  real(8) :: elftau(mg%is(1):mg%ie(1),   &
                    mg%is(2):mg%ie(2),   &
                    mg%is(3):mg%ie(3))
  real(8) :: mrelftau(mg%is(1):mg%ie(1),   &
                      mg%is(2):mg%ie(2),   &
                      mg%is(3):mg%ie(3))
  real(8) :: curden(mg%is(1):mg%ie(1),   &
                    mg%is(2):mg%ie(2),   &
                    mg%is(3):mg%ie(3))
  real(8) :: mrcurden(mg%is(1):mg%ie(1),   &
                      mg%is(2):mg%ie(2),   &
                      mg%is(3):mg%ie(3))
  complex(8) :: gradzpsi(3,mg%is_array(1):mg%ie_array(1),   &
                       mg%is_array(2):mg%ie_array(2),   &
                       mg%is_array(3):mg%ie_array(3))
  real(8) :: gradrho(3,mg%is(1):mg%ie(1),   &
                       mg%is(2):mg%ie(2),   &
                       mg%is(3):mg%ie(3))
  real(8) :: gradrho2(mg%is(1):mg%ie(1),   &
                      mg%is(2):mg%ie(2),   &
                      mg%is(3):mg%ie(3))
  real(8) :: elfc(mg%is(1):mg%ie(1),   &
                  mg%is(2):mg%ie(2),   &
                  mg%is(3):mg%ie(3))
  real(8) :: elfcuni(mg%is(1):mg%ie(1),   &
                     mg%is(2):mg%ie(2),   &
                     mg%is(3):mg%ie(3))
  real(8) :: rho_half(mg%is(1):mg%ie(1),   &
                      mg%is(2):mg%ie(2),   &
                      mg%is(3):mg%ie(3))
  real(8) :: box(mg%is_array(1):mg%ie_array(1), &
                 mg%is_array(2):mg%ie_array(2), &
                 mg%is_array(3):mg%ie_array(3))
  real(8) :: matbox_l(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  complex(8) :: ztmp(mg%is_array(1):mg%ie_array(1), &
                     mg%is_array(2):mg%ie_array(2), &
                     mg%is_array(3):mg%ie_array(3))
  
  if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ calc_elf"

  nspin = system%nspin
  no   = system%no
  nk   = system%nk
  is   = mg%is
  ie   = mg%ie
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  
  rho_half = 0.d0
  ztmp     = 0.d0
  gradzpsi = 0.d0
  mrelftau = 0.d0
  mrcurden = 0.d0
  elftau   = 0.d0
  curden   = 0.d0
  box      = 0.d0
  gradrho  = 0.d0
  gradrho2 = 0.d0
  elfc     = 0.d0
  elfcuni  = 0.d0
  matbox_l = 0.d0
  elf      = 0.d0
  
!$OMP parallel do private(iz,iy,ix)
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
  do ix=is(1),ie(1)
    rho_half(ix,iy,iz)=rho%f(ix,iy,iz)/2.d0
  end do
  end do
  end do
  
  if(info%if_divide_rspace) then
    if(allocated(tpsi%rwf)) then
      call update_overlap_real8(srg, mg, tpsi%rwf)
    else
      call update_overlap_complex8(srg, mg, tpsi%zwf)
    end if
  end if
  
  do ik=ik_s,ik_e
  do io=io_s,io_e
  do ispin=1,nspin
    
    if(allocated(tpsi%rwf)) then
      ztmp = dcmplx(tpsi%rwf(:,:,:,ispin,io,ik,1))
    else
      ztmp = tpsi%zwf(:,:,:,ispin,io,ik,1)
    end if
    
    call calc_gradient_psi(ztmp,gradzpsi,mg%is_array,mg%ie_array,is,ie, &
                           mg%idx,mg%idy,mg%idz,stencil%coef_nab,system%rmatrix_B)
    
!$OMP parallel do private(iz,iy,ix)
    do iz=is(3),ie(3)
    do iy=is(2),ie(2)
    do ix=is(1),ie(1)
      mrelftau(ix,iy,iz)=mrelftau(ix,iy,iz)            &
                         +abs(gradzpsi(1,ix,iy,iz))**2 &
                         +abs(gradzpsi(2,ix,iy,iz))**2 &
                         +abs(gradzpsi(3,ix,iy,iz))**2
      mrcurden(ix,iy,iz)=mrcurden(ix,iy,iz)      &
                        +( abs(conjg(ztmp(ix,iy,iz))*gradzpsi(1,ix,iy,iz)      &
                              -ztmp(ix,iy,iz)*conjg(gradzpsi(1,ix,iy,iz)))**2  &
                          +abs(conjg(ztmp(ix,iy,iz))*gradzpsi(2,ix,iy,iz)      &
                              -ztmp(ix,iy,iz)*conjg(gradzpsi(2,ix,iy,iz)))**2  &
                          +abs(conjg(ztmp(ix,iy,iz))*gradzpsi(3,ix,iy,iz)      &
                              -ztmp(ix,iy,iz)*conjg(gradzpsi(3,ix,iy,iz)))**2 )/2.d0
    end do
    end do
    end do
    
  end do
  end do
  end do
  
  call comm_summation(mrelftau,elftau,mg%num(1)*mg%num(2)*mg%num(3),info%icomm_ko)
  call comm_summation(mrcurden,curden,mg%num(1)*mg%num(2)*mg%num(3),info%icomm_ko)
  
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    box(ix,iy,iz) = rho_half(ix,iy,iz)
  end do
  end do
  end do
  
  if(info%if_divide_rspace) call update_overlap_real8(srg_scalar, mg, box)
  call calc_gradient_field(mg,stencil%coef_nab,system%rmatrix_B,box,gradrho)
  
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    gradrho2(ix,iy,iz)=gradrho(1,ix,iy,iz)**2 &
                      +gradrho(2,ix,iy,iz)**2 &
                      +gradrho(3,ix,iy,iz)**2
    elfc(ix,iy,iz)=elftau(ix,iy,iz)-gradrho2(ix,iy,iz)/rho_half(ix,iy,iz)/4.d0  &
                                   -curden(ix,iy,iz)/rho_half(ix,iy,iz)
  end do
  end do
  end do
  
  ! matbox_l stores ELF
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    elfcuni(ix,iy,iz)=3.d0/5.d0*(6.d0*Pi**2)**(2.d0/3.d0)      &
                     *rho_half(ix,iy,iz)**(5.d0/3.d0)
    matbox_l(ix,iy,iz)=1.d0/(1.d0+elfc(ix,iy,iz)**2/elfcuni(ix,iy,iz)**2)
  end do
  end do
  end do
  
  call comm_summation(matbox_l,elf,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_r)
  
  select case(theory)
  case('dft','dft_band','dft_md') 
    suffix = "elf"
  case('tddft_response','tddft_pulse','single_scale_maxwell_tddft','multi_scale_maxwell_tddft')
    write(filenum, '(i6.6)') itt
    suffix = "elf_"//adjustl(filenum)
  case default
    stop 'invalid theory'
  end select
  
  phys_quantity = "elf"
  if(format_voxel_data=='avs')then
    header_unit = "none"
    call write_avs(lg,103,suffix,header_unit,elf)
  else if(format_voxel_data=='cube')then
    call write_cube(lg,103,suffix,phys_quantity,elf,system)
  else if(format_voxel_data=='vtk')then
    call write_vtk(lg,103,suffix,elf,system%hgs)
  end if
  
end subroutine write_elf

!===================================================================================================================================

subroutine write_estatic(lg,mg,system,stencil,info,Vh,srg_scalar,itt)
  use salmon_global, only: format_voxel_data
  use structures
  use sendrecv_grid, only: update_overlap_real8
  use stencil_sub, only: calc_gradient_field
  use communication, only: comm_summation
  use write_file3d
  implicit none
  type(s_rgrid)   ,  intent(in) :: lg,mg
  type(s_dft_system),intent(in) :: system
  type(s_stencil) ,  intent(in) :: stencil
  type(s_parallel_info),intent(in) :: info
  type(s_scalar)  ,intent(in) :: Vh
  type(s_sendrecv_grid)       :: srg_scalar
  integer,intent(in),optional :: itt
  !
  integer :: ix,iy,iz,jj
  integer,parameter :: Nd=4
  character(60) :: suffix
  character(30) :: phys_quantity
  character(10) :: filenum
  character(20) :: header_unit
  real(8) :: grad_Vh(3,mg%is(1):mg%ie(1),   &
                       mg%is(2):mg%ie(2),   &
                       mg%is(3):mg%ie(3))
  real(8) :: box(mg%is_array(1):mg%ie_array(1), &
  & mg%is_array(2):mg%ie_array(2), &
  & mg%is_array(3):mg%ie_array(3))
  real(8),dimension(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)) :: rmat,rmat2
  
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    box(ix,iy,iz) = Vh%f(ix,iy,iz)
  end do
  end do
  end do

  if(info%if_divide_rspace) call update_overlap_real8(srg_scalar, mg, box)
  call calc_gradient_field(mg,stencil%coef_nab,system%rmatrix_B,box,grad_Vh)

  do jj=1,3
    
    rmat = 0d0
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rmat(ix,iy,iz) = grad_Vh(jj,ix,iy,iz)
    end do
    end do
    end do
 
    call comm_summation(rmat,rmat2,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_r)
  
    write(filenum, '(i6.6)') itt
    if(jj==1)then
      suffix = "Exsta_"//adjustl(filenum)
      phys_quantity = "exsta"
    else if(jj==2)then
      suffix = "Eysta_"//adjustl(filenum)
      phys_quantity = "eysta"
    else if(jj==3)then
      suffix = "Ezsta_"//adjustl(filenum)
      phys_quantity = "ezsta"
    end if

    if(format_voxel_data=='avs')then
      header_unit = "V/A"
      call write_avs(lg,103,suffix,header_unit,rmat2)
    else if(format_voxel_data=='cube')then
      call write_cube(lg,103,suffix,phys_quantity,rmat2,system)
    else if(format_voxel_data=='vtk')then
      call write_vtk(lg,103,suffix,rmat2,system%hgs)
    end if

  end do  
 
end subroutine write_estatic

!===================================================================================================================================

subroutine write_psi(lg,mg,system,info,spsi)
  use inputoutput   ,only: au_length_aa
  use structures
  use salmon_global ,only: format_voxel_data
  use communication ,only: comm_summation
  use write_file3d
  implicit none
  type(s_rgrid)        ,intent(in) :: lg,mg
  type(s_dft_system)   ,intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_orbital)      ,intent(in) :: spsi
  !
  integer :: io,ik,ispin,ix,iy,iz
  complex(8),dimension(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)) :: cmatbox,cmatbox2
  character(60) :: suffix_re, suffix_im
  character(30) :: phys_quantity
  character(10) :: fileid_k, fileid_ob, fileid_spin
  character(20) :: header_unit
  
  cmatbox = 0d0
  
  do ik=1,system%nk
  do io=1,system%no
  do ispin=1,system%nspin
  
!$omp parallel do collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        cmatbox(ix,iy,iz)=0.d0
      end do
      end do
      end do
      if(info%ik_s <= ik .and. ik <= info%ik_e .and.   &
         info%io_s <= io .and. io <= info%io_e) then
        if(allocated(spsi%rwf)) then
!$omp parallel do collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            cmatbox(ix,iy,iz) = spsi%rwf(ix,iy,iz,ispin,io,ik,1)
          end do
          end do
          end do
        else
!$omp parallel do collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            cmatbox(ix,iy,iz) = spsi%zwf(ix,iy,iz,ispin,io,ik,1)
          end do
          end do
          end do
        end if
      end if
      call comm_summation(cmatbox,cmatbox2,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)
      
      phys_quantity = "psi"
      if(allocated(spsi%rwf)) then
        write(fileid_ob, '(i8)') io
        suffix_re = "psi_ob"//trim(adjustl(fileid_ob))
        if(format_voxel_data=='avs')then
          header_unit = "A**(-3/2)"
          call write_avs(lg,103,suffix_re,header_unit,dble(cmatbox2)/sqrt(au_length_aa)**3)
        else if(format_voxel_data=='cube')then
          call write_cube(lg,103,suffix_re,phys_quantity,dble(cmatbox2),system)
        end if
      else
        write(fileid_k,    '(i8)') ik
        write(fileid_ob,   '(i8)') io
        write(fileid_spin, '(i8)') ispin
        suffix_re = "psi_k"//trim(adjustl(fileid_k))   // &
                    "_ob"  //trim(adjustl(fileid_ob))  // &
                    "_spin"//trim(adjustl(fileid_spin))// &
                    "_real" 
        suffix_im = "psi_k"//trim(adjustl(fileid_k))   // &
                    "_ob"  //trim(adjustl(fileid_ob))  // &
                    "_spin"//trim(adjustl(fileid_spin))// &
                    "_imag" 
        if(format_voxel_data=='avs')then
          header_unit = "A**(-3/2)"
          call write_avs(lg,103,suffix_re,header_unit,dble(cmatbox2)/sqrt(au_length_aa)**3)
          call write_avs(lg,103,suffix_im,header_unit,aimag(cmatbox2)/sqrt(au_length_aa)**3)
        else if(format_voxel_data=='cube')then
          call write_cube(lg,103,suffix_re,phys_quantity,dble(cmatbox2),system)
          call write_cube(lg,103,suffix_im,phys_quantity,aimag(cmatbox2),system)
        end if
      endif
      
  end do
  end do
  end do

end subroutine write_psi

!===================================================================================================================================

end module writefield
