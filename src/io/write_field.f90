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
module writefield
  implicit none

contains

!===================================================================================================================================

subroutine write_dns(lg,mg,ng,rho,rmat,rmat2,hgs,iscfrt,rho0,itt)
  use inputoutput, only: format_voxel_data,au_length_aa
  use structures, only: s_rgrid
  use salmon_parallel, only: nproc_group_global
  use salmon_communication, only: comm_summation
  use write_file3d
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  real(8),intent(in) :: rho(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(out) :: rmat(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  real(8),intent(out) :: rmat2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  real(8),intent(in) :: hgs(3)
  integer,intent(in) :: iscfrt
  real(8),intent(in),optional :: rho0(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  integer,intent(in),optional :: itt
  integer :: ix,iy,iz
  character(30) :: suffix
  character(30) :: phys_quantity
  character(10) :: filenum
  character(20) :: header_unit

  !$OMP parallel do collapse(2) private(iz,iy,ix)
  do iz=lg%is(3),lg%ie(3)
  do iy=lg%is(2),lg%ie(2)
  do ix=lg%is(1),lg%ie(1)
    rmat(ix,iy,iz)=0.d0
  end do
  end do
  end do

  !$OMP parallel do collapse(2) private(iz,iy,ix)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    rmat(ix,iy,iz)=rho(ix,iy,iz)
  end do
  end do
  end do

  if(format_voxel_data=='avs')then
    !$OMP parallel do collapse(2) private(iz,iy,ix)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      rmat(ix,iy,iz)=rmat(ix,iy,iz)/(au_length_aa**3)
    end do
    end do
    end do
  end if
  
  call comm_summation(rmat,rmat2,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)

  if(iSCFRT==1)then 
    suffix = "dns"
  else if(iSCFRT==2)then
    write(filenum, '(i6.6)') itt
    suffix = "dns_"//adjustl(filenum)
  end if
  phys_quantity = "dns"
  if(format_voxel_data=='avs')then
    header_unit='A**(-3)'
    call write_avs(lg,103,suffix,header_unit,rmat2)
  else if(format_voxel_data=='cube')then
    call write_cube(lg,103,suffix,phys_quantity,rmat2,hgs)
  else if(format_voxel_data=='vtk')then
    call write_vtk(lg,103,suffix,rmat2,hgs)
  end if

  if(iscfrt==2)then
    !$OMP parallel do collapse(2) private(iz,iy,ix)
    do iz=lg%is(3),lg%ie(3)
    do iy=lg%is(2),lg%ie(2)
    do ix=lg%is(1),lg%ie(1)
      rmat(ix,iy,iz)=0.d0
    end do
    end do
    end do

    !$OMP parallel do collapse(2) private(iz,iy,ix)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      rmat(ix,iy,iz)=rho(ix,iy,iz)-rho0(ix,iy,iz)
    end do
    end do
    end do
  
    if(format_voxel_data=='avs')then
      !$OMP parallel do collapse(2) private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        rmat(ix,iy,iz)=rmat(ix,iy,iz)/(au_length_aa**3)
      end do
      end do
      end do
    end if

    call comm_summation(rmat,rmat2,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)

    write(filenum, '(i6.6)') itt
    suffix = "dnsdiff_"//adjustl(filenum)
    phys_quantity = "dnsdiff"
    if(format_voxel_data=='avs')then
      header_unit='A**(-3)'
      call write_avs(lg,103,suffix,header_unit,rmat2)
    else if(format_voxel_data=='cube')then
      call write_cube(lg,103,suffix,phys_quantity,rmat2,hgs)
    else if(format_voxel_data=='vtk')then
      call write_vtk(lg,103,suffix,rmat2,hgs)
    end if
  end if
 
end subroutine write_dns

!===================================================================================================================================

subroutine write_elf(iscfrt,itt,lg,mg,ng,system,info,stencil,srho,srg,srg_ng,tpsi)
  use salmon_global, only: format_voxel_data
  use structures
  use math_constants, only: pi
  use salmon_communication, only: comm_summation
  use misc_routines, only: get_wtime
  use sendrecv_grid, only: update_overlap_complex8,update_overlap_real8
  use stencil_sub, only: calc_gradient_psi,calc_gradient_field
  use write_file3d
  implicit none
  integer                 ,intent(in) :: iscfrt
  integer                 ,intent(in) :: itt
  type(s_rgrid)           ,intent(in) :: lg,mg,ng
  type(s_dft_system)      ,intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_stencil)         ,intent(in) :: stencil
  type(s_scalar)          ,intent(in) :: srho
  type(s_sendrecv_grid)               :: srg,srg_ng
  type(s_orbital)                     :: tpsi
  !
  real(8):: elf(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  character(30) :: suffix
  character(30) :: phys_quantity
  character(10) :: filenum
  character(20) :: header_unit
  integer :: nspin,no,nk,ik_s,ik_e,io_s,io_e,is(3),ie(3)
  integer :: io,ik,ispin,ix,iy,iz
  integer,parameter :: Nd=4
  complex(8) :: ztmp
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
  real(8) :: gradrho(3,ng%is(1):ng%ie(1),   &
                       ng%is(2):ng%ie(2),   &
                       ng%is(3):ng%ie(3))
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
  real(8) :: box(ng%is(1)-Nd:ng%ie(1)+Nd,   &
                 ng%is(2)-Nd:ng%ie(2)+Nd,   &
                 ng%is(3)-Nd:ng%ie(3)+Nd)
  real(8) :: matbox_l(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))

  if(info%im_s/=1 .or. info%im_e/=1) stop "error: im/=1 @ calc_elf"

  nspin = system%nspin
  no = system%no
  nk = system%nk
  is = mg%is
  ie = mg%ie
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e

  !$OMP parallel do private(iz,iy,ix)
  do iz=is(3),ie(3)
  do iy=is(2),ie(2)
  do ix=is(1),ie(1)
    rho_half(ix,iy,iz)=srho%f(ix,iy,iz)/2.d0
  end do
  end do
  end do
  mrelftau=0.d0
  mrcurden=0.d0

  if(info%if_divide_rspace) then
     call update_overlap_complex8(srg, mg, tpsi%zwf)
  end if

  ik=1    ! --> do loop (future work)
  ispin=1 ! --> do loop (future work)

    do io=1,no

        call calc_gradient_psi(tpsi%zwf(:,:,:,ispin,io,ik,1),gradzpsi,mg%is_array,mg%ie_array,is,ie, &
        & mg%idx,mg%idy,mg%idz,stencil%coef_nab,system%rmatrix_B)

  !$OMP parallel do private(iz,iy,ix)
        do iz=is(3),ie(3)
        do iy=is(2),ie(2)
        do ix=is(1),ie(1)
          ztmp = tpsi%zwf(ix,iy,iz,ispin,io,ik,1)
    
          mrelftau(ix,iy,iz)=mrelftau(ix,iy,iz)+abs(gradzpsi(1,ix,iy,iz))**2      &
                             +abs(gradzpsi(2,ix,iy,iz))**2      &
                             +abs(gradzpsi(3,ix,iy,iz))**2
    
          mrcurden(ix,iy,iz)=mrcurden(ix,iy,iz)      &
               +( abs(conjg(ztmp)*gradzpsi(1,ix,iy,iz)      &
                    -ztmp*conjg(gradzpsi(1,ix,iy,iz)))**2      &
                 +abs(conjg(ztmp)*gradzpsi(2,ix,iy,iz)      &
                    -ztmp*conjg(gradzpsi(2,ix,iy,iz)))**2      &
                 +abs(conjg(ztmp)*gradzpsi(3,ix,iy,iz)      &
                    -ztmp*conjg(gradzpsi(3,ix,iy,iz)))**2 )/2.d0
    
        end do
        end do
        end do
        
    end do

    call comm_summation(mrelftau,elftau,mg%num(1)*mg%num(2)*mg%num(3),info%icomm_o)
    call comm_summation(mrcurden,curden,mg%num(1)*mg%num(2)*mg%num(3),info%icomm_o)
    
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      box(ix,iy,iz) = srho%f(ix,iy,iz)
    end do
    end do
    end do
    
    call update_overlap_real8(srg_ng, ng, box)
    call calc_gradient_field(ng%is,ng%ie,stencil%coef_nab,box,gradrho)

    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      gradrho2(ix,iy,iz)=gradrho(1,ix,iy,iz)**2      &
            +gradrho(2,ix,iy,iz)**2      &
            +gradrho(3,ix,iy,iz)**2
      elfc(ix,iy,iz)=elftau(ix,iy,iz)-gradrho2(ix,iy,iz)/rho_half(ix,iy,iz)/4.d0  &
                                     -curden(ix,iy,iz)/rho_half(ix,iy,iz)
    end do
    end do
    end do

  ! matbox_l stores ELF
  matbox_l=0.d0
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    elfcuni(ix,iy,iz)=3.d0/5.d0*(6.d0*Pi**2)**(2.d0/3.d0)      &
              *rho_half(ix,iy,iz)**(5.d0/3.d0)
    matbox_l(ix,iy,iz)=1.d0/(1.d0+elfc(ix,iy,iz)**2/elfcuni(ix,iy,iz)**2)
  end do
  end do
  end do

  call comm_summation(matbox_l,elf,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)

  if(iSCFRT==1)then 
    suffix = "elf"
  else if(iSCFRT==2)then
    write(filenum, '(i6.6)') itt
    suffix = "elf_"//adjustl(filenum)
  end if
  phys_quantity = "elf"
  if(format_voxel_data=='avs')then
    header_unit = "none"
    call write_avs(lg,103,suffix,header_unit,elf)
  else if(format_voxel_data=='cube')then
    call write_cube(lg,103,suffix,phys_quantity,elf,system%hgs)
  else if(format_voxel_data=='vtk')then
    call write_vtk(lg,103,suffix,elf,system%hgs)
  end if
  
end subroutine write_elf

!===================================================================================================================================

subroutine writeestatic(lg,mg,ng,ex_static,ey_static,ez_static,rmat,rmat2,hgs,itt)
  use inputoutput, only: format_voxel_data
  use structures, only: s_rgrid
  use salmon_parallel, only: nproc_group_global
  use salmon_communication, only: comm_summation
  use write_file3d
  implicit none
  type(s_rgrid),intent(in) :: lg,mg,ng
  real(8),intent(in) :: ex_static(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(in) :: ey_static(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(in) :: ez_static(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8),intent(out) :: rmat(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  real(8),intent(out) :: rmat2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  real(8),intent(in) :: hgs(3)
  integer,intent(in),optional :: itt
  integer :: ix,iy,iz,jj
  character(30) :: suffix
  character(30) :: phys_quantity
  character(10) :: filenum
  character(20) :: header_unit

  do jj=1,3

    !$OMP parallel do collapse(2) private(iz,iy,ix)
    do iz=lg%is(3),lg%ie(3)
    do iy=lg%is(2),lg%ie(2)
    do ix=lg%is(1),lg%ie(1)
      rmat(ix,iy,iz)=0.d0
    end do
    end do
    end do
    
    if(jj==1)then
      !$OMP parallel do collapse(2) private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        rmat(ix,iy,iz)=ex_static(ix,iy,iz)
      end do
      end do
      end do
    else if(jj==2)then
      !$OMP parallel do collapse(2) private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        rmat(ix,iy,iz)=ey_static(ix,iy,iz)
      end do
      end do
      end do
    else if(jj==3)then
      !$OMP parallel do collapse(2) private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        rmat(ix,iy,iz)=ez_static(ix,iy,iz)
      end do
      end do
      end do
    end if
   
    if(format_voxel_data=='avs')then
      !$OMP parallel do collapse(2) private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        rmat(ix,iy,iz)=rmat(ix,iy,iz)*5.14223d1
      end do
      end do
      end do
    end if
 
    call comm_summation(rmat,rmat2,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)
  
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
      call write_cube(lg,103,suffix,phys_quantity,rmat2,hgs)
    else if(format_voxel_data=='vtk')then
      call write_vtk(lg,103,suffix,rmat2,hgs)
    end if

  end do  
 
end subroutine writeestatic

!===================================================================================================================================

subroutine write_psi(lg,mg,system,info,spsi)
  use inputoutput, only: au_length_aa
  use structures
  use salmon_global, only: format_voxel_data
  use salmon_communication, only: comm_summation
  use write_file3d
  implicit none
  type(s_rgrid),intent(in) :: lg,mg
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital),intent(in) :: spsi
  !
  integer :: io,ik,ispin,ix,iy,iz
  complex(8),dimension(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)) :: cmatbox,cmatbox2
  character(30) :: suffix
  character(30) :: phys_quantity
  character(10) :: filenum
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
        
!$omp parallel do collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          cmatbox(ix,iy,iz) = spsi%zwf(ix,iy,iz,ispin,io,ik,1) ! future work: rwf
        end do
        end do
        end do
      end if
      call comm_summation(cmatbox,cmatbox2,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)

      write(filenum, '(i5)') io ! future work: ik,ispin
      suffix = "psi"//trim(adjustl(filenum))
      phys_quantity = "psi"
      if(format_voxel_data=='avs')then
        header_unit = "A**(-3/2)"
      ! future work: real & imaginary part
        call write_avs(lg,103,suffix,header_unit,dble(cmatbox2)/sqrt(au_length_aa)**3)
        call write_avs(lg,103,suffix,header_unit,aimag(cmatbox2)/sqrt(au_length_aa)**3)
      else if(format_voxel_data=='cube')then
      ! future work: real & imaginary part
        call write_cube(lg,103,suffix,phys_quantity,dble(cmatbox2),system%hgs)
        call write_cube(lg,103,suffix,phys_quantity,aimag(cmatbox2),system%hgs)
      end if
      
  end do
  end do
  end do

end subroutine write_psi

!===================================================================================================================================

end module writefield
