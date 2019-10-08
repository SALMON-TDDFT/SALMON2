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
!======================================================================
!======================================================================
subroutine write_psi(lg,mg,system,info,spsi)
  use inputoutput, only: au_length_aa
  use structures
  use salmon_global, only: format_voxel_data
  use salmon_communication, only: comm_summation
  use write_file3d
  use scf_data, only: icoo1d
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
        call write_avs(lg,103,suffix,header_unit,dble(cmatbox2)/sqrt(au_length_aa)**3,icoo1d)
        call write_avs(lg,103,suffix,header_unit,aimag(cmatbox2)/sqrt(au_length_aa)**3,icoo1d)
      else if(format_voxel_data=='cube')then
      ! future work: real & imaginary part
        call write_cube(lg,103,suffix,phys_quantity,dble(cmatbox2),system%hgs)
        call write_cube(lg,103,suffix,phys_quantity,aimag(cmatbox2),system%hgs)
      end if
      
end do
end do
end do

end subroutine write_psi

