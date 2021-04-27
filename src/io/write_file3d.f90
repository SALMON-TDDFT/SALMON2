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
!======================================================================
!======================================================================
module write_file3d
  implicit none

contains

subroutine write_avs(lg,fp,suffix,header_unit,rmat)
  use inputoutput, only: nsplit_voxel_data
  use parallelization, only: nproc_id_global
  use structures, only: s_rgrid
  use communication, only: comm_is_root
  implicit none
  type(s_rgrid),intent(in) :: lg
  integer,intent(in) :: fp
  character(20),intent(in) :: header_unit
  real(8),intent(in) :: rmat(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),  &
                             lg%is(3):lg%ie(3))
  character(60),intent(in):: suffix
  !
  character(60):: filename
  integer::ix,iy,iz
  integer::jj,icount
  integer::jsta,jend
  integer :: icoo1d(3,lg%num(1)*lg%num(2)*lg%num(3))
  character(8)  :: filenumber_data
  
  do iz=lg%is(3),lg%ie(3)
  do iy=lg%is(2),lg%ie(2)
  do ix=lg%is(1),lg%ie(1)
    icount=(iz-lg%is(3))*lg%num(2)*lg%num(1)+(iy-lg%is(2))*lg%num(1)+ix-lg%is(1)+1
    icoo1d(1,icount)=ix
    icoo1d(2,icount)=iy
    icoo1d(3,icount)=iz
  end do
  end do
  end do

  if(nsplit_voxel_data>=2)then
    if(nproc_id_global<nsplit_voxel_data)then
      write(filenumber_data, '(i8)') nproc_id_global
      filename = trim(suffix)//"."//adjustl(filenumber_data)
      open(fp,file=filename)
      if(comm_is_root(nproc_id_global))then
        write(fp,'("# unit is ",a)') header_unit
      end if
      jsta=nproc_id_global*(lg%num(1)*lg%num(2)*lg%num(3))/nsplit_voxel_data+1
      jend=(nproc_id_global+1)*(lg%num(1)*lg%num(2)*lg%num(3))/nsplit_voxel_data
      do jj=jsta,jend
        if(abs(rmat(icoo1d(1,jj),icoo1d(2,jj),icoo1d(3,jj)))>=1.0d-10)then
          write(fp,'(e20.8)') rmat(icoo1d(1,jj),icoo1d(2,jj),icoo1d(3,jj))
        else
          write(fp,'(a1)') "0"
        end if
      enddo
      close(1)
    end if
  else
    if(comm_is_root(nproc_id_global))then
      filename=trim(suffix)
      open(fp,file=filename)
      write(fp,'("# unit is ",a)') header_unit
      do iz=lg%is(3),lg%ie(3)
      do iy=lg%is(2),lg%ie(2)
      do ix=lg%is(1),lg%ie(1)
        if(abs(rmat(ix,iy,iz))>=1.0d-10)then
          write(fp,'(e20.8)') rmat(ix,iy,iz)
        else
          write(fp,'(a1)') "0"
        end if
      enddo
      enddo
      enddo
      close(fp)
    endif
  endif
  
end subroutine write_avs

!======================================================================

subroutine write_cube(lg,fp,suffix,phys_quantity,rmat,system)
  use salmon_global, only: natom,kion,izatom,yn_jm
  use structures, only: s_rgrid, s_dft_system
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root
  implicit none
  type(s_rgrid)     ,intent(in) :: lg
  type(s_dft_system),intent(in) :: system
  integer           ,intent(in) :: fp
  character(60)     ,intent(in) :: suffix
  character(30)     ,intent(in) :: phys_quantity
  real(8)           ,intent(in) :: rmat(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  !
  character(60):: filename
  integer :: j,iatom
  integer :: ix,iy,iz
  integer :: ik
  real(8) :: daa(3,3)
  
  do j=1,3
     daa(1:3,j) = system%primitive_a(1:3,j)/dble(lg%num(j))
  end do

  if(comm_is_root(nproc_id_global))then
    filename=trim(suffix)//".cube"
    open(fp,file=filename)
    if(phys_quantity=="psi")then
      write(fp,*) "Molecular Orbital"
    else if(phys_quantity=="dns")then
      write(fp,*) "Electron Density"
    else if(phys_quantity=="pbcd")then
      write(fp,*) "Positive Background Charge Density"
    else if(phys_quantity=="elf")then
      write(fp,*) "Electron Localization Function"
    else if(phys_quantity=="dnsdiff")then
      write(fp,*) "Difference of Electron Density"
    else if(phys_quantity=="exsta")then
      write(fp,*) "x Component of Static Electric Field"
    else if(phys_quantity=="eysta")then
      write(fp,*) "y Component of Static Electric Field"
    else if(phys_quantity=="ezsta")then
      write(fp,*) "z Component of Static Electric Field"
    end if
    write(fp,*) "All values here are in a.u."
    write(fp,'(i5,3f12.6)') natom,lg%coordinate(lg%is(1),1),lg%coordinate(lg%is(2),2),lg%coordinate(lg%is(3),3)
    write(fp,'(i5,3f12.6)') lg%num(1),daa(1:3,1)
    write(fp,'(i5,3f12.6)') lg%num(2),daa(1:3,2)
    write(fp,'(i5,3f12.6)') lg%num(3),daa(1:3,3)
    if(yn_jm=='n')then
      do iatom=1,natom
        ik=Kion(iatom)
        write(fp,'(i5,4f12.6)') izatom(ik),dble(izatom(ik)),(system%Rion(j,iatom),j=1,3)
      end do
    else
      write(fp,'(i5,4f12.6)') 1,1.0d0,0.0d0,0.0d0,0.0d0
    end if

    do ix=lg%is(1),lg%ie(1)
    do iy=lg%is(2),lg%ie(2)
      write(fp,'(6(1X,E23.15E3))', advance="yes") (rmat(ix,iy,iz),iz=lg%is(3),lg%ie(3))
!    do iz=lg%is(3),lg%ie(3)
!      if(mod(iz+1-lg%is(3),6)==0)then
!        write(fp,'(e13.5)', advance="yes") abs(rmat(ix,iy,iz))
!      else
!        write(fp,'(e13.5)', advance="no") abs(rmat(ix,iy,iz))
!      endif
!    end do
!    write(fp,*)
    end do
    end do
    close(fp)
  end if

end subroutine write_cube

!======================================================================

subroutine write_vtk(lg,fp,suffix,rmat,hgs)
  use inputoutput, only: au_length_aa
  use structures, only: s_rgrid
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root
  implicit none
  type(s_rgrid),intent(in) :: lg
  integer,intent(in) :: fp
  character(60),intent(in):: suffix
  real(8),intent(in) :: rmat(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),  &
                             lg%is(3):lg%ie(3))
  real(8),intent(in) :: hgs(3)
  character(60):: filename
  integer :: ix,iy,iz

  if(comm_is_root(nproc_id_global))then
    filename=trim(suffix)//".vtk"
    open(fp,file=filename)

    write(fp, '(A)') "# vtk DataFile Version 3.0"
    write(fp, '(A)') "vtk output"
    write(fp, '(A)') "ASCII"
    write(fp, '(A)') "DATASET STRUCTURED_POINTS"
    write(fp, '(A,3(1X,I2))') "DIMENSIONS", lg%num(1), lg%num(2), lg%num(3)
    write(fp, '(A,3(1X,F10.5))') "ORIGIN",lg%coordinate(lg%is(1),1)*au_length_aa,  &
                                          lg%coordinate(lg%is(2),2)*au_length_aa,  &
                                          lg%coordinate(lg%is(3),3)*au_length_aa
    write(fp, '(A,3(1X,F10.5))') "SPACING", hgs(1)*au_length_aa,  &
                                            hgs(2)*au_length_aa,  &
                                            hgs(3)*au_length_aa
    write(fp, '(A,1X,I6)') "POINT_DATA", lg%num(1) * lg%num(2) * lg%num(3)
    write(fp, '(A)') "SCALARS scalars float"
    write(fp, '(A)') "LOOKUP_TABLE default"


    do iz=lg%is(3),lg%ie(3)
    do iy=lg%is(2),lg%ie(2)
    do ix=lg%is(1),lg%ie(1)
      write(fp, '(ES12.5)') rmat(ix,iy,iz)
    end do
    end do
    end do
  
    close(fp)
  end if

end subroutine write_vtk

!======================================================================

end module write_file3d
