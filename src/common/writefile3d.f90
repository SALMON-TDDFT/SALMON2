!
!  Copyright 2017 SALMON developers
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
module writefile3d
  implicit none

contains

subroutine writevtk(lg,fp,suffix,rmat,hgs,igc_is,igc_ie,gridcoo)
  use inputoutput, only: au_length_aa
  use structures, only: s_rgrid
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root
  implicit none
  type(s_rgrid),intent(in) :: lg
  integer,intent(in) :: fp
  character(30),intent(in):: suffix
  real(8),intent(in) :: rmat(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),  &
                             lg%is(3):lg%ie(3))
  real(8),intent(in) :: hgs(3)
  integer,intent(in) :: igc_is,igc_ie
  real(8),intent(in) :: gridcoo(igc_is:igc_ie,3)
  character(30):: filename
  integer :: ix,iy,iz

  if(comm_is_root(nproc_id_global))then
    filename=trim(suffix)//".vtk"
    open(fp,file=filename)

    write(fp, '(A)') "# vtk DataFile Version 3.0"
    write(fp, '(A)') "vtk output"
    write(fp, '(A)') "ASCII"
    write(fp, '(A)') "DATASET STRUCTURED_POINTS"
    write(fp, '(A,3(1X,I2))') "DIMENSIONS", lg%num(1), lg%num(2), lg%num(3)
    write(fp, '(A,3(1X,F10.5))') "ORIGIN",gridcoo(lg%is(1),1)*au_length_aa,  &
                                          gridcoo(lg%is(2),2)*au_length_aa,  &
                                          gridcoo(lg%is(3),3)*au_length_aa
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

end subroutine writevtk

!======================================================================

subroutine writecube(lg,fp,suffix,phys_quantity,rmat,hgs,igc_is,igc_ie,gridcoo)
  use inputoutput, only: natom,kion,rion,izatom
  use structures, only: s_rgrid
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root
  implicit none
  type(s_rgrid),intent(in) :: lg
  integer, intent(IN) :: fp
  character(30),intent(in):: suffix
  character(30),intent(in):: phys_quantity
  real(8),intent(IN) :: rmat(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),  &
                              lg%is(3):lg%ie(3))
  real(8),intent(in) :: hgs(3)
  integer,intent(in) :: igc_is,igc_ie
  real(8),intent(in) :: gridcoo(igc_is:igc_ie,3)
  character(30):: filename
  integer :: j,iatom
  integer :: ix,iy,iz
  integer :: ik

  if(comm_is_root(nproc_id_global))then
    filename=trim(suffix)//".cube"
    open(fp,file=filename)
    if(phys_quantity=="psi")then
      write(fp,*) "Molecular Orbital"
    else if(phys_quantity=="dns")then
      write(fp,*) "Electron Density"
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
    write(fp,'(i5,3f12.6)') natom,gridcoo(lg%is(1),1),gridcoo(lg%is(2),2),gridcoo(lg%is(3),3)
    write(fp,'(i5,3f12.6)') lg%num(1),hgs(1),0.d0,0.d0
    write(fp,'(i5,3f12.6)') lg%num(2),0.d0,hgs(2),0.d0
    write(fp,'(i5,3f12.6)') lg%num(3),0.d0,0.d0,hgs(3)
    do iatom=1,natom
      ik=Kion(iatom)
      write(fp,'(i5,4f12.6)') izatom(ik),dble(izatom(ik)),(rion(j,iatom),j=1,3)
    end do

    do ix=lg%is(1),lg%ie(1)
    do iy=lg%is(2),lg%ie(2)
      write(fp,'(6e13.5)', advance="yes") (rmat(ix,iy,iz),iz=lg%is(3),lg%ie(3))
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

end subroutine writecube

end module writefile3d
