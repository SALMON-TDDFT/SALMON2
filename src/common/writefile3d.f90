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
    write(fp, '(A,3(1X,F10.5))') "SPACING", Hgs(1)*au_length_aa,  &
                                            Hgs(2)*au_length_aa,  &
                                            Hgs(3)*au_length_aa
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

