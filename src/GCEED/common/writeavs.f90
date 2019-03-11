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
subroutine writeavs(lg,fp,suffix,header_unit,rmat,icoo1d)
  use inputoutput, only: numfiles_out_3d
  use salmon_parallel, only: nproc_id_global
  use structures, only: s_rgrid
  use salmon_communication, only: comm_is_root
  implicit none
  type(s_rgrid),intent(in) :: lg
  integer,intent(in) :: fp
  character(20),intent(in) :: header_unit
  real(8),intent(in) :: rmat(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),  &
                                   lg%is(3):lg%ie(3))
  integer,intent(in) :: icoo1d(3,lg%num(1)*lg%num(2)*lg%num(3))
  character(30),intent(in):: suffix
  character(30):: filename
  integer::ix,iy,iz
  integer::jj
  integer::jsta,jend
  character(8)  :: filenumber_data

 
  if(numfiles_out_3d>=2)then
    if(nproc_id_global<numfiles_out_3d)then
      write(filenumber_data, '(i8)') nproc_id_global
      filename = trim(suffix)//"."//adjustl(filenumber_data)
      open(fp,file=filename)
      if(comm_is_root(nproc_id_global))then
        write(fp,'("# unit is ",a)') header_unit
      end if
      jsta=nproc_id_global*(lg%num(1)*lg%num(2)*lg%num(3))/numfiles_out_3d+1
      jend=(nproc_id_global+1)*(lg%num(1)*lg%num(2)*lg%num(3))/numfiles_out_3d
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
  
end subroutine writeavs

