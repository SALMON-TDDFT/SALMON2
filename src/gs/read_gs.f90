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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module read_gs
  implicit none

contains

subroutine read_dns(lg,mg,rho)
  use structures
  use salmon_global, only: natom
  use salmon_parallel
  use salmon_communication
  implicit none
  type(s_rgrid),intent(in) :: lg,mg
  real(8) :: rho(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  !
  character(8),parameter :: filename="dns.cube"
  integer,parameter :: fp=103
  integer :: num(3),iatom,ix,iy,iz
  real(8),allocatable :: tmp(:,:,:)
  allocate(tmp(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  if(comm_is_root(nproc_id_global))then
    write(*,*) "read GS density dns.cube"
    open(fp,file=filename)
    read(fp,*)
    read(fp,*)
    read(fp,*) iatom
    if(iatom/=natom) stop "error @ read_rho: natom"
    read(fp,*) num(1)
    read(fp,*) num(2)
    read(fp,*) num(3)
    if(num(1)/=lg%num(1) .or. num(2)/=lg%num(2) .or. num(3)/=lg%num(3)) stop "error @ read_rho: lg"
    do iatom=1,natom
      read(fp,'(i5,4f12.6)')
    end do
    do ix=lg%is(1),lg%ie(1)
    do iy=lg%is(2),lg%ie(2)
      read(fp,*) (tmp(ix,iy,iz),iz=lg%is(3),lg%ie(3))
    end do
    end do
    close(fp)
  end if
  call comm_bcast(tmp,nproc_group_global)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    rho(ix,iy,iz) = tmp(ix,iy,iz)
  enddo
  enddo
  enddo
  deallocate(tmp)
  return
end subroutine read_dns

subroutine read_wfn(lg,mg,psi,info,system,k_rd)
  use structures
  use salmon_parallel
  use salmon_communication
  implicit none
  type(s_rgrid),intent(in) :: lg,mg
  type(s_system) ,intent(in) :: system
  type(s_wf_info),intent(in) :: info
  real(8),intent(in) :: k_rd(3,system%nk)
  type(s_wavefunction) :: psi
  !
  character(7),parameter :: filename="wfn.bin"
  integer,parameter :: fp=104
  integer :: im,ik,io,ispin,ix,iy,iz,ik_s,ik_e,io_s,io_e,no,nk,no0,nk0,ik0,i
  real(8) :: kk(3),dk(3),r,rmin
  real(8),allocatable :: k_tmp(:,:)
  complex(8),allocatable :: tmp1(:,:,:,:,:)

  if(info%im_s/=1 .or. info%im_e/=1 .or. system%nspin/=1) stop "error @ write_wfn: im_s, im_e, nspin"
  im = 1
  ispin = 1
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  no = system%no
  nk = system%nk

  if(comm_is_root(nproc_id_global))then
    write(*,*) "read GS wavefunction wfn.bin"
    open(fp,file=filename,form='unformatted')
    read(fp) ix,iy,iz,io,ik
    if(ix/=lg%num(1) .or. iy/=lg%num(2) .or. iz/=lg%num(3)) stop "error @ read_wfn: array size"
    if(io < no) stop "error @ read_wfn: # of orbitals"
    if(ik > nk) stop "error @ read_wfn: # of k-vectors"
    no0 = io
    nk0 = ik
    allocate(tmp1(lg%num(1),lg%num(2),lg%num(3),no0,nk0),k_tmp(3,nk0))
    read(fp) k_tmp
    read(fp) tmp1
    close(fp)
  end if
  call comm_bcast(no0,info%icomm_rko)
  call comm_bcast(nk0,info%icomm_rko)
  if(.not. comm_is_root(nproc_id_global)) allocate(tmp1(lg%num(1),lg%num(2),lg%num(3),no0,nk0),k_tmp(3,nk0))
  call comm_bcast(k_tmp,info%icomm_rko)
  call comm_bcast(tmp1,info%icomm_rko)

  do ik=ik_s,ik_e
    ik0 = 1
    rmin = sum(system%brl**2)
    do i=1,nk0
      kk = k_rd(:,ik) - k_tmp(:,i)
      do ix=-2,2
      do iy=-2,2
      do iz=-2,2
        dk(1) = ix*system%brl(1,1) + iy*system%brl(1,2) + iz*system%brl(1,3)
        dk(2) = ix*system%brl(2,1) + iy*system%brl(2,2) + iz*system%brl(2,3)
        dk(3) = ix*system%brl(3,1) + iy*system%brl(3,2) + iz*system%brl(3,3)
        dk = dk + kk
        r = dk(1)**2 + dk(2)**2 + dk(3)**2
        if(rmin > r) then
          rmin = r
          ik0 = i
        end if
      enddo
      enddo
      enddo
    end do
    do io=io_s,io_e
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      psi%zwf(ix,iy,iz,ispin,io,ik,im) = tmp1(ix,iy,iz,io,ik0)
    enddo
    enddo
    enddo
    enddo
  enddo

  deallocate(tmp1,k_tmp)
  return
end subroutine read_wfn

subroutine write_wfn(lg,mg,psi,info,system,k_rd)
  use structures
  use salmon_parallel
  use salmon_communication
  implicit none
  type(s_rgrid),intent(in) :: lg,mg
  type(s_system) ,intent(in) :: system
  type(s_wf_info),intent(in) :: info
  real(8),intent(in) :: k_rd(3,system%nk)
  type(s_wavefunction),intent(in) :: psi
  !
  character(7),parameter :: filename="wfn.bin"
  integer,parameter :: fp=104
  integer :: im,ik,io,ispin,ix,iy,iz,ik_s,ik_e,io_s,io_e,no,nk
  complex(8),allocatable :: tmp1(:,:,:,:,:),tmp2(:,:,:,:,:)

  if(info%im_s/=1 .or. info%im_e/=1 .or. system%nspin/=1) stop "error @ write_wfn: im_s, im_e, nspin"
  im = 1
  ispin = 1
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  no = system%no
  nk = system%nk

  allocate(tmp1(lg%num(1),lg%num(2),lg%num(3),no,nk), &
           tmp2(lg%num(1),lg%num(2),lg%num(3),no,nk))
  tmp1 = 0d0
  do ik=ik_s,ik_e
  do io=io_s,io_e
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    tmp1(ix,iy,iz,io,ik) = psi%zwf(ix,iy,iz,ispin,io,ik,im)
  enddo
  enddo
  enddo
  enddo
  enddo
  call comm_summation(tmp1,tmp2,system%ngrid*no*nk,info%icomm_rko)

  if(comm_is_root(nproc_id_global))then
    write(*,*) "write GS wavefunction wfn.bin"
    open(fp,file=filename,form='unformatted')
    write(fp) lg%num(1),lg%num(2),lg%num(3),no,nk
    write(fp) k_rd
    write(fp) tmp2
    close(fp)
  end if

  deallocate(tmp1,tmp2)
  return
end subroutine write_wfn

end module read_gs
