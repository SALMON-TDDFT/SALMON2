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
module salmon_initialization
  implicit none

contains

!===================================================================================================================================

subroutine init_dft(lg,system,stencil)
  use structures
  use lattice
  use salmon_global, only: al_vec1,al_vec2,al_vec3,al,ispin,natom,iperiodic,num_kgrid,num_rgrid,dl
  implicit none
  type(s_rgrid)      :: lg
  type(s_dft_system) :: system
  type(s_stencil)    :: stencil
  !
  real(8) :: rsize(3)

  if(al_vec1(2)==0d0 .and. al_vec1(3)==0d0 .and. al_vec2(1)==0d0 .and. &
     al_vec2(3)==0d0 .and. al_vec3(1)==0d0 .and. al_vec3(2)==0d0) then
    stencil%if_orthogonal = .true.
    system%primitive_a = 0d0
    system%primitive_a(1,1) = al(1)
    system%primitive_a(2,2) = al(2)
    system%primitive_a(3,3) = al(3)
    rsize = al
  else
    stencil%if_orthogonal = .false.
    system%primitive_a(1:3,1) = al_vec1
    system%primitive_a(1:3,2) = al_vec2
    system%primitive_a(1:3,3) = al_vec3
    rsize(1) = sqrt(sum(al_vec1**2))
    rsize(2) = sqrt(sum(al_vec2**2))
    rsize(3) = sqrt(sum(al_vec3**2))
  end if

  if(sum(abs(dl)) == 0d0) then
    system%hgs(1:3) = rsize(1:3) / dble(num_rgrid(1:3))
  else
    system%hgs(1:3) = dl(1:3)
  end if
  call init_grid_whole(rsize,system%hgs,lg)
  system%ngrid = lg%num(1) * lg%num(2) * lg%num(3)

  call init_lattice(system,stencil)
  call init_kvector(num_kgrid,system)

  system%iperiodic = iperiodic
!  system%no    = itotMST ! (future work)
  system%nion  = natom

  if(ispin==0)then
    system%nspin=1
  else
    system%nspin=2
  end if

!  allocate(system%Rion(3,system%nion) & ! (future work)
!          ,system%wtk(system%nk) &
!          ,system%rocc(system%no,system%nk,system%nspin))

  return
end subroutine init_dft

subroutine init_grid_whole(rsize,hgs,lg)
  use structures
  use salmon_global, only: nproc_domain,iperiodic,dl,num_rgrid
  implicit none
  real(8),intent(in) :: rsize(3),hgs(3)
  type(s_rgrid) :: lg
  !
  integer,parameter :: Nd=4
  real(8),parameter :: epsilon=1.d-10
  !
  integer :: lg_sta(3),lg_end(3),lg_num(3)
  integer :: j

  lg%ndir = 3 ! high symmetry nonorthogonal lattice is not implemented
  lg%nd = Nd

  select case(iperiodic)
    case(0)
      lg_end(:)=int((rsize(:)+epsilon)/2.d0/Hgs(:))
      do j=1,3
        if(mod(int(rsize(j)/Hgs(j)+1.d-12),2)==1)then
          lg_sta(j)=-(int((rsize(j)+epsilon)/2.d0/Hgs(j)))
        else
          lg_sta(j)=-(int((rsize(j)+epsilon)/2.d0/Hgs(j)))+1
        end if
      end do
    case(3)
      lg_sta(:)=1
      lg_end(:)=int((rsize(:)+epsilon)/Hgs(:))
  end select
  lg_num(:)=lg_end(:)-lg_sta(:)+1

  lg%is(1:3) = lg_sta(1:3)
  lg%ie(1:3) = lg_end(1:3)
  lg%num(1:3) = lg_num(1:3)
  lg%is_overlap(1:3) = lg_sta(1:3)-nd
  lg%ie_overlap(1:3) = lg_end(1:3)+nd

  allocate(lg%idx(lg%is_overlap(1):lg%ie_overlap(1)) &
          ,lg%idy(lg%is_overlap(2):lg%ie_overlap(2)) &
          ,lg%idz(lg%is_overlap(3):lg%ie_overlap(3)))

  if(iperiodic==3 .and. nproc_domain(1)*nproc_domain(2)*nproc_domain(3)==1) then
    lg%is_array(1:3) = lg_sta(1:3)
    lg%ie_array(1:3) = lg_end(1:3)
    do j=lg%is_overlap(1),lg%ie_overlap(1)
      lg%idx(j) = mod(j+lg%num(1)-1,lg%num(1))+1
    end do
    do j=lg%is_overlap(2),lg%ie_overlap(2)
      lg%idy(j) = mod(j+lg%num(2)-1,lg%num(2))+1
    end do
    do j=lg%is_overlap(3),lg%ie_overlap(3)
      lg%idz(j) = mod(j+lg%num(3)-1,lg%num(3))+1
    end do
  else
    lg%is_array(1:3)=lg_sta(1:3)-nd
    lg%ie_array(1:3)=lg_end(1:3)+nd
    do j=lg%is_overlap(1),lg%ie_overlap(1)
      lg%idx(j) = j
    end do
    do j=lg%is_overlap(2),lg%ie_overlap(2)
      lg%idy(j) = j
    end do
    do j=lg%is_overlap(3),lg%ie_overlap(3)
      lg%idz(j) = j
    end do
  end if

  if(sum(abs(dl)) == 0d0) then
    if( maxval(abs(num_rgrid-lg%num)) > 0) stop "error: num_rgrid /= lg%num"
  else
    if( maxval(abs((rsize/dl)-dble(lg%num))) > 1d-15 ) stop "error: abs((rsize/dl)-dble(lg%num)) is too large"
  end if

  return
end subroutine init_grid_whole

end module salmon_initialization
