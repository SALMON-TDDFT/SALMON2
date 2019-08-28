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
  integer,parameter,private :: Nd=4

contains

!===================================================================================================================================

subroutine init_dft(lg,system,stencil)
  use structures
  use lattice
  use salmon_global, only: al_vec1,al_vec2,al_vec3,al,ispin,natom,nstate &
  & ,iperiodic,num_kgrid,num_rgrid,dl,nproc_domain_orbital,rion
  implicit none
  type(s_rgrid)      :: lg
  type(s_dft_system) :: system
  type(s_stencil)    :: stencil
  !
  integer :: ii,jj
  real(8) :: rsize(3),hgs(3),cnmat(0:12,12),bnmat(4,4)

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
    hgs(1:3) = rsize(1:3) / dble(num_rgrid(1:3))
  else
    hgs(1:3) = dl(1:3)
  end if
  call init_grid_whole(rsize,hgs,lg)
  system%hgs = hgs
  system%ngrid = lg%num(1) * lg%num(2) * lg%num(3)

  call init_lattice(system,stencil)
  call init_kvector(num_kgrid,system)

  system%iperiodic = iperiodic
  system%no   = nstate
  system%nion = natom

  if(ispin==0)then
    system%nspin=1
  else
    system%nspin=2
  end if

  allocate(system%Rion(3,system%nion),system%rocc(system%no,system%nk,system%nspin))
  allocate(system%Velocity(3,system%nion),system%Force(3,system%nion))
  system%rion = rion
  system%rocc = 0d0 ! initial value

  call setbn(bnmat)
  call setcn(cnmat)
  if(stencil%if_orthogonal) then
    stencil%coef_lap0 = -0.5d0*cNmat(0,Nd)*(1.d0/Hgs(1)**2+1.d0/Hgs(2)**2+1.d0/Hgs(3)**2)
  else
    if(nproc_domain_orbital(1)*nproc_domain_orbital(2)*nproc_domain_orbital(3)/=1) &
                                                    stop "error: nonorthogonal lattice and r-space parallelization"
    stencil%coef_lap0 = -0.5d0*cNmat(0,Nd)*  &
                      & ( stencil%coef_F(1)/Hgs(1)**2 + stencil%coef_F(2)/Hgs(2)**2 + stencil%coef_F(3)/Hgs(3)**2 )
  end if
  do jj=1,3
    do ii=1,4
      stencil%coef_lap(ii,jj) = cnmat(ii,4)/hgs(jj)**2
      stencil%coef_nab(ii,jj) = bnmat(ii,4)/hgs(jj)
    end do
  end do

  return
end subroutine init_dft

!===================================================================================================================================

subroutine init_orbital_parallel_singlecell(system,info)
  use structures
  use salmon_global, only: nproc_k,nproc_ob,nproc_domain_orbital
  use salmon_parallel, only: nproc_group_global
use  calc_iobnum_sub ! remove this line
use calc_iroot_sub! remove this line
use scf_data! remove this line
use calc_allob_sub! remove this line
  implicit none
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel)      :: info
  !
  integer :: io,jspin,ik
  integer :: jj,iob

! for single-cell calculations
  info%im_s = 1
  info%im_e = 1
  info%numm = 1

! # of k points
  info%ik_s = info%id_k * system%nk / nproc_k + 1
  info%ik_e = (info%id_k+1) * system%nk / nproc_k
  info%numk = info%ik_e - info%ik_s + 1

! # of orbitals
  info%io_s = 1
  info%io_e = (info%id_o+1) * system%no / nproc_ob - info%id_o * system%no / nproc_ob
  info%numo = info%io_e - info%io_s + 1
!  info%io_s = info%id_o * system%no / nproc_ob + 1
!  info%io_e = (info%id_o+1) * system%no / nproc_ob
!  info%numo = info%io_e - info%io_s + 1

! for parallelization
  info%if_divide_rspace = nproc_domain_orbital(1)*nproc_domain_orbital(2)*nproc_domain_orbital(3).ne.1
  info%if_divide_orbit  = nproc_ob.ne.1
  info%icomm_rko  = nproc_group_global

  allocate(info%occ(info%io_s:info%io_e, info%ik_s:info%ik_e, 1:system%nspin,1) &
            ,info%io_tbl(info%io_s:info%io_e), info%jo_tbl(1:system%no) &
            ,info%irank_jo(1:system%no))

  info%jo_tbl(:) = 0 !(initial value)
  do iob=info%io_s,info%io_e
    call calc_allob(iob,jj,itotmst,mst,info%numo)
    info%io_tbl(iob) = jj
    info%jo_tbl(jj)  = iob
  end do

  do jspin=1,system%nspin
    do ik=info%ik_s,info%ik_e
      do iob=info%io_s,info%io_e
        jj = info%io_tbl(iob)
        info%occ(iob,ik,jspin,1) = system%rocc(jj,ik,jspin)*system%wtk(ik)
      end do
    end do
  end do

  do jj=1, system%no
    call calc_iroot(jj,info%irank_jo(jj),ilsda,nproc_ob,itotmst,mst)
  end do

!  allocate(info%occ(info%io_s:info%io_e, info%ik_s:info%ik_e, 1:system%nspin,1),info%irank_io(1:system%no))
!
!! process ID corresponding to the orbital index io
!  do io=1, system%no
!    if(mod(io*nproc_ob,system%no)==0)then
!      info%irank_io(io) = io*nproc_ob/system%no - 1
!    else
!      info%irank_io(io) = io*nproc_ob/system%no
!    end if
!  end do
!
!! occupation number (initaial value)
!  do jspin=1,system%nspin
!    do ik=info%ik_s,info%ik_e
!      do io=info%io_s,info%io_e
!        info%occ(io,ik,jspin,1) = system%rocc(io,ik,jspin)*system%wtk(ik) ! initial value
!      end do
!    end do
!  end do

end subroutine init_orbital_parallel_singlecell

!===================================================================================================================================

subroutine init_grid_whole(rsize,hgs,lg)
  use structures
  use salmon_global, only: nproc_domain_orbital,iperiodic,dl,num_rgrid
  implicit none
  real(8),intent(in) :: rsize(3),hgs(3)
  type(s_rgrid) :: lg
  !
  real(8),parameter :: epsilon=1.d-10
  integer :: j

  lg%ndir = 3 ! high symmetry nonorthogonal lattice is not implemented
  lg%nd = Nd

  select case(iperiodic)
    case(0)
      lg%ie(:)=int((rsize(:)+epsilon)/2.d0/Hgs(:))
      do j=1,3
        if(mod(int(rsize(j)/Hgs(j)+1.d-12),2)==1)then
          lg%is(j)=-(int((rsize(j)+epsilon)/2.d0/Hgs(j)))
        else
          lg%is(j)=-(int((rsize(j)+epsilon)/2.d0/Hgs(j)))+1
        end if
      end do
    case(3)
      lg%is(:)=1
      lg%ie(:)=int((rsize(:)+epsilon)/Hgs(:))
  end select
  lg%num(:)=lg%ie(:)-lg%is(:)+1

  lg%is_overlap(1:3) = lg%is(1:3)-nd
  lg%ie_overlap(1:3) = lg%ie(1:3)+nd

  allocate(lg%idx(lg%is_overlap(1):lg%ie_overlap(1)) &
          ,lg%idy(lg%is_overlap(2):lg%ie_overlap(2)) &
          ,lg%idz(lg%is_overlap(3):lg%ie_overlap(3)))

  if(iperiodic==3 .and. nproc_domain_orbital(1)*nproc_domain_orbital(2)*nproc_domain_orbital(3)==1) then
    lg%is_array(1:3) = lg%is(1:3)
    lg%ie_array(1:3) = lg%ie(1:3)
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
    lg%is_array(1:3)=lg%is(1:3)-nd
    lg%ie_array(1:3)=lg%ie(1:3)+nd
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

  if(sum(abs(dl)) <= 1d-12) then
    if( maxval(abs(num_rgrid-lg%num)) > 0) stop "error: num_rgrid /= lg%num"
  else
    if( maxval(abs((rsize/dl)-dble(lg%num))) > 1d-4 ) stop "error: abs((rsize/dl)-dble(lg%num)) is too large"
  end if

  return
end subroutine init_grid_whole

!===================================================================================================================================

subroutine init_grid_parallel(lg,mg,ng)
  use salmon_communication, only: comm_is_root
  use salmon_parallel, only: nproc_id_global,nproc_size_global
  use salmon_global, only: calc_mode,iperiodic,process_allocation,nproc_domain_orbital,nproc_domain_general,nproc_k,nproc_ob
  use structures, only: s_rgrid
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid)            :: mg,ng
  !
  integer :: i1,i2,i3,i4,j1,j2,j3,ibox,j,ii
  integer :: nproc,myrank
  integer :: nproc_domain_orbital_mul,ngo(3),ngo_mul

  myrank = nproc_id_global
  nproc = nproc_size_global

  allocate(mg%is_all(3,0:nproc-1),mg%ie_all(3,0:nproc-1),ng%is_all(3,0:nproc-1),ng%ie_all(3,0:nproc-1))

! +-------------------------------+
! | mg: r-space grid for orbitals |
! +-------------------------------+

  mg%ndir = 3 ! high symmetry nonorthogonal lattice is not implemented
  mg%nd = Nd

  nproc_domain_orbital_mul = nproc_domain_orbital(1)*nproc_domain_orbital(2)*nproc_domain_orbital(3)
  if(process_allocation=='orbital_sequential')then
    do j3=0,nproc_domain_orbital(3)-1
    do j2=0,nproc_domain_orbital(2)-1
    do j1=0,nproc_domain_orbital(1)-1
      do i2=0,nproc_k-1
      do i1=0,nproc_ob-1
        ibox = nproc_ob*i2 + i1 + nproc_k*nproc_ob*( j1 + nproc_domain_orbital(1)*j2 &
        & + nproc_domain_orbital(1)*nproc_domain_orbital(2)*j3 )
        mg%is_all(1,ibox) = j1*lg%num(1)/nproc_domain_orbital(1)+lg%is(1)
        mg%ie_all(1,ibox) = (j1+1)*lg%num(1)/nproc_domain_orbital(1)+lg%is(1)-1
        mg%is_all(2,ibox) = j2*lg%num(2)/nproc_domain_orbital(2)+lg%is(2)
        mg%ie_all(2,ibox) = (j2+1)*lg%num(2)/nproc_domain_orbital(2)+lg%is(2)-1
        mg%is_all(3,ibox) = j3*lg%num(3)/nproc_domain_orbital(3)+lg%is(3)
        mg%ie_all(3,ibox) = (j3+1)*lg%num(3)/nproc_domain_orbital(3)+lg%is(3)-1
      end do
      end do
    end do
    end do
    end do
  else if(process_allocation=='grid_sequential')then
    do i2=0,nproc_k-1
    do i1=0,nproc_ob-1
      do j3=0,nproc_domain_orbital(3)-1
      do j2=0,nproc_domain_orbital(2)-1
      do j1=0,nproc_domain_orbital(1)-1
        ibox = j1 + nproc_domain_orbital(1)*j2 + nproc_domain_orbital(1)*nproc_domain_orbital(2)*j3 &
        & + (nproc_ob*i2 + i1)*nproc_domain_orbital_mul
        mg%is_all(1,ibox) = j1*lg%num(1)/nproc_domain_orbital(1)+lg%is(1)
        mg%ie_all(1,ibox) = (j1+1)*lg%num(1)/nproc_domain_orbital(1)+lg%is(1)-1
        mg%is_all(2,ibox) = j2*lg%num(2)/nproc_domain_orbital(2)+lg%is(2)
        mg%ie_all(2,ibox) = (j2+1)*lg%num(2)/nproc_domain_orbital(2)+lg%is(2)-1
        mg%is_all(3,ibox) = j3*lg%num(3)/nproc_domain_orbital(3)+lg%is(3)
        mg%ie_all(3,ibox) = (j3+1)*lg%num(3)/nproc_domain_orbital(3)+lg%is(3)-1
      end do
      end do
      end do
    end do
    end do
  end if

  mg%is(:) = mg%is_all(:,myrank)
  mg%ie(:) = mg%ie_all(:,myrank)

  mg%num(:) = mg%ie(:)-mg%is(:)+1

  mg%is_overlap(1:3) = mg%is(1:3)-nd
  mg%ie_overlap(1:3) = mg%ie(1:3)+nd

  allocate(mg%idx(mg%is_overlap(1):mg%ie_overlap(1)) &
          ,mg%idy(mg%is_overlap(2):mg%ie_overlap(2)) &
          ,mg%idz(mg%is_overlap(3):mg%ie_overlap(3)))

  if(iperiodic==3 .and. nproc_domain_orbital_mul==1) then
    if(comm_is_root(nproc_id_global)) write(*,*) "r-space parallelization: off"
    mg%is_array(1:3) = mg%is(1:3)
    mg%ie_array(1:3) = mg%ie(1:3)
    do j=mg%is_overlap(1),mg%ie_overlap(1)
      mg%idx(j) = mod(j+mg%num(1)-1,mg%num(1))+1
    end do
    do j=mg%is_overlap(2),mg%ie_overlap(2)
      mg%idy(j) = mod(j+mg%num(2)-1,mg%num(2))+1
    end do
    do j=mg%is_overlap(3),mg%ie_overlap(3)
      mg%idz(j) = mod(j+mg%num(3)-1,mg%num(3))+1
    end do
  else
    mg%is_array(1:3) = mg%is(1:3)-nd
    mg%ie_array(1:3) = mg%ie(1:3)+nd
    do j=mg%is_overlap(1),mg%ie_overlap(1)
      mg%idx(j) = j
    end do
    do j=mg%is_overlap(2),mg%ie_overlap(2)
      mg%idy(j) = j
    end do
    do j=mg%is_overlap(3),mg%ie_overlap(3)
      mg%idz(j) = j
    end do
  end if

  if(calc_mode=='RT')then
#ifdef SALMON_STENCIL_PADDING
    mg%ie_array(2)=mg%ie(2)+nd+1
#endif
  end if

  if(mg%num(1)<nd .or.mg%num(2)<nd .or.mg%num(3)<nd)then
    stop "A system is small. Please use less number of processors."
  end if

! +-----------------------------+
! | ng: r-space grid for fields |
! +-----------------------------+

  ng%ndir = 3 ! high symmetry nonorthogonal lattice is not implemented
  ng%nd = Nd

  ngo = nproc_domain_general/nproc_domain_orbital
  ngo_mul = ngo(1)*ngo(2)*ngo(3)

  if(process_allocation=='orbital_sequential')then
    do ii=0,nproc_domain_orbital_mul-1
      do i4=0,nproc/nproc_domain_orbital_mul/ngo_mul-1
      do i3=0,ngo(3)-1
      do i2=0,ngo(2)-1
      do i1=0,ngo(1)-1
        ibox= i1+i2*ngo(1)+i3*ngo(1)*ngo(2)       &
                +i4*ngo_mul   &
                +ii*nproc/nproc_domain_orbital_mul
        ng%is_all(1,ibox) = i1*(mg%ie_all(1,ibox)-mg%is_all(1,ibox)+1)/ngo(1)+mg%is_all(1,ibox)
        ng%ie_all(1,ibox) = (i1+1)*(mg%ie_all(1,ibox)-mg%is_all(1,ibox)+1)/ngo(1)+mg%is_all(1,ibox)-1
        ng%is_all(2,ibox) = i2*(mg%ie_all(2,ibox)-mg%is_all(2,ibox)+1)/ngo(2)+mg%is_all(2,ibox)
        ng%ie_all(2,ibox) = (i2+1)*(mg%ie_all(2,ibox)-mg%is_all(2,ibox)+1)/ngo(2)+mg%is_all(2,ibox)-1
        ng%is_all(3,ibox) = i3*(mg%ie_all(3,ibox)-mg%is_all(3,ibox)+1)/ngo(3)+mg%is_all(3,ibox)
        ng%ie_all(3,ibox) = (i3+1)*(mg%ie_all(3,ibox)-mg%is_all(3,ibox)+1)/ngo(3)+mg%is_all(3,ibox)-1
      end do
      end do
      end do
      end do
    end do
  else if(process_allocation=='grid_sequential')then
    do i4=0,nproc/nproc_domain_orbital_mul/ngo_mul-1
    do i3=0,ngo(3)-1
    do i2=0,ngo(2)-1
    do i1=0,ngo(1)-1
      do ii=0,nproc_domain_orbital_mul-1
        ibox=ii+(i1+i2*ngo(1)+i3*ngo(1)*ngo(2))*nproc_domain_orbital_mul   &
              +i4*nproc_domain_orbital_mul*ngo_mul
        ng%is_all(1,ibox) = i1*(mg%ie_all(1,ii)-mg%is_all(1,ii)+1)/ngo(1)+mg%is_all(1,ii)
        ng%ie_all(1,ibox) = (i1+1)*(mg%ie_all(1,ii)-mg%is_all(1,ii)+1)/ngo(1)+mg%is_all(1,ii)-1
        ng%is_all(2,ibox) = i2*(mg%ie_all(2,ii)-mg%is_all(2,ii)+1)/ngo(2)+mg%is_all(2,ii)
        ng%ie_all(2,ibox) = (i2+1)*(mg%ie_all(2,ii)-mg%is_all(2,ii)+1)/ngo(2)+mg%is_all(2,ii)-1
        ng%is_all(3,ibox) = i3*(mg%ie_all(3,ii)-mg%is_all(3,ii)+1)/ngo(3)+mg%is_all(3,ii)
        ng%ie_all(3,ibox) = (i3+1)*(mg%ie_all(3,ii)-mg%is_all(3,ii)+1)/ngo(3)+mg%is_all(3,ii)-1
      end do
    end do
    end do
    end do
    end do
  end if

  ng%is(1:3) = ng%is_all(:,myrank)
  ng%ie(1:3) = ng%ie_all(:,myrank)

  ng%num(1:3) = ng%ie(1:3)-ng%is(1:3)+1

  ng%is_overlap(1:3) = ng%is(1:3)-nd
  ng%ie_overlap(1:3) = ng%ie(1:3)+nd
  ng%is_array(1:3) = ng%is(1:3)-nd
  ng%ie_array(1:3) = ng%ie(1:3)+nd

  allocate(ng%idx(ng%is_overlap(1):ng%ie_overlap(1)) &
          ,ng%idy(ng%is_overlap(2):ng%ie_overlap(2)) &
          ,ng%idz(ng%is_overlap(3):ng%ie_overlap(3)))
  do j=ng%is_overlap(1),ng%ie_overlap(1)
    ng%idx(j) = j
  end do
  do j=ng%is_overlap(2),ng%ie_overlap(2)
    ng%idy(j) = j
  end do
  do j=ng%is_overlap(3),ng%ie_overlap(3)
    ng%idz(j) = j
  end do

  if(iperiodic==0.and.(ng%num(1)<nd.or.ng%num(2)<nd.or.ng%num(3)<nd))then
    stop "A system is small. Please use less number of processors."
  end if

end subroutine init_grid_parallel

!===================================================================================================================================

subroutine setbn(bnmat)
  implicit none
  real(8) :: bnmat(4,4)

  bNmat(1,1)=1.d0/2.d0

  bNmat(1,2)=2.d0/3.d0
  bNmat(2,2)=-1.d0/12.d0

  bNmat(1,3)=3.d0/4.d0
  bNmat(2,3)=-3.d0/20.d0
  bNmat(3,3)=1.d0/60.d0

  bNmat(1,4)=4.d0/5.d0
  bNmat(2,4)=-1.d0/5.d0
  bNmat(3,4)=4.d0/105.d0
  bNmat(4,4)=-1.d0/280.d0

end subroutine setbN

subroutine setcn(cnmat)
  implicit none
  real(8) :: cnmat(0:12,12)

  cNmat(0,1)=-2.d0
  cNmat(1,1)=1.d0

  cNmat(0,2)=-5.d0/2.d0
  cNmat(1,2)=4.d0/3.d0
  cNmat(2,2)=-1.d0/12.d0

  cNmat(0,3)=-49.d0/18.d0
  cNmat(1,3)=3.d0/2.d0
  cNmat(2,3)=-3.d0/20.d0
  cNmat(3,3)=1.d0/90.d0

  cNmat(0,5)=-5269.d0/1800.d0
  cNmat(1,5)=5.d0/3.d0
  cNmat(2,5)=-5.d0/21.d0
  cNmat(3,5)=5.d0/126.d0
  cNmat(4,5)=-5.d0/1008.d0
  cNmat(5,5)=1.d0/3150.d0

  cNmat(0,4)=-205.d0/72.d0
  cNmat(1,4)=8.d0/5.d0
  cNmat(2,4)=-1.d0/5.d0
  cNmat(3,4)=8.d0/315.d0
  cNmat(4,4)=-1.d0/560.d0

  cNmat(0,6)=-5369.d0/1800.d0
  cNmat(1,6)=12.d0/7.d0
  cNmat(2,6)=-15.d0/56.d0
  cNmat(3,6)=10.d0/189.d0
  cNmat(4,6)=-1.d0/112.d0
  cNmat(5,6)=2.d0/1925.d0
  cNmat(6,6)=-1.d0/16632.d0

  cNmat(0,7)=-266681.d0/88200.d0
  cNmat(1,7)=7.d0/4.d0
  cNmat(2,7)=-7.d0/24.d0
  cNmat(3,7)=7.d0/108.d0
  cNmat(4,7)=-7.d0/528.d0
  cNmat(5,7)=7.d0/3300.d0
  cNmat(6,7)=-7.d0/30888.d0
  cNmat(7,7)=1.d0/84084.d0

  cNmat(0,8)=-1077749.d0/352800.d0
  cNmat(1,8)=16.d0/9.d0
  cNmat(2,8)=-14.d0/45.d0
  cNmat(3,8)=112.d0/1485.d0
  cNmat(4,8)=-7.d0/396.d0
  cNmat(5,8)=112.d0/32175.d0
  cNmat(6,8)=-2.d0/3861.d0
  cNmat(7,8)=16.d0/315315.d0
  cNmat(8,8)=-1.d0/411840.d0

  cNmat(0,9)=-9778141.d0/3175200.d0
  cNmat(1,9)=9.d0/5.d0
  cNmat(2,9)=-18.d0/55.d0
  cNmat(3,9)=14.d0/165.d0
  cNmat(4,9)=-63.d0/2860.d0
  cNmat(5,9)=18.d0/3575.d0
  cNmat(6,9)=-2.d0/2145.d0
  cNmat(7,9)=9.d0/70070.d0
  cNmat(8,9)=-9.d0/777920.d0
  cNmat(9,9)=1.d0/1969110.d0

  cNmat(0,10)=-1968329.d0/635040.d0
  cNmat(1,10)=20.d0/11.d0
  cNmat(2,10)=-15.d0/44.d0
  cNmat(3,10)=40.d0/429.d0
  cNmat(4,10)=-15.d0/572.d0
  cNmat(5,10)=24.d0/3575.d0
  cNmat(6,10)=-5.d0/3432.d0
  cNmat(7,10)=30.d0/119119.d0
  cNmat(8,10)=-5.d0/155584.d0
  cNmat(9,10)=10.d0/3741309.d0
  cNmat(10,10)=-1.d0/9237800.d0

  cNmat(0,11)=-239437889.d0/76839840.d0
  cNmat(1,11)=11.d0/6.d0
  cNmat(2,11)=-55.d0/156.d0
  cNmat(3,11)=55.d0/546.d0
  cNmat(4,11)=-11.d0/364.d0
  cNmat(5,11)=11.d0/1300.d0
  cNmat(6,11)=-11.d0/5304.d0
  cNmat(7,11)=55.d0/129948.d0
  cNmat(8,11)=-55.d0/806208.d0
  cNmat(9,11)=11.d0/1360476.d0
  cNmat(10,11)=-11.d0/17635800.d0
  cNmat(11,11)=1.d0/42678636.d0

  cNmat(0,12)=-240505109.d0/76839840.d0
  cNmat(1,12)=24.d0/13.d0
  cNmat(2,12)=-33.d0/91.d0
  cNmat(3,12)=88.d0/819.d0
  cNmat(4,12)=-99.d0/2912.d0
  cNmat(5,12)=396.d0/38675.d0
  cNmat(6,12)=-11.d0/3978.d0
  cNmat(7,12)=132.d0/205751.d0
  cNmat(8,12)=-33.d0/268736.d0
  cNmat(9,12)=44.d0/2380833.d0
  cNmat(10,12)=-3.d0/1469650.d0
  cNmat(11,12)=12.d0/81800719.d0
  cNmat(12,12)=-1.d0/194699232.d0

end subroutine setcN

end module salmon_initialization
