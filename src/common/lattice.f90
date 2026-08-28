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
module lattice
  implicit none

contains

!===================================================================================================================================

SUBROUTINE init_lattice(system,stencil)
  use math_constants,only : pi
  use structures
  implicit none
  type(s_dft_system) :: system
  type(s_stencil)    :: stencil
  !
  real(8),dimension(3,3) :: A,B,F,wrk
  real(8) :: a1(3),a2(3),a3(3),detA,normA(3),f_uu,f_vv,f_ww,f_uv,f_uw,f_vw

! al = [ a1, a2, a3 ]
  A = system%primitive_a ! primitive lattice vectors
  a1 = A(1:3,1)
  a2 = A(1:3,2)
  a3 = A(1:3,3)
  call calc_inverse(A,wrk,detA)
  ! detA is the signed scalar triple product a1.(a2 x a3); it is negative when the
  ! primitive vectors are given in a left-handed order (e.g. swapped orthorhombic
  ! axes). The cell volume and grid-volume element must be positive, otherwise
  ! sqrt(rbox2*Hvol) in the wavefunction normalization yields NaN. Use |detA| for
  ! the volume; the (signed) inverse wrk is left untouched so the reciprocal
  ! lattice primitive_b = 2*pi*A^{-1} stays correct regardless of handedness.
  system%det_a = abs(detA)
  system%Hvol = abs(detA)/dble(system%ngrid)
  system%primitive_b = 2d0*pi* transpose(wrk) ! reciprocal primitive lattice vectors
  ! [ b1 b2 b3 ]^{T} = 2*pi* [ a1 a2 a3 ]^{-1}

  normA(1) = sqrt(sum(a1**2))
  normA(2) = sqrt(sum(a2**2))
  normA(3) = sqrt(sum(a3**2))

! cf. A. Natan et al., PRB 78, 075109 (2008).
! A = [ u, v, w ], B = A^{-1}
  A(1:3,1) = a1(1:3) / normA(1) ! u
  A(1:3,2) = a2(1:3) / normA(2) ! v
  A(1:3,3) = a3(1:3) / normA(3) ! w
  call calc_inverse(A,B,detA)

  wrk = transpose(B)
  F = matmul(B,wrk)
  f_uu = F(1,1)
  f_vv = F(2,2)
  f_ww = F(3,3)
  f_uv = F(1,2) + F(2,1)
  f_uw = F(1,3) + F(3,1)
  f_vw = F(2,3) + F(3,2)

  system%rmatrix_A = A
  system%rmatrix_B = B

  stencil%coef_F(1) = f_uu
  stencil%coef_F(2) = f_vv
  stencil%coef_F(3) = f_ww
  stencil%coef_F(4) = f_vw ! yz
  stencil%coef_F(5) = f_uw ! zx
  stencil%coef_F(6) = f_uv ! xy

  return
end SUBROUTINE init_lattice

SUBROUTINE calc_inverse(a,b,detA) ! b = a^{-1}
  implicit none
  real(8),intent(in) :: a(3,3)
  real(8)            :: b(3,3),detA

  detA=a(1,1)*a(2,2)*a(3,3)+a(2,1)*a(3,2)*a(1,3)+a(3,1)*a(1,2)*a(2,3) &
    -a(1,3)*a(2,2)*a(3,1)-a(2,3)*a(3,2)*a(1,1)-a(3,3)*a(1,2)*a(2,1)

  b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
  b(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
  b(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)

  b(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3)
  b(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1)
  b(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2)

  b(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
  b(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
  b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)

  b=b/detA
  return
end SUBROUTINE calc_inverse

SUBROUTINE init_kvector(system,unfold)
  use structures
  use sym_kvector, only: init_sym_kvector
  use parallelization, only: nproc_id_global, nproc_group_global
  use communication, only: comm_bcast, comm_sync_all, comm_is_root
  use salmon_global, only: num_kgrid, file_kw, dk_shift, dm_unfold_option, num_lkgrid, num_skgrid, yn_gamma_centered
  implicit none
  type(s_dft_system) :: system
  type(s_unfold) :: unfold
  integer :: num_sk(3), num_lk(3), num_hk(3)
  integer :: eo_kgrid(3), eo_lk(3), eo_hk(3)
  integer :: i,iu,ik,nk,nsk,nlk,nhk,nk0,nsk0,nlk0
  integer :: ik1,ik2,ik3,ilk,ilk1,ilk2,ilk3,ihk,ihk1,ihk2,ihk3,isk
  real(8) :: B(3,3)
  real(8) :: lkp(3),hkp(3)
  real(8),allocatable :: k(:,:), wtk(:)
  real(8),allocatable :: hks(:,:)

  B = system%primitive_b

  ! Read from file_kw
  if(file_kw /= 'none') then
     if (comm_is_root(nproc_id_global)) then
        iu = 410
        open(iu, file=file_kw, status="old")
        read(iu, *) nk    !, nkxyz_dummy
        allocate( k(3,nk), wtk(nk) )
        do ik=1, nk
           read(iu, *) i, k(1:3,ik), wtk(ik)
        enddo
        close(iu)
        wtk(:) = wtk(:) * (1.0d0 / sum(wtk))
     endif

     call comm_bcast(nk,  nproc_group_global)
     if(.not.allocated(k)) allocate( k(3,nk), wtk(nk) )
     call comm_bcast(k,  nproc_group_global)
     call comm_bcast(wtk,nproc_group_global)
     if( abs(sum(wtk(:))-1d0).ge.1d-5 ) stop "error in wtk (wight of k-point)"

     num_kgrid(:) = -1 

  else

    select case (dm_unfold_option) ! 'no' is the ordinary case

    case('no')
!      do ik=1,nk
!        ix=mod(ik-1,num_kgrid(1))+1
!        iy=mod((ik-1)/num_kgrid(1),num_kgrid(2))+1
!        iz=mod((ik-1)/(num_kgrid(1)*num_kgrid(2)),num_kgrid(3))+1
!        k(1,ik) = (dble(ix)-shift_k(1)+dk_shift(1))/dble(num_kgrid(1))-0.5d0
!        k(2,ik) = (dble(iy)-shift_k(2)+dk_shift(2))/dble(num_kgrid(2))-0.5d0
!        k(3,ik) = (dble(iz)-shift_k(3)+dk_shift(3))/dble(num_kgrid(3))-0.5d0
!      end do

      if( yn_gamma_centered == 'n' ) then

        nk = num_kgrid(1)*num_kgrid(2)*num_kgrid(3)
        allocate( k(3,nk), wtk(nk) )
        wtk(:) = 1d0/dble(nk)

        ik = 0
        do ik3 = 1,num_kgrid(3)
        do ik2 = 1,num_kgrid(2)
        do ik1 = 1,num_kgrid(1)
          ik = ik + 1
          k(1,ik) = (dble(ik1)-0.5d0+dk_shift(1))/dble(num_kgrid(1))-0.5d0
          k(2,ik) = (dble(ik2)-0.5d0+dk_shift(2))/dble(num_kgrid(2))-0.5d0
          k(3,ik) = (dble(ik3)-0.5d0+dk_shift(3))/dble(num_kgrid(3))-0.5d0
        enddo
        enddo
        enddo

      else ! yn_gamma_centered = 'y'

        nk0 = num_kgrid(1)*num_kgrid(2)*num_kgrid(3)
        eo_kgrid(:) = mod(num_kgrid(:),2)
        nk = (num_kgrid(1)+1-eo_kgrid(1))*(num_kgrid(2)+1-eo_kgrid(2))*(num_kgrid(3)+1-eo_kgrid(3))
!         for even num_kgrid, both ends are included for k-grid

        allocate( k(3,nk), wtk(nk) )
        wtk(:)  = 1d0/dble(nk0)

        ik = 0
        do ik3 = eo_kgrid(3),num_kgrid(3)
        do ik2 = eo_kgrid(2),num_kgrid(2)
        do ik1 = eo_kgrid(1),num_kgrid(1)
          ik = ik + 1
          k(1,ik) = (dble(ik1)-0.5d0*eo_kgrid(1))/dble(num_kgrid(1))-0.5d0
          k(2,ik) = (dble(ik2)-0.5d0*eo_kgrid(2))/dble(num_kgrid(2))-0.5d0
          k(3,ik) = (dble(ik3)-0.5d0*eo_kgrid(3))/dble(num_kgrid(3))-0.5d0
          if(eo_kgrid(1) == 0 .and. (ik1 == 0 .or. ik1 == num_kgrid(1))) wtk(ik) = wtk(ik)/2d0
          if(eo_kgrid(2) == 0 .and. (ik2 == 0 .or. ik2 == num_kgrid(2))) wtk(ik) = wtk(ik)/2d0
          if(eo_kgrid(3) == 0 .and. (ik3 == 0 .or. ik3 == num_kgrid(3))) wtk(ik) = wtk(ik)/2d0
        enddo
        enddo
        enddo
        if (comm_is_root(nproc_id_global)) then
          write(*,"(A,2x,i8,2x,A,2x,f18.8)") 'k-point gamma-centered, nk=', nk, 'sum(wtk)=', sum(wtk)
        end if

      end if

    case('primitive') ! preparation for dm_unfold calculation

      num_sk(:) = num_kgrid(:)
      num_lk(:) = num_lkgrid(:)
      if( sum(mod(num_sk(:),num_lk(:))) /= 0) stop 'mod(num_sk,num_lk) /= 0 in dm_unfold, primitive'
      num_hk(:) = num_sk(:)/num_lk(:) 
      nsk0 = num_sk(1)*num_sk(2)*num_sk(3)
      eo_lk(:) = mod(num_lk(:),2)
      eo_hk(:) = mod(num_hk(:),2)
      nlk = (num_lk(1)+1-eo_lk(1))*(num_lk(2)+1-eo_lk(2))*(num_lk(3)+1-eo_lk(3))
      nhk = (num_hk(1)+1-eo_hk(1))*(num_hk(2)+1-eo_hk(2))*(num_hk(3)+1-eo_hk(3))
      nsk = nlk * nhk

      nk = nsk
      allocate( k(3,nsk), wtk(nsk) )
      wtk(:) = 1d0/dble(nsk0)
      isk = 0
      do ilk3 = eo_lk(3),num_lk(3)
      do ilk2 = eo_lk(2),num_lk(2)
      do ilk1 = eo_lk(1),num_lk(1)
      do ihk3 = eo_hk(3),num_hk(3)
      do ihk2 = eo_hk(2),num_hk(2)
      do ihk1 = eo_hk(1),num_hk(1)
        isk = isk + 1
        lkp(1) = (dble(ilk1)-0.5d0*eo_lk(1))/dble(num_lk(1))-0.5d0
        lkp(2) = (dble(ilk2)-0.5d0*eo_lk(2))/dble(num_lk(2))-0.5d0
        lkp(3) = (dble(ilk3)-0.5d0*eo_lk(3))/dble(num_lk(3))-0.5d0
        hkp(1) = (dble(ihk1)-0.5d0*eo_hk(1))/dble(num_hk(1))-0.5d0
        hkp(2) = (dble(ihk2)-0.5d0*eo_hk(2))/dble(num_hk(2))-0.5d0
        hkp(3) = (dble(ihk3)-0.5d0*eo_hk(3))/dble(num_hk(3))-0.5d0
        k(1:3,isk) = lkp(1:3)/dble(num_hk(1:3)) + hkp(1:3)
        do i = 1,3
          if( k(i,isk) < -0.50001d0 ) k(i,isk) = k(i,isk) + 1d0
          if( k(i,isk) > +0.50001d0 ) k(i,isk) = k(i,isk) - 1d0
        end do
        if(eo_lk(1) == 0 .and. (ilk1 == 0 .or. ilk1 == num_lk(1))) wtk(isk) = wtk(isk)/2d0
        if(eo_lk(2) == 0 .and. (ilk2 == 0 .or. ilk2 == num_lk(2))) wtk(isk) = wtk(isk)/2d0
        if(eo_lk(3) == 0 .and. (ilk3 == 0 .or. ilk3 == num_lk(3))) wtk(isk) = wtk(isk)/2d0
        if(eo_hk(1) == 0 .and. (ihk1 == 0 .or. ihk1 == num_hk(1))) wtk(isk) = wtk(isk)/2d0
        if(eo_hk(2) == 0 .and. (ihk2 == 0 .or. ihk2 == num_hk(2))) wtk(isk) = wtk(isk)/2d0
        if(eo_hk(3) == 0 .and. (ihk3 == 0 .or. ihk3 == num_hk(3))) wtk(isk) = wtk(isk)/2d0
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      if (comm_is_root(nproc_id_global)) then
        write(*,"(A,2x,i8,2x,A,2x,f18.8)") 'dm_unfold_option=primitive, nk=', nk, 'sum(wtk)=', sum(wtk)
      end if

    case('super') ! dm_unfold calculation

      num_lk(:) = num_kgrid(:)
      nlk0 = num_lk(1)*num_lk(2)*num_lk(3)
      eo_lk(:) = mod(num_lk(:),2)
      nlk = (num_lk(1)+1-eo_lk(1))*(num_lk(2)+1-eo_lk(2))*(num_lk(3)+1-eo_lk(3))

      nk = nlk

      allocate( k(3,nlk), wtk(nlk) )
      wtk(:)  = 1d0/dble(nlk0)

      ilk = 0
      do ilk3 = eo_lk(3),num_lk(3)
      do ilk2 = eo_lk(2),num_lk(2)
      do ilk1 = eo_lk(1),num_lk(1)
        ilk = ilk + 1
        k(1,ilk) = (dble(ilk1)-0.5d0*eo_lk(1))/dble(num_lk(1))-0.5d0
        k(2,ilk) = (dble(ilk2)-0.5d0*eo_lk(2))/dble(num_lk(2))-0.5d0
        k(3,ilk) = (dble(ilk3)-0.5d0*eo_lk(3))/dble(num_lk(3))-0.5d0
        if(eo_lk(1) == 0 .and. (ilk1 == 0 .or. ilk1 == num_lk(1))) wtk(ilk) = wtk(ilk)/2d0
        if(eo_lk(2) == 0 .and. (ilk2 == 0 .or. ilk2 == num_lk(2))) wtk(ilk) = wtk(ilk)/2d0
        if(eo_lk(3) == 0 .and. (ilk3 == 0 .or. ilk3 == num_lk(3))) wtk(ilk) = wtk(ilk)/2d0
      enddo
      enddo
      enddo
      if (comm_is_root(nproc_id_global)) then
        write(*,"(A,2x,i8,2x,A,2x,f18.8)") 'dm_unfold_option=super, nk=', nk, 'sum(wtk)=', sum(wtk)
      end if

      num_sk(:) = num_skgrid(:)
      if( sum(mod(num_sk(:),num_lk(:))) /= 0) stop 'mod(num_sk,num_lk) /= 0 in dm_unfold, super'
      num_hk(:) = num_sk(:)/num_lk(:)
      eo_hk(:) = mod(num_hk(:),2)
      nhk = (num_hk(1)+1-eo_hk(1))*(num_hk(2)+1-eo_hk(2))*(num_hk(3)+1-eo_hk(3))
      nsk = nlk * nhk

      unfold%nhk = nhk
      unfold%num_hkgrid(1:3) = num_hk(:)
      unfold%nsk = nsk
      allocate( hks(3,nhk) )
      allocate( unfold%vec_hk(3, nhk), unfold%isk_tbl(nlk, nhk), unfold%wtk_pr(nsk) )

      ihk = 0
      do ihk3 = eo_hk(3),num_hk(3)
      do ihk2 = eo_hk(2),num_hk(2)
      do ihk1 = eo_hk(1),num_hk(1)
        ihk = ihk + 1
        hks(1,ihk) = ((dble(ihk1)-0.5d0*eo_hk(1))/dble(num_hk(1))-0.5d0) * dble(num_hk(1))
        hks(2,ihk) = ((dble(ihk2)-0.5d0*eo_hk(2))/dble(num_hk(2))-0.5d0) * dble(num_hk(2))
        hks(3,ihk) = ((dble(ihk3)-0.5d0*eo_hk(3))/dble(num_hk(3))-0.5d0) * dble(num_hk(3))
! since B in supercell is folded, num_hk(:) is multiplied.
      enddo
      enddo
      enddo
      do ihk = 1,nhk
        unfold%vec_hk(1,ihk) = hks(1,ihk)*B(1,1) + hks(2,ihk)*B(1,2) + hks(3,ihk)*B(1,3)
        unfold%vec_hk(2,ihk) = hks(1,ihk)*B(2,1) + hks(2,ihk)*B(2,2) + hks(3,ihk)*B(2,3)
        unfold%vec_hk(3,ihk) = hks(1,ihk)*B(3,1) + hks(2,ihk)*B(3,2) + hks(3,ihk)*B(3,3)
      end do

      nsk0 = num_sk(1)*num_sk(2)*num_sk(3)
      unfold%wtk_pr(:) = 1d0/dble(nsk0)
      isk = 0
      ilk=0
      do ilk3 = eo_lk(3),num_lk(3)
      do ilk2 = eo_lk(2),num_lk(2)
      do ilk1 = eo_lk(1),num_lk(1)
        ilk = ilk + 1
        ihk = 0
      do ihk3 = eo_hk(3),num_hk(3)
      do ihk2 = eo_hk(2),num_hk(2)
      do ihk1 = eo_hk(1),num_hk(1)
        isk = isk + 1
        ihk = ihk + 1
        unfold%isk_tbl(ilk,ihk) = isk
        if(eo_lk(1) == 0 .and. (ilk1 == 0 .or. ilk1 == num_lk(1))) unfold%wtk_pr(isk) = unfold%wtk_pr(isk)/2d0
        if(eo_lk(2) == 0 .and. (ilk2 == 0 .or. ilk2 == num_lk(2))) unfold%wtk_pr(isk) = unfold%wtk_pr(isk)/2d0
        if(eo_lk(3) == 0 .and. (ilk3 == 0 .or. ilk3 == num_lk(3))) unfold%wtk_pr(isk) = unfold%wtk_pr(isk)/2d0
        if(eo_hk(1) == 0 .and. (ihk1 == 0 .or. ihk1 == num_hk(1))) unfold%wtk_pr(isk) = unfold%wtk_pr(isk)/2d0
        if(eo_hk(2) == 0 .and. (ihk2 == 0 .or. ihk2 == num_hk(2))) unfold%wtk_pr(isk) = unfold%wtk_pr(isk)/2d0
        if(eo_hk(3) == 0 .and. (ihk3 == 0 .or. ihk3 == num_hk(3))) unfold%wtk_pr(isk) = unfold%wtk_pr(isk)/2d0
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

    end select
  endif

  system%nk = nk
  if ( allocated(system%vec_k) ) deallocate(system%vec_k)
  if ( allocated(system%wtk)   ) deallocate(system%wtk)
  allocate(system%vec_k(3,nk),system%wtk(nk))
  system%wtk  = wtk

  do ik=1,nk
    system%vec_k(1,ik) = k(1,ik)*B(1,1) + k(2,ik)*B(1,2) + k(3,ik)*B(1,3)
    system%vec_k(2,ik) = k(1,ik)*B(2,1) + k(2,ik)*B(2,2) + k(3,ik)*B(2,3)
    system%vec_k(3,ik) = k(1,ik)*B(3,1) + k(2,ik)*B(3,2) + k(3,ik)*B(3,3)
  end do
  call init_sym_kvector( system%vec_k, system%wtk, system%nk, B ) 

  return
end SUBROUTINE init_kvector

end module lattice
