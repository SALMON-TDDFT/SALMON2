program probe
 use structures,only:s_kpoint_config,s_dft_system,s_unfold
 use parallelization,only:setup_parallel,end_parallel
 use lattice,only:init_kvector_config,input_kpoint_config
 use salmon_global,only:num_kgrid,file_kw,dk_shift,dm_unfold_option,num_lkgrid,num_skgrid,yn_gamma_centered
 use sym_sub,only:use_symmetry,DISPLAY,SymMatB
 implicit none
 type(s_kpoint_config)::a,b
 type(s_dft_system)::reference,actual,other
 type(s_unfold)::ur,ua,uo
 integer::mode,i,u
 call setup_parallel
 reference%primitive_b=0d0
 do i=1,3
  reference%primitive_b(i,i)=2d0*i
 enddo
 actual%primitive_b=reference%primitive_b;other%primitive_b=reference%primitive_b
 do mode=1,7
  num_kgrid=[3,2,1];dk_shift=[.125d0,-.25d0,.375d0]
  num_lkgrid=[1,1,1];num_skgrid=[6,4,2]
  file_kw='none';dm_unfold_option='no';yn_gamma_centered='n'
  use_symmetry=.false.;DISPLAY=.false.
  if(mode==2)yn_gamma_centered='y'
  if(mode==3)dm_unfold_option='primitive'
  if(mode==4)dm_unfold_option='super'
  if(mode==5)then
   file_kw='kpoints.dat'
   open(newunit=u,file=trim(file_kw),status='new')
   write(u,*)2
   write(u,*)1,.125d0,.25d0,.375d0,2d0
   write(u,*)2,-.125d0,-.25d0,-.375d0,3d0
   close(u)
  endif
  if(mode==6)then
   num_kgrid=[3,3,3];dk_shift=0d0;use_symmetry=.true.
   if(allocated(SymMatB))deallocate(SymMatB)
   allocate(SymMatB(3,3,2));SymMatB=0d0
   do i=1,3
    SymMatB(i,i,1)=1d0;SymMatB(i,i,2)=-1d0
   enddo
  endif
  if(mode==7)then
   num_kgrid=[3,3,3];yn_gamma_centered='y'
  endif
  call input_kpoint_config(a)
  call reference_init_kvector(reference,ur)
  b=s_kpoint_config();b%num_kgrid=[2,1,1]
  ! Poison every selection input after capture, including symmetry.
  num_kgrid=[7,7,7];dk_shift=99d0;file_kw='must-not-open'
  num_lkgrid=19;num_skgrid=23;dm_unfold_option='invalid';yn_gamma_centered='y'
  use_symmetry=.not.a%use_symmetry
  if(allocated(SymMatB))SymMatB=99d0
  call init_kvector_config(other,uo,b)
  call init_kvector_config(actual,ua,a)
  call compare()
  call init_kvector_config(other,uo,b)
  if(other%nk/=2.or.any(other%kmesh/=b%num_kgrid))stop 8
  if(any(num_kgrid/=[7,7,7]).or.file_kw/='must-not-open')stop 9
  if(a%file_kw=='none')then
   if(any(actual%kmesh/=a%num_kgrid))stop 10
  else
   if(any(actual%kmesh/=-1))stop 11
  endif
  call clear_unfold(ur);call clear_unfold(ua)
 enddo
 print *,'Explicit k configuration PASS'
 call end_parallel
contains
 subroutine compare()
  if(reference%nk/=actual%nk)stop 1
  if(any(reference%vec_k/=actual%vec_k).or.any(reference%wtk/=actual%wtk))stop 2
  if(allocated(ur%vec_hk))then
   if(.not.allocated(ua%vec_hk))stop 3
   if(any(ur%vec_hk/=ua%vec_hk).or.any(ur%wtk_pr/=ua%wtk_pr))stop 4
   if(any(ur%isk_tbl/=ua%isk_tbl).or.any(ur%num_hkgrid/=ua%num_hkgrid))stop 5
  endif
 end subroutine
 subroutine clear_unfold(value)
  type(s_unfold),intent(inout)::value
  if(allocated(value%vec_hk))deallocate(value%vec_hk)
  if(allocated(value%isk_tbl))deallocate(value%isk_tbl)
  if(allocated(value%wtk_pr))deallocate(value%wtk_pr)
 end subroutine
! Frozen reference implementation from runtime/sources/hse-fast-reference/src/
! common/lattice.f90 and symmetry/symmetry_kvector.f90, SALMON developers.
! Kept separately from the explicit overlay; only routine names/import placement
! change. Original arithmetic, file IO and global selection are retained.
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
SUBROUTINE reference_init_kvector(system,unfold)
  use structures
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
  call reference_init_sym_kvector( system%vec_k, system%wtk, system%nk, B ) 

  return
end SUBROUTINE reference_init_kvector

  subroutine reference_init_sym_kvector( kvec_io, weight_io, nkvec, blatvec )
    use sym_sub, only: SymMatB,use_symmetry,DISPLAY
    use salmon_global, only: dm_unfold_option
    use parallelization, only: end_parallel
    implicit none
    real(8),intent(inout) :: kvec_io(:,:), weight_io(:)
    integer,intent(inout) :: nkvec
    real(8),optional,intent(in) :: blatvec(3,3)
    real(8),allocatable :: kvec(:,:), weight(:)
    real(8) :: blatinv(3,3), transformed(3), delta(3)
    integer :: isym, ik, n_unique, i, representative

    if ( .not.use_symmetry ) return

    if( dm_unfold_option == 'super' .or. dm_unfold_option == 'primitive' ) then
      if ( DISPLAY ) write(*,"(A)") 'symmetry option not supported for dm_unfold calculation'
      call end_parallel
      stop
    end if

    if ( DISPLAY ) write(*,'(a60)') repeat("-",36)//" init_sym_kvector(start)"

    allocate( kvec(3,nkvec) ); kvec=0.0d0
    allocate( weight(nkvec) )

    if ( present(blatvec) ) then
       call get_inverse_lattice( blatvec, blatinv )
       kvec=matmul( blatinv, kvec_io )
    else
       kvec=kvec_io
    end if
    weight = weight_io

    kvec_io=0.0d0
    weight_io=0.0d0
    n_unique=0

    do ik=1,nkvec
       representative=0
       do isym=1,size(SymMatB,3)
          transformed=matmul(SymMatB(:,1:3,isym),kvec(:,ik))
          do i=1,n_unique
             delta=transformed-kvec(:,i)
             ! Reciprocal-lattice translations leave a k point equivalent.
             delta=delta-anint(delta)
             if (sum(delta**2)<1.d-5) then
                representative=i
                exit
             end if
          end do
          if (representative/=0) exit
       end do
       if (representative==0) then
          n_unique=n_unique+1
          representative=n_unique
          kvec(:,representative)=kvec(:,ik)
          call shift_vec_into_bz(kvec(:,representative))
       end if
       weight_io(representative)=weight_io(representative)+weight(ik)
    end do
    nkvec=n_unique
    if (present(blatvec)) then
       kvec_io(:,:nkvec)=matmul(blatvec,kvec(:,:nkvec))
    else
       kvec_io(:,:nkvec)=kvec(:,:nkvec)
    end if

    if ( DISPLAY ) then
       do i=1,nkvec
          write(*,'(i4,f10.5,2x,3f10.5)') i, weight_io(i), kvec_io(:,i)
       end do
       write(*,*) "sum(weight)=",sum(weight_io)
    end if

    deallocate( kvec )
    deallocate( weight )

    if ( DISPLAY ) write(*,'(a60)') repeat("-",38)//" init_sym_kvector(end)"

  end subroutine reference_init_sym_kvector

  subroutine get_inverse_lattice( a, ainv )
    implicit none
    real(8),intent(in)  :: a(3,3)
    real(8),intent(out) :: ainv(3,3)
    real(8) :: b(3,3), v
    b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
    b(2,1) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
    b(3,1) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
    b(1,2) = a(2,3)*a(3,1) - a(3,3)*a(2,1)
    b(2,2) = a(3,3)*a(1,1) - a(1,3)*a(3,1)
    b(3,2) = a(1,3)*a(2,1) - a(2,3)*a(1,1)
    b(1,3) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
    b(2,3) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
    b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
    v=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1) &
     +a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1) &
     -a(1,2)*a(2,1)*a(3,3)-a(1,1)*a(2,3)*a(3,2)
    ainv(:,:) = transpose( b(:,:) )/v
  end subroutine get_inverse_lattice

  subroutine shift_vec_into_bz( vec )
    implicit none
    real(8),intent(inout) :: vec(:)
    integer :: i
    do i=1,size(vec)
       do
          if ( vec(i) > 0.5d0 ) then
             vec(i) = vec(i) - 1.0d0
          else if ( vec(i) <= -0.5d0 ) then
             vec(i) = vec(i) + 1.0d0
          else
             exit
          end if
       end do
    end do
  end subroutine shift_vec_into_bz


end program
