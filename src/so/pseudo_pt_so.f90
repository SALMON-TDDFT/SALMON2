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
module pseudo_pt_so_sub

  implicit none

! WARNING: We must not call these except for hpsi routine.

  private
  public :: pseudo_so
  public :: calc_uVpsi_so
  public :: calc_force_uVpsi_so

contains

  subroutine pseudo_so(tpsi,htpsi,info,nspin,ppg,mg)
    use structures, only: s_parallel_info, s_pp_grid, s_rgrid, s_orbital
    use communication, only: comm_summation
    !use parallelization, only: nproc_id_global
    use timer
    implicit none
    integer,intent(in) :: nspin
    type(s_parallel_info),intent(in) :: info
    type(s_pp_grid),intent(in) :: ppg
    type(s_rgrid),intent(in) :: mg
    type(s_orbital),intent(in) :: tpsi
    type(s_orbital) :: htpsi
    !
    integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e
    integer :: ilma,ia,j,ix,iy,iz,Nlma
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8) :: uVpsi(2),wrk
    complex(8),allocatable :: uVpsibox1(:,:,:,:,:)
    complex(8),allocatable :: uVpsibox2(:,:,:,:,:)
    !real(8) :: tmp(2),tmp1(2),tmp2(2),tmp3(2)
    !real(8) :: wf_check(2)
    !integer :: iu

    !write(*,*) "------------------ pseudo_so"

    call timer_begin(LOG_UHPSI_PSEUDO)

    im_s = info%im_s
    im_e = info%im_e
    ik_s = info%ik_s
    ik_e = info%ik_e
    io_s = info%io_s
    io_e = info%io_e

    Nlma = size(ppg%ia_tbl_so)
    !write(*,*) "Nlma=",Nlma

    !iu=200+nproc_id_global
    !rewind iu
    !do ilma=1,Nlma
    !   write(iu,*) ilma,ppg%rinv_uvu_so(ilma)
    !end do

    !rewind 100

    !io=1; ik=1; im=1
    !iu=400+nproc_id_global
    !rewind iu
    !write(iu,*) minval(mg%is_all(1,:)),minval(mg%is_all(2,:)),minval(mg%is_all(3,:)), &
    !            maxval(mg%ie_all(1,:)),maxval(mg%ie_all(2,:)),maxval(mg%ie_all(3,:))
    !write(iu,*) mg%is,mg%ie
    !do ispin=1,2
    !do iz=mg%is(3),mg%ie(3)
    !do iy=mg%is(2),mg%ie(2)
    !do ix=mg%is(1),mg%ie(1)
    !   write(iu,'(4i6,2g23.15)') ix,iy,iz,ispin, &
    !        real(tpsi%zwf(ix,iy,iz,ispin,io,ik,im)), aimag(tpsi%zwf(ix,iy,iz,ispin,io,ik,im))
    !end do
    !end do
    !end do
    !end do

    !ia=maxloc( ppg%mps,1 )
    !iu=500+nproc_id_global
    !rewind iu
    !write(iu,*) Nlma, ppg%mps(ia)
    !do ispin=1,2
    !   do ilma=1,Nlma
    !      ia=ppg%ia_tbl_so(ilma)
    !      write(iu,*) ia, ppg%mps(ia)
    !      do j=1, ppg%mps(ia)
    !         ix = ppg%jxyz(1,j,ia)
    !         iy = ppg%jxyz(2,j,ia)
    !         iz = ppg%jxyz(3,j,ia)
    !         write(iu,'(6i6,2g23.15)') ix,iy,iz,j,ilma,ispin,ppg%zekr_uV_so(j,ilma,ik,ispin,im)
    !      end do
    !   end do
    !end do

!real(ppg%zekr_uV_so(j,ilma,ik,ispin,1)), aimag(ppg%zekr_uV_so(j,ilma,ik,ispin,1)) &
!              end if
    if ( info%if_divide_rspace ) then

      allocate(uVpsibox1(Nspin,Nlma,io_s:io_e,ik_s:ik_e,im_s:im_e)); uVpsibox1=zero
      allocate(uVpsibox2(Nspin,Nlma,io_s:io_e,ik_s:ik_e,im_s:im_e)); uVpsibox2=zero

#ifdef USE_OPENACC
!$acc kernels
!#acc loop collapse(4) default(shared) private(im,ik,io,ilma,ia,ispin,wrk,j,iz,iy,ix)
#else
!$omp parallel do collapse(4) default(shared) private(im,ik,io,ilma,ia,ispin,wrk,j,iz,iy,ix)
#endif
      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e

        do ilma=1,Nlma

          ia = ppg%ia_tbl_so(ilma)

          !if ( ia==1.and.io==1 ) then
          !  do ispin=1,2
          !  tmp(ispin)=sum(abs(ppg%zekr_uV_so(:,ilma,ik,ispin,1)))
          !  end do
          !  call comm_summation( tmp,tmp1,2,info%icomm_r )
          !  if ( nproc_id_global==0 ) then
          !    !write(*,'("ilma,ia,mps,rank",4i6,2f20.15)') ilma,ia,ppg%mps(ia),nproc_id_global,tmp1
          !    write(100,'(3i6,2f20.10)') ilma,ia,nproc_id_global,tmp1
          !  end if
          !end if


          !wf_check(1:2)=zero
          do ispin=1,2
            wrk=zero
            do j=1,ppg%mps(ia)
              ix = ppg%jxyz(1,j,ia)
              iy = ppg%jxyz(2,j,ia)
              iz = ppg%jxyz(3,j,ia)
              wrk = wrk + conjg( ppg%zekr_uV_so(j,ilma,ik,ispin,1) ) &
                        * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
              !wf_check(ispin) = wf_check(ispin) + abs(tpsi%zwf(ix,iy,iz,ispin,io,ik,im))
            end do
            wrk = wrk * ppg%rinv_uvu_so(ilma)
            uVpsibox1(ispin,ilma,io,ik,im) = wrk
          end do !ispin

          !call comm_summation( uVpsibox1(:,ilma,io,ik,im), uVpsibox2(:,ilma,io,ik,im), 2, info%icomm_r )
          !tmp1 = abs(uVpsibox2(:,ilma,io,ik,im))
          !tmp(1)=sum(abs(ppg%zekr_uV_so(:,ilma,ik,1,1)))
          !tmp(2)=sum(abs(ppg%zekr_uV_so(:,ilma,ik,2,1)))
          !call comm_summation(tmp,tmp2,2,info%icomm_r)
          !call comm_summation(wf_check,tmp3,2,info%icomm_r)
          !if ( io==1 ) then
          !if ( nproc_id_global==0)then
          !write(300,'(1x,3i6,6f20.10)') nproc_id_global,ilma,ia,tmp1,tmp2,tmp3
          !end if
          !end if
        end do !ilma

      end do !io
      end do !ik
      end do !im
#ifdef USE_OPENACC
!$acc end kernels
#endif

      call timer_end(LOG_UHPSI_PSEUDO)

      call timer_begin(LOG_UHPSI_PSEUDO_COMM)
      call comm_summation(uVpsibox1,uVpsibox2,size(uVpsibox2),info%icomm_r)
      call timer_end(LOG_UHPSI_PSEUDO_COMM)

      call timer_begin(LOG_UHPSI_PSEUDO)

      !write(*,*) "sum(abs(uVpsibox2(:,:,1,1,1)))",sum(abs(uVpsibox2(:,:,1,1,1)))
      !call mpi_finalize(ix); stop 'xxx'

#ifdef USE_OPENACC
!$acc kernels
!#acc loop collapse(3) private(im,ik,io,ilma,ia,uVpsi,wrk,j,iz,iy,ix)
#else
!$omp parallel do collapse(3) default(shared) private(im,ik,io,ilma,ia,uVpsi,wrk,j,iz,iy,ix)
#endif
      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e

        do ilma=1,Nlma

          ia = ppg%ia_tbl_so(ilma)

          uVpsi(1) = uVpsibox2(1,ilma,io,ik,im)
          uVpsi(2) = uVpsibox2(2,ilma,io,ik,im)

          do j=1,ppg%mps(ia)

            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)

            wrk = ppg%zekr_uV_so(j,ilma,ik,1,1)*( uVpsi(1) + uVpsi(2) )

            htpsi%zwf(ix,iy,iz,1,io,ik,im) = htpsi%zwf(ix,iy,iz,1,io,ik,im) + wrk

            wrk = ppg%zekr_uV_so(j,ilma,ik,2,1)*( uVpsi(1) + uVpsi(2) )

            htpsi%zwf(ix,iy,iz,2,io,ik,im) = htpsi%zwf(ix,iy,iz,2,io,ik,im) + wrk

          end do !j

        end do !ilma

      end do !io
      end do !ik
      end do !im
#ifdef USE_OPENACC
!$acc end kernels
#endif

      deallocate( uVpsibox2 )
      deallocate( uVpsibox1 )

    else !if ( .not. info%if_divide_rspace ) then

      allocate(uVpsibox1(Nspin,Nlma,io_s:io_e,ik_s:ik_e,im_s:im_e)); uVpsibox1=zero

#ifdef USE_OPENACC
!$acc kernels
!$acc loop collapse(3) private(im,ik,io,ilma,ia,ispin,wrk,j,iz,iy,ix,uVpsi)
#else
!$omp parallel do collapse(3) default(shared) private(im,ik,io,ilma,ia,ispin,wrk,j,iz,iy,ix,uVpsi)
#endif
      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e

        do ilma=1,Nlma

          ia = ppg%ia_tbl_so(ilma)

          !if ( ia==1.and.io==1 ) then
          !  do ispin=1,2
          !    tmp(ispin)=sum(abs(ppg%zekr_uV_so(:,ilma,ik,ispin,1)))
          !  end do
          !  call comm_summation( tmp,tmp1,2,info%icomm_r)
          !  if ( nproc_id_global==0 ) then
          !    !write(*,'("ilma,ia,mps,rank",4i6,2f20.15)') ilma,ia,ppg%mps(ia),nproc_id_global,tmp1
          !    write(100,'(3i6,2f20.10)') ilma,ia,nproc_id_global,tmp1
          !  end if
          !end if

          !wf_check(1:2)=zero
          do ispin=1,2

            wrk=zero
            do j=1,ppg%mps(ia)
              ix = ppg%jxyz(1,j,ia)
              iy = ppg%jxyz(2,j,ia)
              iz = ppg%jxyz(3,j,ia)
              wrk = wrk + conjg( ppg%zekr_uV_so(j,ilma,ik,ispin,1) ) &
                        * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
              !wf_check(ispin) = wf_check(ispin) + abs(tpsi%zwf(ix,iy,iz,ispin,io,ik,im))
            end do

            uVpsi(ispin) = wrk * ppg%rinv_uvu_so(ilma)

            uVpsibox1(ispin,ilma,io,ik,im)=uVpsi(ispin)

          end do !ispin

          !if ( io==1 ) then
          !tmp(1)=sum(abs(ppg%zekr_uV_so(:,ilma,ik,1,1)))
          !tmp(2)=sum(abs(ppg%zekr_uV_so(:,ilma,ik,2,1)))
          !call comm_summation(tmp,tmp1,2,info%icomm_r)
          !call comm_summation(wf_check,tmp3,2,info%icomm_r)
          !if ( nproc_id_global==0)then
          !write(300,'(1x,3i6,6f20.10)') nproc_id_global,ilma,ia,abs(uVpsi),tmp1,tmp3
          !end if
          !end if

          do j=1,ppg%mps(ia)

            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)

            wrk = ppg%zekr_uV_so(j,ilma,ik,1,1)*( uVpsi(1) + uVpsi(2) )

            htpsi%zwf(ix,iy,iz,1,io,ik,im) = htpsi%zwf(ix,iy,iz,1,io,ik,im) + wrk

            wrk = ppg%zekr_uV_so(j,ilma,ik,2,1)*( uVpsi(1) + uVpsi(2) )

            htpsi%zwf(ix,iy,iz,2,io,ik,im) = htpsi%zwf(ix,iy,iz,2,io,ik,im) + wrk

          end do !j

        end do !ilma

      end do !io
      end do !ik
      end do !im
#ifdef USE_OPENACC
!$acc end kernels
#endif

      !write(*,*) "sum(abs(uVpsibox1(:,:,1,1,1)))",sum(abs(uVpsibox1(:,:,1,1,1)))
      !call mpi_finalize(ix); stop 'xxx_serial'

    end if

    call timer_end(LOG_UHPSI_PSEUDO)

    return
  end subroutine pseudo_so


  subroutine calc_uVpsi_so(Nspin,info,ppg,tpsi,uVpsibox2)
    use structures, only: s_parallel_info, s_pp_grid, s_rgrid, s_orbital
    use communication, only: comm_summation
    use timer
    implicit none
    integer,intent(in) :: Nspin
    type(s_parallel_info),intent(in) :: info
    type(s_pp_grid),intent(in) :: ppg
    type(s_orbital),intent(in) :: tpsi
    complex(8),allocatable :: uVpsibox2(:,:,:,:,:)
    !
    integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e
    integer :: ilma,ia,j,ix,iy,iz,Nlma
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8) :: uVpsi(2),wrk
    complex(8),allocatable :: uVpsibox1(:,:,:,:,:)

    !write(*,*) "------------------ calc_uVpsi_so"

    call timer_begin(LOG_UHPSI_PSEUDO)

    im_s = info%im_s
    im_e = info%im_e
    ik_s = info%ik_s
    ik_e = info%ik_e
    io_s = info%io_s
    io_e = info%io_e

    Nlma = size(ppg%ia_tbl_so)

    allocate(uVpsibox2(Nspin,Nlma,io_s:io_e,ik_s:ik_e,im_s:im_e)); uVpsibox2=zero

    if ( info%if_divide_rspace ) then

      allocate(uVpsibox1(Nspin,Nlma,io_s:io_e,ik_s:ik_e,im_s:im_e)); uVpsibox1=zero

#ifdef USE_OPENACC
!$acc kernels
!#acc loop collapse(4) default(shared) private(im,ik,io,ilma,ia,ispin,wrk,j,iz,iy,ix)
#else
!$omp parallel do collapse(4) default(shared) private(im,ik,io,ilma,ia,ispin,wrk,j,iz,iy,ix)
#endif
      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e

        do ilma=1,Nlma

          ia = ppg%ia_tbl_so(ilma)

          do ispin=1,2
            wrk=zero
            do j=1,ppg%mps(ia)
              ix = ppg%jxyz(1,j,ia)
              iy = ppg%jxyz(2,j,ia)
              iz = ppg%jxyz(3,j,ia)
              wrk = wrk + conjg( ppg%zekr_uV_so(j,ilma,ik,ispin,1) ) &
                        * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
            end do
            wrk = wrk * ppg%rinv_uvu_so(ilma)
            uVpsibox1(ispin,ilma,io,ik,im) = wrk
          end do !ispin

        end do !ilma

      end do !io
      end do !ik
      end do !im
#ifdef USE_OPENACC
!$acc end kernels
#endif

      call timer_end(LOG_UHPSI_PSEUDO)

      call timer_begin(LOG_UHPSI_PSEUDO_COMM)
      call comm_summation(uVpsibox1,uVpsibox2,size(uVpsibox2),info%icomm_r)
      call timer_end(LOG_UHPSI_PSEUDO_COMM)

      call timer_begin(LOG_UHPSI_PSEUDO)

      deallocate( uVpsibox1 )

    else !if ( .not. info%if_divide_rspace ) then

#ifdef USE_OPENACC
!$acc kernels
!#acc loop collapse(4) private(im,ik,io,ilma,ia,ispin,wrk,j,iz,iy,ix,uVpsi)
#else
!$omp parallel do collapse(4) default(shared) private(im,ik,io,ilma,ia,ispin,wrk,j,iz,iy,ix,uVpsi)
#endif
      do im=im_s,im_e
      do ik=ik_s,ik_e
      do io=io_s,io_e

        do ilma=1,Nlma

          ia = ppg%ia_tbl_so(ilma)

          do ispin=1,2

            wrk=zero
            do j=1,ppg%mps(ia)
              ix = ppg%jxyz(1,j,ia)
              iy = ppg%jxyz(2,j,ia)
              iz = ppg%jxyz(3,j,ia)
              wrk = wrk + conjg( ppg%zekr_uV_so(j,ilma,ik,ispin,1) ) &
                        * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
            end do

            uVpsi(ispin) = wrk * ppg%rinv_uvu_so(ilma)

            uVpsibox2(ispin,ilma,io,ik,im)=uVpsi(ispin)

          end do !ispin

        end do !ilma

      end do !io
      end do !ik
      end do !im
#ifdef USE_OPENACC
!$acc end kernels
#endif

    end if

    call timer_end(LOG_UHPSI_PSEUDO)

    return
  end subroutine calc_uVpsi_so


  subroutine calc_force_uVpsi_so( force, ppg, mg, ik, ispin, gtpsi, uVpsibox2 )
    use structures, only: s_pp_grid, s_rgrid
    use communication, only: comm_summation
    use timer
    implicit none
    real(8),intent(inout) :: force(:,:)
    type(s_pp_grid),intent(in) :: ppg
    type(s_rgrid),intent(in) ::mg
    integer,intent(in) :: ik, ispin
    complex(8),intent(in) :: gtpsi(:,mg%is_array(1):,mg%is_array(2):,mg%is_array(3):)
    complex(8),intent(in) :: uVpsibox2(:,:)
    !
    integer :: ilma,ia,j,ix,iy,iz,Nlma
    complex(8),parameter :: zero=(0.0d0,0.0d0)
    complex(8) :: uVpsi(2),wrk1,wrk2,wrk3,conjg_zekr_uV_so

    !write(*,*) "------------------ calc_force_uVpsi_so"

    call timer_begin(LOG_UHPSI_PSEUDO)

    Nlma = size(ppg%ia_tbl_so)

    do ilma=1,Nlma

      ia = ppg%ia_tbl_so(ilma)

      uVpsi(1) = uVpsibox2(1,ilma)
      uVpsi(2) = uVpsibox2(2,ilma)

      wrk1=zero
      wrk2=zero
      wrk3=zero
      do j=1,ppg%mps(ia)
        ix = ppg%jxyz(1,j,ia)
        iy = ppg%jxyz(2,j,ia)
        iz = ppg%jxyz(3,j,ia)
        conjg_zekr_uV_so = conjg( ppg%zekr_uV_so(j,ilma,ik,ispin,1) )
        wrk1 = wrk1 + conjg_zekr_uV_so*gtpsi(1,ix,iy,iz)
        wrk2 = wrk2 + conjg_zekr_uV_so*gtpsi(2,ix,iy,iz)
        wrk3 = wrk3 + conjg_zekr_uV_so*gtpsi(3,ix,iy,iz)
      end do

      force(1,ia) = force(1,ia) + dble( wrk1*(uVpsi(1)+uVpsi(2)) )
      force(2,ia) = force(2,ia) + dble( wrk2*(uVpsi(1)+uVpsi(2)) )
      force(3,ia) = force(3,ia) + dble( wrk3*(uVpsi(1)+uVpsi(2)) )

    end do !ilma

    call timer_end(LOG_UHPSI_PSEUDO)

    return
  end subroutine calc_force_uVpsi_so

end module pseudo_pt_so_sub
