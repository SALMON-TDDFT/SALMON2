!
!  Copyright 2018 SALMON developers
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
!=======================================================================
!======================================= Conjugate-Gradient minimization

subroutine gscg_periodic(mg,psi_in,iflag)
  use structures, only: s_rgrid
  use salmon_parallel, only: nproc_group_kgrid, nproc_group_korbital, nproc_id_korbital, nproc_group_k
  use salmon_communication, only: comm_bcast, comm_summation
  use misc_routines, only: get_wtime
  use scf_data
  use new_world_sub
  use allocate_mat_sub
  use hpsi2_sub
  !$ use omp_lib
  implicit none
  
  type(s_rgrid),intent(in) :: mg
  complex(8) :: psi_in(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),  &
                 1:iobnum,k_sta:k_end)
  integer :: iter,iob,job,iflag
  integer :: ik
  integer :: ix,iy,iz
  integer :: is,iobsta(2),iobend(2)
  complex(8) :: sum0,sum1
  complex(8) :: sum_ob1(itotMST)
  complex(8) :: sum_obmat0(itotMST,itotMST),sum_obmat1(itotMST,itotMST)
  complex(8) :: xkxk_ob(itotMST),xkHxk_ob(itotMST),xkHpk_ob(itotMST),pkHpk_ob(itotMST),gkgk_ob(itotMST)
  complex(8) :: uk
  real(8) :: ev
  complex(8) :: cx,cp
  complex(8) :: zs_ob(itotMST)
  complex(8) , allocatable :: htpsi(:,:,:)
  real(8) :: elp2(2000)
  complex(8):: tpsi(mg%is_array(1):mg%ie_array(1),  &
                    mg%is_array(2):mg%ie_array(2),  &
                    mg%is_array(3):mg%ie_array(3))
  integer :: iob_myob,job_myob
  integer :: iob_allob
  integer :: icorr_iob,icorr_job
  integer :: iroot
  integer :: is_sta,is_end
  integer :: iter_bak_ob(itotMST)
  
  allocate (htpsi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  
  call set_isstaend(is_sta,is_end,ilsda,nproc_ob,nproc_ob_spin)
  
  iwk_size=2
  call make_iwksta_iwkend
  
  !$OMP parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is_array(3),mg%ie_array(3)
  do iy=mg%is_array(2),mg%ie_array(2)
  do ix=mg%is_array(1),mg%ie_array(1)
    tpsi(ix,iy,iz)=0.d0
  end do
  end do
  end do
  
  elp2(:)=0d0
  elp2(1)=get_wtime()
  
  if(ilsda == 0)then
    iobsta(1)=1
    iobend(1)=itotMST
  else if(ilsda == 1)then
    iobsta(1)=1
    iobend(1)=MST(1)
    iobsta(2)=MST(1)+1
    iobend(2)=itotMST
  end if
  
  do ik=k_sta,k_end
  do is=is_sta,is_end

    iter_bak_ob(:)=0 
  
    do iob_myob=1,iobnum
      call calc_allob(iob_myob,iob_allob)
    
      elp2(2)=get_wtime()
      
    !$omp parallel do private(iz,iy,ix) collapse(2) 
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        zxk_ob(ix,iy,iz,iob_myob)=psi_in(ix,iy,iz,iob_myob,ik)
        tpsi(ix,iy,iz)=zxk_ob(ix,iy,iz,iob_myob)
      end do
      end do
      end do
     
      call hpsi2(tpsi,htpsi,iob_allob,ik,0,0)
      
    !$omp parallel do private(iz,iy,ix) collapse(2) 
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        zhxk_ob(ix,iy,iz,iob_myob)=htpsi(ix,iy,iz)
      end do
      end do
      end do
  
    end do
    call inner_product5(mg,itotmst,iobnum,zxk_ob,zhxk_ob,xkHxk_ob,hvol)

    Iteration : do iter=1,Ncg
      do iob_myob=1,iobnum
        call calc_allob(iob_myob,iob_allob)
    
    !$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          zgk_ob(ix,iy,iz,iob_myob) = zhxk_ob(ix,iy,iz,iob_myob) - xkHxk_ob(iob_allob)*zxk_ob(ix,iy,iz,iob_myob) 
        end do
        end do
        end do
      end do

      if(nproc_ob==1)then
        sum_obmat0(:,:)=0.d0
        elp3(1506)=get_wtime()
        do iob=iobsta(is),iobend(is)
          do job=iobsta(is),iob-1
            sum0=0.d0
    !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+ : sum0)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              sum0=sum0+conjg(psi_in(ix,iy,iz,job,ik))*zgk_ob(ix,iy,iz,iob)
            end do
            end do
            end do
            sum_obmat0(iob,job)=sum0*Hvol
          end do
        end do
        elp3(1507)=get_wtime()
        elp3(1557)=elp3(1557)+elp3(1507)-elp3(1506)
        call comm_summation(sum_obmat0,sum_obmat1,itotMST*itotMST,nproc_group_k)
        elp3(1508)=get_wtime()
        elp3(1558)=elp3(1558)+elp3(1508)-elp3(1507)
        do iob=iobsta(is),iobend(is)
          do job=iobsta(is),iob-1
    !$omp parallel do collapse(2)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              zgk_ob(ix,iy,iz,iob)=zgk_ob(ix,iy,iz,iob)-sum_obmat1(iob,job)*psi_in(ix,iy,iz,job,ik)
            end do
            end do
            end do
          end do
        end do
        elp3(1507)=get_wtime()
        elp3(1557)=elp3(1557)+elp3(1507)-elp3(1508)
      else
        do iob=iobsta(is),iobend(is)
          call calc_myob(iob,iob_myob,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
          call check_corrkob(iob,ik,icorr_iob,ilsda,nproc_ob,iparaway_ob,itotmst,k_sta,k_end,nproc_ob_spin,mst)
          do job=iobsta(is),iob-1
            call calc_myob(job,job_myob,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
            call check_corrkob(job,ik,icorr_job,ilsda,nproc_ob,iparaway_ob,itotmst,k_sta,k_end,nproc_ob_spin,mst)
            if(icorr_job==1)then
    !$omp parallel do private(iz,iy,ix) collapse(2)
              do iz=mg%is(3),mg%ie(3)
              do iy=mg%is(2),mg%ie(2)
              do ix=mg%is(1),mg%ie(1)
                cmatbox_m(ix,iy,iz)=psi_in(ix,iy,iz,job_myob,ik)
              end do
              end do
              end do
            end if
            call calc_iroot(job,iroot,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
            call comm_bcast(cmatbox_m,nproc_group_kgrid,iroot)
            sum0=0.d0
    !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+ : sum0)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              sum0=sum0+conjg(cmatbox_m(ix,iy,iz))*zgk_ob(ix,iy,iz,iob_myob)
            end do
            end do
            end do
            sum0=sum0*Hvol
            call comm_summation(sum0,sum1,nproc_group_korbital)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              zgk_ob(ix,iy,iz,iob_myob)=zgk_ob(ix,iy,iz,iob_myob)-sum1*cmatbox_m(ix,iy,iz)
            end do
            end do
            end do
          end do
        end do
      end if
      call inner_product5(mg,itotmst,iobnum,zgk_ob,zgk_ob,sum_ob1,hvol)
        
      do iob_myob=1,iobnum
        call calc_allob(iob_myob,iob_allob)
    
        if(iter==1)then
    !$OMP parallel do private(iz,iy,ix) collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            zpk_ob(ix,iy,iz,iob_myob)=zgk_ob(ix,iy,iz,iob_myob)
          end do
          end do
          end do
        else
          uk=sum_ob1(iob_allob)/gkgk_ob(iob_allob)
    !$OMP parallel do private(iz,iy,ix) collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            zpk_ob(ix,iy,iz,iob_myob)=zgk_ob(ix,iy,iz,iob_myob)+uk*zpk_ob(ix,iy,iz,iob_myob)
          end do
          end do
          end do
        end if
        gkgk_ob(iob_allob)=sum_ob1(iob_allob)
      end do

      call inner_product5(mg,itotmst,iobnum,zxk_ob,zpk_ob,zs_ob,hvol)

      do iob_myob=1,iobnum
        call calc_allob(iob_myob,iob_allob)
    !$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          zpko_ob(ix,iy,iz,iob_myob)=zpk_ob(ix,iy,iz,iob_myob)-zs_ob(iob_allob)*zxk_ob(ix,iy,iz,iob_myob)
        end do
        end do
        end do
      end do
      call inner_product5(mg,itotmst,iobnum,zpko_ob,zpko_ob,sum_ob1,hvol)

      do iob_myob=1,iobnum
        call calc_allob(iob_myob,iob_allob)
    !$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          zpko_ob(ix,iy,iz,iob_myob)=zpko_ob(ix,iy,iz,iob_myob)/sqrt(sum_ob1(iob_allob))
          tpsi(ix,iy,iz)=zpko_ob(ix,iy,iz,iob_myob)
        end do
        end do
        end do

        call hpsi2(tpsi,htpsi,iob_allob,ik,0,0)

    !$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          zhtpsi_ob(ix,iy,iz,iob_myob)=htpsi(ix,iy,iz)
        end do
        end do
        end do
      end do
      call inner_product5(mg,itotmst,iobnum,zxk_ob,zhtpsi_ob,xkHpk_ob,hvol)
      call inner_product5(mg,itotmst,iobnum,zpko_ob,zhtpsi_ob,pkHpk_ob,hvol)
        
    
      do iob_myob=1,iobnum
        call calc_allob(iob_myob,iob_allob)

        ev=0.5d0*((xkHxk_ob(iob_allob)+pkHpk_ob(iob_allob))   &
                 -sqrt((xkHxk_ob(iob_allob)-pkHpk_ob(iob_allob))**2+4.d0*abs(xkHpk_ob(iob_allob))**2))
        cx=xkHpk_ob(iob_allob)/(ev-xkHxk_ob(iob_allob))
        cp=1.d0/sqrt(1.d0+abs(cx)**2)
        cx=cx*cp
        
    !$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          zxk_ob(ix,iy,iz,iob_myob)=cx*zxk_ob(ix,iy,iz,iob_myob)+cp*zpko_ob(ix,iy,iz,iob_myob)
          zhxk_ob(ix,iy,iz,iob_myob)=cx*zhxk_ob(ix,iy,iz,iob_myob)+cp*zhtpsi_ob(ix,iy,iz,iob_myob)
        end do
        end do
        end do
      end do 
    
      call inner_product5(mg,itotmst,iobnum,zxk_ob,zhxk_ob,xkHxk_ob,hvol)
      call inner_product5(mg,itotmst,iobnum,zxk_ob,zxk_ob,xkxk_ob,hvol)

      do iob_myob=1,iobnum
        call calc_allob(iob_myob,iob_allob)

        if(abs(xkxk_ob(iob_allob))<=1.d30)then
    !$OMP parallel do private(iz,iy,ix) collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            psi_in(ix,iy,iz,iob_myob,ik)=zxk_ob(ix,iy,iz,iob_myob)/sqrt(xkxk_ob(iob_allob))
          end do
          end do
          end do
          iter_bak_ob(iob_allob)=iter
        else
          if(nproc_id_korbital==0)then
            write(*,'("CG termination. orbital ",i6," iter ",i3)') iob_allob,iter_bak_ob(iob_allob)+1
          end if
        end if
    
      end do
    end do Iteration
  
  end do
  end do
  
  
  if(iflag.eq.1) then
    iflag=0
  end if
  
  deallocate (htpsi)
  
  return
  
end subroutine gscg_periodic
