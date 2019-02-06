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
!=======================================================================
!======================================= Conjugate-Gradient minimization

subroutine dtcg_periodic(mg,info,spsi,iflag,itotmst,mst,hvol,ilsda,nproc_ob,nproc_ob_spin,iparaway_ob)
  use inputoutput, only: ncg
  use structures, only: s_rgrid,s_wf_info,s_wavefunction
  use salmon_parallel, only: nproc_group_kgrid, nproc_group_korbital
  use salmon_communication, only: comm_bcast, comm_summation
  use misc_routines, only: get_wtime
  !$ use omp_lib
  implicit none
  type(s_rgrid),intent(in)           :: mg
  type(s_wf_info),intent(in)         :: info
  type(s_wavefunction),intent(inout) :: spsi
  integer,intent(inout) :: iflag
  integer,intent(in)    :: itotmst
  integer,intent(in)    :: mst(2)
  real(8),intent(in)    :: hvol
  integer,intent(in)    :: ilsda
  integer,intent(in)    :: nproc_ob
  integer,intent(in)    :: nproc_ob_spin
  integer,intent(in)    :: iparaway_ob
  integer :: iter,p,q
  integer :: ik
  integer :: ix,iy,iz
  integer :: is,pstart(2),pend(2)
  complex(8) :: sum0,sum1,xkhxk,xkxk,rk,gkgk,pkhpk
  complex(8) :: uk
  real(8) :: ev
  complex(8) :: cx,cp,xkhpk,zs,xktxk
  !real(8) :: xk(ML),hxk(ML),gk(ML),pk(ML)
  complex(8) , allocatable :: xk(:,:,:),hxk(:,:,:),gk(:,:,:),pk(:,:,:)
  complex(8) , allocatable :: txk(:,:,:),htpsi(:,:,:),pko(:,:,:)
  complex(8) , allocatable :: gk2(:,:,:)
  complex(8) , allocatable :: ttpsi(:,:,:)
  real(8) :: elp2(2000)
  complex(8):: tpsi(mg%is_array(1):mg%ie_array(1),  &
                    mg%is_array(2):mg%ie_array(2),  &
                    mg%is_array(3):mg%ie_array(3))
  complex(8):: zmatbox_m(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  integer :: p_myob,q_myob
  integer :: icorr_p,icorr_q
  integer :: iroot
  integer :: is_sta,is_end
  
  allocate (xk(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (hxk(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (gk(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (gk2(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (pk(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  
  allocate (txk(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (htpsi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate (pko(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  
  allocate (ttpsi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  
  call set_isstaend(is_sta,is_end,ilsda,nproc_ob,nproc_ob_spin)
  
  !$OMP parallel do
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
    pstart(1)=1
    pend(1)=itotMST
  else if(ilsda == 1)then
    pstart(1)=1
    pend(1)=MST(1)
    pstart(2)=MST(1)+1
    pend(2)=itotMST
  end if
  
  do ik=info%ik_s,info%ik_e
  do is=is_sta,is_end
  
  orbital : do p=pstart(is),pend(is)
    call calc_myob(p,p_myob,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
    call check_corrkob(p,ik,icorr_p,ilsda,nproc_ob,iparaway_ob,itotmst,info%ik_s,info%ik_e,nproc_ob_spin,mst)
  
    elp2(2)=get_wtime()
  
    if(nproc_ob==1)then
      do q=pstart(is),p-1
        sum0=0.d0
  !$omp parallel do reduction(+ : sum0)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          sum0=sum0+conjg(spsi%zwf(ix,iy,iz,1,q,ik,1))*spsi%zwf(ix,iy,iz,1,p,ik,1)
        end do
        end do
        end do
        sum0=sum0*hvol
        call comm_summation(sum0,sum1,nproc_group_korbital)
  !$omp parallel do
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          spsi%zwf(ix,iy,iz,1,p,ik,1)=spsi%zwf(ix,iy,iz,1,p,ik,1)-sum1*spsi%zwf(ix,iy,iz,1,q,ik,1)
        end do
        end do
        end do
      end do
    else
      do q=pstart(is),p-1
        call calc_myob(q,q_myob,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
        call check_corrkob(q,ik,icorr_q,ilsda,nproc_ob,iparaway_ob,itotmst,info%ik_s,info%ik_e,nproc_ob_spin,mst)
        if(icorr_q==1)then
  !$omp parallel do
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            zmatbox_m(ix,iy,iz)=spsi%zwf(ix,iy,iz,1,q_myob,ik,1)
          end do
          end do
          end do
        end if
        call calc_iroot(q,iroot,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
        call comm_bcast(zmatbox_m,nproc_group_kgrid,iroot)
        sum0=0.d0
        if(icorr_p==1)then
  !$omp parallel do reduction(+ : sum0)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            sum0=sum0+conjg(zmatbox_m(ix,iy,iz))*spsi%zwf(ix,iy,iz,1,p_myob,ik,1)
          end do
          end do
          end do
        end if
        sum0=sum0*hvol
        call comm_summation(sum0,sum1,nproc_group_korbital)
        if(icorr_p==1)then
  !$omp parallel do
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            spsi%zwf(ix,iy,iz,1,p_myob,ik,1)=spsi%zwf(ix,iy,iz,1,p_myob,ik,1)-sum1*zmatbox_m(ix,iy,iz)
          end do
          end do
          end do
        end if
      end do 
    end if
    sum0=0.d0
    if(icorr_p==1)then
  !$omp parallel do reduction(+ : sum0)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        sum0=sum0+abs(spsi%zwf(ix,iy,iz,1,p_myob,ik,1))**2
      end do
      end do
      end do
    end if
    sum0=sum0*hvol
    call comm_summation(sum0,sum1,nproc_group_korbital)
    if(icorr_p==1)then
  !$omp parallel do 
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        xk(ix,iy,iz)=spsi%zwf(ix,iy,iz,1,p_myob,ik,1)/sqrt(sum1)
        tpsi(ix,iy,iz)=xk(ix,iy,iz)
      end do
      end do
      end do
   
      call c_hpsi2_buf(tpsi,hxk,p,ik,0,0)
    
  !$omp parallel do 
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        txk(ix,iy,iz)=ttpsi(ix,iy,iz)
      end do
      end do
      end do
    end if
  
    call calc_iroot(p,iroot,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
    call comm_bcast(xk,nproc_group_kgrid,iroot)
    call comm_bcast(hxk,nproc_group_kgrid,iroot)
    call comm_bcast(txk,nproc_group_kgrid,iroot)
  
    call inner_product4(mg,xk,hxk,xkhxk,hvol)
    call inner_product4(mg,xk,txk,xktxk,hvol)
    
    iteration : do iter=1,ncg
  
  !$OMP parallel do
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        gk(ix,iy,iz) = hxk(ix,iy,iz) - xkhxk*xk(ix,iy,iz) 
      end do
      end do
      end do
  
      if(nproc_ob==1)then
        do q=pstart(is),p-1
          sum0=0.d0
  !$omp parallel do reduction(+ : sum0)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            sum0=sum0+conjg(spsi%zwf(ix,iy,iz,1,q,ik,1))*gk(ix,iy,iz)
          end do
          end do
          end do
          sum0=sum0*hvol
          call comm_summation(sum0,sum1,nproc_group_korbital)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            gk(ix,iy,iz)=gk(ix,iy,iz)-sum1*spsi%zwf(ix,iy,iz,1,q,ik,1)
          end do
          end do
          end do
        end do
      else
        do q=pstart(is),p-1
          call calc_myob(q,q_myob,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
          call check_corrkob(q,ik,icorr_q,ilsda,nproc_ob,iparaway_ob,itotmst,info%ik_s,info%ik_e,nproc_ob_spin,mst)
          if(icorr_q==1)then
  !$omp parallel do
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              zmatbox_m(ix,iy,iz)=spsi%zwf(ix,iy,iz,1,q_myob,ik,1)
            end do
            end do
            end do
          end if
          call calc_iroot(q,iroot,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
          call comm_bcast(zmatbox_m,nproc_group_kgrid,iroot)
          sum0=0.d0
  !$omp parallel do reduction(+ : sum0)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            sum0=sum0+conjg(zmatbox_m(ix,iy,iz))*gk(ix,iy,iz)
          end do
          end do
          end do
          sum0=sum0*hvol
          call comm_summation(sum0,sum1,nproc_group_korbital)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            gk(ix,iy,iz)=gk(ix,iy,iz)-sum1*zmatbox_m(ix,iy,iz)
          end do
          end do
          end do
        end do
      end if
      call inner_product4(mg,gk,gk,sum1,hvol)
      
      if(iter==1)then
  !$OMP parallel do
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          pk(ix,iy,iz)=gk(ix,iy,iz)
        end do
        end do
        end do
      else
        uk=sum1/gkgk
  !$OMP parallel do
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          pk(ix,iy,iz)=gk(ix,iy,iz)+uk*pk(ix,iy,iz)
        end do
        end do
        end do
      end if
      gkgk=sum1
      call inner_product4(mg,xk,pk,zs,hvol)
  !$OMP parallel do
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        pko(ix,iy,iz)=pk(ix,iy,iz)-zs*xk(ix,iy,iz)
      end do
      end do
      end do
      call inner_product4(mg,pko,pko,sum1,hvol)
  !$OMP parallel do
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        pko(ix,iy,iz)=pko(ix,iy,iz)/sqrt(sum1)
        tpsi(ix,iy,iz)=pko(ix,iy,iz)
      end do
      end do
      end do
      call c_hpsi2_buf(tpsi,htpsi,p,ik,0,0)
      call inner_product4(mg,xk,htpsi,xkhpk,hvol)
      call inner_product4(mg,pko,htpsi,pkhpk,hvol)
      
  
      ev=0.5d0*((xkhxk+pkhpk)-sqrt((xkhxk-pkhpk)**2+4.d0*abs(xkhpk)**2))
      cx=xkhpk/(ev-xkhxk)
      cp=1.d0/sqrt(1.d0+abs(cx)**2)
      cx=cx*cp
  
      
  !$OMP parallel do
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        xk(ix,iy,iz)=cx*xk(ix,iy,iz)+cp*pko(ix,iy,iz)
        hxk(ix,iy,iz)=cx*hxk(ix,iy,iz)+cp*htpsi(ix,iy,iz)
        txk(ix,iy,iz)=cx*txk(ix,iy,iz)+cp*ttpsi(ix,iy,iz)
      end do
      end do
      end do
  
      call inner_product4(mg,xk,hxk,xkhxk,hvol)
      call inner_product4(mg,xk,txk,xktxk,hvol)
      call inner_product4(mg,xk,xk,xkxk,hvol)
      rk=xkhxk/xkxk
  
    end do iteration
  
    call inner_product4(mg,xk(mg%is(1),mg%is(2),mg%is(3)),xk(mg%is(1),mg%is(2),mg%is(3)),sum0,hvol)
    if(icorr_p==1)then
  !$OMP parallel do
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        spsi%zwf(ix,iy,iz,1,p_myob,ik,1)=xk(ix,iy,iz)/sqrt(sum0)
      end do
      end do
      end do
    end if
  end do orbital
  
  end do
  end do
  
  if(iflag.eq.1) then
    iflag=0
  end if
  
  deallocate (xk,hxk,gk,pk,gk2)
  
  return
  
end subroutine dtcg_periodic

