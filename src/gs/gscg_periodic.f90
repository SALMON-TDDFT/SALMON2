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

subroutine gscg_periodic(mg,info,spsi,iflag,itotmst,mst,hvol,ilsda,nproc_ob,nproc_ob_spin,iparaway_ob,elp3,  &
                         zxk_ob,zhxk_ob,zgk_ob,zpk_ob,zpko_ob,zhtpsi_ob,   &
                iup_array,idw_array,jup_array,jdw_array,kup_array,kdw_array,bnmat,cnmat,hgs,ppg,vlocal,  &
                nproc_mxin_mul,num_kpoints_rd,k_rd,ksquare)
  use inputoutput, only: ncg,ispin
  use structures, only: s_rgrid,s_wf_info,s_wavefunction,s_stencil,s_scalar,s_pp_grid
  use salmon_parallel, only: nproc_group_kgrid, nproc_group_korbital, nproc_id_korbital, nproc_group_k
  use salmon_communication, only: comm_bcast, comm_summation
  use misc_routines, only: get_wtime
  use hpsi_sub
  !$ use omp_lib
  implicit none
  type(s_rgrid),intent(inout) :: mg
  type(s_wf_info) :: info
  type(s_wavefunction),intent(inout) :: spsi
  type(s_stencil) :: stencil
  type(s_pp_grid) :: ppg
  integer,intent(inout) :: iflag
  integer,intent(in)    :: itotmst
  integer,intent(in)    :: mst(2)
  real(8),intent(in)    :: hvol
  integer,intent(in)    :: ilsda
  integer,intent(in)    :: nproc_ob
  integer,intent(in)    :: nproc_ob_spin
  integer,intent(in)    :: iparaway_ob
  real(8),intent(out)    :: elp3(3000)
  complex(8),intent(out) :: zxk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:info%numo)
  complex(8),intent(out) :: zhxk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:info%numo)
  complex(8),intent(out) :: zgk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:info%numo)
  complex(8),intent(out) :: zpk_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:info%numo)
  complex(8),intent(out) :: zpko_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:info%numo)
  complex(8),intent(out) :: zhtpsi_ob(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:info%numo)
  integer,intent(in)    :: iup_array(4)
  integer,intent(in)    :: idw_array(4)
  integer,intent(in)    :: jup_array(4)
  integer,intent(in)    :: jdw_array(4)
  integer,intent(in)    :: kup_array(4)
  integer,intent(in)    :: kdw_array(4)
  real(8),intent(in)    :: cnmat(0:12,0:12),bnmat(0:12,0:12)
  real(8),intent(in)    :: hgs(3)
  real(8),intent(in)    :: vlocal(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),ispin+1)
  integer,intent(in)    :: nproc_mxin_mul
  integer,intent(in)    :: num_kpoints_rd
  real(8),intent(in)    :: k_rd(3,num_kpoints_rd),ksquare(num_kpoints_rd)
  integer,parameter :: nd=4
  integer :: j,ind
  integer :: iter,iob,job
  integer :: ik
  integer :: ix,iy,iz
  integer :: is,iobsta(2),iobend(2)
  integer :: nspin
  type(s_wf_info)       :: info_ob
  type(s_wavefunction)  :: stpsi
  type(s_wavefunction)  :: shtpsi
  type(s_scalar),allocatable :: v(:)
  complex(8) :: sum0,sum1
  complex(8) :: sum_ob1(itotmst)
  complex(8) :: sum_obmat0(itotmst,itotmst),sum_obmat1(itotmst,itotmst)
  complex(8) :: xkxk_ob(itotmst),xkHxk_ob(itotmst),xkHpk_ob(itotmst),pkHpk_ob(itotmst),gkgk_ob(itotmst)
  complex(8) :: uk
  real(8) :: ev
  complex(8) :: cx,cp
  complex(8) :: zs_ob(itotmst)
  real(8) :: elp2(2000)
  complex(8):: zmatbox_m(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  integer :: iob_myob,job_myob
  integer :: iob_allob
  integer :: icorr_iob,icorr_job
  integer :: iroot
  integer :: is_sta,is_end
  integer :: iter_bak_ob(itotmst)
  
  allocate(stpsi%zwf(mg%is_array(1):mg%ie_array(1),  &
                     mg%is_array(2):mg%ie_array(2),  &
                     mg%is_array(3):mg%ie_array(3),1,1,1,1))
  allocate(shtpsi%zwf(mg%is_array(1):mg%ie_array(1),  &
                      mg%is_array(2):mg%ie_array(2),  &
                      mg%is_array(3):mg%ie_array(3),1,1,1,1))

  allocate(stencil%kAc(ik:ik,3))

  do j=1,3
    do ind=1,4
      stencil%lapt(ind,j) = cnmat(ind,4)/hgs(j)**2
      stencil%nabt(ind,j) = bnmat(ind,4)/hgs(j)
    end do
  end do

  info_ob%im_s = 1
  info_ob%im_e = 1
  info_ob%numm = 1
  info_ob%ik_s = 1
  info_ob%ik_e = 1
  info_ob%numk = 1
  info_ob%io_s = 1
  info_ob%io_e = 1
  info_ob%numo = 1

  mg%is_overlap = mg%is - 4
  mg%ie_overlap = mg%ie + 4

  allocate(mg%idx(mg%is_overlap(1):mg%ie_overlap(1)) &
          ,mg%idy(mg%is_overlap(2):mg%ie_overlap(2)) &
          ,mg%idz(mg%is_overlap(3):mg%ie_overlap(3)))
  do j=mg%is_overlap(1),mg%ie_overlap(1)
    mg%idx(j) = j
  end do
  do j=mg%is_overlap(2),mg%ie_overlap(2)
    mg%idy(j) = j
  end do
  do j=mg%is_overlap(3),mg%ie_overlap(3)
    mg%idz(j) = j
  end do

  info_ob%if_divide_rspace = nproc_mxin_mul.ne.1
  info_ob%irank_overlap(1) = iup_array(1)
  info_ob%irank_overlap(2) = idw_array(1)
  info_ob%irank_overlap(3) = jup_array(1)
  info_ob%irank_overlap(4) = jdw_array(1)
  info_ob%irank_overlap(5) = kup_array(1)
  info_ob%irank_overlap(6) = kdw_array(1)
  info_ob%icomm_overlap = nproc_group_korbital
  info_ob%icomm_pseudo = nproc_group_korbital

  nspin=1
  allocate(v(1))
  allocate(v(1)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  
  call set_isstaend(is_sta,is_end,ilsda,nproc_ob,nproc_ob_spin)
  
  !$OMP parallel do private(iz,iy,ix) collapse(2)
  do iz=mg%is_array(3),mg%ie_array(3)
  do iy=mg%is_array(2),mg%ie_array(2)
  do ix=mg%is_array(1),mg%ie_array(1)
    stpsi%zwf(ix,iy,iz,1,1,1,1)=0.d0
  end do
  end do
  end do
  
  elp2(:)=0d0
  elp2(1)=get_wtime()
  
  if(ilsda == 0)then
    iobsta(1)=1
    iobend(1)=itotmst
  else if(ilsda == 1)then
    iobsta(1)=1
    iobend(1)=mst(1)
    iobsta(2)=mst(1)+1
    iobend(2)=itotmst
  end if
  
  do ik=info%ik_s,info%ik_e
  do is=is_sta,is_end

    stencil%lap0 = 0.5d0*ksquare(ik) -0.5d0*cNmat(0,Nd)*(1.d0/Hgs(1)**2+1.d0/Hgs(2)**2+1.d0/Hgs(3)**2)
    do j=1,3
      stencil%kAc(1,j) = k_rd(j,ik)
    end do

    iter_bak_ob(:)=0 
  
    do iob_myob=1,info%numo
      call calc_allob(iob_myob,iob_allob)
    
      elp2(2)=get_wtime()
      
    !$omp parallel do private(iz,iy,ix) collapse(2) 
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        zxk_ob(ix,iy,iz,iob_myob)=spsi%zwf(ix,iy,iz,1,iob_myob,ik,1)
        stpsi%zwf(ix,iy,iz,1,1,1,1)=zxk_ob(ix,iy,iz,iob_myob)
      end do
      end do
      end do
     
      if(iob_allob<=mst(1))then
  !$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          v(1)%f(ix,iy,iz) = vlocal(ix,iy,iz,1)
        end do
        end do
        end do
      else
  !$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          v(1)%f(ix,iy,iz) = vlocal(ix,iy,iz,2)
        end do
        end do
        end do
      end if

      call hpsi(stpsi,shtpsi,info_ob,mg,v,nspin,stencil,ppg)
      
    !$omp parallel do private(iz,iy,ix) collapse(2) 
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        zhxk_ob(ix,iy,iz,iob_myob)=shtpsi%zwf(ix,iy,iz,1,1,1,1)
      end do
      end do
      end do
  
    end do
    call inner_product5(mg,itotmst,info%numo,zxk_ob,zhxk_ob,xkHxk_ob,hvol)

    Iteration : do iter=1,Ncg
      do iob_myob=1,info%numo
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
              sum0=sum0+conjg(spsi%zwf(ix,iy,iz,1,job,ik,1))*zgk_ob(ix,iy,iz,iob)
            end do
            end do
            end do
            sum_obmat0(iob,job)=sum0*hvol
          end do
        end do
        elp3(1507)=get_wtime()
        elp3(1557)=elp3(1557)+elp3(1507)-elp3(1506)
        call comm_summation(sum_obmat0,sum_obmat1,itotmst*itotmst,nproc_group_k)
        elp3(1508)=get_wtime()
        elp3(1558)=elp3(1558)+elp3(1508)-elp3(1507)
        do iob=iobsta(is),iobend(is)
          do job=iobsta(is),iob-1
    !$omp parallel do collapse(2)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              zgk_ob(ix,iy,iz,iob)=zgk_ob(ix,iy,iz,iob)-sum_obmat1(iob,job)*spsi%zwf(ix,iy,iz,1,job,ik,1)
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
          call check_corrkob(iob,ik,icorr_iob,ilsda,nproc_ob,iparaway_ob,itotmst,info%ik_s,info%ik_e,nproc_ob_spin,mst)
          do job=iobsta(is),iob-1
            call calc_myob(job,job_myob,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
            call check_corrkob(job,ik,icorr_job,ilsda,nproc_ob,iparaway_ob,itotmst,info%ik_s,info%ik_e,nproc_ob_spin,mst)
            if(icorr_job==1)then
    !$omp parallel do private(iz,iy,ix) collapse(2)
              do iz=mg%is(3),mg%ie(3)
              do iy=mg%is(2),mg%ie(2)
              do ix=mg%is(1),mg%ie(1)
                zmatbox_m(ix,iy,iz)=spsi%zwf(ix,iy,iz,1,job_myob,ik,1)
              end do
              end do
              end do
            end if
            call calc_iroot(job,iroot,ilsda,nproc_ob,iparaway_ob,itotmst,nproc_ob_spin,mst)
            call comm_bcast(zmatbox_m,nproc_group_kgrid,iroot)
            sum0=0.d0
    !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+ : sum0)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              sum0=sum0+conjg(zmatbox_m(ix,iy,iz))*zgk_ob(ix,iy,iz,iob_myob)
            end do
            end do
            end do
            sum0=sum0*hvol
            call comm_summation(sum0,sum1,nproc_group_korbital)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              zgk_ob(ix,iy,iz,iob_myob)=zgk_ob(ix,iy,iz,iob_myob)-sum1*zmatbox_m(ix,iy,iz)
            end do
            end do
            end do
          end do
        end do
      end if
      call inner_product5(mg,itotmst,info%numo,zgk_ob,zgk_ob,sum_ob1,hvol)
        
      do iob_myob=1,info%numo
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

      call inner_product5(mg,itotmst,info%numo,zxk_ob,zpk_ob,zs_ob,hvol)

      do iob_myob=1,info%numo
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
      call inner_product5(mg,itotmst,info%numo,zpko_ob,zpko_ob,sum_ob1,hvol)

      do iob_myob=1,info%numo
        call calc_allob(iob_myob,iob_allob)
    !$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          zpko_ob(ix,iy,iz,iob_myob)=zpko_ob(ix,iy,iz,iob_myob)/sqrt(sum_ob1(iob_allob))
          stpsi%zwf(ix,iy,iz,1,1,1,1)=zpko_ob(ix,iy,iz,iob_myob)
        end do
        end do
        end do

        if(iob_allob<=mst(1))then
  !$OMP parallel do private(iz,iy,ix) collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            v(1)%f(ix,iy,iz) = vlocal(ix,iy,iz,1)
          end do
          end do
          end do
        else
  !$OMP parallel do private(iz,iy,ix) collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            v(1)%f(ix,iy,iz) = vlocal(ix,iy,iz,2)
          end do
          end do
          end do
        end if

        call hpsi(stpsi,shtpsi,info_ob,mg,v,nspin,stencil,ppg)

    !$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          zhtpsi_ob(ix,iy,iz,iob_myob)=shtpsi%zwf(ix,iy,iz,1,1,1,1)
        end do
        end do
        end do
      end do
      call inner_product5(mg,itotmst,info%numo,zxk_ob,zhtpsi_ob,xkHpk_ob,hvol)
      call inner_product5(mg,itotmst,info%numo,zpko_ob,zhtpsi_ob,pkHpk_ob,hvol)
        
    
      do iob_myob=1,info%numo
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
    
      call inner_product5(mg,itotmst,info%numo,zxk_ob,zhxk_ob,xkHxk_ob,hvol)
      call inner_product5(mg,itotmst,info%numo,zxk_ob,zxk_ob,xkxk_ob,hvol)

      do iob_myob=1,info%numo
        call calc_allob(iob_myob,iob_allob)

        if(abs(xkxk_ob(iob_allob))<=1.d30)then
    !$OMP parallel do private(iz,iy,ix) collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            spsi%zwf(ix,iy,iz,1,iob_myob,ik,1)=zxk_ob(ix,iy,iz,iob_myob)/sqrt(xkxk_ob(iob_allob))
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
  
  deallocate(stpsi%zwf,shtpsi%zwf)
  deallocate(stencil%kAc)
  deallocate(mg%idx,mg%idy,mg%idz)
  deallocate(v(1)%f)
  deallocate(v)

  return
  
end subroutine gscg_periodic
