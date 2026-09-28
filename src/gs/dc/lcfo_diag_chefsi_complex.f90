!
!  Copyright 2026 SALMON developers
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

module lcfo_diag_chefsi_complex
  implicit none
  private

  integer, parameter :: block_size = 64
  integer, parameter :: lanczos_max_steps=30
  integer(8), parameter :: workspace_limit = 64_8*1024_8*1024_8
  real(8), parameter :: lock_activation_tolerance=1d-3
  real(8), parameter :: lock_gap_relative_tolerance=1d-6
  real(8), parameter :: filter_growth_margin=1d6

  type :: s_grid
    integer :: context=-1,nprow=0,npcol=0,myrow=-1,mycol=-1
  end type s_grid

  type :: s_layout
    integer :: desc(9)=0,nrow=0,ncol=0,lda=1
  end type s_layout

  public :: diag_chefsi_complex

contains

  subroutine diag_chefsi_complex(dc,ik,nspin,filter_degree,filter_chunk_size, &
      max_cycle,residual_tolerance,n_basis,n_mat,n_halo,halo_src,halo_dst, &
      halo_root_src,halo_dvec,h_diag,h_halo,esp_k,coef_frag,ortho_error, &
      residual_error,status)
    use communication, only: comm_bcast,comm_create_group,comm_free_group, &
      comm_get_max,comm_get_min,comm_irecv,comm_isend,comm_summation,comm_wait_all
    use eigen_subdiag_sub, only: eigen_zheev
    use ieee_arithmetic, only: ieee_is_finite
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ik,nspin,filter_degree,filter_chunk_size,max_cycle
    real(8), intent(in) :: residual_tolerance
    integer, intent(in) :: n_basis(:,:),n_mat(:),n_halo
    integer, intent(in) :: halo_src(:),halo_dst(:),halo_root_src(:),halo_dvec(:,:)
    complex(8), intent(in) :: h_diag(:,:,:),h_halo(:,:,:,:)
    real(8), intent(out) :: esp_k(:,:),ortho_error(:),residual_error(:)
    complex(8), intent(inout) :: coef_frag(:,:,:)
    integer, intent(out) :: status
    type(s_grid) :: natural,dense
    type(s_layout) :: ln,ld,lg
    complex(8), allocatable :: hdiag_sym(:,:,:),hrow(:,:,:,:),hrecv(:,:,:,:)
    complex(8), allocatable :: x(:,:),x_backup(:,:),hx(:,:),xd(:,:),hxd(:,:),xrot(:,:)
    complex(8), allocatable :: projected(:,:),z(:,:),q0(:,:),q1(:,:),q2(:,:)
    complex(8), allocatable :: tau(:),work(:)
    real(8), allocatable :: eval(:),rwork(:),residuals(:)
    integer, allocatable :: reqs(:),reqr(:)
    integer :: icomm_row,max_basis,npadded,ispin,nsub,nbuffer,cycle,flag
    integer :: max_local_col,chunk_width,first,last,nchunk,j,local_col,gcol
    integer :: nlocked,active_first,active_count,max_active,nchunks,filter_attempt
    integer :: info,lwork,lrwork,nrow,ncol,owner,local_state,offset,nb
    real(8) :: lower,upper,gersh_upper,cutoff,center,radius,local_value,global_value
    real(8) :: scaled_bound,log_limit,max_amplitude
    real(8) :: max_resid,orth_local,diag_value,max_in(1),max_out(1)
    complex(8) :: alpha,beta

    status=0
    esp_k=0d0
    ortho_error=huge(1d0)
    residual_error=huge(1d0)
    max_basis=size(h_diag,1)
    npadded=max_basis*dc%n_frag
    if(dc%isize_tot/=dc%n_frag*dc%isize_frag .or. &
       size(n_basis,1)/=dc%n_frag .or. size(n_basis,2)/=nspin .or. &
       size(n_mat)/=nspin .or. size(esp_k,1)<dc%nstate_tot .or. &
       size(esp_k,2)/=nspin .or. size(coef_frag,1)<max_basis .or. &
       size(coef_frag,2)<dc%nstate_tot .or. size(coef_frag,3)/=nspin .or. &
       size(h_diag,2)/=max_basis .or. size(h_diag,3)/=nspin .or. &
       size(h_halo,1)/=max_basis .or. size(h_halo,2)/=max_basis .or. &
       size(h_halo,3)/=nspin .or. size(h_halo,4)/=n_halo .or. &
       size(halo_src)/=n_halo .or. size(halo_dst)/=n_halo .or. &
       size(halo_root_src)/=n_halo .or. size(halo_dvec,1)/=3 .or. &
       size(halo_dvec,2)/=n_halo) status=1
    call sync_status(status)
    if(status/=0) return
    if(any(n_basis<0) .or. any(n_basis>max_basis) .or. &
       any(n_mat<dc%nstate_tot) .or. dc%nstate_tot<1) status=1
    call sync_status(status)
    if(status/=0) return
    if(dc%isize_tot<1 .or. dc%isize_frag<1 .or. dc%n_frag<1) then
      status=1
      return
    end if

    call initialize_grids(dc,natural,dense,status)
    if(status/=0) return
    icomm_row=comm_create_group(dc%icomm_tot,dc%id_frag,dc%i_frag)
    allocate(reqs(max(1,n_halo)),reqr(max(1,n_halo)))
    allocate(hdiag_sym(max_basis,max_basis,nspin))
    allocate(hrow(max_basis,max_basis,nspin,n_halo))
    allocate(hrecv(max_basis,max_basis,nspin,n_halo))
    call prepare_blocks
    if(status/=0) goto 900
    flag=0
    if(any(.not.ieee_is_finite(real(hdiag_sym,8))) .or. &
       any(.not.ieee_is_finite(aimag(hdiag_sym))) .or. &
       any(.not.ieee_is_finite(real(hrow,8))) .or. &
       any(.not.ieee_is_finite(aimag(hrow)))) flag=1
    call sync_status(flag)
    if(flag/=0) then
      status=1
      goto 900
    end if

    do ispin=1,nspin
      nbuffer=min(max(1,ceiling(0.05d0*real(dc%nstate_tot,8))), &
        max(0,n_mat(ispin)-dc%nstate_tot))
      nsub=dc%nstate_tot+nbuffer
      if(nsub<1 .or. nsub>n_mat(ispin)) then
        status=1
        call sync_status(status)
        exit
      end if
      call initialize_layouts(npadded,nsub,max_basis,natural,dense,ln,ld,lg,status)
      call sync_status(status)
      if(status/=0) exit
      nrow=ln%lda
      ncol=max(1,ln%ncol)
      allocate(x(nrow,ncol),x_backup(nrow,ncol),hx(nrow,ncol))
      allocate(xd(ld%lda,max(1,ld%ncol)),hxd(ld%lda,max(1,ld%ncol)))
      allocate(xrot(ld%lda,max(1,ld%ncol)))
      allocate(projected(lg%lda,max(1,lg%ncol)),z(lg%lda,max(1,lg%ncol)))
      allocate(eval(nsub))
      allocate(residuals(dc%nstate_tot))
      nlocked=0
      call initialize_subspace(ispin,nsub,ln,x,eval,status)
      if(status/=0) exit
      call rayleigh_ritz(ispin,nsub,0,ln,ld,lg,x,eval,status)
      if(status/=0) exit
      call calculate_residual(ispin,nsub,ln,x,eval,dc%nstate_tot,hx, &
        residuals,max_resid,status)
      if(status/=0) exit
      call update_locked_count(eval,residuals,nlocked)

      if(nsub<n_mat(ispin)) then
        max_local_col=max(1,(nsub+natural%npcol-1)/natural%npcol)
        if(filter_chunk_size>0) then
          chunk_width=min(max_local_col,filter_chunk_size)
        else
          chunk_width=max(1,int(workspace_limit/max(48_8*int(ln%lda,8)+ &
            16_8*int(max_basis,8)*int(max(1,n_halo),8),1_8)))
          chunk_width=min(max_local_col,chunk_width)
        end if
        allocate(q0(ln%lda,chunk_width),q1(ln%lda,chunk_width), &
                 q2(ln%lda,chunk_width))
        do cycle=1,max_cycle
          if(max_resid<=residual_tolerance) exit
          if(nlocked>=dc%nstate_tot) exit
          flag=0
          cutoff=eval(nsub)
          call estimate_upper_bound(ispin,ln,cutoff,lower,gersh_upper,upper)
          if(upper<=cutoff+epsilon(1d0)*max(1d0,abs(cutoff))) then
            status=1
            call sync_status(status)
            exit
          end if
          x_backup=x
          do filter_attempt=1,2
            center=0.5d0*(upper+cutoff)
            radius=0.5d0*(upper-cutoff)
            if(radius<=epsilon(1d0)*max(1d0,abs(center))) then
              flag=1
              call sync_status(flag)
              exit
            end if
            call active_column_range(nsub,nlocked,ln,active_first,active_count)
            max_active=min(nsub-nlocked,max_local_col)
            nchunks=(max_active+chunk_width-1)/chunk_width
            flag=0
            do j=1,nchunks
              first=active_first+(j-1)*chunk_width
              last=min(ln%ncol,first+chunk_width-1)
              nchunk=max(0,last-first+1)
              if(nchunk>0) then
                q0(:,1:nchunk)=x(:,first:last)
                call apply_hamiltonian(ispin,ln,q0(:,1:nchunk),q1(:,1:nchunk),reqs,reqr)
                q1(:,1:nchunk)=(q1(:,1:nchunk)-center*q0(:,1:nchunk))/radius
                do local_col=2,filter_degree
                  call apply_hamiltonian(ispin,ln,q1(:,1:nchunk),q2(:,1:nchunk),reqs,reqr)
                  q2(:,1:nchunk)=2d0*(q2(:,1:nchunk)-center*q1(:,1:nchunk))/radius-q0(:,1:nchunk)
                  q0(:,1:nchunk)=q1(:,1:nchunk)
                  q1(:,1:nchunk)=q2(:,1:nchunk)
                end do
                if(any(.not.ieee_is_finite(real(q1(:,1:nchunk),8))) .or. &
                   any(.not.ieee_is_finite(aimag(q1(:,1:nchunk))))) then
                  flag=1
                else
                  max_amplitude=maxval(abs(q1(:,1:nchunk)))
                  scaled_bound=max(abs((lower-center)/radius),abs((gersh_upper-center)/radius))
                  log_limit=log(filter_growth_margin)
                  if(scaled_bound>1d0) log_limit=log_limit+real(filter_degree,8)* &
                    log(scaled_bound+sqrt(scaled_bound**2-1d0))
                  log_limit=min(log(huge(1d0))-log(10d0),log_limit)
                  if(max_amplitude>0d0) then
                    if(log(max_amplitude)>log_limit) flag=1
                  end if
                end if
                if(flag==0) x(:,first:last)=q1(:,1:nchunk)
              end if
              call sync_status(flag)
              if(flag/=0) exit
            end do
            if(flag==0) exit
            if(filter_attempt==1 .and. upper<gersh_upper*(1d0-10d0*epsilon(1d0))) then
              x=x_backup
              upper=gersh_upper
              if(dc%id_tot==0) write(*,*) &
                'Complex CheFSI retries filter with Gershgorin upper bound:',upper
              cycle
            end if
            status=1
            call sync_status(status)
            exit
          end do
          if(flag/=0) then
            status=flag
            exit
          end if
          call orthonormalize(nsub,nlocked,ln,ld,lg,x,xd,projected,status)
          if(status/=0) exit
          call rayleigh_ritz(ispin,nsub,nlocked,ln,ld,lg,x,eval,status)
          if(status/=0) exit
          call calculate_residual(ispin,nsub,ln,x,eval,dc%nstate_tot,hx, &
            residuals,max_resid,status)
          if(status/=0) exit
          call update_locked_count(eval,residuals,nlocked)
          if(dc%id_tot==0) write(*,'(a,3i6,1x,es12.4)') &
            'Complex CheFSI k/spin/cycle/residual:',ik,ispin,cycle,max_resid
        end do
        deallocate(q2,q1,q0)
        if(status/=0) exit
      end if

      if(allocated(residuals)) deallocate(residuals)

      if(max_resid>residual_tolerance) then
        if(dc%id_tot==0) write(*,*) 'Complex CheFSI did not converge; residual:',max_resid
        status=1
        call sync_status(status)
        exit
      end if
      call calculate_orthogonality(nsub,ln,ld,lg,x,xd,projected, &
        dc%nstate_tot,orth_local,status)
      if(status/=0) exit
      if(orth_local>1d-10) then
        status=1
        call sync_status(status)
        exit
      end if
      esp_k(1:dc%nstate_tot,ispin)=eval(1:dc%nstate_tot)
      call export_coefficients(ispin,nsub,ln,x,dc%nstate_tot)
      ortho_error(ispin)=orth_local
      residual_error(ispin)=max_resid
      deallocate(eval,z,projected,xrot,hxd,xd,hx,x)
      deallocate(x_backup)
    end do
    call sync_status(status)
    if(status==0) call comm_bcast(esp_k,dc%icomm_tot,0)
900 continue
      if(allocated(eval)) deallocate(eval)
      if(allocated(residuals)) deallocate(residuals)
    if(allocated(z)) deallocate(z)
    if(allocated(projected)) deallocate(projected)
    if(allocated(xrot)) deallocate(xrot)
    if(allocated(hxd)) deallocate(hxd)
    if(allocated(xd)) deallocate(xd)
    if(allocated(hx)) deallocate(hx)
    if(allocated(x)) deallocate(x)
    if(allocated(x_backup)) deallocate(x_backup)
    if(allocated(hrecv)) deallocate(hrecv)
    if(allocated(hrow)) deallocate(hrow)
    if(allocated(hdiag_sym)) deallocate(hdiag_sym)
    if(allocated(reqr)) deallocate(reqr)
    if(allocated(reqs)) deallocate(reqs)
    call comm_free_group(icomm_row)
    call finalize_grids(natural,dense)

  contains

    subroutine sync_status(value)
      integer, intent(inout) :: value
      value=merge(1,0,value/=0)
      call comm_get_max(value,dc%icomm_tot)
    end subroutine sync_status

    subroutine prepare_blocks
      integer :: h,s,i,j,tag_send,tag_recv,nreq
      integer, allocatable :: send_req(:),recv_req(:)
      hdiag_sym=(0d0,0d0)
      hrow=(0d0,0d0)
      hrecv=(0d0,0d0)
      if(dc%id_frag==0) then
        do s=1,nspin
          hdiag_sym(:,:,s)=0.5d0*(h_diag(:,:,s)+conjg(transpose(h_diag(:,:,s))))
          do i=1,max_basis
            hdiag_sym(i,i,s)=cmplx(real(hdiag_sym(i,i,s),8),0d0,8)
          end do
        end do
        allocate(send_req(max(1,n_halo)),recv_req(max(1,n_halo)))
        nreq=0
        do h=1,n_halo
          tag_send=direction_tag(-halo_dvec(:,h))
          tag_recv=direction_tag( halo_dvec(:,h))
          nreq=nreq+1
          recv_req(nreq)=comm_irecv(hrecv(:,:,:,h),halo_root_src(h),tag_recv,dc%icomm_tot)
          send_req(nreq)=comm_isend(h_halo(:,:,:,h),halo_root_src(h),tag_send,dc%icomm_tot)
        end do
        if(nreq>0) then
          call comm_wait_all(recv_req(1:nreq))
          call comm_wait_all(send_req(1:nreq))
        end if
        do h=1,n_halo
          do s=1,nspin
            hrow(:,:,s,h)=0.5d0*(conjg(transpose(h_halo(:,:,s,h)))+hrecv(:,:,s,h))
          end do
        end do
        deallocate(recv_req,send_req)
      end if
      call comm_bcast(hdiag_sym,dc%icomm_frag,0)
      call comm_bcast(hrow,dc%icomm_frag,0)
    end subroutine prepare_blocks

    subroutine initialize_subspace(s,n,layout,vector,eigenvalue,istat)
      integer, intent(in) :: s,n
      type(s_layout), intent(in) :: layout
      complex(8), intent(out) :: vector(:,:)
      real(8), intent(out) :: eigenvalue(:)
      integer, intent(out) :: istat
      integer :: nb0,ios0,i0,j0,off0,k0,col0,state0,owner0
      integer, allocatable :: order(:)
      real(8), allocatable :: vals(:),vals_local(:),vals_global(:)
      complex(8), allocatable :: block(:,:),vec(:,:)
      nb0=n_basis(dc%i_frag,s)
      allocate(block(max(1,nb0),max(1,nb0)),vec(max(1,nb0),max(1,nb0)))
      allocate(vals(max(1,nb0)),vals_local(n_mat(s)),vals_global(n_mat(s)))
      block=(0d0,0d0)
      if(nb0>0) block(1:nb0,1:nb0)=hdiag_sym(1:nb0,1:nb0,s)
      vec=(0d0,0d0)
      istat=0
      if(nb0>0) call eigen_zheev(block,vals(1:nb0),vec(1:nb0,1:nb0),lapack_info=istat)
      if(istat/=0) istat=1
      call sync_status(istat)
      if(istat/=0) then
        deallocate(vals_global,vals_local,vals,vec,block)
        return
      end if
      vals_local=0d0
      off0=sum(n_basis(1:dc%i_frag-1,s))
      if(dc%id_frag==0 .and. nb0>0) vals_local(off0+1:off0+nb0)=vals(1:nb0)
      call comm_summation(vals_local,vals_global,n_mat(s),dc%icomm_tot)
      allocate(order(n_mat(s)))
      order=[(i0,i0=1,n_mat(s))]
      call sort_values(vals_global,order,1,n_mat(s))
      vector=(0d0,0d0)
      do k0=1,n
        if(order(k0)<=off0 .or. order(k0)>off0+nb0) cycle
        owner0=(k0-1)/max(1,(n+natural%npcol-1)/natural%npcol)
        if(owner0/=natural%mycol) cycle
        col0=k0-owner0*max(1,(n+natural%npcol-1)/natural%npcol)
        state0=order(k0)-off0
        vector(1:nb0,col0)=vec(1:nb0,state0)
      end do
      eigenvalue(1:n)=0d0
      istat=0
      deallocate(order,vals_global,vals_local,vals,vec,block)
    end subroutine initialize_subspace

    subroutine apply_hamiltonian(s,layout,vin,vout,send_req,recv_req)
      integer, intent(in) :: s
      type(s_layout), intent(in) :: layout
      complex(8), intent(in) :: vin(:,:)
      complex(8), intent(out) :: vout(:,:)
      integer, intent(out) :: send_req(:),recv_req(:)
      complex(8), allocatable :: xhalo(:,:,:),xsend(:,:,:)
      integer :: h,nlocal,nb0,source0,dest0,tag0
      nlocal=min(layout%ncol,size(vin,2))
      vout=(0d0,0d0)
      if(nlocal==0) return
      allocate(xhalo(max_basis,nlocal,max(1,n_halo)))
      allocate(xsend(max_basis,nlocal,max(1,n_halo)))
      xhalo=(0d0,0d0)
      xsend=(0d0,0d0)
      do h=1,n_halo
        source0=halo_src(h)-1
        dest0=halo_dst(h)-1
        tag0=direction_tag(halo_dvec(:,h))
        xsend(:,1:nlocal,h)=vin(:,1:nlocal)
        send_req(h)=comm_isend(xsend(:,:,h:h),dest0,tag0,icomm_row)
        recv_req(h)=comm_irecv(xhalo(:,:,h:h),source0,tag0,icomm_row)
      end do
      if(n_halo>0) call comm_wait_all(recv_req(1:n_halo))
      nb0=n_basis(dc%i_frag,s)
      alpha=(1d0,0d0)
      beta=(0d0,0d0)
      if(nb0>0) call zgemm('N','N',nb0,nlocal,nb0,alpha,hdiag_sym(:,:,s), &
        max_basis,vin,layout%lda,beta,vout,layout%lda)
      do h=1,n_halo
        if(n_basis(halo_src(h),s)==0 .or. nb0==0) cycle
        beta=(1d0,0d0)
        call zgemm('N','N',nb0,nlocal,n_basis(halo_src(h),s),alpha, &
          hrow(:,:,s,h),max_basis,xhalo(:,:,h),max_basis,beta,vout,layout%lda)
      end do
      if(n_halo>0) call comm_wait_all(send_req(1:n_halo))
      deallocate(xsend,xhalo)
    end subroutine apply_hamiltonian

    subroutine orthonormalize(n,nlocked,layout_n,layout_d,layout_g,vector,vector_d,gram,istat)
      integer, intent(in) :: n,nlocked
      type(s_layout), intent(in) :: layout_n,layout_d,layout_g
      complex(8), intent(inout) :: vector(:,:)
      complex(8), intent(inout) :: vector_d(:,:),gram(:,:)
      integer, intent(out) :: istat
      type(s_layout) :: layout_active,layout_small
      integer :: pass,info0,nactive
      real(8) :: one_r,zero_r
      complex(8) :: one_c,zero_c
      istat=0
      vector_d=(0d0,0d0)
      call redistribute_to_dense(layout_n,vector,layout_d,vector_d,n)
      nactive=n-nlocked
      call initialize_active_layouts(nactive,layout_d,layout_g,layout_active,layout_small,istat)
      if(istat/=0) return
      xrot=(0d0,0d0)
      call redistribute_active_to_dense(layout_n,vector,layout_active,xrot,nlocked,nactive)
      one_r=1d0
      zero_r=0d0
      one_c=(1d0,0d0)
      zero_c=(0d0,0d0)
      do pass=1,2
        if(nlocked>0) then
          gram=(0d0,0d0)
          call pzgemm('C','N',nlocked,nactive,npadded,one_c,vector_d,1,1, &
            layout_d%desc,xrot,1,1,layout_active%desc,zero_c,gram,1,1,layout_g%desc)
          call pzgemm('N','N',npadded,nactive,nlocked,-one_c,vector_d,1,1, &
            layout_d%desc,gram,1,1,layout_g%desc,one_c,xrot,1,1,layout_active%desc)
        end if
        gram=(0d0,0d0)
        call pzherk('L','C',nactive,npadded,one_r,xrot,1,1,layout_active%desc, &
          zero_r,gram,1,1,layout_small%desc)
        call pzpotrf('L',nactive,gram,1,1,layout_small%desc,info0)
        if(info0/=0) then
          istat=1
          call sync_status(istat)
          if(istat/=0) exit
        end if
        call pztrsm('R','L','C','N',npadded,nactive,one_c,gram,1,1,layout_small%desc, &
          xrot,1,1,layout_active%desc)
      end do
      if(istat/=0) then
        call distributed_qr(nactive,layout_active,xrot,istat)
      end if
      call redistribute_active_to_natural(layout_active,xrot,layout_n,vector,nlocked,nactive)
      call clear_padding(layout_n,vector)
    end subroutine orthonormalize

    subroutine distributed_qr(n,layout,vector,istat)
      integer, intent(in) :: n
      type(s_layout), intent(in) :: layout
      complex(8), intent(inout) :: vector(:,:)
      integer, intent(out) :: istat
      integer :: info0,lwork0,local_tau
      complex(8) :: query(1)
      complex(8), allocatable :: tau0(:),work0(:)
      integer, external :: numroc
      local_tau=numroc(n,block_size,dense%mycol,0,dense%npcol)
      allocate(tau0(max(1,local_tau)))
      call pzgeqrf(npadded,n,vector,1,1,layout%desc,tau0,query,-1,info0)
      istat=merge(1,0,info0/=0)
      call sync_status(istat)
      if(istat/=0) then
        deallocate(tau0)
        return
      end if
      lwork0=max(1,ceiling(real(query(1),8)))
      call comm_get_max(lwork0,dc%icomm_tot)
      allocate(work0(lwork0))
      call pzgeqrf(npadded,n,vector,1,1,layout%desc,tau0,work0,lwork0,info0)
      istat=merge(1,0,info0/=0)
      call sync_status(istat)
      if(istat==0) then
        call pzungqr(npadded,n,n,vector,1,1,layout%desc,tau0,query,-1,info0)
        istat=merge(1,0,info0/=0)
        call sync_status(istat)
      end if
      if(istat==0) then
        lwork0=max(1,ceiling(real(query(1),8)))
        call comm_get_max(lwork0,dc%icomm_tot)
        if(size(work0)<lwork0) then
          deallocate(work0)
          allocate(work0(lwork0))
        end if
        call pzungqr(npadded,n,n,vector,1,1,layout%desc,tau0,work0,lwork0,info0)
        istat=merge(1,0,info0/=0)
        call sync_status(istat)
      end if
      deallocate(work0,tau0)
    end subroutine distributed_qr

    subroutine rayleigh_ritz(s,n,nlocked,layout_n,layout_d,layout_g,vector,eigenvalue,istat)
      integer, intent(in) :: s,n,nlocked
      type(s_layout), intent(in) :: layout_n,layout_d,layout_g
      complex(8), intent(inout) :: vector(:,:)
      real(8), intent(inout) :: eigenvalue(:)
      integer, intent(out) :: istat
      type(s_layout) :: layout_active,layout_small
      complex(8) :: work_query(1),one_c,zero_c
      real(8) :: rwork_query(1)
      complex(8), allocatable :: work0(:)
      real(8), allocatable :: rwork0(:)
      integer :: lwork0,lrwork0,info0,nactive
      nactive=n-nlocked
      call initialize_active_layouts(nactive,layout_d,layout_g,layout_active,layout_small,istat)
      if(istat/=0) return
      hxd=(0d0,0d0)
      call apply_hamiltonian(s,layout_n,vector,hx,reqs,reqr)
      call redistribute_active_to_dense(layout_n,vector,layout_active,xd,nlocked,nactive)
      call redistribute_active_to_dense(layout_n,hx,layout_active,hxd,nlocked,nactive)
      projected=(0d0,0d0)
      one_c=(1d0,0d0)
      zero_c=(0d0,0d0)
      call pzgemm('C','N',nactive,nactive,npadded,one_c,xd,1,1,layout_active%desc, &
        hxd,1,1,layout_active%desc,zero_c,projected,1,1,layout_small%desc)
      call pzheev('V','L',nactive,projected,1,1,layout_small%desc, &
        eigenvalue(nlocked+1:n),z,1,1,layout_small%desc,work_query,-1,rwork_query,-1,info0)
      istat=merge(1,0,info0/=0)
      call sync_status(istat)
      if(istat/=0) return
      lwork0=max(1,ceiling(real(work_query(1),8)))
      lrwork0=max(1,ceiling(rwork_query(1)))
      call comm_get_max(lwork0,dc%icomm_tot)
      call comm_get_max(lrwork0,dc%icomm_tot)
      allocate(work0(lwork0),rwork0(lrwork0))
      call pzheev('V','L',nactive,projected,1,1,layout_small%desc, &
        eigenvalue(nlocked+1:n),z,1,1,layout_small%desc,work0,lwork0,rwork0,lrwork0,info0)
      istat=merge(1,0,info0/=0)
      if(any(.not.ieee_is_finite(eigenvalue(nlocked+1:n)))) istat=1
      call sync_status(istat)
      deallocate(work0,rwork0)
      if(istat/=0) return
      call comm_bcast(eigenvalue(nlocked+1:n),dc%icomm_tot,0)
      xrot=(0d0,0d0)
      call pzgemm('N','N',npadded,nactive,nactive,one_c,xd,1,1,layout_active%desc, &
        z,1,1,layout_small%desc,zero_c,xrot,1,1,layout_active%desc)
      call redistribute_active_to_natural(layout_active,xrot,layout_n,vector,nlocked,nactive)
      call clear_padding(layout_n,vector)
    end subroutine rayleigh_ritz

    subroutine calculate_residual(s,n,layout,vector,eigenvalue,ntarget,hvector, &
        residual,maximum,istat)
      integer, intent(in) :: s,n,ntarget
      type(s_layout), intent(in) :: layout
      complex(8), intent(in) :: vector(:,:)
      real(8), intent(in) :: eigenvalue(:)
      complex(8), intent(out) :: hvector(:,:)
      real(8), intent(out) :: residual(:)
      real(8), intent(out) :: maximum
      integer, intent(out) :: istat
      integer :: lc,gc,nb0
      real(8), allocatable :: local(:),global(:)
      call apply_hamiltonian(s,layout,vector,hvector,reqs,reqr)
      allocate(local(ntarget),global(ntarget))
      local=0d0
      nb0=n_basis(dc%i_frag,s)
      do lc=1,layout%ncol
        gc=natural_global_column(lc,n, natural%mycol,natural%npcol)
        if(gc<1 .or. gc>ntarget) cycle
        local(gc)=sum(abs(hvector(1:nb0,lc)-eigenvalue(gc)*vector(1:nb0,lc))**2)
      end do
      call comm_summation(local,global,ntarget,dc%icomm_tot)
      maximum=0d0
      istat=0
      do gc=1,ntarget
        global(gc)=sqrt(max(0d0,global(gc)))/max(1d0,abs(eigenvalue(gc)))
        if(.not.ieee_is_finite(global(gc))) istat=1
        residual(gc)=global(gc)
        maximum=max(maximum,global(gc))
      end do
      if(.not.ieee_is_finite(maximum)) istat=1
      call sync_status(istat)
      deallocate(global,local)
    end subroutine calculate_residual

    subroutine update_locked_count(eigenvalue,residual,nlocked0)
      real(8), intent(in) :: eigenvalue(:),residual(:)
      integer, intent(inout) :: nlocked0
      integer :: candidate,i0
      real(8) :: gap_scale
      candidate=nlocked0
      do i0=nlocked0+1,dc%nstate_tot
        if(residual(i0)>residual_tolerance) exit
        candidate=i0
      end do
      do while(candidate>nlocked0 .and. candidate<size(eigenvalue))
        gap_scale=max(1d0,abs(eigenvalue(candidate)),abs(eigenvalue(candidate+1)))
        if(abs(eigenvalue(candidate+1)-eigenvalue(candidate))> &
           lock_gap_relative_tolerance*gap_scale) exit
        candidate=candidate-1
      end do
      if(maxval(residual)<=lock_activation_tolerance) nlocked0=candidate
    end subroutine update_locked_count

    subroutine active_column_range(n,nlocked0,layout,first_active,count_active)
      integer, intent(in) :: n,nlocked0
      type(s_layout), intent(in) :: layout
      integer, intent(out) :: first_active,count_active
      integer :: global_first,nbcol
      nbcol=max(1,(n+natural%npcol-1)/natural%npcol)
      global_first=natural%mycol*nbcol+1
      first_active=max(1,nlocked0-global_first+2)
      first_active=min(first_active,layout%ncol+1)
      count_active=max(0,layout%ncol-first_active+1)
    end subroutine active_column_range

    subroutine initialize_active_layouts(nactive,layout_d,layout_g,layout_a,layout_small,istat)
      integer, intent(in) :: nactive
      type(s_layout), intent(in) :: layout_d,layout_g
      type(s_layout), intent(out) :: layout_a,layout_small
      integer, intent(out) :: istat
      integer :: info0
      integer, external :: numroc
      layout_a=layout_d
      layout_a%ncol=numroc(nactive,block_size,dense%mycol,0,dense%npcol)
      call descinit(layout_a%desc,npadded,nactive,block_size,block_size,0,0, &
        dense%context,layout_a%lda,info0)
      istat=merge(1,0,info0/=0)
      layout_small=layout_g
      layout_small%nrow=numroc(nactive,block_size,dense%myrow,0,dense%nprow)
      layout_small%ncol=numroc(nactive,block_size,dense%mycol,0,dense%npcol)
      layout_small%lda=max(1,layout_small%nrow)
      call descinit(layout_small%desc,nactive,nactive,block_size,block_size, &
        0,0,dense%context,layout_small%lda,info0)
      if(info0/=0) istat=1
      call sync_status(istat)
    end subroutine initialize_active_layouts

    subroutine redistribute_active_to_dense(layout_n,src,layout_a,dst,nlocked0,nactive)
      type(s_layout), intent(in) :: layout_n,layout_a
      integer, intent(in) :: nlocked0,nactive
      complex(8), intent(in) :: src(:,:)
      complex(8), intent(inout) :: dst(:,:)
      call pzgemr2d(npadded,nactive,src,1,nlocked0+1,layout_n%desc, &
        dst,1,1,layout_a%desc,natural%context)
    end subroutine redistribute_active_to_dense

    subroutine redistribute_active_to_natural(layout_a,src,layout_n,dst,nlocked0,nactive)
      type(s_layout), intent(in) :: layout_a,layout_n
      integer, intent(in) :: nlocked0,nactive
      complex(8), intent(in) :: src(:,:)
      complex(8), intent(inout) :: dst(:,:)
      call pzgemr2d(npadded,nactive,src,1,1,layout_a%desc, &
        dst,1,nlocked0+1,layout_n%desc,natural%context)
    end subroutine redistribute_active_to_natural

    subroutine calculate_orthogonality(n,layout_n,layout_d,layout_g,vector,vector_d,gram,ntarget,error,istat)
      integer, intent(in) :: n,ntarget
      type(s_layout), intent(in) :: layout_n,layout_d,layout_g
      complex(8), intent(in) :: vector(:,:)
      complex(8), intent(inout) :: vector_d(:,:),gram(:,:)
      real(8), intent(out) :: error
      integer, intent(out) :: istat
      integer :: i0,j0,gi,gj
      integer, external :: indxl2g
      complex(8) :: one_c,zero_c
      call redistribute_to_dense(layout_n,vector,layout_d,vector_d,n)
      gram=(0d0,0d0)
      one_c=(1d0,0d0)
      zero_c=(0d0,0d0)
      call pzgemm('C','N',n,n,npadded,one_c,vector_d,1,1,layout_d%desc, &
        vector_d,1,1,layout_d%desc,zero_c,gram,1,1,layout_g%desc)
      error=0d0
      do j0=1,layout_g%ncol
        gj=indxl2g(j0,block_size,dense%mycol,0,dense%npcol)
        if(gj>ntarget) cycle
        do i0=1,layout_g%nrow
          gi=indxl2g(i0,block_size,dense%myrow,0,dense%nprow)
          if(gi>ntarget) cycle
          if(gi==gj) then
            error=max(error,abs(gram(i0,j0)-(1d0,0d0)))
          else
            error=max(error,abs(gram(i0,j0)))
          end if
        end do
      end do
      max_in(1)=error
      call comm_get_max(max_in,max_out,1,dc%icomm_tot)
      error=max_out(1)
      istat=0
      if(.not.ieee_is_finite(error) .or. &
         any(.not.ieee_is_finite(real(gram,8))) .or. &
         any(.not.ieee_is_finite(aimag(gram)))) istat=1
      call sync_status(istat)
    end subroutine calculate_orthogonality

    subroutine export_coefficients(s,n,layout,vector,ntarget)
      integer, intent(in) :: s,n,ntarget
      type(s_layout), intent(in) :: layout
      complex(8), intent(in) :: vector(:,:)
      complex(8), allocatable :: local(:,:),merged(:,:)
      integer :: lc,gc,nb0
      nb0=n_basis(dc%i_frag,s)
      allocate(local(max_basis,ntarget),merged(max_basis,ntarget))
      local=(0d0,0d0)
      do lc=1,layout%ncol
        gc=natural_global_column(lc,n,natural%mycol,natural%npcol)
        if(gc<1 .or. gc>ntarget) cycle
        local(1:nb0,gc)=vector(1:nb0,lc)
      end do
      call comm_summation(local,merged,max_basis*ntarget,dc%icomm_frag)
      if(dc%id_frag==0) coef_frag(:,:,s)=merged
      deallocate(merged,local)
    end subroutine export_coefficients

    subroutine spectral_bounds(s,lo,hi)
      integer, intent(in) :: s
      real(8), intent(out) :: lo,hi
      integer :: h,i,nb0
      real(8) :: row_sum,local_lo,local_hi
      nb0=n_basis(dc%i_frag,s)
      local_lo=huge(1d0)
      local_hi=-huge(1d0)
      do i=1,nb0
        diag_value=real(hdiag_sym(i,i,s),8)
        row_sum=sum(abs(hdiag_sym(i,1:nb0,s)))-abs(hdiag_sym(i,i,s))
        do h=1,n_halo
          row_sum=row_sum+sum(abs(hrow(i,1:n_basis(halo_src(h),s),s,h)))
        end do
        local_lo=min(local_lo,diag_value-row_sum)
        local_hi=max(local_hi,diag_value+row_sum)
      end do
      call comm_get_min(local_lo,dc%icomm_tot)
      max_in(1)=local_hi
      call comm_get_max(max_in,max_out,1,dc%icomm_tot)
      local_hi=max_out(1)
      lo=local_lo
      hi=local_hi+max(1d-10,1d-8*abs(local_hi))
    end subroutine spectral_bounds

    subroutine estimate_upper_bound(s,layout,cutoff0,lo,gersh_hi,hi)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      integer, intent(in) :: s
      type(s_layout), intent(in) :: layout
      real(8), intent(in) :: cutoff0
      real(8), intent(out) :: lo,gersh_hi,hi
      real(8) :: ritz,residual,candidate
      call spectral_bounds(s,lo,gersh_hi)
      call lanczos_upper_bound(s,layout,ritz,residual)
      candidate=ritz+max(10d0*residual,0.1d0*max(1d0,ritz-lo))
      if(.not.ieee_is_finite(candidate) .or. candidate<=cutoff0+ &
         epsilon(1d0)*max(1d0,abs(cutoff0))) then
        hi=gersh_hi
      else
        hi=min(gersh_hi,candidate)
      end if
      if(dc%id_tot==0) write(*,'(a,2i6,4(1x,es12.4))') &
        'Complex CheFSI bounds k/spin/lower/Lanczos/Gershgorin/selected:', &
        ik,s,lo,ritz,gersh_hi,hi
    end subroutine estimate_upper_bound

    subroutine lanczos_upper_bound(s,layout,ritz,residual)
      use eigen_subdiag_sub, only: eigen_dsyev
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      integer, intent(in) :: s
      type(s_layout), intent(in) :: layout
      real(8), intent(out) :: ritz,residual
      type(s_layout) :: vector_layout
      complex(8), allocatable :: q(:,:),qold(:,:),hq(:,:),basis(:,:)
      complex(8) :: alpha_c,overlap_c
      real(8), allocatable :: alpha(:),beta_value(:),tri(:,:),eval_small(:),vec_small(:,:)
      real(8) :: local_norm,global_norm,imag_error,send_values(2),recv_values(2)
      integer :: nb0,offset0,k0,j0,pass0,nsteps0,istat0,bad_hermitian
      nb0=n_basis(dc%i_frag,s)
      offset0=sum(n_basis(1:dc%i_frag-1,s))
      vector_layout=layout
      vector_layout%ncol=1
      allocate(q(layout%lda,1),qold(layout%lda,1),hq(layout%lda,1))
      allocate(basis(layout%lda,lanczos_max_steps))
      allocate(alpha(lanczos_max_steps),beta_value(lanczos_max_steps))
      q=(0d0,0d0)
      qold=(0d0,0d0)
      basis=(0d0,0d0)
      do j0=1,nb0
        local_value=real(offset0+j0,8)
        q(j0,1)=cmplx(sin(0.713d0*local_value)+0.37d0*cos(1.117d0*local_value), &
          cos(0.419d0*local_value)-0.23d0*sin(1.371d0*local_value),8)
      end do
      local_norm=sum(abs(q(1:nb0,1))**2)
      call comm_summation(local_norm,global_norm,icomm_row)
      if(global_norm<=tiny(1d0)) then
        q=(0d0,0d0)
        if(offset0==0 .and. nb0>0) q(1,1)=(1d0,0d0)
        global_norm=1d0
      end if
      q=q/sqrt(global_norm)
      alpha=0d0
      beta_value=0d0
      bad_hermitian=0
      nsteps0=lanczos_max_steps
      do k0=1,lanczos_max_steps
        call apply_hamiltonian(s,vector_layout,q,hq,reqs,reqr)
        if(k0>1) hq=hq-beta_value(k0-1)*qold
        alpha_c=dot_product(q(1:nb0,1),hq(1:nb0,1))
        call comm_summation(alpha_c,overlap_c,icomm_row)
        imag_error=abs(aimag(overlap_c))
        if(imag_error>1d-10*max(1d0,abs(real(overlap_c,8)))) then
          bad_hermitian=1
        end if
        alpha_c=cmplx(real(overlap_c,8),0d0,8)
        alpha(k0)=real(alpha_c,8)
        hq=hq-alpha(k0)*q
        basis(:,k0)=q(:,1)
        do pass0=1,2
          do j0=1,k0
            overlap_c=dot_product(basis(1:nb0,j0),hq(1:nb0,1))
            call comm_summation(overlap_c,alpha_c,icomm_row)
            hq(1:nb0,1)=hq(1:nb0,1)-alpha_c*basis(1:nb0,j0)
          end do
        end do
        local_norm=sum(abs(hq(1:nb0,1))**2)
        call comm_summation(local_norm,global_norm,icomm_row)
        beta_value(k0)=sqrt(max(0d0,global_norm))
        if(k0==lanczos_max_steps .or. beta_value(k0)<=100d0*epsilon(1d0)) then
          nsteps0=k0
          exit
        end if
        qold=q
        q(:,1)=hq(:,1)/beta_value(k0)
      end do
      allocate(tri(nsteps0,nsteps0),eval_small(nsteps0),vec_small(nsteps0,nsteps0))
      tri=0d0
      do k0=1,nsteps0
        tri(k0,k0)=alpha(k0)
      end do
      do k0=1,nsteps0-1
        tri(k0,k0+1)=beta_value(k0)
        tri(k0+1,k0)=beta_value(k0)
      end do
      call eigen_dsyev(tri,eval_small,vec_small)
      ritz=eval_small(nsteps0)
      residual=beta_value(nsteps0)*abs(vec_small(nsteps0,nsteps0))
      istat0=bad_hermitian
      if(.not.ieee_is_finite(ritz) .or. .not.ieee_is_finite(residual)) istat0=1
      call sync_status(istat0)
      if(istat0/=0) then
        ritz=-huge(1d0)
        residual=huge(1d0)
      end if
      send_values=[ritz,residual]
      call comm_get_max(send_values,recv_values,2,dc%icomm_tot)
      ritz=recv_values(1)
      residual=recv_values(2)
      deallocate(vec_small,eval_small,tri,beta_value,alpha,basis,hq,qold,q)
    end subroutine lanczos_upper_bound

    subroutine clear_padding(layout,vector)
      type(s_layout), intent(in) :: layout
      complex(8), intent(inout) :: vector(:,:)
      integer :: nb0
      nb0=n_basis(dc%i_frag,ispin)
      if(nb0<max_basis .and. layout%ncol>0) vector(nb0+1:max_basis,1:layout%ncol)=(0d0,0d0)
    end subroutine clear_padding

    subroutine redistribute_to_dense(layout_src,src,layout_dst,dst,n)
      type(s_layout), intent(in) :: layout_src,layout_dst
      integer, intent(in) :: n
      complex(8), intent(in) :: src(:,:)
      complex(8), intent(inout) :: dst(:,:)
      call pzgemr2d(npadded,n,src,1,1,layout_src%desc,dst,1,1,layout_dst%desc,natural%context)
    end subroutine redistribute_to_dense

    subroutine redistribute_to_natural(layout_src,src,layout_dst,dst,n)
      type(s_layout), intent(in) :: layout_src,layout_dst
      integer, intent(in) :: n
      complex(8), intent(in) :: src(:,:)
      complex(8), intent(inout) :: dst(:,:)
      call pzgemr2d(npadded,n,src,1,1,layout_src%desc,dst,1,1,layout_dst%desc,natural%context)
    end subroutine redistribute_to_natural

  end subroutine diag_chefsi_complex

  subroutine initialize_grids(dc,natural,dense,status)
    use communication, only: comm_summation
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    type(s_grid), intent(out) :: natural,dense
    integer, intent(out) :: status
    integer :: iam,nprocs,i,j
    integer, allocatable :: map_local(:,:),map_global(:,:),rank_local(:),rank_global(:),dmap(:,:)
    call blacs_pinfo(iam,nprocs)
    status=0
    if(nprocs<dc%isize_tot) then
      status=1
      return
    end if
    natural%nprow=dc%n_frag
    natural%npcol=dc%isize_frag
    allocate(map_local(natural%nprow,natural%npcol),map_global(natural%nprow,natural%npcol))
    map_local=0
    map_local(dc%i_frag,dc%id_frag+1)=iam+1
    call comm_summation(map_local,map_global,size(map_local),dc%icomm_tot)
    map_global=map_global-1
    call blacs_get(0,0,natural%context)
    call blacs_gridmap(natural%context,map_global,natural%nprow,natural%nprow,natural%npcol)
    call blacs_gridinfo(natural%context,natural%nprow,natural%npcol,natural%myrow,natural%mycol)
    if(natural%myrow/=dc%i_frag-1 .or. natural%mycol/=dc%id_frag) status=1
    deallocate(map_global,map_local)
    dense%nprow=int(sqrt(real(dc%isize_tot,8)))
    do while(dense%nprow>1)
      if(mod(dc%isize_tot,dense%nprow)==0) exit
      dense%nprow=dense%nprow-1
    end do
    dense%npcol=dc%isize_tot/dense%nprow
    allocate(rank_local(dc%isize_tot),rank_global(dc%isize_tot))
    rank_local=0
    rank_local(dc%id_tot+1)=iam+1
    call comm_summation(rank_local,rank_global,dc%isize_tot,dc%icomm_tot)
    rank_global=rank_global-1
    allocate(dmap(dense%nprow,dense%npcol))
    do j=1,dense%npcol
      do i=1,dense%nprow
        dmap(i,j)=rank_global((j-1)*dense%nprow+i)
      end do
    end do
    call blacs_get(0,0,dense%context)
    call blacs_gridmap(dense%context,dmap,dense%nprow,dense%nprow,dense%npcol)
    call blacs_gridinfo(dense%context,dense%nprow,dense%npcol,dense%myrow,dense%mycol)
    deallocate(dmap,rank_global,rank_local)
    call comm_get_dummy_status
  contains
    subroutine comm_get_dummy_status
      use communication, only: comm_get_max
      call comm_get_max(status,dc%icomm_tot)
    end subroutine comm_get_dummy_status
  end subroutine initialize_grids

  subroutine finalize_grids(natural,dense)
    implicit none
    type(s_grid), intent(inout) :: natural,dense
    if(dense%context>=0) call blacs_gridexit(dense%context)
    if(natural%context>=0) call blacs_gridexit(natural%context)
    dense%context=-1
    natural%context=-1
  end subroutine finalize_grids

  subroutine initialize_layouts(npadded,nsub,max_basis,natural,dense,ln,ld,lg,status)
    implicit none
    integer, intent(in) :: npadded,nsub,max_basis
    type(s_grid), intent(in) :: natural,dense
    type(s_layout), intent(out) :: ln,ld,lg
    integer, intent(out) :: status
    integer :: info,nbcol
    integer, external :: numroc
    nbcol=max(1,(nsub+natural%npcol-1)/natural%npcol)
    ln%nrow=numroc(npadded,max_basis,natural%myrow,0,natural%nprow)
    ln%ncol=numroc(nsub,nbcol,natural%mycol,0,natural%npcol)
    ln%lda=max(1,ln%nrow)
    call descinit(ln%desc,npadded,nsub,max_basis,nbcol,0,0,natural%context,ln%lda,info)
    status=merge(1,0,info/=0)
    ld%nrow=numroc(npadded,block_size,dense%myrow,0,dense%nprow)
    ld%ncol=numroc(nsub,block_size,dense%mycol,0,dense%npcol)
    ld%lda=max(1,ld%nrow)
    call descinit(ld%desc,npadded,nsub,block_size,block_size,0,0,dense%context,ld%lda,info)
    if(info/=0) status=1
    lg%nrow=numroc(nsub,block_size,dense%myrow,0,dense%nprow)
    lg%ncol=numroc(nsub,block_size,dense%mycol,0,dense%npcol)
    lg%lda=max(1,lg%nrow)
    call descinit(lg%desc,nsub,nsub,block_size,block_size,0,0,dense%context,lg%lda,info)
    if(info/=0) status=1
  end subroutine initialize_layouts

  integer function natural_global_column(local_column,nsub,mycol,npcol)
    implicit none
    integer, intent(in) :: local_column,nsub,mycol,npcol
    natural_global_column=mycol*max(1,(nsub+npcol-1)/npcol)+local_column
  end function natural_global_column

  integer function direction_tag(dvec)
    implicit none
    integer, intent(in) :: dvec(3)
    direction_tag=1+(dvec(1)+1)*9+(dvec(2)+1)*3+(dvec(3)+1)
  end function direction_tag

  recursive subroutine sort_values(value,index,left,right)
    implicit none
    real(8), intent(inout) :: value(:)
    integer, intent(inout) :: index(:)
    integer, intent(in) :: left,right
    integer :: i,j,itmp
    real(8) :: pivot,vtmp
    if(left>=right) return
    pivot=value((left+right)/2)
    i=left
    j=right
    do
      do while(i<=right)
        if(value(i)>=pivot) exit
        i=i+1
      end do
      do while(j>=left)
        if(value(j)<=pivot) exit
        j=j-1
      end do
      if(i>j) exit
      vtmp=value(i); value(i)=value(j); value(j)=vtmp
      itmp=index(i); index(i)=index(j); index(j)=itmp
      i=i+1; j=j-1
    end do
    if(left<j) call sort_values(value,index,left,j)
    if(i<right) call sort_values(value,index,i,right)
  end subroutine sort_values

end module lcfo_diag_chefsi_complex
