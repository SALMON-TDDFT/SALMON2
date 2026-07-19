program test_dg_wpw_matrix_free_scf_mpi
  use mpi
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_wpw_matrix_free_scf,only:s_dg_wpw_matrix_free_scf_result,run_dg_wpw_matrix_free_scf,&
    run_dg_wpw_matrix_free_algebra_step,initialize_dg_wpw_projected_occupied,&
    complete_dg_wpw_projected_subspace,solve_dg_wpw_metric_projection,&
    initialize_dg_wpw_metric_projected_occupied
  implicit none
  type::s_fixture_context
    integer::rank=0,nrank=1,first=1,nlocal=0
    logical::bad_fixed_energy=.false.,wide_spectrum=.false.,metric_coupled=.false.
  end type
  type(s_fixture_context)::ctx
  type(s_dg_wpw_matrix_free_scf_result)::result
  integer::ierr,info,nlocal,i,j
  real(8)::density(1)=[0d0],eval(2),occ(1)=[2d0]
  complex(8),allocatable::qw(:,:),qp(:,:)
  complex(8),allocatable::qold(:,:)
  complex(8),allocatable::projected_ref(:,:),projected_rot(:,:),srot(:,:)
  complex(8)::rotation(2,2),overlap_occ(2,2)
  complex(8)::metric_w(1,2),metric_p(1,2),metric_bw(1,2),metric_bp(1,2),&
    metric_cw(1,2),metric_cp(1,2)
  real(8)::metric_residual
  real(8)::projection_orth,projector_defect
  integer::projection_rank
  real(8)::gap,residual,orth,projector
  real(8),parameter::seed1(3)=[1d0,0.3d0,0.1d0],seed2(3)=[0.2d0,1d0,0.4d0]
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,ctx%rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,ctx%nrank,ierr)
  if(ctx%nrank/=2)call MPI_Abort(MPI_COMM_WORLD,2,ierr)
  ctx%metric_coupled=.true.
  metric_w=(0d0,0d0);metric_p=(0d0,0d0)
  metric_w(1,ctx%rank+1)=(1d0,0d0)
  metric_p(1,3-(ctx%rank+1))=(0.5d0,0.25d0)
  call apply_s(ctx,metric_w,metric_p,metric_bw,metric_bp,info)
  if(info/=0)call MPI_Abort(MPI_COMM_WORLD,201,ierr)
  metric_cw=(0d0,0d0);metric_cp=(0d0,0d0)
  call solve_dg_wpw_metric_projection(ctx,MPI_COMM_WORLD,apply_s,global_gram,1,1,2,&
    1d-12,64,metric_bw,metric_bp,metric_cw,metric_cp,metric_residual,info)
  if(info/=0.or.metric_residual>1d-10.or.maxval(abs(metric_cw-metric_w))>1d-9.or.&
    maxval(abs(metric_cp-metric_p))>1d-9)call MPI_Abort(MPI_COMM_WORLD,202,ierr)
  metric_bw(:,2)=metric_bw(:,1);metric_bp(:,2)=metric_bp(:,1)
  call initialize_dg_wpw_metric_projected_occupied(ctx,MPI_COMM_WORLD,apply_s,global_gram,1,1,2,&
    1d-12,64,metric_bw,metric_bp,metric_cw,metric_cp,metric_residual,projection_rank,&
    projection_orth,info)
  if(info==0.or.projection_rank>=2)call MPI_Abort(MPI_COMM_WORLD,203,ierr)
  call apply_s(ctx,metric_w,metric_p,metric_bw,metric_bp,info)
  if(ctx%rank==1)metric_bp(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call initialize_dg_wpw_metric_projected_occupied(ctx,MPI_COMM_WORLD,apply_s,global_gram,1,1,2,&
    1d-12,64,metric_bw,metric_bp,metric_cw,metric_cp,metric_residual,projection_rank,&
    projection_orth,info)
  if(info==0)call MPI_Abort(MPI_COMM_WORLD,204,ierr)
  ctx%metric_coupled=.false.
  if(ctx%rank==0)then;ctx%first=1;ctx%nlocal=2;else;ctx%first=3;ctx%nlocal=1;endif
  nlocal=ctx%nlocal;allocate(qw(0,2),qp(nlocal,2));qp=(0d0,0d0)
  do i=1,nlocal
    j=ctx%first+i-1
    qp(i,1)=cmplx(seed1(j),0d0,8)
    qp(i,2)=cmplx(seed2(j),0d0,8)
  enddo
  allocate(projected_ref(nlocal,2),projected_rot(nlocal,2),srot(nlocal,2))
  projected_ref=qp
  call initialize_dg_wpw_projected_occupied(ctx,MPI_COMM_WORLD,apply_s,global_gram,0,nlocal,&
    2,1d-12,qw,projected_ref,projection_rank,projection_orth,info)
  if(info/=0.or.projection_rank/=2.or.projection_orth>1d-10)then
    write(*,'(a,i0,2(a,i0),a,es12.4)')'projection init rank=',ctx%rank,' info=',info,&
      ' effective_rank=',projection_rank,' orth=',projection_orth
    call MPI_Abort(MPI_COMM_WORLD,3,ierr)
  endif
  rotation=reshape([cmplx(cos(0.37d0),0d0,8),cmplx(-sin(0.37d0),0d0,8),&
    cmplx(sin(0.37d0),0d0,8),cmplx(cos(0.37d0),0d0,8)],[2,2])
  projected_rot=matmul(qp,rotation)
  call initialize_dg_wpw_projected_occupied(ctx,MPI_COMM_WORLD,apply_s,global_gram,0,nlocal,&
    2,1d-12,qw,projected_rot,projection_rank,projection_orth,info)
  if(info/=0.or.projection_rank/=2.or.projection_orth>1d-10)call MPI_Abort(MPI_COMM_WORLD,4,ierr)
  call apply_s(ctx,qw,projected_rot,qw,srot,info);if(info/=0)call MPI_Abort(MPI_COMM_WORLD,5,ierr)
  call global_gram(projected_ref,srot,nlocal,2,2,overlap_occ,info)
  projector_defect=sqrt(max(0d0,1d0-sum(abs(overlap_occ)**2)/2d0))
  if(info/=0.or.projector_defect>1d-10)call MPI_Abort(MPI_COMM_WORLD,6,ierr)
  projected_rot=qp
  call initialize_dg_wpw_projected_occupied(ctx,MPI_COMM_WORLD,apply_s,global_gram,0,nlocal,&
    1,1d-12,qw(:,1:1),projected_rot(:,1:1),projection_rank,projection_orth,info)
  if(info/=0)call MPI_Abort(MPI_COMM_WORLD,7,ierr)
  call complete_dg_wpw_projected_subspace(ctx,MPI_COMM_WORLD,apply_s,global_gram,0,nlocal,&
    1,2,1d-12,qw,projected_rot,projection_rank,projection_orth,info)
  if(info/=0.or.projection_rank/=2.or.projection_orth>1d-10)call MPI_Abort(MPI_COMM_WORLD,8,ierr)
  call apply_s(ctx,qw,projected_rot,qw,srot,info);if(info/=0)call MPI_Abort(MPI_COMM_WORLD,9,ierr)
  call global_gram(projected_rot,srot,nlocal,2,2,overlap_occ,info)
  overlap_occ(1,1)=overlap_occ(1,1)-1;overlap_occ(2,2)=overlap_occ(2,2)-1
  if(info/=0.or.maxval(abs(overlap_occ))>1d-10)call MPI_Abort(MPI_COMM_WORLD,91,ierr)
  deallocate(projected_ref,projected_rot,srot)
  call run_dg_wpw_matrix_free_scf(ctx,MPI_COMM_WORLD,apply_h,apply_s,global_gram,potential_map,0,nlocal,1,1,occ,&
    density,0.5d0,30,1d-12,0.5d0,1d-8,qw,qp,eval,result,info)
  if(info/=0.or..not.result%converged)call MPI_Abort(MPI_COMM_WORLD,10+info,ierr)
  if(result%retained_dimension/=2.or.result%gap<0.5d0)then
    call MPI_Abort(MPI_COMM_WORLD,50,ierr)
  endif
  if(max(result%generalized_residual,result%metric_orthonormality)>1d-8)call MPI_Abort(MPI_COMM_WORLD,51,ierr)
  if(max(result%density_residual,result%potential_residual,result%energy_residual,&
    result%projector_residual,result%fixed_point_residual)>1d-8)call MPI_Abort(MPI_COMM_WORLD,52,ierr)
  if(abs(result%charge_error)>1d-12)call MPI_Abort(MPI_COMM_WORLD,53,ierr)
  ctx%bad_fixed_energy=.true.;density=0d0;qp=(0d0,0d0)
  do i=1,nlocal
    j=ctx%first+i-1;qp(i,1)=cmplx(seed1(j),0d0,8);qp(i,2)=cmplx(seed2(j),0d0,8)
  enddo
  call run_dg_wpw_matrix_free_scf(ctx,MPI_COMM_WORLD,apply_h,apply_s,global_gram,potential_map,0,nlocal,1,1,occ,&
    density,0.5d0,30,1d-12,0.5d0,1d-8,qw,qp,eval,result,info)
  if(info==0.or.result%converged)call MPI_Abort(MPI_COMM_WORLD,531,ierr)
  ctx%bad_fixed_energy=.false.
  call run_dg_wpw_matrix_free_scf(ctx,MPI_COMM_WORLD,apply_h,apply_s,global_gram,potential_map,0,nlocal,1,1,occ,&
    density,0.5d0,merge(30,0,ctx%rank==0),1d-12,0.5d0,1d-8,qw,qp,eval,result,info)
  if(info==0)call MPI_Abort(MPI_COMM_WORLD,532,ierr)
  if(ctx%rank==0)write(*,'(a)')'PASS bounded two-rank matrix-free DG-DC fixed point'
  qp=(0d0,0d0)
  do i=1,nlocal
    j=ctx%first+i-1;qp(i,1)=cmplx(seed1(j),0d0,8);qp(i,2)=cmplx(seed2(j),0d0,8)
  enddo
  allocate(qold(nlocal,1));qold=0
  call run_dg_wpw_matrix_free_algebra_step(ctx,MPI_COMM_WORLD,apply_h,apply_s,global_gram,0,nlocal,1,2,1,&
    1d-12,1d-8,qw,qp,qold,eval,gap,residual,orth,projector,info)
  if(info/=0.or.maxval(abs(eval-[0d0,1d0]))>1d-8.or.abs(gap-1d0)>1d-8.or.&
    max(residual,orth)>1d-8)call MPI_Abort(MPI_COMM_WORLD,54,ierr)
  if(ctx%rank==1)qp(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
  call run_dg_wpw_matrix_free_algebra_step(ctx,MPI_COMM_WORLD,apply_h,apply_s,global_gram,0,nlocal,1,2,2,&
    1d-12,1d-8,qw,qp,qold,eval,gap,residual,orth,projector,info)
  if(info==0)call MPI_Abort(MPI_COMM_WORLD,55,ierr)
  deallocate(qw,qp,qold);ctx%wide_spectrum=.true.
  if(ctx%rank==0)then;ctx%first=1;ctx%nlocal=3;else;ctx%first=4;ctx%nlocal=3;endif
  nlocal=ctx%nlocal;allocate(qw(0,2),qp(nlocal,2),qold(nlocal,1));qold=0
  do i=1,nlocal
    j=ctx%first+i-1
    qp(i,1)=cmplx(sin(dble(17*j))+0.1d0*cos(dble(5*j)),0d0,8)
    qp(i,2)=cmplx(cos(dble(13*j))-0.2d0*sin(dble(3*j)),0d0,8)
  enddo
  call run_dg_wpw_matrix_free_algebra_step(ctx,MPI_COMM_WORLD,apply_h,apply_s,global_gram,0,nlocal,1,2,1,&
    1d-12,1d-8,qw,qp,qold,eval,gap,residual,orth,projector,info)
  if(info/=0.or.maxval(abs(eval-[0d0,1d0]))>1d-8.or.abs(gap-1d0)>1d-8.or.&
    max(residual,orth)>1d-8)call MPI_Abort(MPI_COMM_WORLD,56,ierr)
  call MPI_Finalize(ierr)
contains
  subroutine apply_h(context,xw,xp,yw,yp,apply_info)
    class(*),intent(inout)::context;complex(8),intent(in)::xw(:,:),xp(:,:);complex(8),intent(out)::yw(:,:),yp(:,:)
    integer,intent(out)::apply_info;integer::k
    yw=(0d0,0d0);yp=(0d0,0d0);apply_info=0
    select type(c=>context);type is(s_fixture_context)
      if(c%wide_spectrum)then
        do k=1,c%nlocal
          select case(c%first+k-1)
          case(1);yp(k,:)=0d0*xp(k,:)
          case(2);yp(k,:)=xp(k,:)
          case(3);yp(k,:)=1d2*xp(k,:)
          case(4);yp(k,:)=1d3*xp(k,:)
          case(5);yp(k,:)=1d4*xp(k,:)
          case(6);yp(k,:)=1d5*xp(k,:)
          end select
        enddo
      else
        do k=1,c%nlocal;yp(k,:)=dble(c%first+k-2)*dble(c%first+k)*xp(k,:);enddo
      endif
    class default;apply_info=1;end select
  end subroutine
  subroutine apply_s(context,xw,xp,yw,yp,apply_info)
    class(*),intent(inout)::context;complex(8),intent(in)::xw(:,:),xp(:,:);complex(8),intent(out)::yw(:,:),yp(:,:)
    integer,intent(out)::apply_info;integer::k
    yw=xw;apply_info=0
    select type(c=>context);type is(s_fixture_context)
      if(c%metric_coupled)then
        call apply_coupled_metric(c,xw,xp,yw,yp,apply_info)
      else if(c%wide_spectrum)then
        yp=xp
      else
        do k=1,c%nlocal;yp(k,:)=dble(c%first+k)*xp(k,:);enddo
      endif
    class default;apply_info=1;end select
  end subroutine
  subroutine apply_coupled_metric(c,xw,xp,yw,yp,apply_info)
    type(s_fixture_context),intent(in)::c
    complex(8),intent(in)::xw(:,:),xp(:,:)
    complex(8),intent(out)::yw(:,:),yp(:,:)
    integer,intent(out)::apply_info
    complex(8)::xlocal(4,size(xw,2)),xglobal(4,size(xw,2)),smat(4,4)
    integer::mpi_info
    xlocal=(0d0,0d0);xlocal(c%rank+1,:)=xw(1,:);xlocal(c%rank+3,:)=xp(1,:)
    call MPI_Allreduce(xlocal,xglobal,4*size(xw,2),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mpi_info)
    smat=reshape([(2d0,0d0),(0.2d0,0d0),(0.3d0,0.1d0),(0.1d0,0d0),&
      (0.2d0,0d0),(2.5d0,0d0),(0.05d0,0d0),(0.25d0,-0.1d0),&
      (0.3d0,-0.1d0),(0.05d0,0d0),(1.8d0,0d0),(0.15d0,0d0),&
      (0.1d0,0d0),(0.25d0,0.1d0),(0.15d0,0d0),(2.2d0,0d0)],[4,4])
    xlocal=matmul(smat,xglobal)
    yw(1,:)=xlocal(c%rank+1,:);yp(1,:)=xlocal(c%rank+3,:)
    apply_info=merge(0,1,mpi_info==MPI_SUCCESS)
  end subroutine apply_coupled_metric
  subroutine global_gram(x,y,nrow,nx,ny,g,gram_info)
    integer,intent(in)::nrow,nx,ny;complex(8),intent(in)::x(nrow,nx),y(nrow,ny);complex(8),intent(out)::g(nx,ny)
    integer,intent(out)::gram_info;complex(8)::local(nx,ny);integer::mpi_info
    local=matmul(conjg(transpose(x)),y);call MPI_Allreduce(local,g,nx*ny,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mpi_info)
    gram_info=merge(0,1,mpi_info==MPI_SUCCESS)
  end subroutine
  subroutine potential_map(context,cw,cp,f,nocc,density_in,mix,density_out,pot_res,energy,charge,map_info)
    class(*),intent(inout)::context;integer,intent(in)::nocc;complex(8),intent(in)::cw(:,:),cp(:,:)
    real(8),intent(in)::f(nocc),density_in(:),mix;real(8),intent(out)::density_out(:),pot_res,energy,charge
    integer,intent(out)::map_info;real(8)::local_e;integer::k,mpi_info
    density_out=1d0;pot_res=0d0;charge=0d0;local_e=0d0;map_info=0
    select type(c=>context);type is(s_fixture_context)
      do k=1,c%nlocal
        local_e=local_e+f(1)*dble(c%first+k-2)*dble(c%first+k)*abs(cp(k,1))**2
      enddo
    class default;map_info=1;return;end select
    call MPI_Allreduce(local_e,energy,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpi_info)
    select type(c=>context);type is(s_fixture_context)
      if(c%bad_fixed_energy.and.mix==1d0)energy=energy+1d-3
    end select
    if(mpi_info/=MPI_SUCCESS)map_info=2
  end subroutine
end program
