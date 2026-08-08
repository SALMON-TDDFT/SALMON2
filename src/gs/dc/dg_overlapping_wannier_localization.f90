#include "config.h"
module dg_overlapping_wannier_localization
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use,intrinsic::iso_fortran_env,only:real64
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::evaluate_dg_periodic_localization
  public::optimize_dg_wannier_pair
  public::build_dg_overlapping_pair_graph
  public::localize_dg_overlapping_wannier_basis
  public::validate_dg_global_covariant_gauge
contains
  subroutine validate_dg_global_covariant_gauge(transform,representation,tolerance,ok,message)
    complex(real64),intent(in)::transform(:,:),representation(:,:,:)
    real(real64),intent(in)::tolerance
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::work(:,:),identity(:,:)
    real(real64)::defect
    integer::i,operation,n

    ok=.false.;message='';n=size(transform,1)
    if(n<1.or.size(transform,2)/=n.or.size(representation,1)/=n.or.&
        size(representation,2)/=n.or.size(representation,3)<1.or.tolerance<=0d0.or.&
        .not.all(ieee_is_finite(real(transform))).or.&
        .not.all(ieee_is_finite(aimag(transform))).or.&
        .not.all(ieee_is_finite(real(representation))).or.&
        .not.all(ieee_is_finite(aimag(representation))))then
      message='invalid global-covariant gauge contract';return
    end if
    allocate(work(n,n),identity(n,n));identity=(0d0,0d0)
    do i=1,n;identity(i,i)=1d0;end do
    work=matmul(conjg(transpose(transform)),transform)-identity
    defect=maxval(abs(work))
    if(defect>tolerance)then
      message='global-covariant gauge is not unitary';return
    end if
    do operation=1,size(representation,3)
      work=matmul(transform,representation(:,:,operation))-&
        matmul(representation(:,:,operation),transform)
      defect=maxval(abs(work))
      if(defect>tolerance)then
        message='global-covariant gauge breaks full-system symmetry';return
      end if
    end do
    ok=.true.
  end subroutine validate_dg_global_covariant_gauge

  subroutine localize_dg_overlapping_wannier_basis(comm,values,gradients,weights,phases,&
      representation,product_table,support_tolerance,spread_tolerance,gradient_tolerance,&
      symmetry_tolerance,maximum_iterations,initial_spread,final_spread,maximum_pair_gradient,&
      iterations,converged,total_transform,ok,message,spread_evaluations)
    integer,intent(in)::comm,product_table(:,:),maximum_iterations
    complex(real64),intent(inout)::values(:,:),gradients(:,:,:)
    real(real64),intent(in)::weights(:),support_tolerance,spread_tolerance,&
      gradient_tolerance,symmetry_tolerance
    complex(real64),intent(in)::phases(:,:),representation(:,:,:)
    real(real64),intent(out)::initial_spread,final_spread,maximum_pair_gradient
    integer,intent(out)::iterations
    logical,intent(out)::converged,ok
    complex(real64),allocatable,intent(out)::total_transform(:,:)
    character(*),intent(out)::message
    integer,optional,intent(out)::spread_evaluations
    integer,allocatable::pair_first(:),pair_second(:)
    real(real64),allocatable::pair_support(:)
    real(real64),allocatable::graph_gradient_real(:),graph_gradient_imag(:)
    complex(real64),allocatable::raw_generator(:,:),sweep_generator(:,:),scaled_generator(:,:),&
      projection_work(:,:),search_generator(:,:),previous_gradient(:,:),previous_direction(:,:),&
      transported_gradient(:,:),transported_direction(:,:),transport_rotation(:,:),&
      block_rotation(:,:),backup_values(:,:),backup_gradients(:,:,:)
    real(real64)::gradient_real,gradient_imag,current_gradient,theta,phi,trial_spread,&
      antihermiticity_defect,commutator_defect,generator_scale,descent_measure,beta,numerator,denominator
    real(real64)::best_rejected_spread,representation_unitarity_defect,representation_closure_defect
    complex(real64)::amplitude
    integer::nwannier,first,second,left,right,product,line_search,attempt,maximum_attempts,axis,i,edge,&
      operation,rank,ierr,evaluation_count
    logical::step_ok,line_accepted,have_previous,used_conjugate
    character(256)::detail

    ok=.false.;converged=.false.;message='';iterations=0;evaluation_count=0
    if(present(spread_evaluations))spread_evaluations=0
    initial_spread=huge(1d0);final_spread=huge(1d0);maximum_pair_gradient=huge(1d0)
    nwannier=size(values,1)
    if(nwannier<2.or.size(values,2)<1.or.any(shape(gradients)/=[3,nwannier,size(values,2)]).or.&
        size(weights)/=size(values,2).or.size(phases,2)/=size(values,2).or.&
        size(representation,1)/=nwannier.or.size(representation,2)/=nwannier.or.&
        size(representation,3)<1.or.any(shape(product_table)/=&
        [size(representation,3),size(representation,3)]).or.maximum_iterations<1.or.&
        support_tolerance<0d0.or.support_tolerance>1d0.or.spread_tolerance<0d0.or.&
        gradient_tolerance<=0d0.or.symmetry_tolerance<=0d0.or.&
        .not.all(ieee_is_finite(real(representation))).or.&
        .not.all(ieee_is_finite(aimag(representation))))then
      message='invalid symmetry-constrained localization sweep contract';return
    end if
    if(any(product_table<1).or.any(product_table>size(representation,3)))then
      message='localization group product table is invalid';return
    end if
    allocate(raw_generator(nwannier,nwannier),sweep_generator(nwannier,nwannier),&
      scaled_generator(nwannier,nwannier),projection_work(nwannier,nwannier),&
      search_generator(nwannier,nwannier),previous_gradient(nwannier,nwannier),&
      previous_direction(nwannier,nwannier),transported_gradient(nwannier,nwannier),&
      transported_direction(nwannier,nwannier),transport_rotation(nwannier,nwannier))
    raw_generator=(0d0,0d0);do i=1,nwannier;raw_generator(i,i)=1d0;end do
    representation_unitarity_defect=0d0
    do operation=1,size(representation,3)
      projection_work=matmul(conjg(transpose(representation(:,:,operation))),&
        representation(:,:,operation))-raw_generator
      representation_unitarity_defect=max(representation_unitarity_defect,maxval(abs(projection_work)))
    end do
    representation_closure_defect=0d0
    do left=1,size(representation,3);do right=1,size(representation,3)
      product=product_table(left,right)
      projection_work=matmul(representation(:,:,left),representation(:,:,right))-&
        representation(:,:,product)
      representation_closure_defect=max(representation_closure_defect,maxval(abs(projection_work)))
    end do;end do
    if(representation_unitarity_defect>symmetry_tolerance.or.&
        representation_closure_defect>symmetry_tolerance)then
      message='localization group representation is not exact';return
    end if
    call build_dg_overlapping_pair_graph(comm,values,weights,support_tolerance,&
      pair_first,pair_second,pair_support,step_ok,detail)
    if(.not.step_ok)then;message=trim(detail);return;end if
#ifdef USE_MPI
    call MPI_Comm_rank(comm,rank,ierr)
#else
    rank=0
#endif
    if(rank==0)write(*,'(a,i0,a,es12.4)')'[OW-GS-DIAGNOSTIC] localization_pair_count=',&
      size(pair_first),' support_tolerance=',support_tolerance
    allocate(total_transform(nwannier,nwannier));total_transform=(0d0,0d0)
    do i=1,nwannier;total_transform(i,i)=1d0;end do
    call collective_periodic_spread(comm,values,weights,phases,initial_spread,step_ok,detail)
    call record_spread_evaluation()
    if(.not.step_ok)then;message=trim(detail);return;end if
    final_spread=initial_spread
    if(size(pair_first)<1)then
      maximum_pair_gradient=0d0;iterations=0;converged=.true.;ok=.true.;return
    end if
    have_previous=.false.
    do iterations=1,maximum_iterations
      raw_generator=(0d0,0d0)
      call collective_graph_gradients(comm,values,weights,phases,pair_first,pair_second,&
        graph_gradient_real,graph_gradient_imag,step_ok,detail)
      if(.not.step_ok)then;message=trim(detail);return;end if
      do edge=1,size(pair_first)
        first=pair_first(edge);second=pair_second(edge)
        gradient_real=graph_gradient_real(edge);gradient_imag=graph_gradient_imag(edge)
        current_gradient=sqrt(gradient_real**2+gradient_imag**2)
        if(current_gradient<=tiny(1d0))cycle
        phi=atan2(gradient_imag,gradient_real)+acos(-1d0)
        amplitude=current_gradient*exp(cmplx(0d0,phi,real64))
        raw_generator(first,second)=amplitude
        raw_generator(second,first)=-conjg(amplitude)
      end do
      sweep_generator=(0d0,0d0)
      do operation=1,size(representation,3)
        projection_work=matmul(representation(:,:,operation),raw_generator)
        sweep_generator=sweep_generator+matmul(projection_work,&
          conjg(transpose(representation(:,:,operation))))
      end do
      sweep_generator=sweep_generator/real(size(representation,3),real64)
      generator_scale=maxval(abs(sweep_generator));maximum_pair_gradient=generator_scale
      if(maximum_pair_gradient<=gradient_tolerance)then;converged=.true.;exit;end if
      antihermiticity_defect=maxval(abs(sweep_generator+conjg(transpose(sweep_generator))))/&
        max(1d0,generator_scale)
      commutator_defect=0d0
      do operation=1,size(representation,3)
        scaled_generator=matmul(sweep_generator,representation(:,:,operation))-&
          matmul(representation(:,:,operation),sweep_generator)
        commutator_defect=max(commutator_defect,maxval(abs(scaled_generator))/max(1d0,generator_scale))
      end do
      if(antihermiticity_defect>symmetry_tolerance.or.commutator_defect>symmetry_tolerance)then
        message='batched localization generator violates exact symmetry';return
      end if
      search_generator=sweep_generator;used_conjugate=.false.
      if(have_previous)then
        projection_work=matmul(transport_rotation,previous_gradient)
        transported_gradient=matmul(projection_work,conjg(transpose(transport_rotation)))
        projection_work=matmul(transport_rotation,previous_direction)
        transported_direction=matmul(projection_work,conjg(transpose(transport_rotation)))
        numerator=real(sum(conjg(sweep_generator)*(sweep_generator-transported_gradient)),real64)
        denominator=sum(abs(transported_gradient)**2)
        beta=max(0d0,numerator/max(tiny(1d0),denominator))
        search_generator=sweep_generator+beta*transported_direction
        descent_measure=real(sum(conjg(sweep_generator)*search_generator),real64)
        if(beta>0d0.and.descent_measure>tiny(1d0))then
          used_conjugate=.true.
        else
          search_generator=sweep_generator
        end if
      end if
      backup_values=values;backup_gradients=gradients
      line_accepted=.false.;best_rejected_spread=huge(1d0)
      maximum_attempts=merge(2,1,used_conjugate)
      do attempt=1,maximum_attempts
        if(attempt==2)search_generator=sweep_generator
        generator_scale=maxval(abs(search_generator))
        descent_measure=real(sum(conjg(sweep_generator)*search_generator),real64)
        if(generator_scale<=tiny(1d0).or.descent_measure<=tiny(1d0))cycle
        theta=0.25d0*acos(-1d0)
        do line_search=1,40
          scaled_generator=(theta/generator_scale)*search_generator
          call exponentiate_antihermitian_block(scaled_generator,block_rotation,step_ok,detail)
          if(.not.step_ok)then;message=trim(detail);return;end if
          values=matmul(block_rotation,backup_values)
          do axis=1,3
            gradients(axis,:,:)=matmul(block_rotation,backup_gradients(axis,:,:))
          end do
          call collective_periodic_spread(comm,values,weights,phases,trial_spread,step_ok,detail)
          call record_spread_evaluation()
          if(.not.step_ok)then;message=trim(detail);return;end if
          best_rejected_spread=min(best_rejected_spread,trial_spread)
          if(trial_spread<=final_spread-max(spread_tolerance,&
              1d-4*theta*descent_measure/generator_scale))then
            projection_work=matmul(block_rotation,total_transform)
            call validate_dg_global_covariant_gauge(projection_work,representation,&
              symmetry_tolerance,step_ok,detail)
            if(step_ok)then
              total_transform=projection_work
              final_spread=trial_spread;line_accepted=.true.;exit
            end if
          end if
          values=backup_values;gradients=backup_gradients;theta=0.5d0*theta
        end do
        if(line_accepted)exit
      end do
      if(.not.line_accepted)then
        if(rank==0)write(*,'(a,3(a,es24.16))')'[OW-GS-DIAGNOSTIC] localization_line_search_rejected',&
          ' current_spread=',final_spread,' best_trial_spread=',best_rejected_spread,&
          ' projected_gradient=',maximum_pair_gradient
        message='batched symmetry-constrained localization line search failed';return
      end if
      previous_gradient=sweep_generator;previous_direction=search_generator
      transport_rotation=block_rotation;have_previous=.true.
    end do
    if(.not.converged)then
      iterations=maximum_iterations
      message='symmetry-constrained localization did not converge';return
    end if
    ok=.true.
  contains
    subroutine record_spread_evaluation()
      evaluation_count=evaluation_count+1
      if(present(spread_evaluations))spread_evaluations=evaluation_count
    end subroutine record_spread_evaluation
  end subroutine localize_dg_overlapping_wannier_basis

  subroutine collective_periodic_spread(comm,values,weights,phases,spread,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::values(:,:),phases(:,:)
    real(real64),intent(in)::weights(:)
    real(real64),intent(out)::spread
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64),allocatable::local_norm(:),global_norm(:)
    complex(real64),allocatable::local_moment(:,:),global_moment(:,:)
    integer::point,wannier,axis,ierr
    real(real64)::density
    ok=.false.;message='';spread=0d0
    if(size(values,1)<1.or.size(values,2)<1.or.size(weights)/=size(values,2).or.&
        size(phases,2)/=size(values,2).or.any(weights<=0d0).or.&
        .not.all(ieee_is_finite(weights)).or.&
        maxval(abs(abs(phases)-1d0))>64d0*epsilon(1d0))then
      message='invalid collective periodic localization payload';return
    end if
    allocate(local_norm(size(values,1)),global_norm(size(values,1)),&
      local_moment(size(phases,1),size(values,1)),&
      global_moment(size(phases,1),size(values,1)))
    local_norm=0d0;local_moment=(0d0,0d0)
    do point=1,size(values,2);do wannier=1,size(values,1)
      density=weights(point)*abs(values(wannier,point))**2
      local_norm(wannier)=local_norm(wannier)+density
      do axis=1,size(phases,1)
        local_moment(axis,wannier)=local_moment(axis,wannier)+density*phases(axis,point)
      end do
    end do;end do
#ifdef USE_MPI
    call MPI_Allreduce(local_norm,global_norm,size(local_norm),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_moment,global_moment,size(local_moment),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='collective periodic localization reduction failed';return;end if
#else
    global_norm=local_norm;global_moment=local_moment
#endif
    if(any(global_norm<=tiny(1d0)))then;message='collective Wannier norm is zero';return;end if
    do wannier=1,size(values,1);do axis=1,size(phases,1)
      spread=spread+max(0d0,1d0-abs(global_moment(axis,wannier)/global_norm(wannier))**2)
    end do;end do
    ok=ieee_is_finite(spread)
    if(.not.ok)message='collective periodic spread is not finite'
  end subroutine collective_periodic_spread

  subroutine collective_graph_gradients(comm,values,weights,phases,pair_first,pair_second,&
      gradient_real,gradient_imag,ok,message)
    integer,intent(in)::comm,pair_first(:),pair_second(:)
    complex(real64),intent(in)::values(:,:),phases(:,:)
    real(real64),intent(in)::weights(:)
    real(real64),allocatable,intent(out)::gradient_real(:),gradient_imag(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64),allocatable::local_norm(:,:,:),global_norm(:,:,:),&
      local_norm_derivative(:,:,:),global_norm_derivative(:,:,:)
    complex(real64),allocatable::local_moment(:,:,:,:),global_moment(:,:,:,:),&
      local_moment_derivative(:,:,:,:),global_moment_derivative(:,:,:,:)
    complex(real64)::a,b,direction,phase_factor,quotient,quotient_derivative
    real(real64)::density,density_derivative,derivative,phi
    integer::nedge,naxis,edge,point,wannier,axis,component,first,second,ierr

    ok=.false.;message='';nedge=size(pair_first);naxis=size(phases,1)
    if(nedge<1.or.size(pair_second)/=nedge.or.size(values,2)/=size(weights).or.&
        size(phases,2)/=size(weights).or.any(pair_first<1).or.any(pair_second>size(values,1)))then
      message='invalid collective graph-gradient contract';return
    end if
    allocate(local_norm(2,nedge,2),global_norm(2,nedge,2),&
      local_norm_derivative(2,nedge,2),global_norm_derivative(2,nedge,2),&
      local_moment(naxis,2,nedge,2),global_moment(naxis,2,nedge,2),&
      local_moment_derivative(naxis,2,nedge,2),global_moment_derivative(naxis,2,nedge,2),&
      gradient_real(nedge),gradient_imag(nedge))
    local_norm=0d0;local_norm_derivative=0d0;local_moment=(0d0,0d0)
    local_moment_derivative=(0d0,0d0)
    do component=1,2
      phi=merge(0d0,0.5d0*acos(-1d0),component==1)
      phase_factor=exp(cmplx(0d0,phi,real64))
      do edge=1,nedge
        first=pair_first(edge);second=pair_second(edge)
        do point=1,size(values,2)
          a=values(first,point);b=values(second,point)
          do wannier=1,2
            if(wannier==1)then
              direction=phase_factor*b;density=weights(point)*abs(a)**2
              density_derivative=2d0*weights(point)*real(conjg(a)*direction,real64)
            else
              direction=-conjg(phase_factor)*a;density=weights(point)*abs(b)**2
              density_derivative=2d0*weights(point)*real(conjg(b)*direction,real64)
            end if
            local_norm(wannier,edge,component)=local_norm(wannier,edge,component)+density
            local_norm_derivative(wannier,edge,component)=&
              local_norm_derivative(wannier,edge,component)+density_derivative
            do axis=1,naxis
              local_moment(axis,wannier,edge,component)=&
                local_moment(axis,wannier,edge,component)+density*phases(axis,point)
              local_moment_derivative(axis,wannier,edge,component)=&
                local_moment_derivative(axis,wannier,edge,component)+&
                density_derivative*phases(axis,point)
            end do
          end do
        end do
      end do
    end do
#ifdef USE_MPI
    call MPI_Allreduce(local_norm,global_norm,size(local_norm),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_norm_derivative,global_norm_derivative,size(local_norm_derivative),&
      MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_moment,global_moment,size(local_moment),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_moment_derivative,global_moment_derivative,size(local_moment_derivative),&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='collective graph-gradient reduction failed';return;end if
#else
    global_norm=local_norm;global_norm_derivative=local_norm_derivative
    global_moment=local_moment;global_moment_derivative=local_moment_derivative
#endif
    if(any(global_norm<=tiny(1d0)))then;message='collective graph-gradient norm is zero';return;end if
    do component=1,2;do edge=1,nedge
      derivative=0d0
      do wannier=1,2;do axis=1,naxis
        quotient=global_moment(axis,wannier,edge,component)/global_norm(wannier,edge,component)
        quotient_derivative=(global_moment_derivative(axis,wannier,edge,component)*&
          global_norm(wannier,edge,component)-global_moment(axis,wannier,edge,component)*&
          global_norm_derivative(wannier,edge,component))/global_norm(wannier,edge,component)**2
        derivative=derivative-2d0*real(conjg(quotient)*quotient_derivative,real64)
      end do;end do
      if(component==1)then;gradient_real(edge)=derivative;else;gradient_imag(edge)=derivative;end if
    end do;end do
    ok=all(ieee_is_finite(gradient_real)).and.all(ieee_is_finite(gradient_imag))
    if(.not.ok)message='collective graph gradient is not finite'
  end subroutine collective_graph_gradients

  subroutine collective_pair_gradient(comm,pair,weights,phases,phi,derivative,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::pair(:,:),phases(:,:)
    real(real64),intent(in)::weights(:),phi
    real(real64),intent(out)::derivative
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::direction(:,:),local_moment(:,:),global_moment(:,:),&
      local_moment_derivative(:,:),global_moment_derivative(:,:)
    real(real64)::local_norm(2),global_norm(2),local_norm_derivative(2),&
      global_norm_derivative(2),density,density_derivative
    complex(real64)::phase_factor,quotient,quotient_derivative
    integer::point,wannier,axis,ierr
    ok=.false.;message='';derivative=0d0
    if(size(pair,1)/=2.or.size(pair,2)<1.or.size(weights)/=size(pair,2).or.&
        size(phases,2)/=size(pair,2))then;message='invalid collective pair gradient';return;end if
    allocate(direction(2,size(pair,2)),local_moment(size(phases,1),2),&
      global_moment(size(phases,1),2),local_moment_derivative(size(phases,1),2),&
      global_moment_derivative(size(phases,1),2))
    phase_factor=exp(cmplx(0d0,phi,real64));direction(1,:)=phase_factor*pair(2,:)
    direction(2,:)=-conjg(phase_factor)*pair(1,:)
    local_norm=0d0;local_norm_derivative=0d0;local_moment=(0d0,0d0)
    local_moment_derivative=(0d0,0d0)
    do point=1,size(pair,2);do wannier=1,2
      density=weights(point)*abs(pair(wannier,point))**2
      density_derivative=2d0*weights(point)*real(conjg(pair(wannier,point))*&
        direction(wannier,point),real64)
      local_norm(wannier)=local_norm(wannier)+density
      local_norm_derivative(wannier)=local_norm_derivative(wannier)+density_derivative
      do axis=1,size(phases,1)
        local_moment(axis,wannier)=local_moment(axis,wannier)+density*phases(axis,point)
        local_moment_derivative(axis,wannier)=local_moment_derivative(axis,wannier)+&
          density_derivative*phases(axis,point)
      end do
    end do;end do
#ifdef USE_MPI
    call MPI_Allreduce(local_norm,global_norm,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_norm_derivative,global_norm_derivative,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_moment,global_moment,size(local_moment),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_moment_derivative,global_moment_derivative,size(local_moment_derivative),&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='collective pair-gradient reduction failed';return;end if
#else
    global_norm=local_norm;global_norm_derivative=local_norm_derivative
    global_moment=local_moment;global_moment_derivative=local_moment_derivative
#endif
    if(any(global_norm<=tiny(1d0)))then;message='collective pair norm is zero';return;end if
    do wannier=1,2;do axis=1,size(phases,1)
      quotient=global_moment(axis,wannier)/global_norm(wannier)
      quotient_derivative=(global_moment_derivative(axis,wannier)*global_norm(wannier)-&
        global_moment(axis,wannier)*global_norm_derivative(wannier))/global_norm(wannier)**2
      derivative=derivative-2d0*real(conjg(quotient)*quotient_derivative,real64)
    end do;end do
    ok=ieee_is_finite(derivative)
    if(.not.ok)message='collective pair gradient is not finite'
  end subroutine collective_pair_gradient

  subroutine exponentiate_antihermitian_block(generator,rotation,ok,message)
    complex(real64),intent(in)::generator(:,:)
    complex(real64),allocatable,intent(out)::rotation(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::scaled(:,:),term(:,:),identity(:,:)
    real(real64)::matrix_norm,term_norm
    integer::n,i,k,scaling_steps
    n=size(generator,1);ok=.false.;message=''
    if(n<1.or.size(generator,2)/=n.or.&
        maxval(abs(generator+conjg(transpose(generator))))>1d-10*max(1d0,maxval(abs(generator))))then
      message='invalid anti-Hermitian localization block';return
    end if
    allocate(scaled(n,n),term(n,n),identity(n,n),rotation(n,n))
    identity=(0d0,0d0);do i=1,n;identity(i,i)=1d0;end do
    matrix_norm=maxval(sum(abs(generator),dim=2));scaling_steps=0
    if(matrix_norm>0.25d0)scaling_steps=max(0,ceiling(log(matrix_norm/0.25d0)/log(2d0)))
    scaled=generator/(2d0**scaling_steps);rotation=identity;term=identity
    do k=1,64
      term=matmul(term,scaled)/real(k,8);rotation=rotation+term
      term_norm=maxval(abs(term))
      if(term_norm<=epsilon(1d0)*max(1d0,maxval(abs(rotation))))exit
    end do
    if(k>64)then;message='anti-Hermitian block exponential series did not converge';return;end if
    do k=1,scaling_steps;rotation=matmul(rotation,rotation);end do
    identity=matmul(conjg(transpose(rotation)),rotation)
    do i=1,n;identity(i,i)=identity(i,i)-1d0;end do
    if(maxval(abs(identity))>1d-10)then;message='localization block exponential is not unitary';return;end if
    ok=.true.
  end subroutine exponentiate_antihermitian_block

  subroutine build_dg_overlapping_pair_graph(comm,values,weights,support_tolerance,&
      pair_first,pair_second,pair_support,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::values(:,:)
    real(real64),intent(in)::weights(:),support_tolerance
    integer,allocatable,intent(out)::pair_first(:),pair_second(:)
    real(real64),allocatable,intent(out)::pair_support(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,parameter::maximum_neighbors_per_wannier=8
    real(real64),allocatable::local_diagonal(:),global_diagonal(:),local_row(:),global_row(:),&
      temporary_support(:)
    integer,allocatable::temporary_first(:),temporary_second(:)
    logical,allocatable::selected(:)
    real(real64)::score
    integer::nwannier,npoint,first,second,point,npair,ierr,neighbor,best

    ok=.false.;message='';nwannier=size(values,1);npoint=size(values,2)
    if(nwannier<1.or.npoint<1.or.size(weights)/=npoint.or.support_tolerance<0d0.or.&
        support_tolerance>1d0.or..not.ieee_is_finite(support_tolerance).or.&
        .not.all(ieee_is_finite(weights)).or.any(weights<=0d0).or.&
        .not.all(ieee_is_finite(real(values))).or.&
        .not.all(ieee_is_finite(aimag(values))))then
      message='invalid overlapping Wannier pair-graph contract';return
    end if
    allocate(local_diagonal(nwannier),global_diagonal(nwannier),&
      local_row(nwannier),global_row(nwannier));local_diagonal=0d0
    do point=1,npoint
      local_diagonal=local_diagonal+weights(point)*abs(values(:,point))**4
    end do
#ifdef USE_MPI
    call MPI_Allreduce(local_diagonal,global_diagonal,nwannier,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='pair-graph diagonal collective failed';return;end if
#else
    global_diagonal=local_diagonal
#endif
    if(any(global_diagonal<=tiny(1d0)))then
      message='pair-graph Wannier support norm is zero';return
    end if
    allocate(selected(nwannier),temporary_first(maximum_neighbors_per_wannier*nwannier),&
      temporary_second(maximum_neighbors_per_wannier*nwannier),&
      temporary_support(maximum_neighbors_per_wannier*nwannier))
    npair=0
    do first=1,nwannier
      call pair_support_row(first,local_row)
#ifdef USE_MPI
      call MPI_Allreduce(local_row,global_row,nwannier,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='pair-graph row collective failed';return;end if
#else
      global_row=local_row
#endif
      selected=.false.
      do neighbor=1,min(maximum_neighbors_per_wannier,nwannier-1)
        best=0;score=support_tolerance
        do second=1,nwannier
          if(second==first.or.selected(second))cycle
          if(global_row(second)/sqrt(global_diagonal(first)*global_diagonal(second))>score)then
            best=second
            score=global_row(second)/sqrt(global_diagonal(first)*global_diagonal(second))
          end if
        end do
        if(best==0)exit
        selected(best)=.true.
        call append_symmetric_edge(first,best,score)
      end do
    end do
    allocate(pair_first(npair),pair_second(npair),pair_support(npair))
    pair_first=temporary_first(:npair);pair_second=temporary_second(:npair)
    pair_support=temporary_support(:npair)
    ok=.true.
  contains
    subroutine append_symmetric_edge(source,target,edge_support)
      integer,intent(in)::source,target
      real(real64),intent(in)::edge_support
      integer::low,high,edge
      low=min(source,target);high=max(source,target)
      do edge=1,npair
        if(temporary_first(edge)==low.and.temporary_second(edge)==high)then
          temporary_support(edge)=max(temporary_support(edge),min(1d0,max(0d0,edge_support)))
          return
        end if
      end do
      npair=npair+1
      temporary_first(npair)=low;temporary_second(npair)=high
      temporary_support(npair)=min(1d0,max(0d0,edge_support))
    end subroutine append_symmetric_edge

    subroutine pair_support_row(source,row)
      integer,intent(in)::source
      real(real64),intent(out)::row(:)
      integer::p
      row=0d0
      do p=1,npoint
        row=row+weights(p)*abs(values(source,p))**2*abs(values(:,p))**2
      end do
    end subroutine pair_support_row
  end subroutine build_dg_overlapping_pair_graph

  subroutine optimize_dg_wannier_pair(values,gradients,weights,phases,first,second,&
      spread_tolerance,rotation,before,after,gradient,accepted,ok,message)
    complex(real64),intent(inout)::values(:,:),gradients(:,:,:)
    real(real64),intent(in)::weights(:),spread_tolerance
    complex(real64),intent(in)::phases(:,:)
    integer,intent(in)::first,second
    complex(real64),intent(out)::rotation(2,2)
    real(real64),intent(out)::before,after,gradient
    logical,intent(out)::accepted,ok
    character(*),intent(out)::message
    complex(real64),allocatable::pair(:,:),candidate(:,:)
    complex(real64)::trial_rotation(2,2),best_rotation(2,2)
    real(real64)::trial_spread,best_spread,theta,phi,gradient_real,gradient_imag,&
      current_gradient
    integer::iteration,line_search,point,axis
    logical::trial_ok,line_accepted
    character(256)::detail

    ok=.false.;accepted=.false.;message='';before=huge(1d0);after=huge(1d0);gradient=huge(1d0)
    rotation=(0d0,0d0);rotation(1,1)=1d0;rotation(2,2)=1d0
    if(size(values,1)<2.or.size(values,2)<1.or.first<1.or.second<1.or.&
        first>size(values,1).or.second>size(values,1).or.first==second.or.&
        size(gradients,1)/=3.or.size(gradients,2)/=size(values,1).or.&
        size(gradients,3)/=size(values,2).or.spread_tolerance<0d0.or.&
        .not.ieee_is_finite(spread_tolerance))then
      message='invalid Wannier pair-localization contract';return
    end if
    allocate(pair(2,size(values,2)),candidate(2,size(values,2)))
    pair(1,:)=values(first,:);pair(2,:)=values(second,:)
    call pair_spread(pair,weights,phases,before,trial_ok,detail)
    if(.not.trial_ok)then;message=trim(detail);return;end if

    best_spread=before;best_rotation=rotation
    gradient=0d0
    do iteration=1,32
      call pair_spread_directional_derivative(pair,weights,phases,0d0,gradient_real,trial_ok,detail)
      if(.not.trial_ok)then;message=trim(detail);return;end if
      call pair_spread_directional_derivative(pair,weights,phases,0.5d0*acos(-1d0),&
        gradient_imag,trial_ok,detail)
      if(.not.trial_ok)then;message=trim(detail);return;end if
      current_gradient=sqrt(gradient_real**2+gradient_imag**2)
      if(iteration==1)gradient=current_gradient
      if(current_gradient<=max(sqrt(epsilon(1d0)),spread_tolerance))exit
      phi=atan2(gradient_imag,gradient_real)+acos(-1d0)
      theta=0.25d0*acos(-1d0);line_accepted=.false.
      do line_search=1,40
        call make_pair_rotation(theta,phi,trial_rotation)
        candidate=matmul(trial_rotation,pair)
        call pair_spread(candidate,weights,phases,trial_spread,trial_ok,detail)
        if(.not.trial_ok)then;message=trim(detail);return;end if
        if(trial_spread<=best_spread-max(spread_tolerance,&
            1d-4*theta*current_gradient))then
          pair=candidate;best_spread=trial_spread
          best_rotation=matmul(trial_rotation,best_rotation)
          line_accepted=.true.;exit
        end if
        theta=0.5d0*theta
      end do
      if(.not.line_accepted)exit
    end do
    after=before
    if(best_spread<before-spread_tolerance)then
      do point=1,size(values,2)
        values([first,second],point)=matmul(best_rotation,values([first,second],point))
        do axis=1,3
          gradients(axis,[first,second],point)=matmul(best_rotation,&
            gradients(axis,[first,second],point))
        end do
      end do
      rotation=best_rotation;after=best_spread;accepted=.true.
    end if
    ok=.true.
  end subroutine optimize_dg_wannier_pair

  subroutine pair_spread_directional_derivative(pair,weights,phases,phi,derivative,ok,message)
    complex(real64),intent(in)::pair(:,:),phases(:,:)
    real(real64),intent(in)::weights(:),phi
    real(real64),intent(out)::derivative
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::direction(:,:),moment(:,:),moment_derivative(:,:)
    real(real64)::norm(2),norm_derivative(2),spread,density_derivative
    complex(real64)::phase_factor,quotient,quotient_derivative
    integer::point,wannier,axis

    allocate(direction(2,size(pair,2)),moment(size(phases,1),2),&
      moment_derivative(size(phases,1),2))
    call evaluate_dg_periodic_localization(pair,weights,phases,norm,moment,spread,ok,message)
    if(.not.ok)return
    phase_factor=exp(cmplx(0d0,phi,real64))
    direction(1,:)=phase_factor*pair(2,:)
    direction(2,:)=-conjg(phase_factor)*pair(1,:)
    norm_derivative=0d0;moment_derivative=(0d0,0d0)
    do point=1,size(pair,2);do wannier=1,2
      density_derivative=2d0*weights(point)*real(conjg(pair(wannier,point))*&
        direction(wannier,point),real64)
      norm_derivative(wannier)=norm_derivative(wannier)+density_derivative
      do axis=1,size(phases,1)
        moment_derivative(axis,wannier)=moment_derivative(axis,wannier)+&
          density_derivative*phases(axis,point)
      end do
    end do;end do
    derivative=0d0
    do wannier=1,2;do axis=1,size(phases,1)
      quotient=moment(axis,wannier)/norm(wannier)
      quotient_derivative=(moment_derivative(axis,wannier)*norm(wannier)-&
        moment(axis,wannier)*norm_derivative(wannier))/norm(wannier)**2
      derivative=derivative-2d0*real(conjg(quotient)*quotient_derivative,real64)
    end do;end do
    ok=ieee_is_finite(derivative)
    if(.not.ok)message='periodic localization pair gradient is not finite'
  end subroutine pair_spread_directional_derivative

  subroutine make_pair_rotation(theta,phi,rotation)
    real(real64),intent(in)::theta,phi
    complex(real64),intent(out)::rotation(2,2)
    complex(real64)::phase
    phase=exp(cmplx(0d0,phi,real64))
    rotation(1,1)=cos(theta);rotation(1,2)=phase*sin(theta)
    rotation(2,1)=-conjg(phase)*sin(theta);rotation(2,2)=cos(theta)
  end subroutine make_pair_rotation

  subroutine pair_spread(pair,weights,phases,spread,ok,message)
    complex(real64),intent(in)::pair(:,:),phases(:,:)
    real(real64),intent(in)::weights(:)
    real(real64),intent(out)::spread
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64)::norm(2)
    complex(real64),allocatable::moment(:,:)
    allocate(moment(size(phases,1),2))
    call evaluate_dg_periodic_localization(pair,weights,phases,norm,moment,spread,ok,message)
  end subroutine pair_spread

  subroutine evaluate_dg_periodic_localization(values,weights,phases,norm,moment,spread,ok,message)
    complex(real64),intent(in)::values(:,:),phases(:,:)
    real(real64),intent(in)::weights(:)
    real(real64),intent(out)::norm(:),spread
    complex(real64),intent(out)::moment(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64)::density
    integer::nwannier,npoint,naxis,wannier,point,axis

    ok=.false.;message='';spread=0d0;norm=0d0;moment=(0d0,0d0)
    nwannier=size(values,1);npoint=size(values,2);naxis=size(phases,1)
    if(nwannier<1.or.npoint<1.or.naxis<1.or.size(weights)/=npoint.or.&
        size(phases,2)/=npoint.or.size(norm)/=nwannier.or.&
        any(shape(moment)/=[naxis,nwannier]))then
      message='periodic localization arrays have inconsistent dimensions';return
    end if
    if(.not.all(ieee_is_finite(weights)).or.&
        .not.all(ieee_is_finite(real(values))).or.&
        .not.all(ieee_is_finite(aimag(values))).or.&
        .not.all(ieee_is_finite(real(phases))).or.&
        .not.all(ieee_is_finite(aimag(phases))))then
      message='periodic localization payload is not finite';return
    end if
    if(any(weights<=0d0))then
      message='periodic localization weights must be positive';return
    end if
    if(maxval(abs(abs(phases)-1d0))>64d0*epsilon(1d0))then
      message='periodic localization phase is not unit modulus';return
    end if
    do point=1,npoint;do wannier=1,nwannier
      density=weights(point)*abs(values(wannier,point))**2
      norm(wannier)=norm(wannier)+density
      do axis=1,naxis
        moment(axis,wannier)=moment(axis,wannier)+density*phases(axis,point)
      end do
    end do;end do
    if(any(norm<=tiny(1d0)))then
      message='periodic localization Wannier norm is zero';return
    end if
    do wannier=1,nwannier;do axis=1,naxis
      spread=spread+max(0d0,1d0-abs(moment(axis,wannier)/norm(wannier))**2)
    end do;end do
    if(.not.ieee_is_finite(spread))then
      message='periodic localization spread is not finite';return
    end if
    ok=.true.
  end subroutine evaluate_dg_periodic_localization
end module dg_overlapping_wannier_localization
