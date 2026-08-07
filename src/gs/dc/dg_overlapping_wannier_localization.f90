module dg_overlapping_wannier_localization
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use,intrinsic::iso_fortran_env,only:real64
#ifdef USE_MPI
  use mpi
#endif
  use dg_overlapping_wannier_symmetry,only:build_dg_symmetry_constrained_pair_generator
  implicit none
  private
  public::evaluate_dg_periodic_localization
  public::optimize_dg_wannier_pair
  public::build_dg_overlapping_pair_graph
  public::localize_dg_overlapping_wannier_basis
contains
  subroutine localize_dg_overlapping_wannier_basis(comm,values,gradients,weights,phases,&
      representation,product_table,support_tolerance,spread_tolerance,gradient_tolerance,&
      symmetry_tolerance,maximum_iterations,initial_spread,final_spread,maximum_pair_gradient,&
      iterations,converged,total_transform,ok,message)
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
    integer,allocatable::pair_first(:),pair_second(:),active_indices(:)
    real(real64),allocatable::pair_support(:)
    complex(real64),allocatable::generator(:,:),block_rotation(:,:),backup_values(:,:),&
      backup_gradients(:,:,:),transform_rows(:,:)
    real(real64)::gradient_real,gradient_imag,current_gradient,theta,phi,trial_spread,&
      antihermiticity_defect,commutator_defect
    complex(real64)::amplitude
    integer::nwannier,edge,line_search,axis,i
    logical::step_ok
    character(256)::detail

    ok=.false.;converged=.false.;message='';iterations=0
    initial_spread=huge(1d0);final_spread=huge(1d0);maximum_pair_gradient=huge(1d0)
    nwannier=size(values,1)
    if(nwannier<2.or.size(values,2)<1.or.any(shape(gradients)/=[3,nwannier,size(values,2)]).or.&
        size(weights)/=size(values,2).or.size(phases,2)/=size(values,2).or.&
        size(representation,1)/=nwannier.or.size(representation,2)/=nwannier.or.&
        size(representation,3)<1.or.any(shape(product_table)/=&
        [size(representation,3),size(representation,3)]).or.maximum_iterations<1.or.&
        support_tolerance<0d0.or.support_tolerance>1d0.or.spread_tolerance<0d0.or.&
        gradient_tolerance<=0d0.or.symmetry_tolerance<=0d0)then
      message='invalid symmetry-constrained localization sweep contract';return
    end if
    call build_dg_overlapping_pair_graph(comm,values,weights,support_tolerance,&
      pair_first,pair_second,pair_support,step_ok,detail)
    if(.not.step_ok)then;message=trim(detail);return;end if
    if(size(pair_first)<1)then;message='localization pair graph has no overlapping edge';return;end if
    allocate(total_transform(nwannier,nwannier));total_transform=(0d0,0d0)
    do i=1,nwannier;total_transform(i,i)=1d0;end do
    call collective_periodic_spread(comm,values,weights,phases,initial_spread,step_ok,detail)
    if(.not.step_ok)then;message=trim(detail);return;end if
    final_spread=initial_spread
    do iterations=1,maximum_iterations
      maximum_pair_gradient=0d0
      do edge=1,size(pair_first)
        call collective_pair_gradient(comm,values([pair_first(edge),pair_second(edge)],:),&
          weights,phases,0d0,gradient_real,step_ok,detail)
        if(.not.step_ok)then;message=trim(detail);return;end if
        call collective_pair_gradient(comm,values([pair_first(edge),pair_second(edge)],:),&
          weights,phases,0.5d0*acos(-1d0),gradient_imag,step_ok,detail)
        if(.not.step_ok)then;message=trim(detail);return;end if
        current_gradient=sqrt(gradient_real**2+gradient_imag**2)
        maximum_pair_gradient=max(maximum_pair_gradient,current_gradient)
        if(current_gradient<=gradient_tolerance)cycle
        phi=atan2(gradient_imag,gradient_real)+acos(-1d0)
        theta=0.25d0*acos(-1d0)
        do line_search=1,40
          amplitude=theta*exp(cmplx(0d0,phi,real64))
          call build_dg_symmetry_constrained_pair_generator(pair_first(edge),pair_second(edge),&
            amplitude,representation,product_table,symmetry_tolerance,generator,&
            antihermiticity_defect,commutator_defect,active_indices,step_ok,detail)
          if(.not.step_ok)then;message=trim(detail);return;end if
          call exponentiate_antihermitian_block(generator,block_rotation,step_ok,detail)
          if(.not.step_ok)then;message=trim(detail);return;end if
          backup_values=values(active_indices,:);backup_gradients=gradients(:,active_indices,:)
          values(active_indices,:)=matmul(block_rotation,backup_values)
          do axis=1,3
            gradients(axis,active_indices,:)=matmul(block_rotation,backup_gradients(axis,:,:))
          end do
          call collective_periodic_spread(comm,values,weights,phases,trial_spread,step_ok,detail)
          if(.not.step_ok)then;message=trim(detail);return;end if
          if(trial_spread<=final_spread-max(spread_tolerance,&
              1d-4*theta*current_gradient))then
            transform_rows=matmul(block_rotation,total_transform(active_indices,:))
            total_transform(active_indices,:)=transform_rows
            final_spread=trial_spread;exit
          end if
          values(active_indices,:)=backup_values;gradients(:,active_indices,:)=backup_gradients
          theta=0.5d0*theta
        end do
      end do
      call maximum_graph_gradient(comm,values,weights,phases,pair_first,pair_second,&
        maximum_pair_gradient,step_ok,detail)
      if(.not.step_ok)then;message=trim(detail);return;end if
      if(maximum_pair_gradient<=gradient_tolerance)then;converged=.true.;exit;end if
    end do
    if(.not.converged)then
      iterations=maximum_iterations
      message='symmetry-constrained localization did not converge';return
    end if
    ok=.true.
  end subroutine localize_dg_overlapping_wannier_basis

  subroutine maximum_graph_gradient(comm,values,weights,phases,pair_first,pair_second,&
      maximum_gradient,ok,message)
    integer,intent(in)::comm,pair_first(:),pair_second(:)
    complex(real64),intent(in)::values(:,:),phases(:,:)
    real(real64),intent(in)::weights(:)
    real(real64),intent(out)::maximum_gradient
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64)::g0,g1
    integer::edge
    character(256)::detail
    maximum_gradient=0d0;ok=.false.;message=''
    do edge=1,size(pair_first)
      call collective_pair_gradient(comm,values([pair_first(edge),pair_second(edge)],:),&
        weights,phases,0d0,g0,ok,detail)
      if(.not.ok)then;message=trim(detail);return;end if
      call collective_pair_gradient(comm,values([pair_first(edge),pair_second(edge)],:),&
        weights,phases,0.5d0*acos(-1d0),g1,ok,detail)
      if(.not.ok)then;message=trim(detail);return;end if
      maximum_gradient=max(maximum_gradient,sqrt(g0**2+g1**2))
    end do
    ok=.true.
  end subroutine maximum_graph_gradient

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
    real(real64),allocatable::local_diagonal(:),global_diagonal(:),local_row(:),global_row(:)
    real(real64)::score
    integer::nwannier,npoint,first,second,point,npair,index,ierr

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
    npair=0
    do first=1,nwannier-1
      call pair_support_row(first,local_row)
#ifdef USE_MPI
      call MPI_Allreduce(local_row,global_row,nwannier,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='pair-graph row collective failed';return;end if
#else
      global_row=local_row
#endif
      do second=first+1,nwannier
        score=global_row(second)/sqrt(global_diagonal(first)*global_diagonal(second))
        if(score>support_tolerance)npair=npair+1
      end do
    end do
    allocate(pair_first(npair),pair_second(npair),pair_support(npair));index=0
    do first=1,nwannier-1
      call pair_support_row(first,local_row)
#ifdef USE_MPI
      call MPI_Allreduce(local_row,global_row,nwannier,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='pair-graph row collective failed';return;end if
#else
      global_row=local_row
#endif
      do second=first+1,nwannier
        score=global_row(second)/sqrt(global_diagonal(first)*global_diagonal(second))
        if(score<=support_tolerance)cycle
        index=index+1;pair_first(index)=first;pair_second(index)=second
        pair_support(index)=min(1d0,max(0d0,score))
      end do
    end do
    ok=.true.
  contains
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
