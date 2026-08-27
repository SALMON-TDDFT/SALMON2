module dg_hybrid_continuation_residuals
  use,intrinsic::iso_fortran_env,only:real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private

  type,public::s_dg_hybrid_residuals
    real(real64)::r_h=huge(1d0)
    real(real64)::r_rho=huge(1d0)
    real(real64)::r_t=huge(1d0)
    real(real64)::r_s=huge(1d0)
  end type s_dg_hybrid_residuals

  public::build_dg_hybrid_occupied_algebra,build_dg_hybrid_interface_observables,&
    evaluate_dg_hybrid_residuals,validate_cluster_occupations,evaluate_dg_hybrid_projector_change,&
    dg_hybrid_electron_count
contains
  subroutine evaluate_dg_hybrid_projector_change(previous,current,metric,residual,ok,message)
    complex(real64),intent(in)::previous(:,:),current(:,:),metric(:,:)
    real(real64),intent(out)::residual
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::weighted_difference(:,:),weighted_previous(:,:)
    integer::n
    ok=.false.;message='';residual=huge(1d0);n=size(previous,1)
    if(n<1.or.size(previous,2)/=n.or.any(shape(current)/=[n,n]).or.any(shape(metric)/=[n,n]).or.&
        .not.finite_complex(previous).or..not.finite_complex(current).or..not.finite_complex(metric))then
      message='invalid occupied-projector comparison';return
    endif
    allocate(weighted_difference(n,n),weighted_previous(n,n))
    weighted_difference=matmul(metric,current-previous)
    weighted_previous=matmul(metric,previous)
    residual=frobenius(weighted_difference)/max(1d0,frobenius(weighted_previous));ok=.true.
  end subroutine evaluate_dg_hybrid_projector_change

  real(real64) function dg_hybrid_electron_count(gamma,metric) result(count)
    complex(real64),intent(in)::gamma(:,:),metric(:,:)
    complex(real64),allocatable::product(:,:)
    integer::i,n
    count=huge(1d0);n=size(gamma,1)
    if(n<1.or.size(gamma,2)/=n.or.any(shape(metric)/=[n,n]).or.&
        .not.finite_complex(gamma).or..not.finite_complex(metric))return
    allocate(product(n,n));product=matmul(metric,gamma);count=0d0
    do i=1,n;count=count+real(product(i,i));enddo
  end function dg_hybrid_electron_count

  subroutine build_dg_hybrid_occupied_algebra(coefficients,metric,occupations,projector,gamma,ok,message)
    complex(real64),intent(in)::coefficients(:,:),metric(:,:)
    real(real64),intent(in)::occupations(:)
    complex(real64),intent(out)::projector(:,:),gamma(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::overlap(:,:),weighted(:,:)
    integer::i,n,m
    real(real64)::defect
    ok=.false.;message='';n=size(coefficients,1);m=size(coefficients,2)
    if(n<1.or.m<1.or.size(metric,1)/=n.or.size(metric,2)/=n.or.size(occupations)/=m.or.&
        size(projector,1)/=n.or.size(projector,2)/=n.or.size(gamma,1)/=n.or.size(gamma,2)/=n.or.&
        .not.finite_complex(coefficients).or..not.finite_complex(metric).or.&
        .not.all(ieee_is_finite(occupations)).or.any(occupations<0d0))then
      message='invalid occupied algebra input';return
    endif
    if(.not.positive_definite(metric))then;message='DG metric is not Hermitian positive definite';return;endif
    allocate(overlap(m,m),weighted(n,m));overlap=matmul(conjg(transpose(coefficients)),&
      matmul(metric,coefficients))
    do i=1,m;overlap(i,i)=overlap(i,i)-1d0;enddo
    defect=frobenius(overlap)
    if(defect>1d-10)then;message='occupied coefficients are rank deficient or not S-orthonormal';return;endif
    projector=matmul(coefficients,matmul(conjg(transpose(coefficients)),metric))
    weighted=coefficients
    do i=1,m;weighted(:,i)=occupations(i)*weighted(:,i);enddo
    gamma=matmul(weighted,conjg(transpose(coefficients)))
    ok=.true.
  end subroutine build_dg_hybrid_occupied_algebra

  subroutine build_dg_hybrid_interface_observables(value_trace,normal_trace,gamma,value_density,&
      normal_density,cross_density,ok,message)
    complex(real64),intent(in)::value_trace(:,:),normal_trace(:,:),gamma(:,:)
    complex(real64),intent(out)::value_density(:,:),normal_density(:,:),cross_density(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::ntrace,nbasis
    ok=.false.;message='';ntrace=size(value_trace,1);nbasis=size(value_trace,2)
    if(ntrace<1.or.nbasis<1.or.any(shape(normal_trace)/=shape(value_trace)).or.&
        size(gamma,1)/=nbasis.or.size(gamma,2)/=nbasis.or.&
        any(shape(value_density)/=[ntrace,ntrace]).or.any(shape(normal_density)/=[ntrace,ntrace]).or.&
        any(shape(cross_density)/=[ntrace,ntrace]).or..not.finite_complex(value_trace).or.&
        .not.finite_complex(normal_trace).or..not.finite_complex(gamma))then
      message='invalid occupied interface trace input';return
    endif
    value_density=matmul(value_trace,matmul(gamma,conjg(transpose(value_trace))))
    normal_density=matmul(normal_trace,matmul(gamma,conjg(transpose(normal_trace))))
    cross_density=matmul(value_trace,matmul(gamma,conjg(transpose(normal_trace))))
    ok=.true.
  end subroutine build_dg_hybrid_interface_observables

  subroutine evaluate_dg_hybrid_residuals(hc,sc_epsilon,coefficients,metric,rho_output,rho_input,&
      trace_output,trace_input,residuals,ok,message)
    complex(real64),intent(in)::hc(:,:),sc_epsilon(:,:),coefficients(:,:),metric(:,:),&
      trace_output(:,:),trace_input(:,:)
    real(real64),intent(in)::rho_output(:),rho_input(:)
    type(s_dg_hybrid_residuals),intent(out)::residuals
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::orthogonality(:,:)
    integer::i,m
    ok=.false.;message='';m=size(coefficients,2)
    if(any(shape(hc)/=shape(sc_epsilon)).or.size(hc,1)/=size(coefficients,1).or.&
        size(hc,2)/=m.or.any(shape(metric)/=[size(coefficients,1),size(coefficients,1)]).or.&
        size(rho_output)/=size(rho_input).or.any(shape(trace_output)/=shape(trace_input)).or.&
        .not.finite_complex(hc).or..not.finite_complex(sc_epsilon).or.&
        .not.finite_complex(coefficients).or..not.finite_complex(metric).or.&
        .not.all(ieee_is_finite(rho_output)).or..not.all(ieee_is_finite(rho_input)).or.&
        .not.finite_complex(trace_output).or..not.finite_complex(trace_input))then
      message='invalid continuation residual input';return
    endif
    residuals%r_h=frobenius(hc-sc_epsilon)/max(1d0,frobenius(hc),frobenius(sc_epsilon))
    residuals%r_rho=sqrt(sum((rho_output-rho_input)**2))/max(1d0,sqrt(sum(rho_input**2)))
    residuals%r_t=frobenius(trace_output-trace_input)/max(1d0,frobenius(trace_input))
    allocate(orthogonality(m,m));orthogonality=matmul(conjg(transpose(coefficients)),&
      matmul(metric,coefficients))
    do i=1,m;orthogonality(i,i)=orthogonality(i,i)-1d0;enddo
    residuals%r_s=frobenius(orthogonality);ok=.true.
  end subroutine evaluate_dg_hybrid_residuals

  subroutine validate_cluster_occupations(occupations,cluster_ids,ok,message)
    real(real64),intent(in)::occupations(:)
    integer,intent(in)::cluster_ids(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j
    ok=.false.;message=''
    if(size(occupations)<1.or.size(cluster_ids)/=size(occupations).or.any(cluster_ids<1).or.&
        .not.all(ieee_is_finite(occupations)).or.any(occupations<0d0))then
      message='invalid occupation cluster';return
    endif
    do i=1,size(occupations);do j=i+1,size(occupations)
      if(cluster_ids(i)==cluster_ids(j).and.abs(occupations(i)-occupations(j))>1d-12)then
        message='symmetry-incompatible occupations in a degenerate cluster';return
      endif
    enddo;enddo
    ok=.true.
  end subroutine validate_cluster_occupations

  real(real64) function frobenius(values) result(norm)
    complex(real64),intent(in)::values(:,:)
    norm=sqrt(sum(abs(values)**2))
  end function frobenius

  logical function finite_complex(values) result(ok)
    complex(real64),intent(in)::values(:,:)
    ok=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_complex

  logical function positive_definite(matrix) result(ok)
    complex(real64),intent(in)::matrix(:,:)
    complex(real64),allocatable::lower(:,:)
    complex(real64)::pivot
    real(real64)::scale,tolerance
    integer::i,j,k,n
    n=size(matrix,1);ok=.false.
    if(n<1.or.size(matrix,2)/=n)return
    scale=max(1d0,maxval(abs(matrix)));tolerance=1d-12*scale
    if(maxval(abs(matrix-conjg(transpose(matrix))))>tolerance)return
    allocate(lower(n,n));lower=(0d0,0d0)
    do i=1,n
      do j=1,i
        pivot=matrix(i,j)
        do k=1,j-1;pivot=pivot-lower(i,k)*conjg(lower(j,k));enddo
        if(i==j)then
          if(abs(aimag(pivot))>tolerance.or.real(pivot)<=tolerance)return
          lower(i,j)=sqrt(real(pivot))
        else
          lower(i,j)=pivot/lower(j,j)
        endif
      enddo
    enddo
    ok=.true.
  end function positive_definite
end module dg_hybrid_continuation_residuals
