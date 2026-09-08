program test_dg_buffer_window_projector
  use,intrinsic :: ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf
  use dg_buffer_window_projector, only: s_dg_buffer_projector_diagnostics, &
    build_dg_buffer_window_projector
  implicit none
  real(8), parameter :: tolerance=1d-12
  real(8) :: phi(5,3),buffer(5,2),weights(5),coefficients(3,2),reconstructed(5,2)
  real(8) :: transformed(5,3),transformed_coefficients(3,2),transformed_reconstructed(5,2)
  real(8) :: exact(5,2),outside(5),rotation(3,3),rank_deficient(5,4)
  real(8) :: deficient_coefficients(4,2),deficient_reconstructed(5,2)
  real(8) :: occupied_coefficients(2,2),occupied_reconstructed(5,2)
  real(8) :: cutoff_phi(2,2),cutoff_buffer(2,1),cutoff_coefficients(2,1)
  real(8) :: cutoff_reconstructed(2,1),scaled_phi(5,3),scaled_coefficients(3,2)
  real(8) :: scaled_reconstructed(5,2),clean_buffer(5,1),clean_coefficients(3,1)
  real(8) :: clean_reconstructed(5,1),invalid_phi(5,3),invalid_buffer(5,2)
  real(8) :: extreme_phi(2,1),extreme_buffer(2,1),extreme_weight(2)
  real(8) :: extreme_coefficients(1,1),extreme_reconstructed(2,1)
  real(8),allocatable :: empty_phi(:,:),empty_buffer(:,:),empty_weights(:)
  real(8),allocatable :: empty_coefficients(:,:),empty_reconstructed(:,:)
  real(8) :: expected_residual,expected_escape
  type(s_dg_buffer_projector_diagnostics) :: diagnostics,transformed_diagnostics
  logical :: ok
  character(len=256) :: message

  weights=[1d0,2d0,1.5d0,0.5d0,3d0]
  phi(:,1)=[1d0,0d0,1d0,2d0,-1d0]
  phi(:,2)=[0d0,1d0,1d0,-1d0,2d0]
  phi(:,3)=[1d0,1d0,0d0,1d0,1d0]
  exact(:,1)=matmul(phi,[0.5d0,-1.25d0,2d0])
  outside=[1d0,-2d0,1d0,0.5d0,-1d0]
  outside=outside-matmul(phi,solve_weighted_normal(phi,weights,outside))
  outside=outside/sqrt(sum(weights*outside**2))
  exact(:,2)=matmul(phi,[-0.75d0,0.25d0,1.5d0])
  buffer(:,1)=exact(:,1)
  buffer(:,2)=exact(:,2)+0.2d0*outside

  call build_dg_buffer_window_projector(phi,buffer,weights,1d-12,coefficients, &
    reconstructed,diagnostics,ok,message)
  if(.not.ok) error stop trim(message)
  if(maxval(abs(reconstructed(:,1)-exact(:,1)))>tolerance) error stop 'projection mismatch'
  if(maxval(abs(reconstructed(:,2)-exact(:,2)))>tolerance) error stop 'outside projection mismatch'
  expected_escape=0.2d0
  expected_residual=expected_escape/sqrt(sum(spread(weights,2,2)*buffer**2))
  if(abs(diagnostics%projection_residual-expected_residual)>tolerance) error stop 'bad residual'
  if(abs(diagnostics%escape_norm-expected_escape)>tolerance) error stop 'bad escape norm'
  if(diagnostics%configured_states/=3.or.diagnostics%retained_rank/=3) &
    error stop 'bad full-window diagnostics'

  transformed=phi
  transformed(:,[1,3])=-transformed(:,[1,3])
  call build_dg_buffer_window_projector(transformed,buffer,weights,1d-12, &
    transformed_coefficients,transformed_reconstructed,transformed_diagnostics,ok,message)
  if(.not.ok) error stop trim(message)
  call assert_invariant(reconstructed,diagnostics,transformed_reconstructed, &
    transformed_diagnostics,'independent sign invariance')

  transformed=phi(:,[3,1,2])
  call build_dg_buffer_window_projector(transformed,buffer,weights,1d-12, &
    transformed_coefficients,transformed_reconstructed,transformed_diagnostics,ok,message)
  if(.not.ok) error stop trim(message)
  call assert_invariant(reconstructed,diagnostics,transformed_reconstructed, &
    transformed_diagnostics,'permutation invariance')

  rotation=0d0
  rotation(1,1)=sqrt(0.5d0);rotation(1,2)=-sqrt(0.5d0)
  rotation(2,1)=sqrt(0.5d0);rotation(2,2)=sqrt(0.5d0)
  rotation(3,3)=1d0
  transformed=matmul(phi,rotation)
  call build_dg_buffer_window_projector(transformed,buffer,weights,1d-12, &
    transformed_coefficients,transformed_reconstructed,transformed_diagnostics,ok,message)
  if(.not.ok) error stop trim(message)
  call assert_invariant(reconstructed,diagnostics,transformed_reconstructed, &
    transformed_diagnostics,'orthogonal rotation invariance')

  rank_deficient(:,1:3)=phi
  rank_deficient(:,4)=phi(:,1)+phi(:,2)
  call build_dg_buffer_window_projector(rank_deficient,buffer,weights,1d-12, &
    deficient_coefficients,deficient_reconstructed,transformed_diagnostics,ok,message)
  if(.not.ok) error stop trim(message)
  if(transformed_diagnostics%configured_states/=4.or.transformed_diagnostics%retained_rank/=3) &
    error stop 'rank-deficient metric was not truncated'
  if(maxval(abs(deficient_reconstructed-reconstructed))>2d-12) &
    error stop 'rank-deficient reconstruction changed retained window'

  clean_buffer(:,1)=phi(:,3)
  call build_dg_buffer_window_projector(phi,clean_buffer,weights,1d-12, &
    clean_coefficients,clean_reconstructed,transformed_diagnostics,ok,message)
  if(.not.ok.or.transformed_diagnostics%projection_residual>tolerance) &
    error stop 'full window failed known empty-state target'
  call build_dg_buffer_window_projector(phi(:,1:2),clean_buffer,weights,1d-12, &
    occupied_coefficients(:,1:1),occupied_reconstructed(:,1:1),transformed_diagnostics,ok,message)
  if(.not.ok) error stop trim(message)
  if(transformed_diagnostics%projection_residual<=1d-3) &
    error stop 'occupied-only window did not expose missing empty-state content'

  cutoff_phi=0d0
  cutoff_phi(1,1)=1d0
  cutoff_phi(2,2)=0.5d0
  cutoff_buffer(:,1)=[1d0,0.5d0]
  call build_dg_buffer_window_projector(cutoff_phi,cutoff_buffer,[1d0,1d0],0.25d0, &
    cutoff_coefficients,cutoff_reconstructed,transformed_diagnostics,ok,message)
  if(.not.ok.or.transformed_diagnostics%retained_rank/=2) error stop 'cutoff equality was not retained'
  if(abs(transformed_diagnostics%minimum_retained_singular_value-0.25d0)>tolerance) &
    error stop 'bad minimum retained singular value'
  call build_dg_buffer_window_projector(cutoff_phi,cutoff_buffer,[1d0,1d0],1d0, &
    cutoff_coefficients,cutoff_reconstructed,transformed_diagnostics,ok,message)
  if(.not.ok.or.transformed_diagnostics%retained_rank/=1) error stop 'unit cutoff was rejected'
  cutoff_phi(2,2)=0.4d0
  call build_dg_buffer_window_projector(cutoff_phi,cutoff_buffer,[1d0,1d0],0.25d0, &
    cutoff_coefficients,cutoff_reconstructed,transformed_diagnostics,ok,message)
  if(.not.ok.or.transformed_diagnostics%retained_rank/=1) error stop 'below-cutoff value was retained'

  scaled_phi=phi*1d200
  call build_dg_buffer_window_projector(scaled_phi,buffer,weights,1d-12,scaled_coefficients, &
    scaled_reconstructed,transformed_diagnostics,ok,message)
  if(.not.ok.or.maxval(abs(scaled_reconstructed-reconstructed))>2d-12) &
    error stop 'large basis scaling changed projection'
  scaled_phi=phi*1d-200
  call build_dg_buffer_window_projector(scaled_phi,buffer,weights,1d-12,scaled_coefficients, &
    scaled_reconstructed,transformed_diagnostics,ok,message)
  if(.not.ok.or.maxval(abs(scaled_reconstructed-reconstructed))>2d-12) &
    error stop 'small basis scaling changed projection'
  call build_dg_buffer_window_projector(phi,buffer,weights*1d307,1d-12,scaled_coefficients, &
    scaled_reconstructed,transformed_diagnostics,ok,message)
  if(.not.ok.or.maxval(abs(scaled_reconstructed-reconstructed))>2d-12) &
    error stop 'weight scaling changed projection'
  call build_dg_buffer_window_projector(phi,buffer,[1d0,1d0,1d0,1d0,1d0],1d-12, &
    transformed_coefficients,transformed_reconstructed,transformed_diagnostics,ok,message)
  if(.not.ok) error stop trim(message)
  call build_dg_buffer_window_projector(phi,buffer,[1d308,1d308,1d308,1d308,1d308],1d-12, &
    scaled_coefficients,scaled_reconstructed,transformed_diagnostics,ok,message)
  if(.not.ok.or.maxval(abs(scaled_reconstructed-transformed_reconstructed))>2d-12) &
    error stop 'overflowing weight sum changed projection'
  extreme_phi=1d308
  extreme_buffer=1d308
  extreme_weight=1d308
  call build_dg_buffer_window_projector(extreme_phi,extreme_buffer,extreme_weight,1d-12, &
    extreme_coefficients,extreme_reconstructed,transformed_diagnostics,ok,message)
  if(.not.ok.or.maxval(abs(extreme_reconstructed(:,1)/1d308-1d0))>1d-15.or.&
    abs(transformed_diagnostics%escape_norm)>0d0) &
    error stop 'exact extreme-scale projection failed'

  weights(1)=-1d0
  call expect_failure(phi,buffer,weights,1d-12,coefficients,reconstructed,'weight')
  weights=[1d0,2d0,1.5d0,0.5d0,3d0]
  invalid_phi=phi
  invalid_phi(1,1)=ieee_value(0d0,ieee_quiet_nan)
  call expect_failure(invalid_phi,buffer,weights,1d-12,coefficients,reconstructed,'finite')
  invalid_buffer=buffer
  invalid_buffer(1,1)=ieee_value(0d0,ieee_positive_inf)
  call expect_failure(phi,invalid_buffer,weights,1d-12,coefficients,reconstructed,'finite')
  weights(1)=ieee_value(0d0,ieee_quiet_nan)
  call expect_failure(phi,buffer,weights,1d-12,coefficients,reconstructed,'finite')
  weights=[1d0,2d0,1.5d0,0.5d0,3d0]
  call expect_failure(phi,buffer,weights,0d0,coefficients,reconstructed,'tolerance')
  call expect_failure(phi,buffer,weights,ieee_value(0d0,ieee_quiet_nan),coefficients, &
    reconstructed,'tolerance')
  call expect_failure(phi,buffer,weights,1d-12,coefficients(:,1:1),reconstructed,'dimension')
  allocate(empty_phi(0,3),empty_buffer(0,2),empty_weights(0),empty_coefficients(3,2), &
    empty_reconstructed(0,2))
  call expect_failure(empty_phi,empty_buffer,empty_weights,1d-12,empty_coefficients, &
    empty_reconstructed,'nonempty')
  deallocate(empty_phi,empty_buffer,empty_weights,empty_coefficients,empty_reconstructed)
  allocate(empty_phi(2,0),empty_buffer(2,1),empty_weights(2),empty_coefficients(0,1), &
    empty_reconstructed(2,1))
  empty_weights=1d0
  call expect_failure(empty_phi,empty_buffer,empty_weights,1d-12,empty_coefficients, &
    empty_reconstructed,'nonempty')
  deallocate(empty_phi,empty_buffer,empty_weights,empty_coefficients,empty_reconstructed)
  allocate(empty_phi(2,1),empty_buffer(2,0),empty_weights(2),empty_coefficients(1,0), &
    empty_reconstructed(2,0))
  empty_phi=1d0
  empty_weights=1d0
  call expect_failure(empty_phi,empty_buffer,empty_weights,1d-12,empty_coefficients, &
    empty_reconstructed,'nonempty')

  print '(a)','PASS full-window overlap-metric projector algebra and invariances'
contains
  function solve_weighted_normal(basis,point_weights,rhs) result(solution)
    real(8),intent(in) :: basis(:,:),point_weights(:),rhs(:)
    real(8) :: solution(size(basis,2)),normal(size(basis,2),size(basis,2))
    real(8) :: projected(size(basis,2)),factor
    integer :: i,k,n
    normal=matmul(transpose(basis*spread(point_weights,2,size(basis,2))),basis)
    projected=matmul(transpose(basis*spread(point_weights,2,size(basis,2))),rhs)
    n=size(solution)
    do k=1,n-1
      do i=k+1,n
        factor=normal(i,k)/normal(k,k)
        normal(i,k:n)=normal(i,k:n)-factor*normal(k,k:n)
        projected(i)=projected(i)-factor*projected(k)
      enddo
    enddo
    solution(n)=projected(n)/normal(n,n)
    do i=n-1,1,-1
      solution(i)=(projected(i)-sum(normal(i,i+1:n)*solution(i+1:n)))/normal(i,i)
    enddo
  end function

  real(8) function weighted_relative_residual(values,projection,point_weights)
    real(8),intent(in) :: values(:,:),projection(:,:),point_weights(:)
    weighted_relative_residual=sqrt(sum(spread(point_weights,2,size(values,2))* &
      (values-projection)**2)/sum(spread(point_weights,2,size(values,2))*values**2))
  end function

  real(8) function weighted_escape_norm(values,projection,point_weights)
    real(8),intent(in) :: values(:,:),projection(:,:),point_weights(:)
    weighted_escape_norm=sqrt(sum(spread(point_weights,2,size(values,2))* &
      (values-projection)**2))
  end function

  subroutine assert_invariant(reference,reference_diagnostics,candidate,candidate_diagnostics,label)
    real(8),intent(in) :: reference(:,:),candidate(:,:)
    type(s_dg_buffer_projector_diagnostics),intent(in) :: reference_diagnostics,candidate_diagnostics
    character(*),intent(in) :: label
    if(maxval(abs(reference-candidate))>2d-12) error stop label//': reconstruction'
    if(reference_diagnostics%retained_rank/=candidate_diagnostics%retained_rank) &
      error stop label//': rank'
    if(reference_diagnostics%configured_states/=candidate_diagnostics%configured_states) &
      error stop label//': configured states'
    if(abs(reference_diagnostics%minimum_retained_singular_value-&
      candidate_diagnostics%minimum_retained_singular_value)>2d-12*&
      max(1d0,abs(reference_diagnostics%minimum_retained_singular_value))) &
      error stop label//': minimum singular value'
    if(abs(reference_diagnostics%projection_residual-candidate_diagnostics%projection_residual)>2d-12) &
      error stop label//': residual'
    if(abs(reference_diagnostics%escape_norm-candidate_diagnostics%escape_norm)>2d-12) &
      error stop label//': escape norm'
  end subroutine

  subroutine expect_failure(test_phi,test_buffer,test_weights,test_tolerance,test_coefficients, &
      test_reconstructed,expected_message)
    real(8),intent(in) :: test_phi(:,:),test_buffer(:,:),test_weights(:),test_tolerance
    real(8),intent(out) :: test_coefficients(:,:),test_reconstructed(:,:)
    character(*),intent(in) :: expected_message
    type(s_dg_buffer_projector_diagnostics) :: test_diagnostics
    logical :: test_ok
    character(len=256) :: test_message
    call build_dg_buffer_window_projector(test_phi,test_buffer,test_weights,test_tolerance, &
      test_coefficients,test_reconstructed,test_diagnostics,test_ok,test_message)
    if(test_ok.or.index(test_message,expected_message)==0) error stop 'invalid input was accepted: '//expected_message
  end subroutine
end program
