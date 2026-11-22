program test_projection_capture
  use hsex_projection, only: s_hsex_capture,projection_capture,projection_seed
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  complex(8) :: old(8,2),rt(8,2),grown(8,6),overlap(6,2,1),saved(6,2,1),gram(6,6)
  real(8) :: norms(2,1),occupations(2,1),weights(1)
  type(s_hsex_capture) :: capture
  integer :: i,status
  old=0d0;old(1,1)=1d0;old(2,2)=1d0
  rt=0d0;rt(1,1)=sqrt(.7d0);rt(3,1)=sqrt(.3d0)
  rt(2,2)=sqrt(.9d0);rt(4,2)=sqrt(.1d0)
  norms=1d0;occupations=2d0;weights=1d0
  overlap=0d0;overlap(:2,:,1)=matmul(conjg(transpose(old)),rt);saved=overlap
  call projection_capture(overlap(:2,:,:),norms,occupations,weights,2,1d-8,capture,status)
  if(status/=0)error stop 'valid incomplete projection rejected'
  if(abs(capture%missing_norm-.3d0)>1d-14)error stop 'missing orbital norm incorrect'
  if(abs(capture%lost_electrons-.8d0)>1d-14)error stop 'missing electrons incorrect'
  if(any(overlap/=saved))error stop 'capture check modified occupations'
  call projection_seed(old,rt,grown,status)
  if(status/=0)error stop 'auxiliary expansion rejected'
  if(any(grown(:,:2)/=old))error stop 'existing auxiliaries changed'
  gram=matmul(conjg(transpose(grown)),grown)
  do i=1,6
    gram(i,i)=gram(i,i)-1d0
  enddo
  if(maxval(abs(gram))>1d-12)error stop 'expanded auxiliary orbitals not orthonormal'
  overlap(:,:,1)=matmul(conjg(transpose(grown)),rt);saved=overlap
  call projection_capture(overlap,norms,occupations,weights,2,1d-8,capture,status)
  if(status/=0.or.capture%missing_norm>1d-12)error stop 'expanded space loses RT norm'
  if(abs(capture%excited_electrons-.8d0)>1d-12.or.abs(capture%holes-.8d0)>1d-12) &
    error stop 'nonthermal electron-hole populations changed'
  if(any(overlap/=saved))error stop 'expanded populations renormalized'
  overlap(6,2,1)=ieee_value(0d0,ieee_quiet_nan)
  call projection_capture(overlap,norms,occupations,weights,2,1d-8,capture,status)
  if(status==0)error stop 'NaN capture accepted'
  overlap=saved;norms(2,1)=ieee_value(0d0,ieee_positive_inf)
  call projection_capture(overlap,norms,occupations,weights,2,1d-8,capture,status)
  if(status==0)error stop 'infinite RT norm accepted'
  norms=1d0;overlap=saved;overlap(1,1,1)=2d0
  call projection_capture(overlap,norms,occupations,weights,2,1d-8,capture,status)
  if(status==0)error stop 'overcapture accepted'
  old(8,2)=ieee_value(0d0,ieee_positive_inf)
  call projection_seed(old,rt,grown,status)
  if(status==0)error stop 'infinite seed accepted'
  old=0d0;old(1,:)=1d0
  call projection_seed(old,rt,grown,status)
  if(status==0)error stop 'dependent old auxiliaries accepted'
  print *, 'NONTHERMAL_CAPTURE_EXPANSION_IEEE_PASS'
end program test_projection_capture
