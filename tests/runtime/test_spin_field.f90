program test_spin_field
 use noncollinear_module,only:spin_field_matrix,spin_field_density_energy
 implicit none
 complex(8)::h(2,2),rho(1,1,1,2,2),psi(2),expected(2,2)
 real(8)::b(3),plus(3),minus(3),step,derivative,e
 integer::i,j,k
 psi=[cmplx(.3d0,.4d0,8),cmplx(.5d0,-sqrt(.5d0),8)]
 psi=psi/sqrt(real(dot_product(psi,psi),8))
 do j=1,2
 do i=1,2
  rho(1,1,1,i,j)=psi(i)*conjg(psi(j))
 enddo
 enddo
 b=[.01d0,-.02d0,.03d0]
 call spin_field_matrix(b,h)
 if(maxval(abs(h-transpose(conjg(h))))>1d-15)error stop 'field matrix not Hermitian'
 e=spin_field_density_energy(b,rho,1d0)
 if(abs(e-real(dot_product(psi,matmul(h,psi)),8))>1d-15)error stop 'Hamiltonian energy mismatch'
 step=1d-6
 do k=1,3
  plus=b;minus=b;plus(k)=plus(k)+step;minus(k)=minus(k)-step
  derivative=(spin_field_density_energy(plus,rho,1d0)-spin_field_density_energy(minus,rho,1d0))/(2*step)
  plus=0d0;plus(k)=1d0;call spin_field_matrix(plus,expected)
  if(abs(derivative-real(dot_product(psi,matmul(expected,psi)),8))>1d-11)error stop 'energy derivative mismatch'
 enddo
 call spin_field_matrix([0d0,0d0,0d0],h)
 if(any(h/=(0d0,0d0)))error stop 'zero field changed matrix'
 print *, 'SPIN_FIELD_MATRIX_ENERGY_PASS'
end program
