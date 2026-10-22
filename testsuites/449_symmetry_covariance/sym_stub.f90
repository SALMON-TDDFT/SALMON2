module sym_sub
 real(8)::Amat(3,3),Ainv(3,3)
 real(8),allocatable::SymMatA(:,:,:),SymMatB(:,:,:)
 logical::use_symmetry=.true.
end module
