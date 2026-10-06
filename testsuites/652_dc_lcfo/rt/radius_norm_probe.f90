program radius_norm_probe
 use lcfo_wf_support,only:lcfo_wf_sphere_norm
 implicit none
 real(8) :: p(3,5),c(3,2),length(3),norm(2)
 complex(8) :: w(5,2)
 length=10d0;c=0d0;c(:,2)=5d0
 p(:,1)=[9.5d0,0d0,0d0];p(:,2)=[0d0,9.5d0,0d0];p(:,3)=[0d0,0d0,9.5d0]
 p(:,4)=[0.5d0,0.5d0,0.5d0];p(:,5)=5d0
 w(:,1)=[(1d0,1d0),(2d0,0d0),(0d0,3d0),(4d0,0d0),(5d0,0d0)]
 w(:,2)=2*w(:,1)
 call lcfo_wf_sphere_norm(w,p,c,length,0.5d0,0.5d0,norm)
 if(maxval(abs(norm-[7.5d0,50d0]))>1d-12)error stop '3D periodic or inclusive boundary error'
 call lcfo_wf_sphere_norm(w,p,c,length,0.5d0,0d0,norm)
 if(maxval(abs(norm-[28d0,112d0]))>1d-12)error stop 'Full support error'
 call lcfo_wf_sphere_norm(w(:0,:),p(:,:0),c,length,0.5d0,1d0,norm)
 if(any(norm/=0d0))error stop 'Empty local domain error'
 print *,'Sphere norm: 3D periodic distances, boundary, full and empty core passed'
end program
