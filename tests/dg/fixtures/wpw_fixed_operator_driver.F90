program wpw_fixed_operator_driver
  use dg_wpw_fixed_operator
  implicit none
  integer,parameter::nw=2,np=1,n=nw+np
  type(s_dg_wpw_fixed_operator)::op
  complex(8)::sww(nw,nw),swp(nw,np),spp(np,np)
  complex(8)::hww(nw,nw),hwp(nw,np),hpp(np,np),hfaceww(nw,nw),hfacewp(nw,np)
  complex(8)::v(n,n),h(n,n),href(n,n)
  integer::info,nwgot,npgot,ngot
  sww=0;sww(1,1)=1;sww(2,2)=1;swp=reshape([cmplx(.1d0,.05d0,8),cmplx(-.08d0,.03d0,8)],[nw,np]);spp=1
  hww=0;hww(1,1)=-.7d0;hww(2,2)=.2d0;hww(1,2)=cmplx(.04d0,.02d0,8);hww(2,1)=conjg(hww(1,2))
  hwp=reshape([cmplx(.12d0,-.03d0,8),cmplx(-.05d0,.01d0,8)],[nw,np]);hpp=reshape([cmplx(1.3d0,0d0,8)],[np,np])
  hfaceww=0;hfaceww(1,1)=.03d0;hfaceww(2,2)=.04d0
  hfacewp=reshape([cmplx(.01d0,.02d0,8),cmplx(-.02d0,.01d0,8)],[nw,np])
  call initialize_dg_wpw_fixed_operator(op,sww,swp,spp,hww,hwp,hpp,hfaceww,hfacewp,nw,np,info)
  call dg_wpw_fixed_dimensions(op,nwgot,npgot,ngot);call copy_dg_wpw_metric(op,h,info);call copy_dg_wpw_fixed_hamiltonian(op,href,info)
  call req(info==0.and.nwgot==nw.and.npgot==np.and.ngot==n,'initialize')
  call req(maxval(abs(h(1:nw,nw+1:n)-swp))<1d-14,'WP overlap')
  call req(maxval(abs(href(1:nw,1:nw)-(hww+hfaceww)))<1d-14,'WW fixed face')
  call req(maxval(abs(href(1:nw,nw+1:n)-(hwp+hfacewp)))<1d-14,'WP fixed face')
  call req(maxval(abs(href(nw+1:n,nw+1:n)-hpp))<1d-14,'PP has no face term')
  v=0;v(1,1)=.2d0;v(2,2)=.1d0;v(3,3)=-.05d0
  call compose_dg_wpw_hamiltonian(op,v,h,info)
  call req(info==0.and.maxval(abs(h-(href+v)))<1d-14,'potential composition')
  v=2*v;call compose_dg_wpw_hamiltonian(op,v,h,info)
  call copy_dg_wpw_fixed_hamiltonian(op,h,info);call req(maxval(abs(h-href))<1d-14,'fixed block immutable')
  hfaceww(1,2)=cmplx(.3d0,.2d0,8)
  call initialize_dg_wpw_fixed_operator(op,sww,swp,spp,hww,hwp,hpp,hfaceww,hfacewp,nw,np,info)
  call req(info/=0,'non-Hermitian rejection')
  write(*,'(a)')'PASS WPW fixed operator'
contains
  subroutine req(ok,msg);logical,intent(in)::ok;character(*),intent(in)::msg
    if(.not.ok)then;write(*,*)trim(msg);error stop 1;endif
  end subroutine
end program
