program probe
 use builtin_pz_sp,only:exc_cor_pz_sp
 use builtin_pz,only:exc_cor_pz
 use structures
 use salmon_xc,only:calc_xc,salmon_xctype_pz,ensure_xc_workspace
 use code_optimization,only:set_modulo_tables,configure_openmp_parallelization
 implicit none
 type(s_dft_system)::a,b
 type(s_stencil)::sa,sb
 type(s_xc_functional)::xc
 type(s_pp_info)::pp
 real(8)::ra(2,3,4),rb(3,4,5),va(2,3,4),vb(3,4,5),refa(2,3,4),refb(3,4,5)
 real(8)::e1(24),e2(24),v1(24),vs(24,2),rs(2,3,4,2),result_s(2,3,4,2)
 integer::i
 call configure_openmp_parallelization([8,8,8],sa,'orbital','orbital','orbital')
 call configure_openmp_parallelization([8,8,8],sb,'rgrid','rgrid','rgrid')
 if(sa%stencil_is_parallelized_by_omp.or..not.sa%current_omp_mode.or..not.sa%force_omp_mode)stop 20
 if(.not.sb%stencil_is_parallelized_by_omp.or.sb%current_omp_mode.or.sb%force_omp_mode)stop 21
 call configure_openmp_parallelization([32,32,32],sb,'auto','auto','auto')
 if(.not.sb%stencil_is_parallelized_by_omp.or.sb%current_omp_mode)stop 22
 if(sa%stencil_is_parallelized_by_omp.or..not.sa%current_omp_mode)stop 23
 xc%xctype=[salmon_xctype_pz,0,0];xc%ispin=0
 ra=.2d0;rb=.7d0
 call exc_cor_pz(24,reshape(ra,[24])*.5d0,e1,e2,v1)
 refa=reshape(v1,[2,3,4])
 call calc_xc(xc,pp,rho=rb,vxc=refb)
 call ensure_xc_workspace(a%xc_work,[2,3,4],1)
 call ensure_xc_workspace(b%xc_work,[3,4,5],1)
 a%xc_work%rho_tmp=ra;b%xc_work%rho_tmp=rb
 call calc_xc(xc,pp,rho=a%xc_work%rho_tmp,vxc=va,workspace=a%xc_work%pz)
 call calc_xc(xc,pp,rho=b%xc_work%rho_tmp,vxc=vb,workspace=b%xc_work%pz)
 if(any(va/=refa).or.any(vb/=refb))stop 1
 if(size(a%xc_work%pz%exc_1d)/=24.or.size(b%xc_work%pz%exc_1d)/=60)stop 2
 call set_modulo_tables([8,12,16],sa)
 call set_modulo_tables([12,8,20],sb)
 do i=0,23
  if(sa%modx(i)/=mod(i,8).or.sb%modx(i)/=mod(i,12))stop 3
 enddo
 if(any([lbound(sa%modx,1),lbound(sa%mody,1),lbound(sa%modz,1)]/=0))stop 8
 do i=0,size(sa%mody)-1
  if(sa%mody(i)/=mod(i,12))stop 9
 enddo
 do i=0,size(sb%modz)-1
  if(sb%modz(i)/=mod(i,20))stop 10
 enddo
 ! Changing B's shape/spin must not change A's contents or allocation bounds.
 call ensure_xc_workspace(b%xc_work,[5,3,2],2)
 if(any(shape(b%xc_work%rho_s_tmp)/=[5,3,2,2]))stop 4
 if(any(a%xc_work%rho_tmp/=ra).or.size(a%xc_work%pz%exc_1d)/=24)stop 5
 b%xc_work=s_xc_workspace()
 deallocate(sb%modx,sb%mody,sb%modz)
 call calc_xc(xc,pp,rho=a%xc_work%rho_tmp,vxc=va,workspace=a%xc_work%pz)
 if(any(va/=refa).or.size(sa%modx)/=24)stop 6
 ! Reusing the same PZ workspace at another size must also resize safely.
 call calc_xc(xc,pp,rho=rb,vxc=vb,workspace=a%xc_work%pz)
 if(any(vb/=refb).or.size(a%xc_work%pz%exc_1d)/=60)stop 7
 xc%ispin=1
 rs(:,:,:,1)=.2d0;rs(:,:,:,2)=.4d0
 call exc_cor_pz_sp(24,reshape(rs,[24,2]),e1,e2,vs)
 call calc_xc(xc,pp,rho_s=rs,vxc_s=result_s,workspace=a%xc_work%pz)
 if(any(result_s/=reshape(vs,[2,3,4,2])))stop 11
 if(a%xc_work%pz%spin/=1.or.size(a%xc_work%pz%exc_1d)/=24)stop 12
 print *,'Owned XC/stencil workspace probe PASS'
end program
