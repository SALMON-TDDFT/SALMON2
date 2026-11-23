program ion_initialization
 use initialization_sub,only:s_ion_initial_config,ion_input_config,initialize_ions
 use structures,only:s_dft_system
 use salmon_global,only:natom,nelem,Rion,Rion_red,kion,iflag_atom_coor,ntype_atom_coor_reduced
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
 implicit none
 type(s_ion_initial_config)::a,b,bad
 type(s_dft_system)::x,y
 real(8)::saved(3,2),reference(3,2)
 integer::status,i
 natom=2;nelem=2;iflag_atom_coor=ntype_atom_coor_reduced
 allocate(Rion(3,2),Rion_red(3,2),kion(2))
 Rion=999d0;Rion_red(:,1)=[0.25d0,0.5d0,0.75d0];Rion_red(:,2)=[0.125d0,0.25d0,0.5d0];kion=[1,2]
 call ion_input_config(a)
 b=a;b%reduced=.false.;b%positions=-7d0;b%species=[2,1]
 x%primitive_a=0d0;x%primitive_a(:,1)=[8d0,1d0,0d0]
 x%primitive_a(:,2)=[0d0,12d0,0d0];x%primitive_a(:,3)=[0d0,0d0,16d0]
 y%primitive_a=0d0
 reference(:,1)=[2d0,6.25d0,12d0];reference(:,2)=[1d0,3.125d0,8d0]
 ! Captured data and explicit initialization survive unrelated global inputs.
 natom=99;nelem=99;Rion=-123d0;Rion_red=-456d0;kion=99;iflag_atom_coor=-99
 do i=1,3
  call initialize_ions(x,a,status)
  if(status/=0.or.x%nion/=2.or.any(x%Rion/=reference).or.any(x%kion/=[1,2]))stop 1
  if(size(x%mass)/=2.or.any(x%Velocity/=0d0))stop 2
  saved=x%Rion
  call initialize_ions(y,b,status)
  if(status/=0.or.any(y%Rion/=-7d0).or.any(y%kion/=[2,1]).or.any(x%Rion/=saved))stop 3
 enddo
 if(any(Rion/=-123d0).or.any(Rion_red/=-456d0).or.any(kion/=99).or.natom/=99)stop 4
 bad=a;bad%positions(1,1)=ieee_value(0d0,ieee_quiet_nan)
 call reject(bad)
 bad=a;bad%species(1)=3
 call reject(bad)
 bad=a;bad%nion=3
 call reject(bad)
 bad=a;bad%nelem=-1
 call reject(bad)
 bad=a;deallocate(bad%positions)
 call reject(bad)
 bad=a;deallocate(bad%positions);allocate(bad%positions(0:2,0:1));bad%positions=0d0
 call reject(bad)
 bad=a;deallocate(bad%species);allocate(bad%species(0:1));bad%species=1
 call reject(bad)
 print *,'Owned ion initialization PASS' 
contains
 subroutine reject(value)
 type(s_ion_initial_config),intent(in)::value
 call initialize_ions(x,value,status)
 if(status==0.or.x%nion/=2.or.any(x%Rion/=saved).or.any(x%kion/=[1,2]))stop 5
 if(size(x%mass)/=2.or.any(x%Velocity/=0d0))stop 6
 end subroutine
end program
