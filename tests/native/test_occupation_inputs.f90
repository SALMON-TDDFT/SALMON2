program test_occupation_inputs
 use structures
 use occupation,only:occupation_input_config,initialize_occupations,ne2mu
 use reference_occupation,only:reference_ne2mu=>ne2mu
 use salmon_global
 use parallelization
 implicit none
 type(s_dft_system)::a,b,ref
 type(s_dft_energy)::energy
 real(8)::saved(4,2,1),mu_a,temp_b
 integer::i,scenario
 call setup_parallel
 a%nspin=1;a%nk=2;a%no=4
 allocate(a%rocc(4,2,1),a%wtk(2),energy%esp(4,2,1))
 a%wtk=[0.25d0,0.75d0]
 do i=1,4
  energy%esp(i,:,1)=[0.1d0*i,0.12d0*i]
 enddo
 b=a;ref=a
 do scenario=1,2
  nelec=3;nelec_spin=0;nstate=4;spin='unpolarized';yn_spinorbit='n'
  temperature=0.08d0;nscf_init_redistribution=2
  call occupation_input_config(a%occupation_input);a%occupation_input_ready=.true.
  temp_b=0.2d0
  if(scenario==2)temp_b=0d0
  nelec=4;temperature=temp_b;nscf_init_redistribution=7
  call occupation_input_config(b%occupation_input);b%occupation_input_ready=.true.
  nelec=1;nelec_spin=99;temperature=9d0;nstate=1;spin='polarized';yn_spinorbit='y'
  nscf_init_redistribution=99
  if(a%occupation_input%nstate/=4.or.a%occupation_input%redistribution_delay/=2)stop 9
  if(b%occupation_input%redistribution_delay/=7)stop 10
  call initialize_occupations(a)
  if(any(a%rocc(1,:,1)/=2d0).or.any(a%rocc(2,:,1)/=1d0))stop 1
  if(any(a%rocc(3:4,:,1)/=0d0))stop 2
  call initialize_occupations(b)
  if(any(b%rocc(1:2,:,1)/=2d0).or.any(b%rocc(3:4,:,1)/=0d0))stop 3
  call ne2mu(energy,a,0)
  saved=a%rocc;mu_a=a%mu
  call ne2mu(energy,b,0)
  if(any(a%rocc/=saved).or.a%mu/=mu_a)stop 4
  nelec=3;nelec_spin=0;temperature=0.08d0;yn_spinorbit='n'
  call reference_ne2mu(energy,ref,0)
  if(any(a%rocc/=ref%rocc).or.a%mu/=ref%mu)stop 5
  nelec=4;temperature=temp_b
  call reference_ne2mu(energy,ref,0)
  if(any(b%rocc/=ref%rocc).or.b%mu/=ref%mu)stop 6
  if(abs(sum(a%rocc(:,1,1))*0.25d0+sum(a%rocc(:,2,1))*0.75d0-3d0)>1d-9)stop 7
  if(abs(sum(b%rocc(:,1,1))*0.25d0+sum(b%rocc(:,2,1))*0.75d0-4d0)>1d-9)stop 8
 enddo
 call end_parallel
 print *,'Occupation input isolation PASS'
end program
