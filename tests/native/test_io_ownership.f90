! Synthetic fixtures only: no GS wavefunction/density data or SCF calculation.
program test_io_ownership
 use structures
 use salmon_global
 use inputoutput,only:initialize_inputoutput_units,uenergy_from_au,ulength_from_au, &
   utime_from_au,t_unit_energy,au_energy_ev,au_length_aa,au_time_fs,nsplit_voxel_data
 use salmonpy_model_inputs,only:s_model_inputs
 use parallelization,only:setup_parallel,end_parallel,nproc_id_global,nproc_group_global
 use math_constants,only:set_math_constants
 use write_sub,only:write_dos,write_pdos
 use band_dft_sub,only:write_band,calc_band_write
 use write_file3d,only:write_avs
 implicit none
 type(s_model_inputs)::a,b,bad
 type(s_dft_system)::system
 type(s_dft_energy)::energy
 type(s_rgrid)::lg,mg
 type(s_parallel_info)::info
 type(s_pp_info)::pp
 type(s_orbital)::psi
 type(s_band_dft)::band
 integer::status,j,unit,i,host(4)=[1,77,100,101]
 character(32)::mode,variant
 character(80)::names(4),actual
 character(60)::suffix
 character(20)::header
 logical::opened,reference
 real(8)::field(2,1,1)
 call setup_parallel
 call set_math_constants
 call get_command_argument(1,mode)
 call get_command_argument(2,variant)
 reference=variant=='reference'
 unit_system='au';unit_time='au';unit_length='au';unit_energy='au';unit_charge='au'
 call initialize_inputoutput_units
 call a%capture()
 unit_time='fs';unit_length='AA';unit_energy='eV'
 call initialize_inputoutput_units
 call b%capture()
 if(mode=='units')then
   call a%activate(status)
   if(status/=0)stop 11
   if(uenergy_from_au/=1d0.or.ulength_from_au/=1d0.or.utime_from_au/=1d0)stop 12
   if(t_unit_energy%name/='a.u.')stop 13
   call b%activate(status)
   if(status/=0)stop 14
   if(uenergy_from_au/=au_energy_ev.or.ulength_from_au/=au_length_aa.or.utime_from_au/=au_time_fs)stop 15
   if(t_unit_energy%name/='eV')stop 16
   call a%activate(status)
   if(status/=0.or.uenergy_from_au/=1d0)stop 17
   unit_energy='invalid';call bad%capture()
   unit_energy='au'
   call bad%activate(status)
   if(status==0.or.unit_energy/='au'.or.uenergy_from_au/=1d0)stop 19
 else
   ! Exercise byte-preserving output, separately from model-switch conversion.
   unit_time='au';unit_length='au';unit_energy='au';unit_charge='au'
   call initialize_inputoutput_units
   if(.not.reference)then
     do j=1,4
       write(names(j),'(a,i0,a,i0)')'host-',nproc_id_global,'-',host(j)
       open(unit=host(j),file=trim(names(j)),status='replace')
       write(host(j),'(a)')'before'
     enddo
   endif
   system%nspin=1;system%nk=1;system%no=2;system%Hvol=.5d0
   allocate(system%wtk(1),system%vec_k(3,1),system%Rion(3,1))
   system%wtk=1d0;system%vec_k=0d0;system%Rion=0d0
   system%primitive_b=0d0
   do j=1,3
     system%primitive_b(j,j)=1d0
   enddo
   allocate(energy%esp(2,1,1));energy%esp(:,1,1)=[.25d0,.75d0]
   base_directory='./';sysname='synthetic';natom=1;nelec=2;nstate=2;yn_spinorbit='n'
   temperature=-1d0;yn_out_dos_set_fe_origin='n';out_dos_function='gaussian'
   out_dos_start=0d0;out_dos_end=1d0;out_dos_width=.1d0;out_dos_nenergy=5
   select case(trim(mode))
   case('dos')
     call write_dos(system,energy)
   case('band')
     allocate(band%band_kpt(3,1),band%check_conv_esp(2,1,1));band%band_kpt=0d0
     info%ik_s=1;info%ik_e=1
     call calc_band_write(1,system,band,info)
     call write_band(system,energy)
   case('avs','avs-unsplit')
     nsplit_voxel_data=2
     if(mode=='avs-unsplit')nsplit_voxel_data=1
     lg%is=[1,1,1];lg%ie=[2,1,1];lg%num=[2,1,1]
     field(:,1,1)=[.125d0,.25d0];suffix='synthetic-field';header='fixture'
     call write_avs(lg,77,suffix,header,field)
   case('pdos')
     allocate(kion(1));kion=1
     lg%is=1;lg%ie=1;lg%num=1;mg%is=1;mg%ie=1;mg%num=1
     allocate(lg%coordinate(1,3));lg%coordinate(1,:)=[1d0,0d0,0d0]
     pp%nrmax=3;allocate(pp%mlps(1),pp%upp_f(3,0:0,1),pp%rad(4,1))
     pp%mlps=0;pp%upp_f=1d0;pp%rad(:,1)=[0d0,2d0,4d0,6d0]
     allocate(psi%rwf(1,1,1,1,2,1,1));psi%rwf=1d0
     info%ik_s=1;info%ik_e=1;info%io_s=1;info%io_e=2
     info%icomm_r=nproc_group_global;info%icomm_ko=nproc_group_global
     call write_pdos(lg,mg,system,info,pp,energy,psi)
   case default
     stop 18
   end select
   if(.not.reference)then
     do j=1,4
       inquire(unit=host(j),opened=opened,name=actual)
       if(.not.opened)stop 21
       if(trim(actual)/=trim(names(j)))stop 22
       write(host(j),'(a)')'after'
       close(host(j))
     enddo
   endif
 endif
 call end_parallel
 print *, 'IO ownership PASS ',trim(mode)
end program
