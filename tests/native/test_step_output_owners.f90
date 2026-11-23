! Copyright 2026 SALMONpy contributors. Apache-2.0.
program step_output_owners
 use structures,only:s_ofile
 use filesystem,only:open_filehandle
 use salmonpy_step_output
 implicit none
 type(s_ofile)::a,b
 type(step_output_state)::saved
 integer::status
 logical::connected
 character(256)::actual
 a%file_rt_data='a_rt.data';a%file_rt_energy_data='a_rt_energy.data'
 a%fh_rt=open_filehandle(a%file_rt_data,'new')
 a%fh_rt_energy=open_filehandle(a%file_rt_energy_data,'new')
 call step_output_capture(saved,a,status);if(status/=0)stop 1
 call step_output_close(saved,status);if(status/=0)stop 2
 b%file_rt_data='b_rt.data';b%file_rt_energy_data='b_rt_energy.data'
 b%fh_rt=open_filehandle(b%file_rt_data,'new')
 b%fh_rt_energy=open_filehandle(b%file_rt_energy_data,'new')
 if(b%fh_rt/=saved%units(1))stop 3 ! Exercise the actual positive-unit reuse.
 call step_output_reopen(saved,status,a)
 if(status/=0)stop 4
 if(any(saved%units/=[a%fh_rt,a%fh_rt_energy]))stop 8
 inquire(unit=b%fh_rt,opened=connected,name=actual)
 if(.not.connected.or.trim(actual)/=trim(b%file_rt_data))stop 5
 inquire(unit=b%fh_rt_energy,opened=connected,name=actual)
 if(.not.connected.or.trim(actual)/=trim(b%file_rt_energy_data))stop 6
 call step_output_close(saved,status);if(status/=0)stop 7
 write(b%fh_rt,*)'B still owns its RT stream'
 write(b%fh_rt_energy,*)'B still owns its energy stream'
 close(b%fh_rt);close(b%fh_rt_energy)
 print *,'Step output unit reuse PASS'
end program
