program probe
 use salmon_global, only: xc,hse_omega,exx_mlwf_radius,theory,yn_conventional_from_dcdft,yn_restart
 use exx_functional, only: lcfo_write_functional,lcfo_check_functional
 implicit none
 integer :: status
 xc='hse06';hse_omega=.11d0;exx_mlwf_radius=0d0
 theory='dft';yn_conventional_from_dcdft='n';yn_restart='n'
 call lcfo_write_functional('metadata-probe.txt','probe',status)
 if(status/=0)error stop 'write'
 theory='tddft_pulse';yn_conventional_from_dcdft='y';exx_mlwf_radius=10d0
 call lcfo_check_functional('metadata-probe.txt','probe',status)
 if(status/=0)error stop 'fresh RT radius rejected'
 hse_omega=.12d0
 call lcfo_check_functional('metadata-probe.txt','probe',status)
 if(status==0)error stop 'changed kernel accepted'
 hse_omega=.11d0;yn_restart='y'
 call lcfo_check_functional('metadata-probe.txt','probe',status)
 if(status==0)error stop 'changed restart support accepted'
 yn_restart='n';theory='dft'
 call lcfo_check_functional('metadata-probe.txt','probe',status)
 if(status==0)error stop 'changed GS support accepted'
 exx_mlwf_radius=2d0
 call lcfo_write_functional('metadata-probe.txt','probe',status)
 theory='tddft_pulse';exx_mlwf_radius=10d0
 call lcfo_check_functional('metadata-probe.txt','probe',status)
 if(status==0)error stop 'previously truncated GS accepted'
 print *, 'PASS fresh RT radius and unchanged metadata rejection guards'
end program
