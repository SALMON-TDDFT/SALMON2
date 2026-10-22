! Exercise the production reader, normalization, masking and MPI broadcasts.
program ky_projector_probe
  use structures, only: s_pp_info
  use salmon_pp, only: init_pp
  use input_pp_sub, only: input_pp
  use salmon_global, only: nelem,izatom,file_pseudo,lloc_ps,lmax_ps,yn_psmask,yn_spinorbit, &
                          alpha_mask,gamma_mask,eta_mask,method_init_density,quiet,base_directory
  use parallelization, only: setup_parallel,end_parallel,nproc_id_global
  use math_constants, only: set_math_constants
  implicit none
  type(s_pp_info) :: pp
  character(256) :: filename,mask_mode,soc_mode,output
  integer :: i,l,unit
  call setup_parallel
  call set_math_constants
  call get_command_argument(1,filename)
  call get_command_argument(2,mask_mode)
  call get_command_argument(3,soc_mode)
  nelem=1;izatom=14;file_pseudo(1)=filename;lloc_ps=2;lmax_ps=-1
  yn_psmask=mask_mode(1:1);yn_spinorbit=soc_mode(1:1)
  alpha_mask=0.8d0;gamma_mask=1.8d0;eta_mask=15d0
  method_init_density='wf';quiet=.true.;base_directory='./'
  call init_pp(pp,3000,8,.false.)
  call input_pp(pp,0.5d0,0.5d0,0.5d0)
  write(output,'(a,i0)')'tables.rank',nproc_id_global
  open(newunit=unit,file=trim(output),status='replace')
  write(unit,*)pp%nrps(1),pp%mlps(1),pp%rps(1)
  do l=0,pp%mlps(1)
    write(unit,*)l,pp%inorm(l,1),pp%inorm_so(l,1),pp%anorm(l,1),pp%anorm_so(l,1)
    do i=1,pp%nrps(1)
      write(unit,'(5es25.16)')pp%radnl(i,1),pp%udvtbl(i,l,1),pp%udvtbl_so(i,l,1), &
                            pp%dudvtbl(i,l,1),pp%dudvtbl_so(i,l,1)
    enddo
  enddo
  close(unit)
  call end_parallel
end program ky_projector_probe
