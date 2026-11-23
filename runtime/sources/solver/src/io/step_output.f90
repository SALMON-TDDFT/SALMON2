! Copyright 2026 SALMONpy contributors. Apache-2.0.
! Metadata for the audited live RT streams (plus spin for SOC). No duplicated owning handles.
module salmonpy_step_output
  use filesystem,only:resolve_path
  use structures,only:s_ofile
  implicit none
  private
  public::step_output_state,step_output_capture,step_output_close,step_output_reopen,step_output_flush
  type step_output_state
    integer::streams=2
    integer::units(3)=0
    character(256)::names(3)=''
    logical::opened(3)=.false.,closed(3)=.false.
  end type
contains
  subroutine step_output_capture(value,ofl,status)
    use salmon_global,only:yn_spinorbit
    type(step_output_state),intent(out)::value
    type(s_ofile),intent(in)::ofl
    integer,intent(out)::status
    integer::i,ios,n
    character(256)::name,actual
    value%units(1:2)=[ofl%fh_rt,ofl%fh_rt_energy]
    value%names(1:2)=[ofl%file_rt_data,ofl%file_rt_energy_data]
    status=0
    if(yn_spinorbit=='y')then
      n=len_trim(ofl%file_rt_data)
      if(n<8)then
        status=15;return
      endif
      value%streams=3
      value%units(3)=ofl%fh_rt_spin
      value%names(3)=ofl%file_rt_data(:n-5)//'_spin.data'
    endif
    do i=1,value%streams
      name=trim(value%names(i))
      if(name(1:2)=='./')name=name(3:)
      if(len_trim(name)<6.or.index(name,'/')/=0.or.index(name,achar(92))/=0)status=max(status,15)
      if(len_trim(name)>=5)then
        if(name(len_trim(name)-4:len_trim(name))/='.data')status=max(status,15)
      endif
      actual=''
      inquire(unit=value%units(i),opened=value%opened(i),name=actual,iostat=ios)
      if(ios/=0.or..not.value%opened(i))then
        status=90
      else if(trim(actual)/=trim(value%names(i)).and. &
              trim(actual)/=trim(resolve_path(value%names(i))))then
        status=max(status,15)
      endif
    enddo
    if(status==0)call step_output_flush(value,status)
  end subroutine
  subroutine step_output_flush(value,status)
    type(step_output_state),intent(in)::value
    integer,intent(out)::status
    integer::i,ios
    status=0
    do i=1,value%streams
      if(value%opened(i).and..not.value%closed(i))then
        flush(value%units(i),iostat=ios)
        if(ios/=0)status=90
      endif
    enddo
  end subroutine
  subroutine step_output_close(value,status)
    type(step_output_state),intent(inout)::value
    integer,intent(out)::status
    integer::i,ios
    status=0
    do i=1,value%streams
      if(value%opened(i).and..not.value%closed(i))then
        close(value%units(i),iostat=ios)
        if(ios/=0)then
          status=90
        else
          value%closed(i)=.true.
        endif
      endif
    enddo
  end subroutine
  subroutine step_output_reopen(value,status,ofl)
    use filesystem,only:get_filehandle
    type(step_output_state),intent(inout)::value
    type(s_ofile),intent(inout)::ofl
    integer,intent(out)::status
    integer::i,ios,unit
    status=0
    do i=1,value%streams
      if(value%closed(i))then
        ! Another context may have reused the old positive unit while closed.
        ! Select an available unit and update the owning ofile after reconnection.
        unit=get_filehandle()
        open(unit=unit,file=resolve_path(trim(value%names(i))),status='old', &
          form='formatted',access='sequential',action='readwrite',position='append',iostat=ios)
        if(ios/=0)then
          status=90
        else
          value%units(i)=unit
          if(i==1)ofl%fh_rt=unit
          if(i==2)ofl%fh_rt_energy=unit
          if(i==3)ofl%fh_rt_spin=unit
          value%closed(i)=.false.
        endif
      endif
    enddo
  end subroutine
end module
