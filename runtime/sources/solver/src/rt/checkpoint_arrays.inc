! Shared bounded numerical-array stream codec. Caller supplies checkpoint_io_limit().
subroutine checkpoint_io_s_scalar(unit,reading,value,ios)
  integer,intent(in)::unit
  logical,intent(in)::reading
  type(s_scalar),intent(inout)::value
  integer,intent(inout)::ios
  logical::present_array
  integer::lo(7),hi(7),i1,i2,k,allocation_status
  integer(8)::elements,extent
  if(ios/=0)return
  present_array=allocated(value%f)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%f))deallocate(value%f)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%f))then
        if(any(lo(1:3)/=lbound(value%f)).or.any(hi(1:3)/=ubound(value%f)))then
          ios=1;return
        endif
      else
        allocate(value%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%f);hi(1:3)=ubound(value%f)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%f
    else
      write(unit,iostat=ios) value%f
    endif
    if(ios/=0)return
  endif
end subroutine
subroutine checkpoint_io_s_vector(unit,reading,value,ios)
  integer,intent(in)::unit
  logical,intent(in)::reading
  type(s_vector),intent(inout)::value
  integer,intent(inout)::ios
  logical::present_array
  integer::lo(7),hi(7),i1,i2,k,allocation_status
  integer(8)::elements,extent
  if(ios/=0)return
  present_array=allocated(value%v)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%v))deallocate(value%v)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:4),hi(1:4)
      if(ios/=0)return
      elements=1_8
      do k=1,4
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%v))then
        if(any(lo(1:4)/=lbound(value%v)).or.any(hi(1:4)/=ubound(value%v)))then
          ios=1;return
        endif
      else
        allocate(value%v(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),lo(4):hi(4)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:4)=lbound(value%v);hi(1:4)=ubound(value%v)
      write(unit,iostat=ios)lo(1:4),hi(1:4)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%v
    else
      write(unit,iostat=ios) value%v
    endif
    if(ios/=0)return
  endif
end subroutine
subroutine checkpoint_io_s_orbital(unit,reading,value,ios)
  integer,intent(in)::unit
  logical,intent(in)::reading
  type(s_orbital),intent(inout)::value
  integer,intent(inout)::ios
  logical::present_array
  integer::lo(7),hi(7),i1,i2,k,allocation_status
  integer(8)::elements,extent
  if(ios/=0)return
  present_array=allocated(value%rwf)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%rwf))deallocate(value%rwf)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:7),hi(1:7)
      if(ios/=0)return
      elements=1_8
      do k=1,7
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%rwf))then
        if(any(lo(1:7)/=lbound(value%rwf)).or.any(hi(1:7)/=ubound(value%rwf)))then
          ios=1;return
        endif
      else
        allocate(value%rwf(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),lo(4):hi(4), &
          lo(5):hi(5),lo(6):hi(6),lo(7):hi(7)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:7)=lbound(value%rwf);hi(1:7)=ubound(value%rwf)
      write(unit,iostat=ios)lo(1:7),hi(1:7)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%rwf
    else
      write(unit,iostat=ios) value%rwf
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%zwf)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%zwf))deallocate(value%zwf)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:7),hi(1:7)
      if(ios/=0)return
      elements=1_8
      do k=1,7
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%zwf))then
        if(any(lo(1:7)/=lbound(value%zwf)).or.any(hi(1:7)/=ubound(value%zwf)))then
          ios=1;return
        endif
      else
        allocate(value%zwf(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),lo(4):hi(4), &
          lo(5):hi(5),lo(6):hi(6),lo(7):hi(7)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:7)=lbound(value%zwf);hi(1:7)=ubound(value%zwf)
      write(unit,iostat=ios)lo(1:7),hi(1:7)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%zwf
    else
      write(unit,iostat=ios) value%zwf
    endif
    if(ios/=0)return
  endif
  if(reading)then
    read(unit,iostat=ios) value%update_zwf_overlap
  else
    write(unit,iostat=ios) value%update_zwf_overlap
  endif
  if(ios/=0)return
end subroutine
subroutine checkpoint_io_s_dft_energy(unit,reading,value,ios)
  integer,intent(in)::unit
  logical,intent(in)::reading
  type(s_dft_energy),intent(inout)::value
  integer,intent(inout)::ios
  logical::present_array
  integer::lo(7),hi(7),i1,i2,k,allocation_status
  integer(8)::elements,extent
  if(ios/=0)return
  present_array=allocated(value%esp)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%esp))deallocate(value%esp)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%esp))then
        if(any(lo(1:3)/=lbound(value%esp)).or.any(hi(1:3)/=ubound(value%esp)))then
          ios=1;return
        endif
      else
        allocate(value%esp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%esp);hi(1:3)=ubound(value%esp)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%esp
    else
      write(unit,iostat=ios) value%esp
    endif
    if(ios/=0)return
  endif
  if(reading)then
    read(unit,iostat=ios) value%E_tot
  else
    write(unit,iostat=ios) value%E_tot
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_kin
  else
    write(unit,iostat=ios) value%E_kin
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_h
  else
    write(unit,iostat=ios) value%E_h
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_xc
  else
    write(unit,iostat=ios) value%E_xc
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_ion_ion
  else
    write(unit,iostat=ios) value%E_ion_ion
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_ion_loc
  else
    write(unit,iostat=ios) value%E_ion_loc
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_ion_nloc
  else
    write(unit,iostat=ios) value%E_ion_nloc
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_U
  else
    write(unit,iostat=ios) value%E_U
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_tot0
  else
    write(unit,iostat=ios) value%E_tot0
  endif
  if(ios/=0)return
end subroutine
subroutine checkpoint_io_s_poisson(unit,reading,value,ios)
  integer,intent(in)::unit
  logical,intent(in)::reading
  type(s_poisson),intent(inout)::value
  integer,intent(inout)::ios
  logical::present_array
  integer::lo(7),hi(7),i1,i2,k,allocation_status
  integer(8)::elements,extent
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%iterVh
  else
    write(unit,iostat=ios) value%iterVh
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%npole_partial
  else
    write(unit,iostat=ios) value%npole_partial
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%npole_total
  else
    write(unit,iostat=ios) value%npole_total
  endif
  if(ios/=0)return
  present_array=allocated(value%ipole_tbl)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ipole_tbl))deallocate(value%ipole_tbl)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ipole_tbl))then
        if(any(lo(1:1)/=lbound(value%ipole_tbl)).or.any(hi(1:1)/=ubound(value%ipole_tbl)))then
          ios=1;return
        endif
      else
        allocate(value%ipole_tbl(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%ipole_tbl);hi(1:1)=ubound(value%ipole_tbl)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ipole_tbl
    else
      write(unit,iostat=ios) value%ipole_tbl
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ig_num)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ig_num))deallocate(value%ig_num)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ig_num))then
        if(any(lo(1:1)/=lbound(value%ig_num)).or.any(hi(1:1)/=ubound(value%ig_num)))then
          ios=1;return
        endif
      else
        allocate(value%ig_num(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%ig_num);hi(1:1)=ubound(value%ig_num)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ig_num
    else
      write(unit,iostat=ios) value%ig_num
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ig)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ig))deallocate(value%ig)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ig))then
        if(any(lo(1:3)/=lbound(value%ig)).or.any(hi(1:3)/=ubound(value%ig)))then
          ios=1;return
        endif
      else
        allocate(value%ig(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ig);hi(1:3)=ubound(value%ig)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ig
    else
      write(unit,iostat=ios) value%ig
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ig_bound)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ig_bound))deallocate(value%ig_bound)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ig_bound))then
        if(any(lo(1:3)/=lbound(value%ig_bound)).or.any(hi(1:3)/=ubound(value%ig_bound)))then
          ios=1;return
        endif
      else
        allocate(value%ig_bound(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ig_bound);hi(1:3)=ubound(value%ig_bound)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ig_bound
    else
      write(unit,iostat=ios) value%ig_bound
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%wkbound)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%wkbound))deallocate(value%wkbound)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%wkbound))then
        if(any(lo(1:1)/=lbound(value%wkbound)).or.any(hi(1:1)/=ubound(value%wkbound)))then
          ios=1;return
        endif
      else
        allocate(value%wkbound(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%wkbound);hi(1:1)=ubound(value%wkbound)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%wkbound
    else
      write(unit,iostat=ios) value%wkbound
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%wkbound2)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%wkbound2))deallocate(value%wkbound2)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%wkbound2))then
        if(any(lo(1:1)/=lbound(value%wkbound2)).or.any(hi(1:1)/=ubound(value%wkbound2)))then
          ios=1;return
        endif
      else
        allocate(value%wkbound2(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%wkbound2);hi(1:1)=ubound(value%wkbound2)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%wkbound2
    else
      write(unit,iostat=ios) value%wkbound2
    endif
    if(ios/=0)return
  endif
  if(reading)then
    read(unit,iostat=ios) value%n_multipole_xyz
  else
    write(unit,iostat=ios) value%n_multipole_xyz
  endif
  if(ios/=0)return
  present_array=allocated(value%zrhoG_ele)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%zrhoG_ele))deallocate(value%zrhoG_ele)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%zrhoG_ele))then
        if(any(lo(1:3)/=lbound(value%zrhoG_ele)).or.any(hi(1:3)/=ubound(value%zrhoG_ele)))then
          ios=1;return
        endif
      else
        allocate(value%zrhoG_ele(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%zrhoG_ele);hi(1:3)=ubound(value%zrhoG_ele)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%zrhoG_ele
    else
      write(unit,iostat=ios) value%zrhoG_ele
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff1x)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff1x))deallocate(value%ff1x)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff1x))then
        if(any(lo(1:3)/=lbound(value%ff1x)).or.any(hi(1:3)/=ubound(value%ff1x)))then
          ios=1;return
        endif
      else
        allocate(value%ff1x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff1x);hi(1:3)=ubound(value%ff1x)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff1x
    else
      write(unit,iostat=ios) value%ff1x
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff1y)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff1y))deallocate(value%ff1y)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff1y))then
        if(any(lo(1:3)/=lbound(value%ff1y)).or.any(hi(1:3)/=ubound(value%ff1y)))then
          ios=1;return
        endif
      else
        allocate(value%ff1y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff1y);hi(1:3)=ubound(value%ff1y)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff1y
    else
      write(unit,iostat=ios) value%ff1y
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff1z)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff1z))deallocate(value%ff1z)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff1z))then
        if(any(lo(1:3)/=lbound(value%ff1z)).or.any(hi(1:3)/=ubound(value%ff1z)))then
          ios=1;return
        endif
      else
        allocate(value%ff1z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff1z);hi(1:3)=ubound(value%ff1z)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff1z
    else
      write(unit,iostat=ios) value%ff1z
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff2x)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff2x))deallocate(value%ff2x)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff2x))then
        if(any(lo(1:3)/=lbound(value%ff2x)).or.any(hi(1:3)/=ubound(value%ff2x)))then
          ios=1;return
        endif
      else
        allocate(value%ff2x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff2x);hi(1:3)=ubound(value%ff2x)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff2x
    else
      write(unit,iostat=ios) value%ff2x
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff2y)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff2y))deallocate(value%ff2y)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff2y))then
        if(any(lo(1:3)/=lbound(value%ff2y)).or.any(hi(1:3)/=ubound(value%ff2y)))then
          ios=1;return
        endif
      else
        allocate(value%ff2y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff2y);hi(1:3)=ubound(value%ff2y)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff2y
    else
      write(unit,iostat=ios) value%ff2y
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff2z)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff2z))deallocate(value%ff2z)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff2z))then
        if(any(lo(1:3)/=lbound(value%ff2z)).or.any(hi(1:3)/=ubound(value%ff2z)))then
          ios=1;return
        endif
      else
        allocate(value%ff2z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff2z);hi(1:3)=ubound(value%ff2z)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff2z
    else
      write(unit,iostat=ios) value%ff2z
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff1)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff1))deallocate(value%ff1)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff1))then
        if(any(lo(1:3)/=lbound(value%ff1)).or.any(hi(1:3)/=ubound(value%ff1)))then
          ios=1;return
        endif
      else
        allocate(value%ff1(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff1);hi(1:3)=ubound(value%ff1)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff1
    else
      write(unit,iostat=ios) value%ff1
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff2)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff2))deallocate(value%ff2)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff2))then
        if(any(lo(1:3)/=lbound(value%ff2)).or.any(hi(1:3)/=ubound(value%ff2)))then
          ios=1;return
        endif
      else
        allocate(value%ff2(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff2);hi(1:3)=ubound(value%ff2)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff2
    else
      write(unit,iostat=ios) value%ff2
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff3x)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff3x))deallocate(value%ff3x)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff3x))then
        if(any(lo(1:3)/=lbound(value%ff3x)).or.any(hi(1:3)/=ubound(value%ff3x)))then
          ios=1;return
        endif
      else
        allocate(value%ff3x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff3x);hi(1:3)=ubound(value%ff3x)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff3x
    else
      write(unit,iostat=ios) value%ff3x
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff3y)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff3y))deallocate(value%ff3y)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff3y))then
        if(any(lo(1:3)/=lbound(value%ff3y)).or.any(hi(1:3)/=ubound(value%ff3y)))then
          ios=1;return
        endif
      else
        allocate(value%ff3y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff3y);hi(1:3)=ubound(value%ff3y)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff3y
    else
      write(unit,iostat=ios) value%ff3y
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff3z)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff3z))deallocate(value%ff3z)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff3z))then
        if(any(lo(1:3)/=lbound(value%ff3z)).or.any(hi(1:3)/=ubound(value%ff3z)))then
          ios=1;return
        endif
      else
        allocate(value%ff3z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff3z);hi(1:3)=ubound(value%ff3z)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff3z
    else
      write(unit,iostat=ios) value%ff3z
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff4x)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff4x))deallocate(value%ff4x)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff4x))then
        if(any(lo(1:3)/=lbound(value%ff4x)).or.any(hi(1:3)/=ubound(value%ff4x)))then
          ios=1;return
        endif
      else
        allocate(value%ff4x(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff4x);hi(1:3)=ubound(value%ff4x)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff4x
    else
      write(unit,iostat=ios) value%ff4x
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff4y)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff4y))deallocate(value%ff4y)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff4y))then
        if(any(lo(1:3)/=lbound(value%ff4y)).or.any(hi(1:3)/=ubound(value%ff4y)))then
          ios=1;return
        endif
      else
        allocate(value%ff4y(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff4y);hi(1:3)=ubound(value%ff4y)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff4y
    else
      write(unit,iostat=ios) value%ff4y
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%ff4z)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%ff4z))deallocate(value%ff4z)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%ff4z))then
        if(any(lo(1:3)/=lbound(value%ff4z)).or.any(hi(1:3)/=ubound(value%ff4z)))then
          ios=1;return
        endif
      else
        allocate(value%ff4z(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%ff4z);hi(1:3)=ubound(value%ff4z)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%ff4z
    else
      write(unit,iostat=ios) value%ff4z
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%dgf)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%dgf))deallocate(value%dgf)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%dgf))then
        if(any(lo(1:3)/=lbound(value%dgf)).or.any(hi(1:3)/=ubound(value%dgf)))then
          ios=1;return
        endif
      else
        allocate(value%dgf(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%dgf);hi(1:3)=ubound(value%dgf)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%dgf
    else
      write(unit,iostat=ios) value%dgf
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%a_ffte)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%a_ffte))deallocate(value%a_ffte)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%a_ffte))then
        if(any(lo(1:3)/=lbound(value%a_ffte)).or.any(hi(1:3)/=ubound(value%a_ffte)))then
          ios=1;return
        endif
      else
        allocate(value%a_ffte(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%a_ffte);hi(1:3)=ubound(value%a_ffte)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%a_ffte
    else
      write(unit,iostat=ios) value%a_ffte
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%b_ffte)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%b_ffte))deallocate(value%b_ffte)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%b_ffte))then
        if(any(lo(1:3)/=lbound(value%b_ffte)).or.any(hi(1:3)/=ubound(value%b_ffte)))then
          ios=1;return
        endif
      else
        allocate(value%b_ffte(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%b_ffte);hi(1:3)=ubound(value%b_ffte)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%b_ffte
    else
      write(unit,iostat=ios) value%b_ffte
    endif
    if(ios/=0)return
  endif
end subroutine
subroutine checkpoint_io_s_rt(unit,reading,value,ios)
  integer,intent(in)::unit
  logical,intent(in)::reading
  type(s_rt),intent(inout)::value
  integer,intent(inout)::ios
  logical::present_array
  integer::lo(7),hi(7),i1,i2,k,allocation_status
  integer(8)::elements,extent
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%Dp0_e
  else
    write(unit,iostat=ios) value%Dp0_e
  endif
  if(ios/=0)return
  present_array=allocated(value%dDp_e)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%dDp_e))deallocate(value%dDp_e)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%dDp_e))then
        if(any(lo(1:2)/=lbound(value%dDp_e)).or.any(hi(1:2)/=ubound(value%dDp_e)))then
          ios=1;return
        endif
      else
        allocate(value%dDp_e(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%dDp_e);hi(1:2)=ubound(value%dDp_e)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%dDp_e
    else
      write(unit,iostat=ios) value%dDp_e
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%Dp_e)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%Dp_e))deallocate(value%Dp_e)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%Dp_e))then
        if(any(lo(1:2)/=lbound(value%Dp_e)).or.any(hi(1:2)/=ubound(value%Dp_e)))then
          ios=1;return
        endif
      else
        allocate(value%Dp_e(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%Dp_e);hi(1:2)=ubound(value%Dp_e)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%Dp_e
    else
      write(unit,iostat=ios) value%Dp_e
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%Dp_i)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%Dp_i))deallocate(value%Dp_i)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%Dp_i))then
        if(any(lo(1:2)/=lbound(value%Dp_i)).or.any(hi(1:2)/=ubound(value%Dp_i)))then
          ios=1;return
        endif
      else
        allocate(value%Dp_i(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%Dp_i);hi(1:2)=ubound(value%Dp_i)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%Dp_i
    else
      write(unit,iostat=ios) value%Dp_i
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%Qp_e)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%Qp_e))deallocate(value%Qp_e)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%Qp_e))then
        if(any(lo(1:3)/=lbound(value%Qp_e)).or.any(hi(1:3)/=ubound(value%Qp_e)))then
          ios=1;return
        endif
      else
        allocate(value%Qp_e(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%Qp_e);hi(1:3)=ubound(value%Qp_e)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%Qp_e
    else
      write(unit,iostat=ios) value%Qp_e
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%rIe)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%rIe))deallocate(value%rIe)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%rIe))then
        if(any(lo(1:1)/=lbound(value%rIe)).or.any(hi(1:1)/=ubound(value%rIe)))then
          ios=1;return
        endif
      else
        allocate(value%rIe(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%rIe);hi(1:1)=ubound(value%rIe)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%rIe
    else
      write(unit,iostat=ios) value%rIe
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%curr)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%curr))deallocate(value%curr)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%curr))then
        if(any(lo(1:2)/=lbound(value%curr)).or.any(hi(1:2)/=ubound(value%curr)))then
          ios=1;return
        endif
      else
        allocate(value%curr(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%curr);hi(1:2)=ubound(value%curr)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%curr
    else
      write(unit,iostat=ios) value%curr
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%E_ext)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%E_ext))deallocate(value%E_ext)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%E_ext))then
        if(any(lo(1:2)/=lbound(value%E_ext)).or.any(hi(1:2)/=ubound(value%E_ext)))then
          ios=1;return
        endif
      else
        allocate(value%E_ext(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%E_ext);hi(1:2)=ubound(value%E_ext)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%E_ext
    else
      write(unit,iostat=ios) value%E_ext
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%E_ind)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%E_ind))deallocate(value%E_ind)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%E_ind))then
        if(any(lo(1:2)/=lbound(value%E_ind)).or.any(hi(1:2)/=ubound(value%E_ind)))then
          ios=1;return
        endif
      else
        allocate(value%E_ind(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%E_ind);hi(1:2)=ubound(value%E_ind)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%E_ind
    else
      write(unit,iostat=ios) value%E_ind
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%E_tot)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%E_tot))deallocate(value%E_tot)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%E_tot))then
        if(any(lo(1:2)/=lbound(value%E_tot)).or.any(hi(1:2)/=ubound(value%E_tot)))then
          ios=1;return
        endif
      else
        allocate(value%E_tot(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%E_tot);hi(1:2)=ubound(value%E_tot)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%E_tot
    else
      write(unit,iostat=ios) value%E_tot
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%Ac_ext)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%Ac_ext))deallocate(value%Ac_ext)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%Ac_ext))then
        if(any(lo(1:2)/=lbound(value%Ac_ext)).or.any(hi(1:2)/=ubound(value%Ac_ext)))then
          ios=1;return
        endif
      else
        allocate(value%Ac_ext(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%Ac_ext);hi(1:2)=ubound(value%Ac_ext)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%Ac_ext
    else
      write(unit,iostat=ios) value%Ac_ext
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%Ac_ind)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%Ac_ind))deallocate(value%Ac_ind)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%Ac_ind))then
        if(any(lo(1:2)/=lbound(value%Ac_ind)).or.any(hi(1:2)/=ubound(value%Ac_ind)))then
          ios=1;return
        endif
      else
        allocate(value%Ac_ind(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%Ac_ind);hi(1:2)=ubound(value%Ac_ind)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%Ac_ind
    else
      write(unit,iostat=ios) value%Ac_ind
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%Ac_tot)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%Ac_tot))deallocate(value%Ac_tot)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%Ac_tot))then
        if(any(lo(1:2)/=lbound(value%Ac_tot)).or.any(hi(1:2)/=ubound(value%Ac_tot)))then
          ios=1;return
        endif
      else
        allocate(value%Ac_tot(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%Ac_tot);hi(1:2)=ubound(value%Ac_tot)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%Ac_tot
    else
      write(unit,iostat=ios) value%Ac_tot
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%zc)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%zc))deallocate(value%zc)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%zc))then
        if(any(lo(1:1)/=lbound(value%zc)).or.any(hi(1:1)/=ubound(value%zc)))then
          ios=1;return
        endif
      else
        allocate(value%zc(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%zc);hi(1:1)=ubound(value%zc)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%zc
    else
      write(unit,iostat=ios) value%zc
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%vloc_t)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%vloc_t))deallocate(value%vloc_t)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%vloc_t))then
        if(any(lo(1:1)/=lbound(value%vloc_t)).or.any(hi(1:1)/=ubound(value%vloc_t)))then
          ios=1;return
        endif
      else
        allocate(value%vloc_t(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%vloc_t);hi(1:1)=ubound(value%vloc_t)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    do i1=lbound(value%vloc_t,1),ubound(value%vloc_t,1)
    call checkpoint_io_s_scalar(unit,reading,value%vloc_t(i1),ios)
    if(ios/=0)return
    enddo
  endif
  present_array=allocated(value%vloc_new)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%vloc_new))deallocate(value%vloc_new)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%vloc_new))then
        if(any(lo(1:1)/=lbound(value%vloc_new)).or.any(hi(1:1)/=ubound(value%vloc_new)))then
          ios=1;return
        endif
      else
        allocate(value%vloc_new(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%vloc_new);hi(1:1)=ubound(value%vloc_new)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    do i1=lbound(value%vloc_new,1),ubound(value%vloc_new,1)
    call checkpoint_io_s_scalar(unit,reading,value%vloc_new(i1),ios)
    if(ios/=0)return
    enddo
  endif
  present_array=allocated(value%vloc_old)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%vloc_old))deallocate(value%vloc_old)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%vloc_old))then
        if(any(lo(1:2)/=lbound(value%vloc_old)).or.any(hi(1:2)/=ubound(value%vloc_old)))then
          ios=1;return
        endif
      else
        allocate(value%vloc_old(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%vloc_old);hi(1:2)=ubound(value%vloc_old)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    do i1=lbound(value%vloc_old,1),ubound(value%vloc_old,1)
    do i2=lbound(value%vloc_old,2),ubound(value%vloc_old,2)
    call checkpoint_io_s_scalar(unit,reading,value%vloc_old(i1,i2),ios)
    if(ios/=0)return
    enddo
    enddo
  endif
  present_array=allocated(value%rho0_s)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%rho0_s))deallocate(value%rho0_s)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
      elements=1_8
      do k=1,1
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%rho0_s))then
        if(any(lo(1:1)/=lbound(value%rho0_s)).or.any(hi(1:1)/=ubound(value%rho0_s)))then
          ios=1;return
        endif
      else
        allocate(value%rho0_s(lo(1):hi(1)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:1)=lbound(value%rho0_s);hi(1:1)=ubound(value%rho0_s)
      write(unit,iostat=ios)lo(1:1),hi(1:1)
      if(ios/=0)return
    endif
    do i1=lbound(value%rho0_s,1),ubound(value%rho0_s,1)
    call checkpoint_io_s_scalar(unit,reading,value%rho0_s(i1),ios)
    if(ios/=0)return
    enddo
  endif
  if(reading)then
    read(unit,iostat=ios)value%solid_lmax,value%solid_center,value%solid_radius
  else
    write(unit,iostat=ios)value%solid_lmax,value%solid_center,value%solid_radius
  endif
  if(ios/=0)return
  present_array=allocated(value%solid_coeff)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%solid_coeff))deallocate(value%solid_coeff)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
      elements=1_8
      do k=1,2
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%solid_coeff))then
        if(any(lo(1:2)/=lbound(value%solid_coeff)).or.any(hi(1:2)/=ubound(value%solid_coeff)))then
          ios=1;return
        endif
      else
        allocate(value%solid_coeff(lo(1):hi(1),lo(2):hi(2)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:2)=lbound(value%solid_coeff);hi(1:2)=ubound(value%solid_coeff)
      write(unit,iostat=ios)lo(1:2),hi(1:2)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%solid_coeff
    else
      write(unit,iostat=ios) value%solid_coeff
    endif
    if(ios/=0)return
  endif
  if(reading)then
    read(unit,iostat=ios)value%near_lmax,value%near_index,value%near_geometry,value%near_magnetic
  else
    write(unit,iostat=ios)value%near_lmax,value%near_index,value%near_geometry,value%near_magnetic
  endif
  if(ios/=0)return
  present_array=allocated(value%near_Ecoeff)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%near_Ecoeff))deallocate(value%near_Ecoeff)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%near_Ecoeff))then
        if(any(lo(1:3)/=lbound(value%near_Ecoeff)).or.any(hi(1:3)/=ubound(value%near_Ecoeff)))then
          ios=1;return
        endif
      else
        allocate(value%near_Ecoeff(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%near_Ecoeff);hi(1:3)=ubound(value%near_Ecoeff)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%near_Ecoeff
    else
      write(unit,iostat=ios) value%near_Ecoeff
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%near_Hcoeff)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%near_Hcoeff))deallocate(value%near_Hcoeff)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%near_Hcoeff))then
        if(any(lo(1:3)/=lbound(value%near_Hcoeff)).or.any(hi(1:3)/=ubound(value%near_Hcoeff)))then
          ios=1;return
        endif
      else
        allocate(value%near_Hcoeff(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%near_Hcoeff);hi(1:3)=ubound(value%near_Hcoeff)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%near_Hcoeff
    else
      write(unit,iostat=ios) value%near_Hcoeff
    endif
    if(ios/=0)return
  endif
  if(reading)then
    read(unit,iostat=ios)value%near_grid_counts,value%near_grid_origins,value%near_grid_spacing,value%near_grid_center
  else
    write(unit,iostat=ios)value%near_grid_counts,value%near_grid_origins,value%near_grid_spacing,value%near_grid_center
  endif
  if(ios/=0)return
  present_array=allocated(value%near_grid_values)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%near_grid_values))deallocate(value%near_grid_values)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:6),hi(1:6)
      if(ios/=0)return
      elements=1_8
      do k=1,6
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(any(lo(1:5)/=1).or.lo(6)/=0.or.hi(1)/=3.or.hi(5)/=2.or. &
         any(hi(2:4)/=value%near_grid_counts).or.hi(6)/=value%max_steps+1)then
        ios=1;return
      endif
      if(allocated(value%near_grid_values))then
        if(any(lo(1:6)/=lbound(value%near_grid_values)).or.any(hi(1:6)/=ubound(value%near_grid_values)))then
          ios=1;return
        endif
      else
        allocate(value%near_grid_values(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3), &
          lo(4):hi(4),lo(5):hi(5),lo(6):hi(6)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:6)=lbound(value%near_grid_values);hi(1:6)=ubound(value%near_grid_values)
      write(unit,iostat=ios)lo(1:6),hi(1:6)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%near_grid_values
    else
      write(unit,iostat=ios) value%near_grid_values
    endif
    if(ios/=0)return
  endif
  present_array=allocated(value%near_projectors)
  if(reading)then
    read(unit,iostat=ios) present_array
  else
    write(unit,iostat=ios) present_array
  endif
  if(ios/=0)return
  if(reading.and..not.present_array)then
    if(allocated(value%near_projectors))deallocate(value%near_projectors)
  endif
  if(present_array)then
    if(reading)then
      read(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
      elements=1_8
      do k=1,3
        extent=int(hi(k),8)-int(lo(k),8)+1_8
        if(extent<0_8.or.extent>checkpoint_io_limit())then
          ios=1;return
        endif
        if(extent>0_8)then
          if(elements>checkpoint_io_limit()/extent)then
            ios=1;return
          endif
        endif
        elements=elements*extent
      enddo
      if(allocated(value%near_projectors))then
        if(any(lo(1:3)/=lbound(value%near_projectors)).or.any(hi(1:3)/=ubound(value%near_projectors)))then
          ios=1;return
        endif
      else
        allocate(value%near_projectors(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),stat=allocation_status)
        if(allocation_status/=0)then
          ios=1;return
        endif
      endif
    else
      lo(1:3)=lbound(value%near_projectors);hi(1:3)=ubound(value%near_projectors)
      write(unit,iostat=ios)lo(1:3),hi(1:3)
      if(ios/=0)return
    endif
    if(reading)then
      read(unit,iostat=ios) value%near_projectors
    else
      write(unit,iostat=ios) value%near_projectors
    endif
    if(ios/=0)return
  endif
  call checkpoint_io_s_vector(unit,reading,value%near_Ac,ios)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%near_div,ios)
  if(ios/=0)return
  call checkpoint_io_s_scalar(unit,reading,value%vonf,ios)
  if(ios/=0)return
  call checkpoint_io_s_vector(unit,reading,value%j_e,ios)
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%E_old
  else
    write(unit,iostat=ios) value%E_old
  endif
  if(ios/=0)return
  if(reading)then
    read(unit,iostat=ios) value%rho0_loaded
  else
    write(unit,iostat=ios) value%rho0_loaded
  endif
  if(ios/=0)return
end subroutine
