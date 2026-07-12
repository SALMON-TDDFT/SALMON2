module lcfo_wannier_sawf_seed
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: write_sawf_local_eig_amn_mmn

contains

  subroutine write_sawf_local_eig_amn_mmn(directory,seedname,energy_ev,amn,mmn,ok,message)
    character(*),intent(in)::directory,seedname
    real(8),intent(in)::energy_ev(:)
    complex(8),intent(in)::amn(:,:),mmn(:,:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,parameter::gvec(3,3)=reshape([1,0,0,0,1,0,0,0,1],[3,3])
    integer::unit,ios,iband,iproj,ineighbor,jband
    character(1024)::filename

    ok=.false.;message=''
    if(len_trim(directory)==0.or.len_trim(seedname)==0.or.size(energy_ev)<=0.or. &
        size(amn,1)/=size(energy_ev).or.size(amn,2)<=0.or. &
        size(mmn,1)/=size(energy_ev).or.size(mmn,2)/=size(energy_ev).or.size(mmn,3)/=3)then
      message='SAWF local seed dimensions are inconsistent';return
    end if
    if(.not.all(ieee_is_finite(energy_ev)).or. &
        .not.all(ieee_is_finite(real(amn))).or..not.all(ieee_is_finite(aimag(amn))).or. &
        .not.all(ieee_is_finite(real(mmn))).or..not.all(ieee_is_finite(aimag(mmn))))then
      message='SAWF local seed contains a non-finite value';return
    end if

    filename=trim(directory)//'/'//trim(seedname)//'.eig'
    open(newunit=unit,file=trim(filename),status='replace',action='write',iostat=ios)
    if(ios/=0)then;message='SAWF local .eig open failed';return;end if
    do iband=1,size(energy_ev);write(unit,'(2i8,1x,es23.15)',iostat=ios)iband,1,energy_ev(iband);if(ios/=0)exit;end do
    close(unit);if(ios/=0)then;message='SAWF local .eig write failed';return;end if

    filename=trim(directory)//'/'//trim(seedname)//'.amn'
    open(newunit=unit,file=trim(filename),status='replace',action='write',iostat=ios)
    if(ios/=0)then;message='SAWF local .amn open failed';return;end if
    write(unit,'(a)',iostat=ios)'SALMON local SAWF projections'
    if(ios==0)write(unit,'(3i10)',iostat=ios)size(amn,1),1,size(amn,2)
    do iproj=1,size(amn,2);do iband=1,size(amn,1)
      if(ios==0)write(unit,'(3i8,2(1x,es23.15))',iostat=ios)iband,iproj,1,amn(iband,iproj)
    end do;end do
    close(unit);if(ios/=0)then;message='SAWF local .amn write failed';return;end if

    filename=trim(directory)//'/'//trim(seedname)//'.mmn'
    open(newunit=unit,file=trim(filename),status='replace',action='write',iostat=ios)
    if(ios/=0)then;message='SAWF local .mmn open failed';return;end if
    write(unit,'(a)',iostat=ios)'SALMON local SAWF overlaps'
    if(ios==0)write(unit,'(3i10)',iostat=ios)size(mmn,1),1,3
    do ineighbor=1,3
      if(ios==0)write(unit,'(5i8)',iostat=ios)1,1,gvec(:,ineighbor)
      do jband=1,size(mmn,2);do iband=1,size(mmn,1)
        if(ios==0)write(unit,'(2(1x,es23.15))',iostat=ios)mmn(iband,jband,ineighbor)
      end do;end do
    end do
    close(unit);if(ios/=0)then;message='SAWF local .mmn write failed';return;end if
    ok=.true.
  end subroutine write_sawf_local_eig_amn_mmn

end module lcfo_wannier_sawf_seed
