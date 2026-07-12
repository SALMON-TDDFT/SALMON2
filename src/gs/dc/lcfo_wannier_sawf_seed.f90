module lcfo_wannier_sawf_seed
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: write_sawf_local_eig_amn_mmn
  public :: build_sawf_local_seed_matrices
  public :: select_sawf_local_complete_shells
  public :: solve_sawf_local_generalized_eigensystem
  public :: read_sawf_nnkp_neighbors
  public :: write_sawf_local_eig_amn
  public :: restrict_sawf_stabilizer_representation

contains

  subroutine restrict_sawf_stabilizer_representation(global_representation,selected_channel, &
      stabilizer_operation,tolerance,local_representation,ok,message)
    complex(8),intent(in)::global_representation(:,:,:)
    integer,intent(in)::selected_channel(:),stabilizer_operation(:)
    real(8),intent(in)::tolerance
    complex(8),intent(out)::local_representation(:,:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,allocatable::selected(:)
    complex(8),allocatable::identity(:,:)
    real(8)::leakage,unitarity
    integer::operation,local_operation,row,column,i,j,nchannel

    ok=.false.;message='';local_representation=(0d0,0d0);nchannel=size(global_representation,1)
    if(nchannel<=0.or.size(global_representation,2)/=nchannel.or.size(selected_channel)<=0.or. &
        size(stabilizer_operation)<=0.or.size(local_representation,1)/=size(selected_channel).or. &
        size(local_representation,2)/=size(selected_channel).or. &
        size(local_representation,3)/=size(stabilizer_operation).or. &
        any(selected_channel<1).or.any(selected_channel>nchannel).or. &
        any(stabilizer_operation<1).or.any(stabilizer_operation>size(global_representation,3)).or. &
        .not.ieee_is_finite(tolerance).or.tolerance<=0d0)then
      message='SAWF local stabilizer representation dimensions are invalid';return
    end if
    allocate(selected(nchannel),identity(size(selected_channel),size(selected_channel)))
    selected=.false.;do i=1,size(selected_channel)
      if(selected(selected_channel(i)))then;message='SAWF selected local channel is duplicated';return;end if
      selected(selected_channel(i))=.true.
    end do
    identity=(0d0,0d0);do i=1,size(selected_channel);identity(i,i)=(1d0,0d0);end do
    do local_operation=1,size(stabilizer_operation)
      operation=stabilizer_operation(local_operation);leakage=0d0
      do column=1,nchannel;do row=1,nchannel
        if(selected(row).neqv.selected(column))leakage=max(leakage,abs(global_representation(row,column,operation)))
      end do;end do
      if(leakage>tolerance)then
        message='SAWF global representation leaks outside the local complete-shell subspace';return
      end if
      do j=1,size(selected_channel);do i=1,size(selected_channel)
        local_representation(i,j,local_operation)= &
          global_representation(selected_channel(i),selected_channel(j),operation)
      end do;end do
      unitarity=maxval(abs(matmul(conjg(transpose(local_representation(:,:,local_operation))), &
        local_representation(:,:,local_operation))-identity))
      if(unitarity>tolerance)then;message='SAWF restricted local representation is not unitary';return;end if
    end do
    ok=.true.
  end subroutine restrict_sawf_stabilizer_representation

  subroutine read_sawf_nnkp_neighbors(filename,neighbor_gvec,ok,message)
    character(*),intent(in)::filename
    integer,allocatable,intent(out)::neighbor_gvec(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    character(512)::line
    integer::unit,ios,nneighbor,neighbor,ikpoint,jkpoint
    logical::found

    ok=.false.;message='';found=.false.
    open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
    if(ios/=0)then;message='SAWF local .nnkp is missing';return;end if
    do
      read(unit,'(a)',iostat=ios)line
      if(ios/=0)exit
      if(trim(adjustl(line))=='begin nnkpts')then;found=.true.;exit;end if
    end do
    if(.not.found)then;close(unit);message='SAWF local .nnkp has no nnkpts block';return;end if
    read(unit,*,iostat=ios)nneighbor
    if(ios/=0.or.nneighbor<=0)then;close(unit);message='SAWF local .nnkp neighbor count is invalid';return;end if
    allocate(neighbor_gvec(3,nneighbor))
    do neighbor=1,nneighbor
      read(unit,*,iostat=ios)ikpoint,jkpoint,neighbor_gvec(:,neighbor)
      if(ios/=0.or.ikpoint/=1.or.jkpoint/=1)then
        close(unit);message='SAWF local .nnkp is not a Gamma-only neighbor list';return
      end if
    end do
    read(unit,'(a)',iostat=ios)line;close(unit)
    if(ios/=0.or.trim(adjustl(line))/='end nnkpts')then
      message='SAWF local .nnkp nnkpts block is truncated';return
    end if
    ok=.true.
  end subroutine read_sawf_nnkp_neighbors

  subroutine solve_sawf_local_generalized_eigensystem(buffer_basis,weight,h_basis,rank_tolerance, &
      states,energies,ok,message)
    real(8),intent(in)::buffer_basis(:,:),weight,h_basis(:,:),rank_tolerance
    complex(8),allocatable,intent(out)::states(:,:)
    real(8),allocatable,intent(out)::energies(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8),allocatable::overlap(:,:),overlap_eval(:),whitener(:,:),h_orth(:,:),h_eval(:),coeff(:,:), &
      diagonal(:,:)
    complex(8),allocatable::identity(:,:)
    real(8)::cutoff,orthogonality_residual,diagonalization_residual
    integer::nbasis,nkeep,index,mode

    ok=.false.;message='';nbasis=size(buffer_basis,2)
    if(size(buffer_basis,1)<=0.or.nbasis<=0.or.size(h_basis,1)/=nbasis.or. &
        size(h_basis,2)/=nbasis.or..not.ieee_is_finite(weight).or.weight<=0d0.or. &
        .not.ieee_is_finite(rank_tolerance).or.rank_tolerance<=0d0.or. &
        .not.all(ieee_is_finite(buffer_basis)).or..not.all(ieee_is_finite(h_basis)))then
      message='SAWF local generalized eigensystem inputs are invalid';return
    end if
    if(maxval(abs(h_basis-transpose(h_basis)))>rank_tolerance*max(1d0,maxval(abs(h_basis))))then
      message='SAWF local Hamiltonian is not symmetric';return
    end if
    allocate(overlap(nbasis,nbasis),overlap_eval(nbasis))
    overlap=weight*matmul(transpose(buffer_basis),buffer_basis)
    call diagonalize_sawf_real_symmetric(overlap,overlap_eval,ok,message)
    if(.not.ok)return
    cutoff=rank_tolerance*max(1d0,maxval(overlap_eval));nkeep=count(overlap_eval>cutoff)
    if(nkeep<=0)then;message='SAWF buffered overlap is rank deficient at all modes';ok=.false.;return;end if
    allocate(whitener(nbasis,nkeep));index=0
    do mode=1,nbasis
      if(overlap_eval(mode)<=cutoff)cycle
      index=index+1;whitener(:,index)=overlap(:,mode)/sqrt(overlap_eval(mode))
    end do
    allocate(h_orth(nkeep,nkeep),h_eval(nkeep))
    h_orth=matmul(transpose(whitener),matmul(h_basis,whitener))
    h_orth=0.5d0*(h_orth+transpose(h_orth))
    call diagonalize_sawf_real_symmetric(h_orth,h_eval,ok,message)
    if(.not.ok)return
    allocate(coeff(nbasis,nkeep),states(size(buffer_basis,1),nkeep),energies(nkeep))
    coeff=matmul(whitener,h_orth)
    states=cmplx(matmul(buffer_basis,coeff),0d0,kind=8);energies=h_eval
    allocate(identity(nkeep,nkeep),diagonal(nkeep,nkeep));identity=(0d0,0d0);diagonal=0d0
    do mode=1,nkeep;identity(mode,mode)=(1d0,0d0);diagonal(mode,mode)=energies(mode);end do
    orthogonality_residual=maxval(abs(weight*matmul(conjg(transpose(states)),states)-identity))
    diagonalization_residual=maxval(abs(matmul(transpose(coeff),matmul(h_basis,coeff))-diagonal))
    if(orthogonality_residual>100d0*rank_tolerance.or. &
        diagonalization_residual>100d0*rank_tolerance*max(1d0,maxval(abs(energies))))then
      message='SAWF local generalized eigensystem residual exceeds tolerance';ok=.false.;return
    end if
    ok=.true.
  end subroutine solve_sawf_local_generalized_eigensystem

  subroutine diagonalize_sawf_real_symmetric(matrix,eigenvalue,ok,message)
    real(8),intent(inout)::matrix(:,:)
    real(8),intent(out)::eigenvalue(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8),allocatable::work(:)
    integer::info,lwork,n
    external::dsyev
    ok=.false.;n=size(matrix,1)
    if(size(matrix,2)/=n.or.size(eigenvalue)/=n.or.n<=0)then
      message='SAWF symmetric eigensolver dimensions are invalid';return
    end if
    allocate(work(1));lwork=-1
    call dsyev('V','U',n,matrix,n,eigenvalue,work,lwork,info)
    if(info/=0.or..not.ieee_is_finite(work(1)))then
      message='SAWF symmetric eigensolver workspace query failed';return
    end if
    lwork=max(1,int(work(1)));deallocate(work);allocate(work(lwork))
    call dsyev('V','U',n,matrix,n,eigenvalue,work,lwork,info)
    if(info/=0.or..not.all(ieee_is_finite(eigenvalue)).or..not.all(ieee_is_finite(matrix)))then
      message='SAWF symmetric eigensolver failed';return
    end if
    ok=.true.;message=''
  end subroutine diagonalize_sawf_real_symmetric

  subroutine select_sawf_local_complete_shells(channel_atom,expected_per_atom,inside_atom, &
      selected_channel,ok,message)
    integer,intent(in)::channel_atom(:),expected_per_atom(:)
    logical,intent(in)::inside_atom(:)
    integer,allocatable,intent(out)::selected_channel(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::channel_count(:)
    integer::channel,nselected

    ok=.false.;message=''
    if(size(channel_atom)<=0.or.size(expected_per_atom)<=0.or. &
        size(inside_atom)/=size(expected_per_atom).or.any(expected_per_atom<=0).or. &
        any(channel_atom<1).or.any(channel_atom>size(expected_per_atom)))then
      message='SAWF local projection-shell dimensions are invalid';return
    end if
    do channel=2,size(channel_atom)
      if(channel_atom(channel)<channel_atom(channel-1))then
        message='SAWF projection channels are not atom-major';return
      end if
    end do
    allocate(channel_count(size(expected_per_atom)));channel_count=0
    do channel=1,size(channel_atom)
      channel_count(channel_atom(channel))=channel_count(channel_atom(channel))+1
    end do
    if(any(channel_count/=expected_per_atom))then
      message='SAWF projection channels do not contain complete atomic shells';return
    end if
    nselected=count(inside_atom(channel_atom));allocate(selected_channel(nselected));nselected=0
    do channel=1,size(channel_atom)
      if(.not.inside_atom(channel_atom(channel)))cycle
      nselected=nselected+1;selected_channel(nselected)=channel
    end do
    if(nselected<=0)then
      message='SAWF local core and buffer contain no complete projection shell';return
    end if
    ok=.true.
  end subroutine select_sawf_local_complete_shells

  subroutine build_sawf_local_seed_matrices(states,projections,phase_factor,weight,amn,mmn,ok,message)
    complex(8),intent(in)::states(:,:),projections(:,:),phase_factor(:,:)
    real(8),intent(in)::weight
    complex(8),intent(out)::amn(:,:),mmn(:,:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::normalized_projection(:,:),phased_states(:,:)
    real(8)::norm2
    integer::projection,neighbor

    ok=.false.;message='';amn=(0d0,0d0);mmn=(0d0,0d0)
    if(size(states,1)<=0.or.size(states,2)<=0.or.size(projections,1)/=size(states,1).or. &
        size(phase_factor,1)/=size(states,1).or.size(phase_factor,2)<=0.or. &
        size(amn,1)/=size(states,2).or.size(amn,2)/=size(projections,2).or. &
        size(mmn,1)/=size(states,2).or.size(mmn,2)/=size(states,2).or. &
        size(mmn,3)/=size(phase_factor,2).or. &
        .not.ieee_is_finite(weight).or.weight<=0d0)then
      message='SAWF local seed matrix dimensions or integration weight are invalid';return
    end if
    if(.not.all(ieee_is_finite(real(states))).or..not.all(ieee_is_finite(aimag(states))).or. &
        .not.all(ieee_is_finite(real(projections))).or..not.all(ieee_is_finite(aimag(projections))).or. &
        .not.all(ieee_is_finite(real(phase_factor))).or..not.all(ieee_is_finite(aimag(phase_factor))))then
      message='SAWF local seed matrix input contains a non-finite value';return
    end if
    allocate(normalized_projection(size(projections,1),size(projections,2)), &
      phased_states(size(states,1),size(states,2)))
    do projection=1,size(projections,2)
      norm2=weight*sum(abs(projections(:,projection))**2)
      if(norm2<=1d-300)then
        message='SAWF local projection has zero norm';return
      end if
      normalized_projection(:,projection)=projections(:,projection)/sqrt(norm2)
    end do
    amn=weight*matmul(conjg(transpose(states)),normalized_projection)
    do neighbor=1,size(phase_factor,2)
      phased_states=spread(phase_factor(:,neighbor),2,size(states,2))*states
      mmn(:,:,neighbor)=weight*matmul(conjg(transpose(states)),phased_states)
    end do
    ok=.true.
  end subroutine build_sawf_local_seed_matrices

  subroutine write_sawf_local_eig_amn_mmn(directory,seedname,energy_ev,amn,mmn,neighbor_gvec,ok,message)
    character(*),intent(in)::directory,seedname
    real(8),intent(in)::energy_ev(:)
    complex(8),intent(in)::amn(:,:),mmn(:,:,:)
    integer,intent(in)::neighbor_gvec(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::unit,ios,iband,iproj,ineighbor,jband
    character(1024)::filename

    ok=.false.;message=''
    if(len_trim(directory)==0.or.len_trim(seedname)==0.or.size(energy_ev)<=0.or. &
        size(amn,1)/=size(energy_ev).or.size(amn,2)<=0.or. &
        size(mmn,1)/=size(energy_ev).or.size(mmn,2)/=size(energy_ev).or. &
        size(neighbor_gvec,1)/=3.or.size(neighbor_gvec,2)/=size(mmn,3).or.size(mmn,3)<=0)then
      message='SAWF local seed dimensions are inconsistent';return
    end if
    if(.not.all(ieee_is_finite(energy_ev)).or. &
        .not.all(ieee_is_finite(real(amn))).or..not.all(ieee_is_finite(aimag(amn))).or. &
        .not.all(ieee_is_finite(real(mmn))).or..not.all(ieee_is_finite(aimag(mmn))))then
      message='SAWF local seed contains a non-finite value';return
    end if

    call write_sawf_local_eig_amn(directory,seedname,energy_ev,amn,ok,message)
    if(.not.ok)return

    filename=trim(directory)//'/'//trim(seedname)//'.mmn'
    open(newunit=unit,file=trim(filename),status='replace',action='write',iostat=ios)
    if(ios/=0)then;message='SAWF local .mmn open failed';ok=.false.;return;end if
    write(unit,'(a)',iostat=ios)'SALMON local SAWF overlaps'
    if(ios==0)write(unit,'(3i10)',iostat=ios)size(mmn,1),1,size(mmn,3)
    do ineighbor=1,size(mmn,3)
      if(ios==0)write(unit,'(5i8)',iostat=ios)1,1,neighbor_gvec(:,ineighbor)
      do jband=1,size(mmn,2);do iband=1,size(mmn,1)
        if(ios==0)write(unit,'(2(1x,es23.15))',iostat=ios)mmn(iband,jband,ineighbor)
      end do;end do
    end do
    close(unit);if(ios/=0)then;message='SAWF local .mmn write failed';ok=.false.;return;end if
    ok=.true.
  end subroutine write_sawf_local_eig_amn_mmn

  subroutine write_sawf_local_eig_amn(directory,seedname,energy_ev,amn,ok,message)
    character(*),intent(in)::directory,seedname
    real(8),intent(in)::energy_ev(:)
    complex(8),intent(in)::amn(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    character(1024)::filename
    integer::unit,ios,iband,iproj
    ok=.false.;message=''
    if(len_trim(directory)==0.or.len_trim(seedname)==0.or.size(energy_ev)<=0.or. &
        size(amn,1)/=size(energy_ev).or.size(amn,2)<=0.or. &
        .not.all(ieee_is_finite(energy_ev)).or..not.all(ieee_is_finite(real(amn))).or. &
        .not.all(ieee_is_finite(aimag(amn))))then
      message='SAWF local eig/amn payload is invalid';return
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

    ok=.true.
  end subroutine write_sawf_local_eig_amn

end module lcfo_wannier_sawf_seed
