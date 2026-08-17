#include "config.h"
module dg_overlapping_wannier_w90
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use lcfo_wannier_sawf_seed,only:write_sawf_local_eig_amn_mmn
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::estimate_dg_w90_coordinator_bytes,validate_dg_w90_result
  public::setup_dg_w90_gamma_library,run_dg_w90_gamma_library
  public::assemble_dg_w90_gamma_matrices
  public::apply_dg_w90_gamma_transform
  public::inherit_dg_w90_affine_receipts
  public::validate_dg_w90_generator_covariance
  public::validate_dg_w90_convergence_log
  public::align_dg_w90_character_sector_gauge
  public::align_dg_w90_cross_character_sector_gauge
  public::align_dg_w90_character_sectors_by_periodic_phase
  public::sew_dg_w90_periodic_phase_conjugate_sector
  public::anchor_dg_w90_reference_character_sector
  public::project_dg_w90_reference_sector_operators
  public::validate_dg_w90_localization_cluster
  public::build_dg_sector_periodic_position_tuple
  public::canonicalize_dg_sector_periodic_position_gauge
  public::jointly_canonicalize_dg_sector_periodic_position_gauge
  public::build_dg_orbital_major_periodic_position_tuple
  public::apply_dg_orbital_rotation_tiled
  public::export_dg_w90_replay_bundle
  public::convert_dg_w90_library_geometry
contains

  subroutine convert_dg_w90_library_geometry(real_lattice_au,reciprocal_lattice_au,atoms_cart_au,&
      real_lattice_angstrom,reciprocal_lattice_inv_angstrom,atoms_cart_angstrom,ok,message)
    real(real64),parameter::bohr_to_angstrom=0.52917721067_real64
    real(real64),intent(in)::real_lattice_au(3,3),reciprocal_lattice_au(3,3),atoms_cart_au(:,:)
    real(real64),intent(out)::real_lattice_angstrom(3,3),reciprocal_lattice_inv_angstrom(3,3),&
      atoms_cart_angstrom(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    ok=.false.;message='';real_lattice_angstrom=0d0;reciprocal_lattice_inv_angstrom=0d0
    atoms_cart_angstrom=0d0
    if(any(shape(atoms_cart_au)/=shape(atoms_cart_angstrom)).or.size(atoms_cart_au,1)/=3.or.&
        .not.all(ieee_is_finite(real_lattice_au)).or.&
        .not.all(ieee_is_finite(reciprocal_lattice_au)).or.&
        .not.all(ieee_is_finite(atoms_cart_au)))then
      message='invalid atomic-unit Wannier90 library geometry';return
    endif
    real_lattice_angstrom=bohr_to_angstrom*real_lattice_au
    reciprocal_lattice_inv_angstrom=reciprocal_lattice_au/bohr_to_angstrom
    atoms_cart_angstrom=bohr_to_angstrom*atoms_cart_au
    if(.not.all(ieee_is_finite(real_lattice_angstrom)).or.&
        .not.all(ieee_is_finite(reciprocal_lattice_inv_angstrom)).or.&
        .not.all(ieee_is_finite(atoms_cart_angstrom)))then
      message='Wannier90 library geometry unit conversion overflowed';return
    endif
    ok=.true.
  end subroutine convert_dg_w90_library_geometry

  subroutine export_dg_w90_replay_bundle(comm,source_directory,source_seed,output_directory,&
      output_seed,energy_ev,amn,mmn,neighbor_gvec,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::source_directory,source_seed,output_directory,output_seed
    real(real64),intent(in)::energy_ev(:)
    complex(real64),intent(in)::amn(:,:),mmn(:,:,:)
    integer,intent(in)::neighbor_gvec(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,ierr,bad,global_bad,nband,nproj,nneighbor
    character(len(source_directory))::agreed_source_directory
    character(len(source_seed))::agreed_source_seed
    character(len(output_directory))::agreed_output_directory
    character(len(output_seed))::agreed_output_seed
    logical::writer_ok
    character(len(message))::writer_message

    ok=.false.;message='';bad=0
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier replay communicator query failed';return;endif
    agreed_source_directory=source_directory;agreed_source_seed=source_seed
    agreed_output_directory=output_directory;agreed_output_seed=output_seed
    call MPI_Bcast(agreed_source_directory,len(agreed_source_directory),MPI_CHARACTER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier replay source-directory agreement failed';return;endif
    call MPI_Bcast(agreed_source_seed,len(agreed_source_seed),MPI_CHARACTER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier replay source-seed agreement failed';return;endif
    call MPI_Bcast(agreed_output_directory,len(agreed_output_directory),MPI_CHARACTER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier replay output-directory agreement failed';return;endif
    call MPI_Bcast(agreed_output_seed,len(agreed_output_seed),MPI_CHARACTER,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier replay output-seed agreement failed';return;endif
    if(source_directory/=agreed_source_directory.or.source_seed/=agreed_source_seed.or.&
        output_directory/=agreed_output_directory.or.output_seed/=agreed_output_seed.or.&
        len_trim(source_directory)==0.or.len_trim(source_seed)==0.or.&
        len_trim(output_directory)==0.or.len_trim(output_seed)==0)bad=1
    nband=0;nproj=0;nneighbor=0
    if(rank==0)then
      nband=size(energy_ev);nproj=size(amn,2);nneighbor=size(mmn,3)
      if(nband<1.or.size(amn,1)/=nband.or.nproj<1.or.size(mmn,1)/=nband.or.&
          size(mmn,2)/=nband.or.nneighbor<1.or.any(shape(neighbor_gvec)/=[3,nneighbor]).or.&
          .not.all(ieee_is_finite(energy_ev)).or..not.all(ieee_is_finite(real(amn))).or.&
          .not.all(ieee_is_finite(aimag(amn))).or..not.all(ieee_is_finite(real(mmn))).or.&
          .not.all(ieee_is_finite(aimag(mmn))))bad=1
    endif
    call MPI_Bcast(nband,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)bad=1
    call MPI_Bcast(nproj,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)bad=1
    call MPI_Bcast(nneighbor,1,MPI_INTEGER,0,comm,ierr);if(ierr/=MPI_SUCCESS)bad=1
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid Wannier replay export contract';return;endif
    writer_ok=.false.;writer_message=''
    if(rank==0)then
      call write_sawf_local_eig_amn_mmn(trim(output_directory),trim(output_seed),energy_ev,amn,mmn,&
        neighbor_gvec,writer_ok,writer_message)
      if(writer_ok)call copy_text_file(trim(source_directory)//'/'//trim(source_seed)//'.win',&
        trim(output_directory)//'/'//trim(output_seed)//'.win',writer_ok,writer_message)
      if(writer_ok)call copy_text_file(trim(source_directory)//'/'//trim(source_seed)//'.dmn',&
        trim(output_directory)//'/'//trim(output_seed)//'.dmn',writer_ok,writer_message)
    endif
    call MPI_Bcast(writer_ok,1,MPI_LOGICAL,0,comm,ierr)
    call MPI_Bcast(writer_message,len(writer_message),MPI_CHARACTER,0,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.writer_ok
    if(.not.ok)then
      if(len_trim(writer_message)>0)then;message=trim(writer_message)
      else;message='Wannier replay export failed collectively';endif
    endif
#else
    ok=.false.;message='Wannier replay export requires MPI'
#endif
  contains
    subroutine copy_text_file(source,target,copy_ok,detail)
      character(*),intent(in)::source,target
      logical,intent(out)::copy_ok
      character(*),intent(out)::detail
      integer::input_unit,output_unit,ios,write_ios
      character(4096)::line
      copy_ok=.false.;detail=''
      open(newunit=input_unit,file=source,status='old',action='read',iostat=ios)
      if(ios/=0)then;detail='Wannier replay source file is unreadable: '//trim(source);return;endif
      open(newunit=output_unit,file=target,status='replace',action='write',iostat=ios)
      if(ios/=0)then;close(input_unit);detail='Wannier replay output file cannot be opened: '//trim(target);return;endif
      do
        read(input_unit,'(a)',iostat=ios)line
        if(ios<0)exit
        if(ios/=0)then;detail='Wannier replay source file read failed: '//trim(source);exit;endif
        write(output_unit,'(a)',iostat=write_ios)trim(line)
        if(write_ios/=0)then;detail='Wannier replay output file write failed: '//trim(target);ios=write_ios;exit;endif
      enddo
      close(input_unit)
      close(output_unit,iostat=write_ios)
      if(ios<0.and.write_ios==0)then;copy_ok=.true.;detail='';endif
    end subroutine copy_text_file
  end subroutine export_dg_w90_replay_bundle

  subroutine validate_dg_w90_generator_covariance(comm,transform,d_band,d_wann,tolerance,&
      covariance_defect,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::transform(:,:),d_band(:,:),d_wann(:,:)
    real(real64),intent(in)::tolerance
    real(real64),intent(out)::covariance_defect
    integer(int64),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,ierr,status,nstate,allocation_status
    integer(int64)::elements
    real(real64)::minimum_tolerance,maximum_tolerance
    complex(real64),allocatable::image(:,:),rotated(:,:)
    ok=.false.;message='';covariance_defect=huge(1d0);workspace_peak_bytes=0_int64;status=0
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier covariance communicator query failed';return;endif
    call MPI_Allreduce(tolerance,minimum_tolerance,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier covariance tolerance MIN reduction failed';return;endif
    call MPI_Allreduce(tolerance,maximum_tolerance,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.transfer(minimum_tolerance,0_int64)/=transfer(maximum_tolerance,0_int64).or.&
        .not.ieee_is_finite(tolerance).or.tolerance<=0d0)status=1
    nstate=0
    if(rank==0)then
      nstate=size(transform,1)
      if(nstate<1.or.size(transform,2)/=nstate.or.any(shape(d_band)/=[nstate,nstate]).or.&
          any(shape(d_wann)/=[nstate,nstate]).or..not.all(ieee_is_finite(real(transform))).or.&
          .not.all(ieee_is_finite(aimag(transform))).or..not.all(ieee_is_finite(real(d_band))).or.&
          .not.all(ieee_is_finite(aimag(d_band))).or..not.all(ieee_is_finite(real(d_wann))).or.&
          .not.all(ieee_is_finite(aimag(d_wann))))status=1
    endif
    call MPI_Bcast(nstate,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.status/=0)then;message='invalid coordinator Wannier covariance contract';return;endif
    if(int(nstate,int64)>huge(0_int64)/int(nstate,int64)/2_int64)status=1
    elements=2_int64*int(nstate,int64)*int(nstate,int64)
    if(elements>huge(0_int64)/16_int64)status=1
    allocation_status=0
    if(rank==0.and.status==0)allocate(image(nstate,nstate),rotated(nstate,nstate),stat=allocation_status)
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,allocation_status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.status/=0.or.allocation_status/=0)then
      if(allocated(image))deallocate(image)
      if(allocated(rotated))deallocate(rotated)
      message='Wannier covariance workspace allocation or extent failed';return
    endif
    if(rank==0)then
      image=matmul(d_band,transform)
      rotated=matmul(conjg(transpose(transform)),image)
      covariance_defect=maxval(abs(rotated-d_wann))
      workspace_peak_bytes=16_int64*elements
      if(.not.ieee_is_finite(covariance_defect))status=1
      deallocate(image,rotated)
    endif
    call MPI_Bcast(covariance_defect,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(workspace_peak_bytes,1,MPI_INTEGER8,0,comm,ierr)
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.status==0.and.covariance_defect<=tolerance
    if(.not.ok)message='Wannier transform violates the supplied generator covariance'
#else
    ok=.false.;message='Wannier covariance validation requires MPI'
    covariance_defect=huge(1d0);workspace_peak_bytes=0_int64
#endif
  end subroutine validate_dg_w90_generator_covariance

  subroutine build_dg_orbital_major_periodic_position_tuple(comm,orbital_values,weights,&
      periodic_phases,tolerance,provenance_fingerprint,position_tuple,gram_defect,fingerprint,&
      workspace_peak_bytes,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::orbital_values(:,:),periodic_phases(:,:)
    real(real64),intent(in)::weights(:),tolerance
    integer(int64),intent(in)::provenance_fingerprint
    complex(real64),allocatable,intent(out)::position_tuple(:,:,:)
    real(real64),intent(out)::gram_defect
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,parameter::tile_size=256
    integer::m,nlocal,first,count,i,j,axis,ierr,bad,gbad,status,minm,maxm
    complex(real64),allocatable::local_tuple(:,:,:),gram(:,:),left_tile(:,:)
    real(real64)::mintol,maxtol,scale,global_scale,local_weight,global_weight,safe_scale,quantum,value
    integer(int64)::minhash,maxhash,global_points,local_points,elements,term,bytes,quantized
    logical::receipt_valid
    m=size(orbital_values,1);nlocal=size(orbital_values,2)
    ok=.false.;message='';gram_defect=huge(1d0);fingerprint=0_int64;workspace_peak_bytes=0_int64
    bad=merge(0,1,m>0.and.nlocal>0.and.size(weights)==nlocal.and.&
      all(shape(periodic_phases)==[3,nlocal]).and.provenance_fingerprint/=0_int64.and.&
      tolerance>=1d-15.and.tolerance<=1d-2.and.ieee_is_finite(tolerance).and.&
      all(weights>0d0).and.all(ieee_is_finite(weights)).and.&
      all(ieee_is_finite(real(orbital_values))).and.all(ieee_is_finite(aimag(orbital_values))).and.&
      all(ieee_is_finite(real(periodic_phases))).and.all(ieee_is_finite(aimag(periodic_phases))).and.&
      maxval(abs(abs(periodic_phases)-1d0))<=10d0*tolerance)
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='invalid orbital-major periodic-position contract';return;endif
    call MPI_Allreduce(m,minm,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(m,maxm,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minm/=maxm)then;message='orbital-major periodic-position rank disagrees';return;endif
    call MPI_Allreduce(tolerance,mintol,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tolerance,maxtol,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.mintol/=maxtol)then;message='orbital-major periodic-position tolerance disagrees';return;endif
    call MPI_Allreduce(provenance_fingerprint,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(provenance_fingerprint,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)then;message='orbital-major periodic-position provenance disagrees';return;endif
    local_points=int(nlocal,int64);call MPI_Allreduce(local_points,global_points,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    scale=maxval(abs(orbital_values));local_weight=maxval(weights)
    call MPI_Allreduce(scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_weight,global_weight,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    safe_scale=sqrt((huge(1d0)/real(max(1_int64,global_points),real64))/global_weight)/4d0
    bad=merge(0,1,ierr==MPI_SUCCESS.and.global_scale<=safe_scale)
    elements=0_int64;receipt_valid=.true.
    call checked_product([5_int64,int(m,int64),int(m,int64)],term,receipt_valid);call checked_add(elements,term,receipt_valid)
    call checked_product([int(m,int64),int(min(tile_size,nlocal),int64)],term,receipt_valid);call checked_add(elements,term,receipt_valid)
    call checked_product([elements,16_int64],bytes,receipt_valid)
    if(.not.receipt_valid)bad=1
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='orbital-major periodic-position magnitude or workspace overflows';return;endif
    workspace_peak_bytes=bytes
    allocate(position_tuple(m,m,3),local_tuple(m,m,3),gram(m,m),left_tile(m,min(tile_size,nlocal)),stat=status)
    call MPI_Allreduce(status,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;call cleanup();message='orbital-major periodic-position allocation failed';return;endif
    gram=(0d0,0d0);local_tuple=(0d0,0d0)
    do first=1,nlocal,tile_size
      count=min(tile_size,nlocal-first+1)
      do j=1,count
        left_tile(:,j)=weights(first+j-1)*conjg(orbital_values(:,first+j-1))
      enddo
      call zgemm('N','T',m,m,count,(1d0,0d0),left_tile,m,orbital_values(:,first:first+count-1),m,&
        (1d0,0d0),gram,m)
      do axis=1,3
        do j=1,count
          left_tile(:,j)=weights(first+j-1)*periodic_phases(axis,first+j-1)*&
            conjg(orbital_values(:,first+j-1))
        enddo
        call zgemm('N','T',m,m,count,(1d0,0d0),left_tile,m,orbital_values(:,first:first+count-1),m,&
          (1d0,0d0),local_tuple(:,:,axis),m)
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    call MPI_Allreduce(local_tuple,position_tuple,3*m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    gram_defect=0d0
    do j=1,m;do i=1,m
      gram_defect=max(gram_defect,abs(gram(i,j)-merge((1d0,0d0),(0d0,0d0),i==j)))
    enddo;enddo
    call MPI_Allreduce(MPI_IN_PLACE,gram_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gram_defect>10d0*tolerance.or.&
        .not.all(ieee_is_finite(real(position_tuple))).or..not.all(ieee_is_finite(aimag(position_tuple))))then
      call cleanup();message='orbital-major periodic-position tuple or Gram is invalid';return
    endif
    fingerprint=ieor(provenance_fingerprint,int(m,int64));quantum=100d0*tolerance
    do axis=1,3;do j=1,m;do i=1,m
      value=real(position_tuple(i,j,axis),real64)
      if(abs(value)/quantum>0.25d0*real(huge(0_int64),real64))then;call cleanup();message='orbital-major tuple fingerprint overflows';return;endif
      quantized=nint(value/quantum,int64);fingerprint=ieor(ishftc(fingerprint,11),quantized)
      value=aimag(position_tuple(i,j,axis))
      if(abs(value)/quantum>0.25d0*real(huge(0_int64),real64))then;call cleanup();message='orbital-major tuple fingerprint overflows';return;endif
      quantized=nint(value/quantum,int64);fingerprint=ieor(ishftc(fingerprint,7),quantized)
    enddo;enddo;enddo
    if(fingerprint==0_int64)fingerprint=ieor(provenance_fingerprint,7927_int64)
    ok=.true.;message='';call cleanup(.false.)
  contains
    subroutine cleanup(remove_output)
      logical,intent(in),optional::remove_output
      logical::drop
      drop=.true.;if(present(remove_output))drop=remove_output
      if(drop.and.allocated(position_tuple))deallocate(position_tuple)
      if(allocated(local_tuple))deallocate(local_tuple)
      if(allocated(gram))deallocate(gram)
      if(allocated(left_tile))deallocate(left_tile)
    end subroutine
#else
    ok=.false.;message='orbital-major periodic-position tuple requires MPI';gram_defect=huge(1d0)
    fingerprint=0_int64;workspace_peak_bytes=0_int64
#endif
  end subroutine build_dg_orbital_major_periodic_position_tuple

  subroutine apply_dg_orbital_rotation_tiled(comm,orbital_values,rotation,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(inout)::orbital_values(:,:)
    complex(real64),intent(in)::rotation(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,parameter::tile_size=256
    complex(real64),allocatable::tile(:,:)
    integer::m,nlocal,first,count,status,ierr,bad,gbad
    m=size(orbital_values,1);nlocal=size(orbital_values,2);ok=.false.;message=''
    bad=merge(0,1,m>=1.and.nlocal>=1.and.all(shape(rotation)==[m,m]).and.&
      all(ieee_is_finite(real(rotation))).and.all(ieee_is_finite(aimag(rotation))))
#ifdef USE_MPI
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='invalid orbital rotation contract';return;endif
#else
    if(bad/=0)then;message='invalid orbital rotation contract';return;endif
#endif
    allocate(tile(m,min(tile_size,nlocal)),stat=status)
#ifdef USE_MPI
    call MPI_Allreduce(status,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then
      if(allocated(tile))deallocate(tile);message='orbital rotation tile allocation failed';return
    endif
#else
    if(status/=0)then;message='orbital rotation tile allocation failed';return;endif
#endif
    do first=1,nlocal,tile_size
      count=min(tile_size,nlocal-first+1)
      call zgemm('T','N',m,count,m,(1d0,0d0),rotation,m,orbital_values(:,first:first+count-1),m,&
        (0d0,0d0),tile,m)
      orbital_values(:,first:first+count-1)=tile(:,1:count)
    enddo
    ok=all(ieee_is_finite(real(orbital_values))).and.all(ieee_is_finite(aimag(orbital_values)))
#ifdef USE_MPI
    bad=merge(0,1,ok);call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr);ok=ierr==MPI_SUCCESS.and.gbad==0
#endif
    if(.not.ok)message='orbital rotation produced nonfinite values'
    deallocate(tile)
  end subroutine apply_dg_orbital_rotation_tiled

  subroutine build_dg_sector_periodic_position_tuple(comm,row_ids,global_row_count,sector_rows,&
      periodic_phases,tolerance,provenance_fingerprint,position_tuple,gram_defect,fingerprint,&
      workspace_peak_bytes,ok,message,integration_weights)
    integer,intent(in)::comm,global_row_count
    integer(int64),intent(in)::row_ids(:),provenance_fingerprint
    complex(real64),intent(in)::sector_rows(:,:),periodic_phases(:,:)
    real(real64),intent(in)::tolerance
    real(real64),intent(in),optional::integration_weights(:)
    complex(real64),allocatable,intent(out)::position_tuple(:,:,:)
    real(real64),intent(out)::gram_defect
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nlocal,m,i,j,axis,rank,ierr,bad,gbad,status,minint,maxint
    integer,allocatable::owner(:),position(:),count(:)
    complex(real64),allocatable::local_tuple(:,:,:),gram(:,:),phased(:,:)
    real(real64)::mintol,maxtol,local_scale,global_scale,safe_scale,quantum,value,local_weight,global_weight
    integer(int64)::minhash,maxhash,elements,term,bytes,bits,quantized
    logical::receipt_valid
    nlocal=size(row_ids);m=size(sector_rows,2)
    ok=.false.;message='';gram_defect=huge(1d0);fingerprint=0_int64;workspace_peak_bytes=0_int64
    bad=merge(0,1,global_row_count>=1.and.nlocal>=1.and.m>=1.and.size(sector_rows,1)==nlocal.and.&
      all(shape(periodic_phases)==[nlocal,3]).and.provenance_fingerprint/=0_int64.and.&
      tolerance>=1d-15.and.tolerance<=1d-2.and.ieee_is_finite(tolerance).and.&
      all(row_ids>=1_int64).and.all(row_ids<=int(global_row_count,int64)).and.&
      all(ieee_is_finite(real(sector_rows))).and.all(ieee_is_finite(aimag(sector_rows))).and.&
      all(ieee_is_finite(real(periodic_phases))).and.all(ieee_is_finite(aimag(periodic_phases))).and.&
      maxval(abs(abs(periodic_phases)-1d0))<=10d0*tolerance)
    if(present(integration_weights))then
      if(size(integration_weights)/=nlocal)bad=1
      if(size(integration_weights)==nlocal)then
        if(.not.all(ieee_is_finite(integration_weights)).or.any(integration_weights<=0d0))bad=1
      endif
    endif
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='invalid sector periodic-position tuple contract';return;endif
    do i=1,2
      j=merge(global_row_count,m,i==1)
      call MPI_Allreduce(j,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(j,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='periodic-position tuple dimensions disagree';return;endif
    enddo
    call MPI_Allreduce(tolerance,mintol,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tolerance,maxtol,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.mintol/=maxtol)then;message='periodic-position tuple tolerance disagrees';return;endif
    call MPI_Allreduce(provenance_fingerprint,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(provenance_fingerprint,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)then;message='periodic-position provenance disagrees';return;endif
    local_scale=maxval(abs(sector_rows));local_weight=1d0
    if(present(integration_weights))local_weight=maxval(integration_weights)
    call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(local_weight,global_weight,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    safe_scale=sqrt((huge(1d0)/real(global_row_count,real64))/global_weight)/4d0
    bad=merge(0,1,ierr==MPI_SUCCESS.and.global_scale<=safe_scale)
    elements=0_int64;bytes=0_int64;receipt_valid=.true.
    call checked_product([7_int64,int(m,int64),int(m,int64)],term,receipt_valid);call checked_add(elements,term,receipt_valid)
    call checked_product([int(nlocal,int64),int(m,int64)],term,receipt_valid);call checked_add(elements,term,receipt_valid)
    call checked_product([elements,16_int64],bytes,receipt_valid)
    call checked_product([3_int64,int(global_row_count,int64),4_int64],term,receipt_valid);call checked_add(bytes,term,receipt_valid)
    if(.not.receipt_valid)bad=1
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='sector periodic-position magnitude or workspace overflows';return;endif
    workspace_peak_bytes=bytes
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(position_tuple(m,m,3),local_tuple(m,m,3),gram(m,m),phased(nlocal,m),&
      owner(global_row_count),position(global_row_count),count(global_row_count),stat=status)
    call MPI_Allreduce(status,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then
      call cleanup();message='sector periodic-position allocation failed';return
    endif
    owner=0;position=0;count=0
    do i=1,nlocal
      owner(int(row_ids(i)))=rank+1;position(int(row_ids(i)))=i
      count(int(row_ids(i)))=count(int(row_ids(i)))+1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,owner,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)then;call cleanup();return;endif
    call MPI_Allreduce(MPI_IN_PLACE,position,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)then;call cleanup();return;endif
    call MPI_Allreduce(MPI_IN_PLACE,count,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    bad=merge(0,1,ierr==MPI_SUCCESS.and.all(count==1))
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;call cleanup();message='sector periodic-position row ownership is incomplete';return;endif
    if(present(integration_weights))then
      gram=matmul(conjg(transpose(sector_rows)),spread(integration_weights,2,m)*sector_rows)
    else
      gram=matmul(conjg(transpose(sector_rows)),sector_rows)
    endif
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    gram_defect=0d0
    do j=1,m;do i=1,m
      gram_defect=max(gram_defect,abs(gram(i,j)-merge((1d0,0d0),(0d0,0d0),i==j)))
    enddo;enddo
    call MPI_Allreduce(MPI_IN_PLACE,gram_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gram_defect>10d0*tolerance)then;call cleanup();message='sector periodic-position frame is not orthonormal';return;endif
    do axis=1,3
      do j=1,m
        phased(:,j)=periodic_phases(:,axis)*sector_rows(:,j)
        if(present(integration_weights))phased(:,j)=integration_weights*phased(:,j)
      enddo
      local_tuple(:,:,axis)=matmul(conjg(transpose(sector_rows)),phased)
    enddo
    call MPI_Allreduce(local_tuple,position_tuple,3*m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.all(ieee_is_finite(real(position_tuple))).or.&
        .not.all(ieee_is_finite(aimag(position_tuple))))then
      call cleanup();message='sector periodic-position tuple is nonfinite';return
    endif
    fingerprint=ieor(provenance_fingerprint,int(m,int64));quantum=100d0*tolerance
    do axis=1,3;do j=1,m;do i=1,m
      value=real(position_tuple(i,j,axis),real64)
      if(abs(value)/quantum>0.25d0*real(huge(0_int64),real64))then;call cleanup();message='periodic-position fingerprint range overflows';return;endif
      quantized=nint(value/quantum,int64);fingerprint=ieor(ishftc(fingerprint,11),quantized)
      value=aimag(position_tuple(i,j,axis))
      if(abs(value)/quantum>0.25d0*real(huge(0_int64),real64))then;call cleanup();message='periodic-position fingerprint range overflows';return;endif
      quantized=nint(value/quantum,int64);fingerprint=ieor(ishftc(fingerprint,7),quantized)
    enddo;enddo;enddo
    if(fingerprint==0_int64)fingerprint=ieor(provenance_fingerprint,7919_int64)
    call MPI_Allreduce(workspace_peak_bytes,term,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)workspace_peak_bytes=term
    ok=ierr==MPI_SUCCESS
    if(ok)message=''
    call cleanup(.false.)
  contains
    subroutine cleanup(remove_output)
      logical,intent(in),optional::remove_output
      logical::drop
      drop=.true.;if(present(remove_output))drop=remove_output
      if(drop.and.allocated(position_tuple))deallocate(position_tuple)
      if(allocated(local_tuple))deallocate(local_tuple)
      if(allocated(gram))deallocate(gram)
      if(allocated(phased))deallocate(phased)
      if(allocated(owner))deallocate(owner)
      if(allocated(position))deallocate(position)
      if(allocated(count))deallocate(count)
    end subroutine
#else
    ok=.false.;message='sector periodic-position tuple requires MPI';gram_defect=huge(1d0)
    fingerprint=0_int64;workspace_peak_bytes=0_int64
#endif
  end subroutine build_dg_sector_periodic_position_tuple

  subroutine canonicalize_dg_sector_periodic_position_gauge(comm,row_ids,sector_rows,position_tuple,&
      lcfo_operator,tolerance,tuple_fingerprint,aligned_rows,rotation,canonical_defect,fingerprint,&
      workspace_peak_bytes,ok,message)
    integer,intent(in)::comm
    integer(int64),intent(in)::row_ids(:),tuple_fingerprint
    complex(real64),intent(in)::sector_rows(:,:),position_tuple(:,:,:),lcfo_operator(:,:)
    real(real64),intent(in)::tolerance
    complex(real64),allocatable,intent(out)::aligned_rows(:,:)
    complex(real64),intent(out)::rotation(:,:)
    real(real64),intent(out)::canonical_defect
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::discriminator(:,:),local_rotation(:,:)
    real(real64),parameter::alpha(3)=[sqrt(2d0),sqrt(3d0),sqrt(5d0)]
    real(real64),parameter::beta(3)=[sqrt(7d0),sqrt(11d0),sqrt(13d0)]
    real(real64)::mintol,maxtol,invariant,quantum
    integer::nlocal,m,axis,i,j,ierr,bad,gbad,status,local_count,global_count,minint,maxint
    integer(int64)::minhash,maxhash,invariant_hash,bits,extra_bytes,peak
    logical::receipt_valid
    nlocal=size(row_ids);m=size(sector_rows,2)
    ok=.false.;message='';canonical_defect=huge(1d0);fingerprint=0_int64;workspace_peak_bytes=0_int64
    rotation=(0d0,0d0)
    bad=merge(0,1,m>=1.and.size(sector_rows,1)==nlocal.and.all(shape(position_tuple)==[m,m,3]).and.&
      all(shape(lcfo_operator)==[m,m]).and.all(shape(rotation)==[m,m]).and.tuple_fingerprint/=0_int64.and.&
      tolerance>=1d-15.and.tolerance<=1d-2.and.ieee_is_finite(tolerance).and.&
      all(row_ids>=1_int64).and.all(ieee_is_finite(real(sector_rows))).and.&
      all(ieee_is_finite(aimag(sector_rows))).and.all(ieee_is_finite(real(position_tuple))).and.&
      all(ieee_is_finite(aimag(position_tuple))).and.all(ieee_is_finite(real(lcfo_operator))).and.&
      all(ieee_is_finite(aimag(lcfo_operator))).and.&
      maxval(abs(lcfo_operator-conjg(transpose(lcfo_operator))))<=10d0*tolerance)
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='invalid periodic-position canonical gauge contract';return;endif
    call MPI_Allreduce(m,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(m,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='periodic-position canonical rank disagrees';return;endif
    call MPI_Allreduce(tolerance,mintol,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tolerance,maxtol,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.mintol/=maxtol)then;message='periodic-position canonical tolerance disagrees';return;endif
    call MPI_Allreduce(tuple_fingerprint,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tuple_fingerprint,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)then;message='periodic-position tuple receipt disagrees';return;endif
    do axis=1,3;do j=1,m;do i=1,m
      bits=transfer(real(position_tuple(i,j,axis),real64),bits)
      call MPI_Allreduce(bits,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(bits,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)return
      bits=transfer(aimag(position_tuple(i,j,axis)),bits)
      call MPI_Allreduce(bits,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(bits,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)return
    enddo;enddo;enddo
    local_count=0;if(nlocal>0)local_count=int(maxval(row_ids))
    call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_count<1)then;message='periodic-position canonical rows are empty';return;endif
    allocate(discriminator(m,m),local_rotation(m,m),stat=status)
    call MPI_Allreduce(status,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then
      if(allocated(discriminator))deallocate(discriminator)
      if(allocated(local_rotation))deallocate(local_rotation)
      message='periodic-position canonical allocation failed';return
    endif
    discriminator=(0d0,0d0)
    do axis=1,3
      discriminator=discriminator+0.5d0*alpha(axis)*(position_tuple(:,:,axis)+&
        conjg(transpose(position_tuple(:,:,axis))))+cmplx(0d0,-0.5d0*beta(axis),real64)*&
        (position_tuple(:,:,axis)-conjg(transpose(position_tuple(:,:,axis))))
    enddo
    quantum=100d0*tolerance;invariant_hash=ieor(int(z'6A09E667F3BCC909',int64),int(m,int64))
    do axis=1,3
      invariant=real(sum([(position_tuple(i,i,axis),i=1,m)]),real64)
      if(abs(invariant)/quantum>0.25d0*real(huge(0_int64),real64))then
        deallocate(discriminator,local_rotation);message='periodic-position invariant quantization overflows';return
      endif
      bits=nint(invariant/quantum,int64);invariant_hash=ieor(ishftc(invariant_hash,11),bits)
      invariant=sum(abs(position_tuple(:,:,axis))**2)
      if(abs(invariant)/quantum>0.25d0*real(huge(0_int64),real64))then
        deallocate(discriminator,local_rotation);message='periodic-position invariant quantization overflows';return
      endif
      bits=nint(invariant/quantum,int64);invariant_hash=ieor(ishftc(invariant_hash,11),bits)
    enddo
    if(invariant_hash==0_int64)invariant_hash=1_int64
    call anchor_dg_w90_reference_character_sector(comm,row_ids,sector_rows,discriminator,lcfo_operator,&
      global_count,invariant_hash,invariant_hash,0d0,0d0,tolerance,aligned_rows,canonical_defect,&
      fingerprint,workspace_peak_bytes,ok,message,projector_diagonal_only=.true.)
    if(ok)then
      local_rotation=matmul(conjg(transpose(sector_rows)),aligned_rows)
      call MPI_Allreduce(local_rotation,rotation,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      ok=ierr==MPI_SUCCESS
      if(.not.ok)message='periodic-position canonical rotation reduction failed'
    endif
    receipt_valid=.true.;call checked_product([2_int64,int(m,int64),int(m,int64),16_int64],extra_bytes,receipt_valid)
    if(receipt_valid.and.ok)then
      if(workspace_peak_bytes<=huge(0_int64)-extra_bytes)workspace_peak_bytes=workspace_peak_bytes+extra_bytes
      call MPI_Allreduce(workspace_peak_bytes,peak,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
      if(ierr==MPI_SUCCESS)workspace_peak_bytes=peak
    endif
    deallocate(discriminator,local_rotation)
#else
    ok=.false.;message='periodic-position canonical gauge requires MPI';canonical_defect=huge(1d0)
    fingerprint=0_int64;workspace_peak_bytes=0_int64;rotation=(0d0,0d0)
#endif
  end subroutine canonicalize_dg_sector_periodic_position_gauge

  subroutine jointly_canonicalize_dg_sector_periodic_position_gauge(comm,row_ids,sector_rows,position_tuple,&
      lcfo_operator,tolerance,tuple_fingerprint,aligned_rows,rotation,centers,final_objective,maximum_update,&
      sweep_count,canonical_defect,fingerprint,workspace_peak_bytes,ok,message,point_representations)
    integer,intent(in)::comm
    integer(int64),intent(in)::row_ids(:),tuple_fingerprint
    complex(real64),intent(in)::sector_rows(:,:),position_tuple(:,:,:),lcfo_operator(:,:)
    complex(real64),intent(in),optional::point_representations(:,:,:)
    real(real64),intent(in)::tolerance
    complex(real64),allocatable,intent(out)::aligned_rows(:,:)
    complex(real64),intent(out)::rotation(:,:)
    real(real64),allocatable,intent(out)::centers(:,:)
    real(real64),intent(out)::final_objective,maximum_update,canonical_defect
    integer,intent(out)::sweep_count
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::hmats(:,:,:),unitary(:,:),gram(:,:),stream(:),block(:,:),zwork(:),&
      orbit_vectors(:,:),orbit_basis(:,:),orbit_residual(:),cover_basis(:,:)
    integer,allocatable::owner(:),position(:),count(:),order(:),orbit_labels(:),cluster_sequence(:),basis_labels(:)
    integer(int64),allocatable::payload_bits(:),payload_minimum(:),payload_maximum(:)
    real(real64),allocatable::block_eval(:),zrwork(:),orbit_centers(:,:),cluster_centers(:,:),inverse_sqrt(:)
    real(real64)::gmat(3,3),geval(3),gwork(9),gvec(3),two_theta,c,sabs,phase_angle,pivot,&
      local_objective,global_objective,local_update,global_update,quantum
    real(real64)::previous_objective,objective_scale,objective_threshold,objective_change
    complex(real64)::sphase,jacobi(2,2),left_pair(2),right_pair(2),tmp,phase_fix,probe
    integer::nlocal,m,global_count,local_count,rank,ierr,bad,gbad,status,axis,q,pair_i,pair_j,&
      i,j,k,l,r,block_size,sweep,info,minint,maxint,payload_count,npoint,point,hmatrix_count,&
      maximum_sweeps,sweep_bad
    integer::worst_point,worst_column
    real(real64)::minimum_tolerance,maximum_tolerance,safe_position_magnitude,safe_lcfo_magnitude
    real(real64)::gram_defect,point_leakage
    integer(int64)::bits,quantized,term,elements,bytes,peak,minimum_fingerprint,maximum_fingerprint,&
      payload_elements
    logical::receipt_valid,swapped,point_orbit_built
    interface
      subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
        character,intent(in)::jobz,uplo
        integer,intent(in)::n,lda,lwork
        real(8),intent(inout)::a(lda,*),work(*)
        real(8),intent(out)::w(*)
        integer,intent(out)::info
      end subroutine
      subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
        character,intent(in)::jobz,uplo
        integer,intent(in)::n,lda,lwork
        complex(8),intent(inout)::a(lda,*),work(*)
        real(8),intent(out)::w(*),rwork(*)
        integer,intent(out)::info
      end subroutine
    end interface
    nlocal=size(row_ids);m=size(sector_rows,2);npoint=1
    if(present(point_representations))npoint=size(point_representations,3)
    ok=.false.;message='';final_objective=huge(1d0);maximum_update=huge(1d0);canonical_defect=huge(1d0)
    fingerprint=0_int64;workspace_peak_bytes=0_int64;sweep_count=0;rotation=(0d0,0d0)
    bad=merge(0,1,nlocal>=1.and.m>=1.and.size(sector_rows,1)==nlocal.and.&
      all(shape(position_tuple)==[m,m,3]).and.all(shape(lcfo_operator)==[m,m]).and.&
      all(shape(rotation)==[m,m]).and.tuple_fingerprint/=0_int64.and.&
      tolerance>=1d-15.and.tolerance<=1d-2.and.ieee_is_finite(tolerance).and.&
      all(row_ids>=1_int64).and.all(ieee_is_finite(real(sector_rows))).and.&
      all(ieee_is_finite(aimag(sector_rows))).and.all(ieee_is_finite(real(position_tuple))).and.&
      all(ieee_is_finite(aimag(position_tuple))).and.all(ieee_is_finite(real(lcfo_operator))).and.&
      all(ieee_is_finite(aimag(lcfo_operator))))
    if(present(point_representations))then
      if(any(shape(point_representations,kind=int64)/=[int(m,int64),int(m,int64),int(npoint,int64)]).or.npoint<1)bad=1
      if(bad==0)then
        if(.not.all(ieee_is_finite(real(point_representations))).or.&
            .not.all(ieee_is_finite(aimag(point_representations))))bad=1
        if(maxval(abs(point_representations))>2d0)bad=1
      endif
    endif
    if(npoint>huge(0)/6)bad=1
    if(npoint>0)then
      if(m>huge(0)/npoint)bad=1
    endif
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='invalid joint periodic-center gauge contract';return;endif
    call MPI_Allreduce(m,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(m,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='joint periodic-center rank disagrees';return;endif
    call MPI_Allreduce(npoint,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(npoint,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='joint periodic-center point count disagrees';return;endif
    call MPI_Allreduce(tolerance,minimum_tolerance,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tolerance,maximum_tolerance,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_tolerance/=maximum_tolerance)then
      message='joint periodic-center tolerance disagrees';return
    endif
    call MPI_Allreduce(tuple_fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tuple_fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_fingerprint/=maximum_fingerprint)then
      message='joint periodic-center provenance disagrees';return
    endif
    safe_position_magnitude=sqrt(huge(1d0))/32d0
    safe_lcfo_magnitude=(huge(1d0)/32d0)/real(m,real64)
    bad=merge(0,1,maxval(abs(position_tuple))<=safe_position_magnitude.and.&
      maxval(abs(lcfo_operator))<=safe_lcfo_magnitude)
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then
      message='joint periodic-center input magnitude is unsafe';return
    endif
    local_count=int(maxval(row_ids));call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    elements=0_int64;bytes=0_int64;receipt_valid=ierr==MPI_SUCCESS
    if(receipt_valid)call checked_product([8_int64+2_int64*int(npoint,int64),int(m,int64),int(m,int64)],&
      payload_elements,receipt_valid)
    if(receipt_valid)receipt_valid=payload_elements<=int(huge(0),int64)
    if(receipt_valid)call checked_product([3_int64+6_int64*int(npoint,int64),int(m,int64),int(m,int64)],&
      term,receipt_valid)
    call checked_add(elements,term,receipt_valid)
    if(receipt_valid)call checked_product([int(m,int64),int(npoint,int64)],term,receipt_valid)
    call checked_add(elements,term,receipt_valid)
    if(receipt_valid)call checked_product([int(m,int64),int(m,int64)],term,receipt_valid)
    call checked_add(elements,term,receipt_valid)
    if(receipt_valid)call checked_product([int(m,int64),int(m,int64)],term,receipt_valid)
    call checked_add(elements,term,receipt_valid)
    call checked_add(elements,int(m,int64),receipt_valid)
    if(receipt_valid)call checked_product([int(nlocal,int64),int(m,int64)],term,receipt_valid)
    call checked_add(elements,term,receipt_valid)
    if(receipt_valid)call checked_product([3_int64,int(m,int64)],term,receipt_valid)
    call checked_add(elements,term,receipt_valid)
    if(receipt_valid)call checked_product([elements,16_int64],bytes,receipt_valid)
    if(receipt_valid)call checked_product([7_int64,int(m,int64),8_int64],term,receipt_valid)
    call checked_add(bytes,term,receipt_valid)
    if(receipt_valid)call checked_product([3_int64,1_int64+int(m,int64),int(npoint,int64),8_int64],term,receipt_valid)
    call checked_add(bytes,term,receipt_valid)
    if(receipt_valid)call checked_product([int(m,int64),8_int64],term,receipt_valid)
    call checked_add(bytes,term,receipt_valid)
    if(receipt_valid)call checked_product([4_int64,int(global_count,int64),4_int64],term,receipt_valid)
    call checked_add(bytes,term,receipt_valid)
    if(receipt_valid)call checked_product([1_int64+int(m,int64),int(npoint,int64),4_int64],term,receipt_valid)
    call checked_add(bytes,term,receipt_valid)
    if(receipt_valid)call checked_product([int(m,int64),4_int64],term,receipt_valid)
    call checked_add(bytes,term,receipt_valid)
    if(receipt_valid)call checked_product([3_int64,payload_elements,8_int64],term,receipt_valid)
    call checked_add(bytes,term,receipt_valid)
    bad=merge(0,1,receipt_valid)
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='joint periodic-center workspace overflows';return;endif
    workspace_peak_bytes=bytes
    payload_count=int(payload_elements)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    hmatrix_count=6*npoint
    allocate(hmats(m,m,hmatrix_count),unitary(m,m),gram(m,m),stream(m),block(m,m),zwork(max(1,2*m)),&
      block_eval(m),zrwork(max(1,3*m)),aligned_rows(nlocal,m),centers(3,m),&
      owner(global_count),position(global_count),count(global_count),order(m),payload_bits(payload_count),&
      payload_minimum(payload_count),payload_maximum(payload_count),orbit_vectors(m,npoint),orbit_basis(m,m),&
      orbit_residual(m),cover_basis(m,m),orbit_centers(3,npoint),cluster_centers(3,m*npoint),inverse_sqrt(m),&
      orbit_labels(npoint),cluster_sequence(m*npoint),basis_labels(m),stat=status)
    call MPI_Allreduce(status,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;call cleanup();message='joint periodic-center allocation failed';return;endif
    owner=0;position=0;count=0
    do i=1,nlocal
      if(row_ids(i)>int(global_count,int64))then;bad=1;cycle;endif
      owner(int(row_ids(i)))=rank+1;position(int(row_ids(i)))=i;count(int(row_ids(i)))=count(int(row_ids(i)))+1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,owner,global_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)then;call cleanup();return;endif
    call MPI_Allreduce(MPI_IN_PLACE,position,global_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)then;call cleanup();return;endif
    call MPI_Allreduce(MPI_IN_PLACE,count,global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(count/=1))then;call cleanup();message='joint periodic-center rows are not uniquely owned';return;endif
    k=0
    do axis=1,3;do j=1,m;do i=1,m
      k=k+1;payload_bits(k)=transfer(real(position_tuple(i,j,axis),real64),bits)
      k=k+1;payload_bits(k)=transfer(aimag(position_tuple(i,j,axis)),bits)
    enddo;enddo;enddo
    do j=1,m;do i=1,m
      k=k+1;payload_bits(k)=transfer(real(lcfo_operator(i,j),real64),bits)
      k=k+1;payload_bits(k)=transfer(aimag(lcfo_operator(i,j)),bits)
    enddo;enddo
    if(present(point_representations))then
      do point=1,npoint;do j=1,m;do i=1,m
        k=k+1;payload_bits(k)=transfer(real(point_representations(i,j,point),real64),bits)
        k=k+1;payload_bits(k)=transfer(aimag(point_representations(i,j,point)),bits)
      enddo;enddo;enddo
    else
      do j=1,m;do i=1,m
        k=k+1;payload_bits(k)=transfer(merge(1d0,0d0,i==j),bits)
        k=k+1;payload_bits(k)=transfer(0d0,bits)
      enddo;enddo
    endif
    call MPI_Allreduce(payload_bits,payload_minimum,payload_count,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();return;endif
    call MPI_Allreduce(payload_bits,payload_maximum,payload_count,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(payload_minimum/=payload_maximum))then
      call cleanup();message='joint periodic-center payload disagrees';return
    endif
    do axis=1,3
      gram=0.5d0*(position_tuple(:,:,axis)+conjg(transpose(position_tuple(:,:,axis))))
      do point=1,npoint
        q=6*(point-1)+2*axis-1
        if(present(point_representations))then
          hmats(:,:,q)=matmul(conjg(transpose(point_representations(:,:,point))),&
            matmul(gram,point_representations(:,:,point)))
        else
          hmats(:,:,q)=gram
        endif
      enddo
      gram=cmplx(0d0,-0.5d0,real64)*(position_tuple(:,:,axis)-conjg(transpose(position_tuple(:,:,axis))))
      do point=1,npoint
        q=6*(point-1)+2*axis
        if(present(point_representations))then
          hmats(:,:,q)=matmul(conjg(transpose(point_representations(:,:,point))),&
            matmul(gram,point_representations(:,:,point)))
        else
          hmats(:,:,q)=gram
        endif
      enddo
    enddo
    if(present(point_representations))then
      do point=1,npoint
        gram=matmul(conjg(transpose(point_representations(:,:,point))),point_representations(:,:,point))
        do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
        if(maxval(abs(gram))>10d0*tolerance)bad=1
      enddo
      gram=point_representations(:,:,1);do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
      if(maxval(abs(gram))>10d0*tolerance)bad=1
      call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.gbad/=0)then
        call cleanup();message='joint periodic-center point representations are not unitary with identity first';return
      endif
    endif
    unitary=(0d0,0d0);do i=1,m;unitary(i,i)=1d0;enddo
    objective_scale=sum(abs(hmats)**2)
    previous_objective=0d0
    do q=1,hmatrix_count;do j=1,m;do i=1,m
      if(i/=j)previous_objective=previous_objective+abs(hmats(i,j,q))**2
    enddo;enddo;enddo
    objective_threshold=max(100d0*epsilon(1d0),tolerance*tolerance)*max(1d0,objective_scale)
    maximum_update=0d0
    maximum_sweeps=max(100,m)
    do sweep=1,maximum_sweeps
      local_update=0d0
      sweep_bad=0
      do pair_j=2,m;do pair_i=1,pair_j-1
        gmat=0d0
        do q=1,hmatrix_count
          gvec=[real(hmats(pair_i,pair_i,q)-hmats(pair_j,pair_j,q),real64),&
            2d0*real(hmats(pair_i,pair_j,q),real64),-2d0*aimag(hmats(pair_i,pair_j,q))]
          do j=1,3;do i=1,3;gmat(i,j)=gmat(i,j)+gvec(i)*gvec(j);enddo;enddo
        enddo
        if(maxval(abs(gmat))<=100d0*tolerance*tolerance)cycle
        call dsyev('V','U',3,gmat,3,geval,gwork,9,info)
        if(info/=0.or..not.all(ieee_is_finite(geval)))then
          sweep_bad=1;cycle
        endif
        gvec=gmat(:,3);if(gvec(1)<0d0)gvec=-gvec
        c=sqrt(max(0d0,0.5d0*(1d0+min(1d0,gvec(1)))))
        if(c<=epsilon(1d0))cycle
        sphase=cmplx(gvec(2),gvec(3),real64)/(2d0*c)
        sabs=abs(sphase);if(sabs<=10d0*tolerance)cycle
        jacobi=reshape([cmplx(c,0d0,real64),sphase,-conjg(sphase),cmplx(c,0d0,real64)],[2,2])
        local_update=max(local_update,sabs)
        do q=1,hmatrix_count
          do k=1,m
            left_pair=[hmats(k,pair_i,q),hmats(k,pair_j,q)]
            hmats(k,pair_i,q)=left_pair(1)*jacobi(1,1)+left_pair(2)*jacobi(2,1)
            hmats(k,pair_j,q)=left_pair(1)*jacobi(1,2)+left_pair(2)*jacobi(2,2)
          enddo
          do k=1,m
            right_pair=[hmats(pair_i,k,q),hmats(pair_j,k,q)]
            hmats(pair_i,k,q)=conjg(jacobi(1,1))*right_pair(1)+conjg(jacobi(2,1))*right_pair(2)
            hmats(pair_j,k,q)=conjg(jacobi(1,2))*right_pair(1)+conjg(jacobi(2,2))*right_pair(2)
          enddo
        enddo
        do k=1,m
          left_pair=[unitary(k,pair_i),unitary(k,pair_j)]
          unitary(k,pair_i)=left_pair(1)*jacobi(1,1)+left_pair(2)*jacobi(2,1)
          unitary(k,pair_j)=left_pair(1)*jacobi(1,2)+left_pair(2)*jacobi(2,2)
        enddo
      enddo;enddo
      call MPI_Allreduce(sweep_bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.gbad/=0)then
        call cleanup();message='joint periodic-center local diagonalization failed';return
      endif
      call MPI_Allreduce(local_update,global_update,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      maximum_update=global_update;sweep_count=sweep
      local_objective=0d0
      do q=1,hmatrix_count;do j=1,m;do i=1,m
        if(i/=j)local_objective=local_objective+abs(hmats(i,j,q))**2
      enddo;enddo;enddo
      objective_change=previous_objective-local_objective
      bad=merge(0,1,ieee_is_finite(local_objective).and.&
        objective_change>=-objective_threshold)
      call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.gbad/=0)then
        call cleanup();message='joint periodic-center objective became unstable';return
      endif
      if(global_update<=10d0*tolerance.or.abs(objective_change)<=objective_threshold)exit
      previous_objective=local_objective
    enddo
    local_objective=0d0
    do q=1,hmatrix_count;do j=1,m;do i=1,m
      if(i/=j)local_objective=local_objective+abs(hmats(i,j,q))**2
    enddo;enddo;enddo
    call MPI_Allreduce(local_objective,global_objective,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    final_objective=global_objective
    if(ierr/=MPI_SUCCESS.or.sweep_count>=maximum_sweeps.and.maximum_update>10d0*tolerance.and.&
        abs(objective_change)>objective_threshold)then
      write(message,'(a,i0,4(a,es16.8))')'joint periodic-center sweeps did not converge count=',sweep_count,&
        ' objective=',final_objective,' objective_change=',objective_change,&
        ' maximum_update=',maximum_update,' threshold=',objective_threshold
      call cleanup();return
    endif
    point_orbit_built=.false.
    if(present(point_representations))then
      call build_point_orbit_blocks(point_orbit_built)
      if(.not.point_orbit_built)then
        call cleanup()
        if(len_trim(message)==0)message='joint periodic-center could not construct complete point-orbit blocks'
        return
      endif
    else
      do j=1,m;do axis=1,3
        phase_angle=atan2(real(hmats(j,j,2*axis),real64),real(hmats(j,j,2*axis-1),real64))
        centers(axis,j)=modulo(phase_angle/(2d0*acos(-1d0)),1d0)
      enddo;enddo
    endif
    order=[(i,i=1,m)]
    do i=2,m;k=i
      do while(k>1)
        swapped=.false.
        do axis=1,3
          if(centers(axis,order(k))<centers(axis,order(k-1))-10d0*tolerance)then;swapped=.true.;exit;endif
          if(centers(axis,order(k))>centers(axis,order(k-1))+10d0*tolerance)exit
        enddo
        if(.not.swapped)exit
        j=order(k);order(k)=order(k-1);order(k-1)=j;k=k-1
      enddo
    enddo
    gram=unitary;do j=1,m;unitary(:,j)=gram(:,order(j));enddo
    centers=centers(:,order)
    l=1
    do while(l<=m)
      r=l
      do while(r<m)
        if(any(min(abs(centers(:,r+1)-centers(:,l)),&
          1d0-abs(centers(:,r+1)-centers(:,l)))>10d0*tolerance))exit
        r=r+1
      enddo
      block_size=r-l+1
      if(block_size>1)then
        block=(0d0,0d0)
        block(1:block_size,1:block_size)=matmul(conjg(transpose(unitary(:,l:r))),&
          matmul(lcfo_operator,unitary(:,l:r)))
        call zheev('V','U',block_size,block,m,block_eval,zwork,max(1,2*m),zrwork,info)
        bad=merge(0,1,info==0.and.all(ieee_is_finite(block_eval(1:block_size))))
        call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
        if(ierr/=MPI_SUCCESS.or.gbad/=0)then
          call cleanup();message='joint periodic-center LCFO block diagonalization failed';return
        endif
        gram(:,1:block_size)=matmul(unitary(:,l:r),block(1:block_size,1:block_size))
        unitary(:,l:r)=gram(:,1:block_size)
      endif
      l=r+1
    enddo
    aligned_rows=matmul(sector_rows,unitary)
    do j=1,m
      pivot=maxval(abs(aligned_rows(:,j)));call MPI_Allreduce(MPI_IN_PLACE,pivot,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      do k=1,global_count
        tmp=(0d0,0d0);if(rank==owner(k)-1)tmp=aligned_rows(position(k),j)
        call MPI_Bcast(tmp,1,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr);if(ierr/=MPI_SUCCESS)then;call cleanup();return;endif
        if(abs(tmp)>=pivot-10d0*tolerance)then
          if(abs(tmp)>tolerance)then
            phase_fix=conjg(tmp)/abs(tmp);aligned_rows(:,j)=aligned_rows(:,j)*phase_fix;unitary(:,j)=unitary(:,j)*phase_fix
          endif
          exit
        endif
      enddo
    enddo
    rotation=unitary
    gram=matmul(conjg(transpose(aligned_rows)),aligned_rows)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    gram_defect=maxval(abs(gram));point_leakage=0d0;worst_point=0;worst_column=0
    if(present(point_representations))then
      do point=1,npoint
        block=matmul(conjg(transpose(unitary)),matmul(point_representations(:,:,point),unitary))
        do j=1,m
          local_update=0d0;l=1
          do while(l<=m)
            r=l
            do while(r<m)
              if(any(min(abs(centers(:,r+1)-centers(:,l)),1d0-abs(centers(:,r+1)-centers(:,l)))>&
                  tolerance**0.25d0))exit
              r=r+1
            enddo
            local_update=max(local_update,sum(abs(block(l:r,j))**2));l=r+1
          enddo
          if(max(0d0,1d0-local_update)>point_leakage)then
            point_leakage=max(0d0,1d0-local_update);worst_point=point;worst_column=j
          endif
        enddo
      enddo
    endif
    call MPI_Allreduce(MPI_IN_PLACE,gram_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,point_leakage,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    canonical_defect=max(gram_defect,point_leakage)
    if(ierr/=MPI_SUCCESS.or.gram_defect>10d0*tolerance)then
      call cleanup();message='joint periodic-center frame is not orthonormal';return
    endif
    if(point_leakage>tolerance**0.25d0)then
      write(message,'(a,es12.4,a,i0,a,i0,a,es12.4,a,es12.4,a,i0)')&
        'joint periodic-center point action leakage=',point_leakage,' point=',worst_point,&
        ' column=',worst_column,' objective=',final_objective,' update=',maximum_update,' sweeps=',sweep_count
      call cleanup();return
    endif
    quantum=100d0*tolerance;fingerprint=int(z'510E527FADE682D1',int64)
    fingerprint=ieor(ishftc(fingerprint,9),int(global_count,int64))
    fingerprint=ieor(ishftc(fingerprint,9),int(m,int64))
    do q=1,2
      stream=(0d0,0d0)
      do i=1,nlocal
        phase_angle=2d0*acos(-1d0)*modulo(real(row_ids(i),real64)*&
          merge(0.6180339887498948482d0,0.4142135623730950488d0,q==1),1d0)
        probe=cmplx(cos(phase_angle),sin(phase_angle),real64)/sqrt(real(global_count,real64))
        stream=stream+conjg(aligned_rows(i,:))*probe
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,stream,m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();return;endif
      do k=1,global_count
        tmp=(0d0,0d0)
        if(rank==owner(k)-1)tmp=sum(aligned_rows(position(k),:)*stream)
        call MPI_Bcast(tmp,1,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;call cleanup();return;endif
        quantized=nint(real(tmp,real64)/quantum,int64);fingerprint=ieor(ishftc(fingerprint,9),quantized)
        quantized=nint(aimag(tmp)/quantum,int64);fingerprint=ieor(ishftc(fingerprint,9),quantized)
      enddo
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    call MPI_Allreduce(workspace_peak_bytes,peak,1,MPI_INTEGER8,MPI_MAX,comm,ierr);workspace_peak_bytes=peak
    ok=ierr==MPI_SUCCESS
    if(ok)message=''
    call cleanup(.false.)
  contains
    subroutine build_point_orbit_blocks(built)
      logical,intent(out)::built
      integer::seed,pidx,cidx,ncluster,ncandidate,ncovered,first_column,last_column,destination,&
        source_start,source_end,norm_info
      real(real64)::center_tolerance,rank_tolerance,residual_norm,cover_norm,&
        block_weight,best_weight,leakage
      complex(real64)::expectation
      logical::insert_before
      built=.false.;center_tolerance=tolerance**0.25d0;rank_tolerance=sqrt(tolerance)
      ncluster=0;ncandidate=0;ncovered=0;norm_info=0;orbit_labels=0;basis_labels=0
      cluster_centers=0d0;orbit_basis=(0d0,0d0);cover_basis=(0d0,0d0);centers=0d0
      do seed=1,m
        orbit_residual=unitary(:,seed)
        if(ncovered>0)orbit_residual=orbit_residual-matmul(cover_basis(:,1:ncovered),&
          matmul(conjg(transpose(cover_basis(:,1:ncovered))),orbit_residual))
        residual_norm=sqrt(max(0d0,real(dot_product(orbit_residual,orbit_residual),real64)))
        if(residual_norm<=rank_tolerance)cycle
        orbit_residual=orbit_residual/residual_norm
        do pidx=1,npoint
          orbit_vectors(:,pidx)=matmul(point_representations(:,:,pidx),orbit_residual)
        enddo
        do pidx=1,npoint
          do axis=1,3
            expectation=dot_product(orbit_vectors(:,pidx),&
              matmul(position_tuple(:,:,axis),orbit_vectors(:,pidx)))
            if(abs(expectation)>rank_tolerance)then
              orbit_centers(axis,pidx)=modulo(atan2(aimag(expectation),real(expectation,real64))/&
                (2d0*acos(-1d0)),1d0)
            else
              orbit_centers(axis,pidx)=0d0
            endif
          enddo
          orbit_labels(pidx)=0
          do cidx=1,ncluster
            if(all(min(abs(orbit_centers(:,pidx)-cluster_centers(:,cidx)),&
                1d0-abs(orbit_centers(:,pidx)-cluster_centers(:,cidx)))<=center_tolerance))then
              orbit_labels(pidx)=cidx;exit
            endif
          enddo
          if(orbit_labels(pidx)==0)then
            ncluster=ncluster+1;orbit_labels(pidx)=ncluster
            cluster_centers(:,ncluster)=orbit_centers(:,pidx)
          endif
          orbit_residual=orbit_vectors(:,pidx)
          do cidx=1,ncandidate
            if(basis_labels(cidx)==orbit_labels(pidx))orbit_residual=orbit_residual-&
              orbit_basis(:,cidx)*dot_product(orbit_basis(:,cidx),orbit_residual)
          enddo
          residual_norm=sqrt(max(0d0,real(dot_product(orbit_residual,orbit_residual),real64)))
          if(residual_norm<=rank_tolerance)cycle
          stream=orbit_residual/residual_norm
          orbit_residual=stream
          if(ncovered>0)orbit_residual=orbit_residual-matmul(cover_basis(:,1:ncovered),&
            matmul(conjg(transpose(cover_basis(:,1:ncovered))),orbit_residual))
          cover_norm=sqrt(max(0d0,real(dot_product(orbit_residual,orbit_residual),real64)))
          if(cover_norm<=rank_tolerance)cycle
          if(ncandidate>=m.or.ncovered>=m)then;norm_info=1;exit;endif
          ncandidate=ncandidate+1;orbit_basis(:,ncandidate)=stream;basis_labels(ncandidate)=orbit_labels(pidx)
          ncovered=ncovered+1;cover_basis(:,ncovered)=orbit_residual/cover_norm
        enddo
        if(norm_info/=0.or.ncovered>=m)exit
      enddo
      if(norm_info/=0.or.ncandidate/=m.or.ncovered/=m)then
        write(message,'(a,i0,a,i0,a,i0)')'point-orbit cover is incomplete: candidates=',ncandidate,&
          ' covered=',ncovered,' required=',m
        return
      endif
      cluster_sequence(1:ncluster)=[(cidx,cidx=1,ncluster)]
      do cidx=2,ncluster
        destination=cidx
        do while(destination>1)
          insert_before=.false.
          do axis=1,3
            if(cluster_centers(axis,cluster_sequence(destination))<&
                cluster_centers(axis,cluster_sequence(destination-1))-center_tolerance)then
              insert_before=.true.;exit
            endif
            if(cluster_centers(axis,cluster_sequence(destination))>&
                cluster_centers(axis,cluster_sequence(destination-1))+center_tolerance)exit
          enddo
          if(.not.insert_before)exit
          pidx=cluster_sequence(destination);cluster_sequence(destination)=cluster_sequence(destination-1)
          cluster_sequence(destination-1)=pidx;destination=destination-1
        enddo
      enddo
      gram=orbit_basis;destination=0
      do cidx=1,ncluster
        do pidx=1,m
          if(basis_labels(pidx)/=cluster_sequence(cidx))cycle
          destination=destination+1;orbit_basis(:,destination)=gram(:,pidx)
          centers(:,destination)=cluster_centers(:,cluster_sequence(cidx))
        enddo
      enddo
      gram=matmul(conjg(transpose(orbit_basis)),orbit_basis);block=gram
      call zheev('V','U',m,block,m,block_eval,zwork,max(1,2*m),zrwork,norm_info)
      bad=merge(0,1,norm_info==0.and.all(ieee_is_finite(block_eval)).and.&
        minval(block_eval)>rank_tolerance*rank_tolerance)
      call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.gbad/=0)then
        if(ierr==MPI_SUCCESS)write(message,'(a,es14.6,a,es14.6)')&
          'point-orbit Gram is rank deficient: minimum_eigenvalue=',minval(block_eval),&
          ' threshold=',rank_tolerance*rank_tolerance
        return
      endif
      inverse_sqrt=1d0/sqrt(block_eval)
      gram=matmul(block*spread(inverse_sqrt,1,m),conjg(transpose(block)))
      orbit_basis=matmul(orbit_basis,gram);leakage=0d0
      do point=1,npoint
        block=matmul(conjg(transpose(orbit_basis)),matmul(point_representations(:,:,point),orbit_basis))
        source_start=1
        do while(source_start<=m)
          source_end=source_start
          do while(source_end<m)
            if(any(min(abs(centers(:,source_end+1)-centers(:,source_start)),&
                1d0-abs(centers(:,source_end+1)-centers(:,source_start)))>center_tolerance))exit
            source_end=source_end+1
          enddo
          best_weight=0d0;first_column=1
          do while(first_column<=m)
            last_column=first_column
            do while(last_column<m)
              if(any(min(abs(centers(:,last_column+1)-centers(:,first_column)),&
                  1d0-abs(centers(:,last_column+1)-centers(:,first_column)))>center_tolerance))exit
              last_column=last_column+1
            enddo
            block_weight=sum(abs(block(first_column:last_column,source_start:source_end))**2)
            best_weight=max(best_weight,block_weight);first_column=last_column+1
          enddo
          leakage=max(leakage,max(0d0,1d0-best_weight/real(source_end-source_start+1,real64)))
          source_start=source_end+1
        enddo
      enddo
      if(leakage>center_tolerance)then
        write(message,'(a,es14.6,a,es14.6)')'point-orbit cluster leakage is excessive: leakage=',leakage,&
          ' threshold=',center_tolerance
        return
      endif
      unitary=orbit_basis;built=.true.
    end subroutine build_point_orbit_blocks

    subroutine cleanup(drop_outputs)
      logical,intent(in),optional::drop_outputs
      logical::drop
      drop=.true.;if(present(drop_outputs))drop=drop_outputs
      if(drop.and.allocated(aligned_rows))deallocate(aligned_rows)
      if(drop.and.allocated(centers))deallocate(centers)
      if(allocated(hmats))deallocate(hmats)
      if(allocated(unitary))deallocate(unitary)
      if(allocated(gram))deallocate(gram)
      if(allocated(stream))deallocate(stream)
      if(allocated(block))deallocate(block)
      if(allocated(zwork))deallocate(zwork)
      if(allocated(block_eval))deallocate(block_eval)
      if(allocated(zrwork))deallocate(zrwork)
      if(allocated(owner))deallocate(owner)
      if(allocated(position))deallocate(position)
      if(allocated(count))deallocate(count)
      if(allocated(order))deallocate(order)
      if(allocated(payload_bits))deallocate(payload_bits)
      if(allocated(payload_minimum))deallocate(payload_minimum)
      if(allocated(payload_maximum))deallocate(payload_maximum)
      if(allocated(orbit_vectors))deallocate(orbit_vectors)
      if(allocated(orbit_basis))deallocate(orbit_basis)
      if(allocated(orbit_residual))deallocate(orbit_residual)
      if(allocated(cover_basis))deallocate(cover_basis)
      if(allocated(orbit_centers))deallocate(orbit_centers)
      if(allocated(cluster_centers))deallocate(cluster_centers)
      if(allocated(inverse_sqrt))deallocate(inverse_sqrt)
      if(allocated(orbit_labels))deallocate(orbit_labels)
      if(allocated(cluster_sequence))deallocate(cluster_sequence)
      if(allocated(basis_labels))deallocate(basis_labels)
    end subroutine
#else
    ok=.false.;message='joint periodic-center gauge requires MPI';final_objective=huge(1d0)
    maximum_update=huge(1d0);canonical_defect=huge(1d0);fingerprint=0_int64;workspace_peak_bytes=0_int64
    sweep_count=0;rotation=(0d0,0d0)
#endif
  end subroutine jointly_canonicalize_dg_sector_periodic_position_gauge

  subroutine project_dg_w90_reference_sector_operators(comm,row_ids,sector_rows,w90_rows,w90_values,&
      lcfo_rows,lcfo_values,global_row_count,w90_fingerprint,lcfo_fingerprint,w90_frame_defect,lcfo_source_defect,tolerance,&
      w90_operator,lcfo_operator,projection_defect,fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,global_row_count
    integer(int64),intent(in)::row_ids(:),w90_fingerprint,lcfo_fingerprint
    complex(real64),intent(in)::sector_rows(:,:),w90_rows(:,:),lcfo_rows(:,:)
    real(real64),intent(in)::w90_values(:),lcfo_values(:),w90_frame_defect,lcfo_source_defect,tolerance
    complex(real64),allocatable,intent(out)::w90_operator(:,:),lcfo_operator(:,:)
    real(real64),intent(out)::projection_defect
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::wp(:,:),lp(:,:),gram(:,:),stream_w(:),stream_l(:)
    integer,allocatable::owner(:),position(:),count(:)
    integer::nlocal,m,nw,i,j,rank,ierr,bad,gbad,status,minint,maxint
    integer(int64)::bits,minhash,maxhash,elements,term,bytes,frame_hash
    real(real64)::mintol,maxtol,scale,local_scale,global_scale,sector_scale,w90_scale,lcfo_scale,safe_limit,bound
    logical::receipt_valid
    nlocal=size(row_ids);m=size(sector_rows,2);nw=size(w90_rows,2)
    ok=.false.;message='';projection_defect=huge(1d0);fingerprint=0_int64;workspace_peak_bytes=0_int64
    bad=merge(0,1,global_row_count>=1.and.m>=1.and.nw>=m.and.size(sector_rows,1)==nlocal.and.&
      all(shape(w90_rows)==[nlocal,nw]).and.all(shape(lcfo_rows)==[nlocal,nw]).and.&
      size(w90_values)==nw.and.size(lcfo_values)==nw.and.w90_fingerprint/=0_int64.and.lcfo_fingerprint/=0_int64.and.&
      w90_frame_defect>=0d0.and.w90_frame_defect<=tolerance.and.&
      lcfo_source_defect>=0d0.and.lcfo_source_defect<=tolerance.and.&
      all(row_ids>=1_int64).and.all(row_ids<=int(global_row_count,int64)).and.&
      tolerance>=1d-15.and.tolerance<=1d-2.and.ieee_is_finite(tolerance).and.&
      all(ieee_is_finite(w90_values)).and.all(ieee_is_finite(lcfo_values)).and.&
      all(ieee_is_finite(real(sector_rows))).and.all(ieee_is_finite(aimag(sector_rows))).and.&
      all(ieee_is_finite(real(w90_rows))).and.all(ieee_is_finite(aimag(w90_rows))).and.&
      all(ieee_is_finite(real(lcfo_rows))).and.all(ieee_is_finite(aimag(lcfo_rows))))
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='invalid reference-sector operator projection contract';return;endif
    do i=1,3
      if(i==1)j=global_row_count;if(i==2)j=m;if(i==3)j=nw
      call MPI_Allreduce(j,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(j,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='reference operator dimensions disagree';return;endif
    enddo
    call MPI_Allreduce(tolerance,mintol,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tolerance,maxtol,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.mintol/=maxtol)then;message='reference operator tolerance disagrees';return;endif
    do i=1,2
      scale=merge(w90_frame_defect,lcfo_source_defect,i==1)
      call MPI_Allreduce(scale,mintol,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(scale,maxtol,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.mintol/=maxtol)then;message='reference frame defects disagree';return;endif
      bits=merge(w90_fingerprint,lcfo_fingerprint,i==1)
      call MPI_Allreduce(bits,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(bits,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)then;message='reference frame fingerprints disagree';return;endif
    enddo
    do i=1,nw
      bits=transfer(w90_values(i),bits);call MPI_Allreduce(bits,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)return;call MPI_Allreduce(bits,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)then;message='W90 physical operator values disagree';return;endif
      bits=transfer(lcfo_values(i),bits);call MPI_Allreduce(bits,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)return;call MPI_Allreduce(bits,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)then;message='LCFO physical operator values disagree';return;endif
    enddo
    sector_scale=0d0;w90_scale=0d0;lcfo_scale=0d0
    if(nlocal>0)then
      sector_scale=maxval(abs(sector_rows));w90_scale=maxval(abs(w90_rows));lcfo_scale=maxval(abs(lcfo_rows))
    endif
    call MPI_Allreduce(MPI_IN_PLACE,sector_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,w90_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,lcfo_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS)return
    scale=max(maxval(abs(w90_values)),maxval(abs(lcfo_values)))
    safe_limit=huge(1d0)/64d0;bad=0
    do i=1,2
      local_scale=merge(w90_scale,lcfo_scale,i==1);global_scale=merge(maxval(abs(w90_values)),maxval(abs(lcfo_values)),i==1)
      if(global_scale>0d0)then
        bound=sqrt((safe_limit/real(nw,real64))/global_scale)/real(global_row_count,real64)
      else
        bound=huge(1d0)
      endif
      if(local_scale>0d0)then
        if(sector_scale>bound/local_scale)bad=1
      endif
    enddo
    if(m>huge(0)/m.or.m>huge(0)/nw)bad=1
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='reference operator projection magnitude or extent overflows';return;endif
    elements=0_int64;bytes=0_int64;receipt_valid=.true.
    call checked_product([2_int64,int(m,int64),int(nw,int64)],term,receipt_valid);call checked_add(elements,term,receipt_valid)
    call checked_product([3_int64,int(m,int64),int(m,int64)],term,receipt_valid);call checked_add(elements,term,receipt_valid)
    call checked_product([2_int64,int(nw,int64)],term,receipt_valid);call checked_add(elements,term,receipt_valid)
    call checked_product([elements,16_int64],bytes,receipt_valid)
    call checked_product([3_int64,int(global_row_count,int64),4_int64],term,receipt_valid);call checked_add(bytes,term,receipt_valid)
    bad=merge(0,1,receipt_valid);call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='reference operator workspace overflows';return;endif
    workspace_peak_bytes=bytes;call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(wp(m,nw),lp(m,nw),gram(m,m),stream_w(nw),stream_l(nw),w90_operator(m,m),lcfo_operator(m,m),&
      owner(global_row_count),position(global_row_count),count(global_row_count),stat=status)
    call MPI_Allreduce(status,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then
      if(allocated(w90_operator))deallocate(w90_operator);if(allocated(lcfo_operator))deallocate(lcfo_operator)
      message='reference operator projection allocation failed';return
    endif
    owner=0;position=0;count=0
    do i=1,nlocal
      owner(int(row_ids(i)))=rank+1;position(int(row_ids(i)))=i;count(int(row_ids(i)))=count(int(row_ids(i)))+1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,owner,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,position,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,count,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(count/=1))then;message='reference operator rows are not uniquely owned';return;endif
    frame_hash=int(z'510E527FADE682D1',int64)
    do i=1,global_row_count
      stream_w=(0d0,0d0);stream_l=(0d0,0d0)
      if(rank==owner(i)-1)then;stream_w=w90_rows(position(i),:);stream_l=lcfo_rows(position(i),:);endif
      call MPI_Bcast(stream_w,nw,MPI_DOUBLE_COMPLEX,owner(i)-1,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Bcast(stream_l,nw,MPI_DOUBLE_COMPLEX,owner(i)-1,comm,ierr);if(ierr/=MPI_SUCCESS)return
      frame_hash=ieor(ishftc(frame_hash,7),int(i,int64))
      do j=1,nw
        bits=transfer(real(stream_w(j),real64),bits);frame_hash=ieor(ishftc(frame_hash,7),bits)
        bits=transfer(aimag(stream_w(j)),bits);frame_hash=ieor(ishftc(frame_hash,7),bits)
        bits=transfer(real(stream_l(j),real64),bits);frame_hash=ieor(ishftc(frame_hash,7),bits)
        bits=transfer(aimag(stream_l(j)),bits);frame_hash=ieor(ishftc(frame_hash,7),bits)
      enddo
    enddo
    if(frame_hash==0_int64)frame_hash=1_int64
    gram=matmul(conjg(transpose(sector_rows)),sector_rows)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo;projection_defect=maxval(abs(gram))
    wp=matmul(conjg(transpose(sector_rows)),w90_rows);lp=matmul(conjg(transpose(sector_rows)),lcfo_rows)
    call MPI_Allreduce(MPI_IN_PLACE,wp,m*nw,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,lp,m*nw,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    w90_operator=matmul(wp*spread(w90_values,1,m),conjg(transpose(wp)))
    lcfo_operator=matmul(lp*spread(lcfo_values,1,m),conjg(transpose(lp)))
    bad=merge(0,1,all(ieee_is_finite(real(w90_operator))).and.all(ieee_is_finite(aimag(w90_operator))).and.&
      all(ieee_is_finite(real(lcfo_operator))).and.all(ieee_is_finite(aimag(lcfo_operator))))
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='projected reference operators are nonfinite';return;endif
    projection_defect=max(projection_defect,maxval(abs(w90_operator-conjg(transpose(w90_operator)))),&
      maxval(abs(lcfo_operator-conjg(transpose(lcfo_operator)))))
    call MPI_Allreduce(MPI_IN_PLACE,projection_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.projection_defect>10d0*tolerance)then;message='reference operator projection defect too large';return;endif
    fingerprint=ieor(ieor(w90_fingerprint,lcfo_fingerprint),frame_hash)
    do i=1,m;do j=1,m
      if(max(abs(real(w90_operator(i,j),real64)),abs(aimag(w90_operator(i,j))),&
          abs(real(lcfo_operator(i,j),real64)),abs(aimag(lcfo_operator(i,j))))/(100d0*tolerance)>&
          0.25d0*real(huge(0_int64),real64))then;message='reference operator fingerprint overflows';return;endif
      bits=nint(real(w90_operator(i,j),real64)/(100d0*tolerance),int64);fingerprint=ieor(ishftc(fingerprint,7),bits)
      bits=nint(aimag(w90_operator(i,j))/(100d0*tolerance),int64);fingerprint=ieor(ishftc(fingerprint,7),bits)
      bits=nint(real(lcfo_operator(i,j),real64)/(100d0*tolerance),int64);fingerprint=ieor(ishftc(fingerprint,7),bits)
      bits=nint(aimag(lcfo_operator(i,j))/(100d0*tolerance),int64);fingerprint=ieor(ishftc(fingerprint,7),bits)
    enddo;enddo
    if(fingerprint==0_int64)fingerprint=1_int64;ok=.true.
  end subroutine project_dg_w90_reference_sector_operators

  subroutine anchor_dg_w90_reference_character_sector(comm,row_ids,sector_rows,w90_operator,lcfo_operator,&
      global_row_count,w90_fingerprint,lcfo_fingerprint,w90_frame_defect,lcfo_source_defect,tolerance,&
      anchored_rows,anchor_defect,fingerprint,workspace_peak_bytes,ok,message,projector_diagonal_only)
    integer,intent(in)::comm,global_row_count
    integer(int64),intent(in)::row_ids(:),w90_fingerprint,lcfo_fingerprint
    complex(real64),intent(in)::sector_rows(:,:)
    complex(real64),intent(in)::w90_operator(:,:),lcfo_operator(:,:)
    real(real64),intent(in)::tolerance
    complex(real64),allocatable,intent(out)::anchored_rows(:,:)
    real(real64),intent(out)::anchor_defect
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,intent(in),optional::projector_diagonal_only
    complex(real64),allocatable::hmat(:,:),kmat(:,:),unitary(:,:),block(:,:),work(:),stream(:),remote_stream(:),gram(:,:)
    real(real64),allocatable::eval(:),block_eval(:),rwork(:)
    integer,allocatable::owner(:),position(:),count(:),same_cluster(:)
    integer::nlocal,m,i,j,k,l,r,rank,ierr,bad,gbad,status,lwork,info,minint,maxint,diagonal_flag
    integer(int64)::bits,minhash,maxhash,complex_elements,real_elements,integer_elements,byte_term,operator_hash
    real(real64),intent(in)::w90_frame_defect,lcfo_source_defect
    real(real64)::mintol,maxtol,pivot,scale,operator_scale,local_defect,global_defect
    complex(real64)::projector_value
    logical::receipt_valid
    interface
      subroutine zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
        character,intent(in)::jobz,uplo;integer,intent(in)::n,lda,lwork
        complex(8),intent(inout)::a(lda,*),work(*);real(8),intent(out)::w(*),rwork(*)
        integer,intent(out)::info
      end subroutine
    end interface
    nlocal=size(row_ids);m=size(sector_rows,2)
    ok=.false.;message='';anchor_defect=huge(1d0);fingerprint=0_int64;workspace_peak_bytes=0_int64
    bad=merge(0,1,global_row_count>=1.and.m>=1.and.size(sector_rows,1)==nlocal.and.&
      all(shape(w90_operator)==[m,m]).and.all(shape(lcfo_operator)==[m,m]).and.&
      w90_fingerprint/=0_int64.and.lcfo_fingerprint/=0_int64.and.&
      w90_frame_defect>=0d0.and.w90_frame_defect<=tolerance.and.&
      lcfo_source_defect>=0d0.and.lcfo_source_defect<=tolerance.and.&
      all(row_ids>=1_int64).and.all(row_ids<=int(global_row_count,int64)).and.&
      tolerance>=1d-15.and.tolerance<=1d-2.and.ieee_is_finite(tolerance).and.&
      all(ieee_is_finite(real(w90_operator))).and.all(ieee_is_finite(aimag(w90_operator))).and.&
      all(ieee_is_finite(real(lcfo_operator))).and.all(ieee_is_finite(aimag(lcfo_operator))).and.&
      all(ieee_is_finite(real(sector_rows))).and.all(ieee_is_finite(aimag(sector_rows))).and.&
      maxval(abs(w90_operator-conjg(transpose(w90_operator))))<=10d0*tolerance.and.&
      maxval(abs(lcfo_operator-conjg(transpose(lcfo_operator))))<=10d0*tolerance)
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='invalid W90 reference-sector anchor contract';return;endif
    do i=1,2
      select case(i);case(1);j=global_row_count;case default;j=m;endselect
      call MPI_Allreduce(j,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(j,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='W90 anchor metadata disagree';return;endif
    enddo
    call MPI_Allreduce(tolerance,mintol,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tolerance,maxtol,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.mintol/=maxtol)then;message='W90 anchor tolerance disagrees';return;endif
    diagonal_flag=0;if(present(projector_diagonal_only))diagonal_flag=merge(1,0,projector_diagonal_only)
    call MPI_Allreduce(diagonal_flag,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(diagonal_flag,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='W90 anchor fingerprint mode disagrees';return;endif
    do k=1,2
      bits=merge(w90_fingerprint,lcfo_fingerprint,k==1)
      call MPI_Allreduce(bits,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(bits,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)then;message='W90 anchor provenance disagrees';return;endif
    enddo
    do k=1,2
      pivot=merge(w90_frame_defect,lcfo_source_defect,k==1)
      call MPI_Allreduce(pivot,mintol,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(pivot,maxtol,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.mintol/=maxtol)then;message='W90 anchor frame receipts disagree';return;endif
    enddo
    do j=1,m;do i=1,m
      bits=transfer(real(w90_operator(i,j),real64),bits)
      call MPI_Allreduce(bits,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(bits,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)return
      bits=transfer(aimag(w90_operator(i,j)),bits)
      call MPI_Allreduce(bits,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(bits,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)return
      bits=transfer(real(lcfo_operator(i,j),real64),bits)
      call MPI_Allreduce(bits,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(bits,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr);if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)return
      bits=transfer(aimag(lcfo_operator(i,j)),bits)
      call MPI_Allreduce(bits,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(bits,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)then;message='W90 anchor operator payloads disagree';return;endif
    enddo;enddo
    operator_scale=max(maxval(abs(w90_operator)),maxval(abs(lcfo_operator)))
    bad=merge(0,1,operator_scale<=sqrt(huge(1d0))/(16d0*real(max(1,m),real64)))
    if(m>huge(0)/m.or.m>huge(0)/5)bad=1
    complex_elements=0_int64;real_elements=0_int64;integer_elements=0_int64;receipt_valid=bad==0
    if(receipt_valid)then
      call checked_product([int(m,int64),int(m,int64),5_int64],byte_term,receipt_valid)
      call checked_add(complex_elements,byte_term,receipt_valid)
      call checked_product([int(nlocal,int64),int(m,int64)],byte_term,receipt_valid)
      call checked_add(complex_elements,byte_term,receipt_valid)
      call checked_product([4_int64,int(m,int64)],byte_term,receipt_valid)
      call checked_add(complex_elements,byte_term,receipt_valid)
      call checked_product([7_int64,int(m,int64)],byte_term,receipt_valid)
      call checked_add(real_elements,byte_term,receipt_valid)
      call checked_product([3_int64,int(global_row_count,int64)],byte_term,receipt_valid)
      call checked_add(integer_elements,byte_term,receipt_valid)
      call checked_product([complex_elements,16_int64],workspace_peak_bytes,receipt_valid)
      call checked_product([real_elements,8_int64],byte_term,receipt_valid)
      call checked_add(workspace_peak_bytes,byte_term,receipt_valid)
      call checked_product([integer_elements,4_int64],byte_term,receipt_valid)
      call checked_add(workspace_peak_bytes,byte_term,receipt_valid)
    endif
    bad=merge(0,1,receipt_valid)
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='W90 anchor extent, magnitude, or workspace overflows';return;endif
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(hmat(m,m),kmat(m,m),unitary(m,m),gram(m,m),eval(m),block_eval(m),block(m,m),&
      rwork(max(1,3*m)),work(max(1,2*m)),anchored_rows(nlocal,m),owner(global_row_count),&
      position(global_row_count),count(global_row_count),same_cluster(max(1,m-1)),stream(m),remote_stream(m),stat=status)
    call MPI_Allreduce(status,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then
      if(allocated(anchored_rows))deallocate(anchored_rows)
      message='W90 anchor allocation failed';return
    endif
    owner=0;position=0;count=0
    do i=1,nlocal
      owner(int(row_ids(i)))=rank+1;position(int(row_ids(i)))=i
      count(int(row_ids(i)))=count(int(row_ids(i)))+1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,owner,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,position,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,count,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(count/=1))then;message='W90 anchor rows are not uniquely owned';return;endif
    gram=matmul(conjg(transpose(sector_rows)),sector_rows)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    bad=merge(0,1,maxval(abs(gram))<=10d0*tolerance)
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='W90 anchor sector frame is not orthonormal';return;endif
    hmat=w90_operator;kmat=lcfo_operator
    unitary=hmat;lwork=size(work);call zheev('V','U',m,unitary,m,eval,work,lwork,rwork,info)
    bad=merge(0,1,info==0.and.all(ieee_is_finite(eval)))
    call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.gbad/=0)then;message='W90 anchor position diagonalization failed';return;endif
    local_defect=maxval(abs(matmul(hmat,unitary)-unitary*spread(eval,1,m)))
    same_cluster=0
    do i=1,m-1
      bad=merge(1,0,abs(eval(i+1)-eval(i))<=10d0*tolerance*max(1d0,maxval(abs(eval))))
      call MPI_Allreduce(bad,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(bad,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='W90 anchor H cluster boundaries disagree';return;endif
      same_cluster(i)=minint
    enddo
    l=1
    do while(l<=m)
      r=l
      do while(r<m)
        if(same_cluster(r)==0)exit
        r=r+1
      enddo
      if(r>l)then
        block(1:r-l+1,1:r-l+1)=matmul(conjg(transpose(unitary(:,l:r))),matmul(kmat,unitary(:,l:r)))
        call zheev('V','U',r-l+1,block,r-l+1,block_eval,work,lwork,rwork,info)
        bad=merge(0,1,info==0.and.all(ieee_is_finite(block_eval(1:r-l+1))))
        call MPI_Allreduce(bad,gbad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
        if(ierr/=MPI_SUCCESS.or.gbad/=0)then
          message='LCFO anchor block diagonalization failed';return
        endif
        scale=max(1d0,maxval(abs(block_eval(1:r-l+1))))
        ! Exact repeated eigenvalues are complete symmetry multiplets, not a
        ! numerical rank cut.  ZHEEV resolves every nondegenerate direction and
        ! leaves only an arbitrary unitary gauge inside each exact multiplet;
        ! downstream character alignment transports that whole block gauge.
        unitary(:,l:r)=matmul(unitary(:,l:r),block(1:r-l+1,1:r-l+1))
        hmat(:,1:r-l+1)=matmul(kmat,unitary(:,l:r))
        hmat(:,1:r-l+1)=matmul(unitary(:,l:r),matmul(conjg(transpose(unitary(:,l:r))),hmat(:,1:r-l+1)))
        local_defect=max(local_defect,maxval(abs(hmat(:,1:r-l+1)-&
          unitary(:,l:r)*spread(block_eval(1:r-l+1),1,m))))
      endif
      l=r+1
    enddo
    anchored_rows=matmul(sector_rows,unitary)
    gram=matmul(conjg(transpose(anchored_rows)),anchored_rows)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    local_defect=max(local_defect,maxval(abs(gram)))
    call MPI_Allreduce(local_defect,global_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    anchor_defect=global_defect
    if(anchor_defect>10d0*tolerance)then;message='anchored W90 reference sector is not orthonormal';return;endif
    do j=1,m
      pivot=maxval(abs(anchored_rows(:,j)));call MPI_Allreduce(MPI_IN_PLACE,pivot,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      do k=1,global_row_count
        stream=(0d0,0d0);if(rank==owner(k)-1)stream(1)=anchored_rows(position(k),j)
        call MPI_Bcast(stream,1,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr);if(ierr/=MPI_SUCCESS)return
        if(abs(stream(1))>=pivot-10d0*tolerance)then
          if(abs(stream(1))>tolerance)anchored_rows(:,j)=anchored_rows(:,j)*conjg(stream(1))/abs(stream(1))
          exit
        endif
      enddo
    enddo
    operator_hash=ieor(w90_fingerprint,lcfo_fingerprint)
    bits=transfer(w90_frame_defect,bits);operator_hash=ieor(ishftc(operator_hash,7),bits)
    bits=transfer(lcfo_source_defect,bits);operator_hash=ieor(ishftc(operator_hash,7),bits)
    do i=1,m
      if(abs(eval(i))/(100d0*tolerance)>0.25d0*real(huge(0_int64),real64))then
        message='W90 anchor eigenvalue quantization overflows';return
      endif
      bits=nint(eval(i)/(100d0*tolerance),int64);operator_hash=ieor(ishftc(operator_hash,7),bits)
    enddo
    fingerprint=operator_hash
    do k=1,global_row_count
      stream=(0d0,0d0);if(rank==owner(k)-1)stream=anchored_rows(position(k),:)
      call MPI_Bcast(stream,m,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr);if(ierr/=MPI_SUCCESS)return
      fingerprint=ieor(ishftc(fingerprint,11),int(k,int64))
      if(diagonal_flag==1)then
        projector_value=cmplx(sum(abs(stream)**2),0d0,real64)
        if(abs(real(projector_value,real64))/(100d0*tolerance)>&
            0.25d0*real(huge(0_int64),real64))then
          message='W90 anchor diagonal fingerprint quantization overflows';return
        endif
        bits=nint(real(projector_value,real64)/(100d0*tolerance),int64)
        fingerprint=ieor(ishftc(fingerprint,11),bits)
      else
        do i=1,global_row_count
          remote_stream=(0d0,0d0)
          if(rank==owner(i)-1)remote_stream=anchored_rows(position(i),:)
          call MPI_Bcast(remote_stream,m,MPI_DOUBLE_COMPLEX,owner(i)-1,comm,ierr);if(ierr/=MPI_SUCCESS)return
          projector_value=sum(stream*conjg(remote_stream))
          if(max(abs(real(projector_value,real64)),abs(aimag(projector_value)))/(100d0*tolerance)>&
              0.25d0*real(huge(0_int64),real64))then
            message='W90 anchor projector fingerprint quantization overflows';return
          endif
          bits=nint(real(projector_value,real64)/(100d0*tolerance),int64)
          fingerprint=ieor(ishftc(fingerprint,11),ieor(int(i,int64),bits))
          bits=nint(aimag(projector_value)/(100d0*tolerance),int64)
          fingerprint=ieor(ishftc(fingerprint,11),bits)
        enddo
      endif
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    ok=.true.
  end subroutine anchor_dg_w90_reference_character_sector
  subroutine validate_dg_w90_localization_cluster(eigenvalues,selected_count,tolerance,ok,message)
    ! The localization eigensolver must supply this spectrum in ascending order.
    real(real64),intent(in)::eigenvalues(:),tolerance
    integer,intent(in)::selected_count
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i
    real(real64)::scale
    ok=.false.;message=''
    if(size(eigenvalues)<1.or.selected_count<1.or.selected_count>size(eigenvalues).or.&
        tolerance<=0d0.or..not.ieee_is_finite(tolerance).or..not.all(ieee_is_finite(eigenvalues)))then
      message='invalid Wannier90 localization-cluster contract';return
    endif
    do i=2,size(eigenvalues)
      if(eigenvalues(i)<eigenvalues(i-1))then;message='Wannier90 localization spectrum is not ordered';return;endif
    enddo
    if(selected_count<size(eigenvalues))then
      scale=max(1d0,maxval(abs(eigenvalues)))
      if(abs(eigenvalues(selected_count+1)-eigenvalues(selected_count))<=tolerance*scale)then
        message='Wannier90 selection splits a degenerate localization cluster';return
      endif
    endif
    ok=.true.
  end subroutine validate_dg_w90_localization_cluster

#ifdef USE_MPI
  subroutine align_dg_w90_character_sectors_by_periodic_phase(comm,row_ids,reference_rows,target_rows,&
      periodic_phase,global_row_count,phase_fingerprint,phase_payload_fingerprint,tolerance,&
      aligned_rows,singular_values,polar_defect,fingerprint,&
      workspace_peak_bytes,ok,message,integration_weights)
    integer,intent(in)::comm,global_row_count
    integer(int64),intent(in)::phase_fingerprint,phase_payload_fingerprint
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::reference_rows(:,:),target_rows(:,:),periodic_phase(:)
    real(real64),intent(in)::tolerance
    real(real64),intent(in),optional::integration_weights(:)
    complex(real64),allocatable,intent(out)::aligned_rows(:,:)
    real(real64),allocatable,intent(out)::singular_values(:)
    real(real64),intent(out)::polar_defect
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::link(:,:),left(:,:),right(:,:),polar(:,:),gram(:,:),work(:),sketch(:)
    real(real64),allocatable::rwork(:),metric_weights(:)
    integer,allocatable::owner(:),position(:),ownership_count(:)
    integer::nlocal,m,i,k,q,rank,ierr,local_bad,global_bad,allocation_status,lwork,info
    integer::minint,maxint
    real(real64)::mintol,maxtol,scale,global_defect,phase_angle,quantum,safe_quantized
    complex(real64)::probe,projector_value
    integer(int64)::phase_bits,elements,term,min_fingerprint,max_fingerprint,recomputed_phase_fingerprint,&
      row_hash,local_sketch_hash,global_sketch_hash
    interface
      subroutine zgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)
        character,intent(in)::jobu,jobvt
        integer,intent(in)::m,n,lda,ldu,ldvt,lwork
        complex(8),intent(inout)::a(lda,*),work(*)
        real(8),intent(out)::s(*),rwork(*)
        complex(8),intent(out)::u(ldu,*),vt(ldvt,*)
        integer,intent(out)::info
      end subroutine
    end interface
    nlocal=size(row_ids);m=size(reference_rows,2);ok=.false.;message=''
    polar_defect=huge(1d0);fingerprint=0_int64;workspace_peak_bytes=0_int64
    local_bad=merge(0,1,global_row_count>=1.and.m>=1.and.size(reference_rows,1)==nlocal.and.&
      all(shape(target_rows)==[nlocal,m]).and.size(periodic_phase)==nlocal.and.&
      all(row_ids>=1_int64).and.all(row_ids<=int(global_row_count,int64)).and.&
      phase_fingerprint/=0_int64.and.phase_payload_fingerprint/=0_int64.and.&
      tolerance>=1d-15.and.tolerance<=1d-2.and.ieee_is_finite(tolerance).and.&
      all(ieee_is_finite(real(reference_rows))).and.all(ieee_is_finite(aimag(reference_rows))).and.&
      all(ieee_is_finite(real(target_rows))).and.all(ieee_is_finite(aimag(target_rows))).and.&
      all(ieee_is_finite(real(periodic_phase))).and.all(ieee_is_finite(aimag(periodic_phase))).and.&
      maxval(abs(abs(periodic_phase)-1d0))<=10d0*tolerance)
    if(present(integration_weights))then
      if(size(integration_weights)/=nlocal)local_bad=1
      if(local_bad==0)then
        if(.not.all(ieee_is_finite(integration_weights)).or.any(integration_weights<=0d0))local_bad=1
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid periodic-phase sector alignment contract';return;endif
    call MPI_Allreduce(global_row_count,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(global_row_count,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='periodic-phase row extent disagrees';return;endif
    call MPI_Allreduce(m,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(m,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='periodic-phase sector rank disagrees';return;endif
    call MPI_Allreduce(tolerance,mintol,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tolerance,maxtol,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.mintol/=maxtol)then;message='periodic-phase tolerance disagrees';return;endif
    call MPI_Allreduce(phase_fingerprint,min_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(phase_fingerprint,max_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.min_fingerprint/=max_fingerprint)then;message='periodic-phase provenance disagrees';return;endif
    call MPI_Allreduce(phase_payload_fingerprint,min_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(phase_payload_fingerprint,max_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.min_fingerprint/=max_fingerprint)then;message='periodic-phase payload receipt disagrees';return;endif
    local_bad=merge(1,0,m>huge(0)/m.or.m>huge(0)/5)
    if(local_bad==0)then
      if(int(nlocal,int64)>huge(0_int64)/int(m,int64).or.&
          int(m,int64)>huge(0_int64)/int(m,int64))local_bad=1
    endif
    if(local_bad==0)then
      elements=int(nlocal,int64)*int(m,int64);term=int(m,int64)*int(m,int64)
      if(term>huge(0_int64)/5_int64)then
        local_bad=1
      elseif(elements>huge(0_int64)-5_int64*term-int(m,int64))then
        local_bad=1
      endif
    endif
    if(local_bad==0)then
      elements=elements+5_int64*term+int(m,int64)
      if(elements>huge(0_int64)/16_int64)then
        local_bad=1
      else
        workspace_peak_bytes=16_int64*elements
      endif
    endif
    if(local_bad==0)then
      term=6_int64*int(m,int64)
      if(term>huge(0_int64)/8_int64.or.workspace_peak_bytes>huge(0_int64)-8_int64*term)then
        local_bad=1
      else
        workspace_peak_bytes=workspace_peak_bytes+8_int64*term
      endif
    endif
    if(local_bad==0)then
      term=3_int64*int(global_row_count,int64)
      if(term>huge(0_int64)/4_int64.or.workspace_peak_bytes>huge(0_int64)-4_int64*term)then
        local_bad=1
      else
        workspace_peak_bytes=workspace_peak_bytes+4_int64*term
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='periodic-phase workspace overflows';return;endif
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(link(m,m),left(m,m),right(m,m),polar(m,m),gram(m,m),singular_values(m),&
      aligned_rows(nlocal,m),rwork(max(1,5*m)),metric_weights(nlocal),work(1),owner(global_row_count),position(global_row_count),&
      ownership_count(global_row_count),sketch(m),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(aligned_rows))deallocate(aligned_rows)
      if(allocated(singular_values))deallocate(singular_values)
      message='periodic-phase alignment allocation failed';return
    endif
    metric_weights=1d0
    if(present(integration_weights))metric_weights=integration_weights
    owner=0;position=0;ownership_count=0
    do i=1,nlocal
      owner(int(row_ids(i)))=rank+1;position(int(row_ids(i)))=i;ownership_count(int(row_ids(i)))=1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,owner,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,position,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership_count/=1))then;message='periodic-phase rows are not uniquely owned';return;endif
    recomputed_phase_fingerprint=int(z'243F6A8885A308D3',int64)
    do k=1,global_row_count
      phase_bits=0_int64;if(rank==owner(k)-1)phase_bits=transfer(real(periodic_phase(position(k)),real64),phase_bits)
      call MPI_Bcast(phase_bits,1,MPI_INTEGER8,owner(k)-1,comm,ierr);if(ierr/=MPI_SUCCESS)return
      recomputed_phase_fingerprint=ieor(ishftc(recomputed_phase_fingerprint,11),phase_bits)
      if(rank==owner(k)-1)phase_bits=transfer(aimag(periodic_phase(position(k))),phase_bits)
      call MPI_Bcast(phase_bits,1,MPI_INTEGER8,owner(k)-1,comm,ierr);if(ierr/=MPI_SUCCESS)return
      recomputed_phase_fingerprint=ieor(ishftc(recomputed_phase_fingerprint,11),phase_bits)
    enddo
    if(recomputed_phase_fingerprint==0_int64)recomputed_phase_fingerprint=1_int64
    local_bad=merge(0,1,recomputed_phase_fingerprint==phase_payload_fingerprint)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='periodic-phase payload does not match its receipt';return;endif
    gram=matmul(conjg(transpose(reference_rows)),spread(metric_weights,2,m)*reference_rows)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    local_bad=merge(0,1,maxval(abs(gram))<=10d0*tolerance)
    gram=matmul(conjg(transpose(target_rows)),spread(metric_weights,2,m)*target_rows)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    if(maxval(abs(gram))>10d0*tolerance)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='periodic-phase sector frames are not orthonormal';return;endif
    link=matmul(conjg(transpose(target_rows)),spread(metric_weights*periodic_phase,2,m)*reference_rows)
    call MPI_Allreduce(MPI_IN_PLACE,link,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='periodic-phase localization link reduction failed';return;endif
    left=link;lwork=-1
    call zgesvd('A','A',m,m,left,m,singular_values,polar,m,right,m,work,lwork,rwork,info)
    local_bad=merge(0,1,info==0.and.ieee_is_finite(real(work(1))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='periodic-phase SVD query failed';return;endif
    local_bad=merge(0,1,real(work(1),real64)<=real(huge(0),real64))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='periodic-phase SVD workspace overflows';return;endif
    lwork=max(1,ceiling(real(work(1),real64)))
    term=16_int64*(int(lwork,int64)-1_int64)
    local_bad=merge(0,1,term>=0_int64.and.workspace_peak_bytes<=huge(0_int64)-term)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='periodic-phase workspace receipt overflows';return;endif
    workspace_peak_bytes=workspace_peak_bytes+term
    deallocate(work);allocate(work(lwork),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(aligned_rows))deallocate(aligned_rows)
      if(allocated(singular_values))deallocate(singular_values)
      message='periodic-phase SVD allocation failed';return
    endif
    left=link;call zgesvd('A','A',m,m,left,m,singular_values,polar,m,right,m,work,lwork,rwork,info)
    local_bad=merge(0,1,info==0.and.all(ieee_is_finite(singular_values)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='periodic-phase SVD failed';return;endif
    scale=max(1d0,maxval(singular_values))
    local_bad=merge(0,1,minval(singular_values)>tolerance*scale)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='singular periodic-phase sector link';return;endif
    polar=matmul(polar,right);aligned_rows=matmul(target_rows,polar)
    gram=matmul(conjg(transpose(polar)),polar);do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    polar_defect=maxval(abs(gram));call MPI_Allreduce(polar_defect,global_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    polar_defect=global_defect
    if(ierr/=MPI_SUCCESS.or.polar_defect>10d0*tolerance)then;message='periodic-phase polar is not unitary';return;endif
    fingerprint=ieor(int(z'A54FF53A5F1D36F1',int64),phase_fingerprint)
    fingerprint=ieor(ishftc(fingerprint,13),int(global_row_count,int64))
    fingerprint=ieor(ishftc(fingerprint,13),int(m,int64))
    quantum=100d0*tolerance;safe_quantized=0.25d0*real(huge(0_int64),real64)
    do q=1,2
      sketch=(0d0,0d0)
      do i=1,nlocal
        phase_angle=2d0*acos(-1d0)*modulo(real(row_ids(i),real64)*&
          merge(0.6180339887498948482d0,0.4142135623730950488d0,q==1),1d0)
        probe=cmplx(cos(phase_angle),sin(phase_angle),real64)/sqrt(real(global_row_count,real64))
        sketch=sketch+conjg(sqrt(metric_weights(i))*aligned_rows(i,:))*probe
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,sketch,m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      local_bad=0
      do i=1,nlocal
        projector_value=sqrt(metric_weights(i))*sum(aligned_rows(i,:)*sketch)
        if(.not.ieee_is_finite(real(projector_value)).or..not.ieee_is_finite(aimag(projector_value)).or.&
            abs(real(projector_value,real64)/quantum)>safe_quantized.or.&
            abs(aimag(projector_value)/quantum)>safe_quantized)then
          local_bad=1
        endif
      enddo
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
        message='periodic-phase projector sketch is not quantizable';return
      endif
      local_sketch_hash=0_int64
      do i=1,nlocal
        projector_value=sqrt(metric_weights(i))*sum(aligned_rows(i,:)*sketch)
        row_hash=ieor(int(z'9E3779B97F4A7C15',int64),row_ids(i))
        row_hash=ieor(ishftc(row_hash,13),int(q,int64))
        phase_bits=nint(real(projector_value,real64)/quantum,int64)
        row_hash=ieor(ishftc(row_hash,13),phase_bits)
        phase_bits=nint(aimag(projector_value)/quantum,int64)
        row_hash=ieor(ishftc(row_hash,13),phase_bits)
        local_sketch_hash=ieor(local_sketch_hash,row_hash)
      enddo
      call MPI_Allreduce(local_sketch_hash,global_sketch_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      fingerprint=ieor(ishftc(fingerprint,13),global_sketch_hash)
    enddo
    ok=.true.
  end subroutine align_dg_w90_character_sectors_by_periodic_phase

  subroutine sew_dg_w90_periodic_phase_conjugate_sector(comm,row_ids,aligned_rows,gamma_rows,&
      conjugate_rows,global_row_count,gamma_fingerprint,self_conjugate,gamma_sewing_defect,tolerance,aligned_conjugate_rows,&
      gamma_defect,workspace_peak_bytes,ok,message,implicit_identity,integration_weights)
    integer,intent(in)::comm,global_row_count
    integer(int64),intent(in)::gamma_fingerprint
    logical,intent(in)::self_conjugate
    logical,intent(in),optional::implicit_identity
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::aligned_rows(:,:),gamma_rows(:,:),conjugate_rows(:,:)
    real(real64),intent(in)::gamma_sewing_defect,tolerance
    real(real64),intent(in),optional::integration_weights(:)
    complex(real64),allocatable,intent(out)::aligned_conjugate_rows(:,:)
    real(real64),intent(out)::gamma_defect
    integer(int64),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::generated(:,:),remote_row(:),remote_gamma(:),local_operator_vector(:),&
      gram(:,:),sewing(:,:),candidates(:,:),gauge(:,:)
    integer,allocatable::owner(:),position(:),ownership_count(:)
    integer::nlocal,m,i,j,k,nfixed,rank,ierr,local_bad,global_bad,status,minint,maxint,flagint,identity_flag
    real(real64)::minimum,maximum,norm_value,operator_defect,local_operator_defect
    real(real64),allocatable::metric_weights(:)
    complex(real64)::overlap_value
    integer(int64)::elements,bits,recomputed_gamma_fingerprint,minhash,maxhash,term
    logical::receipt_valid
    logical::use_identity
    nlocal=size(row_ids);m=size(aligned_rows,2);use_identity=.false.
    if(present(implicit_identity))use_identity=implicit_identity
    ok=.false.;message=''
    gamma_defect=huge(1d0);workspace_peak_bytes=0_int64
    local_bad=merge(0,1,global_row_count>=1.and.m>=1.and.size(aligned_rows,1)==nlocal.and.&
      all(shape(conjugate_rows)==[nlocal,m]).and.&
      (use_identity.or.all(shape(gamma_rows)==[nlocal,global_row_count])).and.&
      all(row_ids>=1_int64).and.all(row_ids<=int(global_row_count,int64)).and.&
      gamma_fingerprint/=0_int64.and.gamma_sewing_defect>=0d0.and.gamma_sewing_defect<=tolerance.and.&
      ieee_is_finite(gamma_sewing_defect).and.&
      tolerance>=1d-15.and.tolerance<=1d-2.and.ieee_is_finite(tolerance).and.&
      all(ieee_is_finite(real(aligned_rows))).and.all(ieee_is_finite(aimag(aligned_rows))).and.&
      (use_identity.or.(all(ieee_is_finite(real(gamma_rows))).and.all(ieee_is_finite(aimag(gamma_rows))))) .and.&
      all(ieee_is_finite(real(conjugate_rows))).and.all(ieee_is_finite(aimag(conjugate_rows))))
    if(present(integration_weights))then
      if(size(integration_weights)/=nlocal)local_bad=1
      if(local_bad==0)then
        if(.not.all(ieee_is_finite(integration_weights)).or.any(integration_weights<=0d0))local_bad=1
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid periodic-phase Gamma sewing contract';return;endif
    call MPI_Allreduce(global_row_count,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(global_row_count,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='Gamma sewing row extent disagrees';return;endif
    call MPI_Allreduce(m,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(m,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='Gamma sewing sector rank disagrees';return;endif
    call MPI_Allreduce(tolerance,minimum,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(tolerance,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum/=maximum)then;message='Gamma sewing tolerance disagrees';return;endif
    call MPI_Allreduce(gamma_sewing_defect,minimum,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(gamma_sewing_defect,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum/=maximum)then;message='Gamma sewing receipt disagrees';return;endif
    call MPI_Allreduce(gamma_fingerprint,minhash,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(gamma_fingerprint,maxhash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minhash/=maxhash)then;message='Gamma sewing fingerprint disagrees';return;endif
    flagint=merge(1,0,self_conjugate)
    call MPI_Allreduce(flagint,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(flagint,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='Gamma conjugacy branch disagrees';return;endif
    identity_flag=merge(1,0,use_identity)
    call MPI_Allreduce(identity_flag,minint,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(identity_flag,maxint,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minint/=maxint)then;message='Gamma operator representation branch disagrees';return;endif
    local_bad=merge(1,0,int(nlocal,int64)>huge(0_int64)/int(m,int64).or.&
      int(m,int64)>huge(0_int64)/int(m,int64))
    if(local_bad==0)then
      if(int(nlocal,int64)>huge(0_int64)/(2_int64*int(m,int64)))local_bad=1
      if(int(m,int64)>huge(0_int64)/(4_int64*int(m,int64)))local_bad=1
    endif
    elements=0_int64;receipt_valid=local_bad==0
    if(receipt_valid)then
      term=2_int64*int(nlocal,int64)*int(m,int64);call checked_add(elements,term,receipt_valid)
      term=4_int64*int(m,int64)*int(m,int64);call checked_add(elements,term,receipt_valid)
      call checked_add(elements,int(m,int64),receipt_valid)
      call checked_add(elements,int(nlocal,int64),receipt_valid)
      call checked_add(elements,int(global_row_count,int64),receipt_valid)
      if(elements>huge(0_int64)/16_int64)receipt_valid=.false.
      if(int(global_row_count,int64)>huge(0_int64)/12_int64)receipt_valid=.false.
      if(receipt_valid)then
        term=16_int64*elements
        if(term>huge(0_int64)-12_int64*int(global_row_count,int64))receipt_valid=.false.
      endif
    endif
    local_bad=merge(0,1,receipt_valid)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='Gamma sewing workspace overflows';return;endif
    workspace_peak_bytes=16_int64*elements+12_int64*int(global_row_count,int64)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    allocate(aligned_conjugate_rows(nlocal,m),generated(nlocal,m),remote_row(m),&
      remote_gamma(merge(0,global_row_count,use_identity)),&
      local_operator_vector(nlocal),gram(m,m),sewing(m,m),candidates(m,2*m),gauge(m,m),metric_weights(nlocal),&
      owner(merge(0,global_row_count,use_identity)),position(merge(0,global_row_count,use_identity)),&
      ownership_count(merge(0,global_row_count,use_identity)),stat=status)
    call MPI_Allreduce(status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(aligned_conjugate_rows))deallocate(aligned_conjugate_rows)
      message='Gamma sewing allocation failed';return
    endif
    metric_weights=1d0
    if(present(integration_weights))metric_weights=integration_weights
    if(.not.use_identity)then
      owner=0;position=0;ownership_count=0
      do i=1,nlocal
        owner(int(row_ids(i)))=rank+1;position(int(row_ids(i)))=i
        ownership_count(int(row_ids(i)))=ownership_count(int(row_ids(i)))+1
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,owner,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(MPI_IN_PLACE,position,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
      call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.any(ownership_count/=1))then;message='Gamma sewing rows are not uniquely owned';return;endif
    endif
    recomputed_gamma_fingerprint=int(z'6A09E667F3BCC909',int64)
    if(use_identity)then
      recomputed_gamma_fingerprint=ieor(ishftc(recomputed_gamma_fingerprint,9),int(global_row_count,int64))
    else;do k=1,global_row_count
      remote_gamma=(0d0,0d0)
      if(use_identity)then
        remote_gamma(k)=1d0
      else
        if(rank==owner(k)-1)remote_gamma=gamma_rows(position(k),:)
        call MPI_Bcast(remote_gamma,global_row_count,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr);if(ierr/=MPI_SUCCESS)return
      endif
      do j=1,global_row_count
        bits=transfer(real(remote_gamma(j),real64),bits)
        recomputed_gamma_fingerprint=ieor(ishftc(recomputed_gamma_fingerprint,9),bits)
        bits=transfer(aimag(remote_gamma(j)),bits)
        recomputed_gamma_fingerprint=ieor(ishftc(recomputed_gamma_fingerprint,9),bits)
      enddo
    enddo;endif
    if(recomputed_gamma_fingerprint==0_int64)recomputed_gamma_fingerprint=1_int64
    local_bad=merge(0,1,recomputed_gamma_fingerprint==gamma_fingerprint)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='Gamma sewing payload does not match its receipt';return;endif
    local_operator_defect=0d0
    if(.not.use_identity)then
    do j=1,global_row_count;do k=1,global_row_count
      overlap_value=sum(conjg(gamma_rows(:,j))*gamma_rows(:,k))
      call MPI_Allreduce(MPI_IN_PLACE,overlap_value,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
      if(j==k)overlap_value=overlap_value-1d0
      local_operator_defect=max(local_operator_defect,abs(overlap_value))
    enddo;enddo
    do j=1,global_row_count
      local_operator_vector=(0d0,0d0)
      do k=1,global_row_count
        overlap_value=(0d0,0d0)
        if(rank==owner(k)-1)overlap_value=gamma_rows(position(k),j)
        call MPI_Bcast(overlap_value,1,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr);if(ierr/=MPI_SUCCESS)return
        local_operator_vector=local_operator_vector+gamma_rows(:,k)*conjg(overlap_value)
      enddo
      do i=1,nlocal
        overlap_value=local_operator_vector(i)
        if(int(row_ids(i))==j)overlap_value=overlap_value-1d0
        local_operator_defect=max(local_operator_defect,abs(overlap_value))
      enddo
    enddo
    endif
    call MPI_Allreduce(local_operator_defect,operator_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.operator_defect>10d0*tolerance.or.operator_defect>10d0*gamma_sewing_defect+tolerance)then
      message='Gamma sewing operator is not unitary and involutory';return
    endif
    if(use_identity)then
      generated=conjg(aligned_rows)
    else
      generated=(0d0,0d0)
      do k=1,global_row_count
        remote_row=(0d0,0d0)
        if(rank==owner(k)-1)remote_row=conjg(aligned_rows(position(k),:))
        call MPI_Bcast(remote_row,m,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr);if(ierr/=MPI_SUCCESS)return
        do i=1,nlocal;generated(i,:)=generated(i,:)+gamma_rows(i,k)*remote_row;enddo
      enddo
    endif
    gram=matmul(conjg(transpose(conjugate_rows)),spread(metric_weights,2,m)*conjugate_rows)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    local_bad=merge(0,1,maxval(abs(gram))<=10d0*tolerance)
    gram=matmul(conjg(transpose(generated)),spread(metric_weights,2,m)*generated)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    if(maxval(abs(gram))>10d0*tolerance)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='periodic-phase Gamma frames are not orthonormal';return;endif
    gram=matmul(conjg(transpose(conjugate_rows)),spread(metric_weights,2,m)*generated)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
    sewing=gram;gram=matmul(conjg(transpose(gram)),gram)
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    gamma_defect=maxval(abs(gram))
    call MPI_Allreduce(gamma_defect,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);gamma_defect=maximum
    if(ierr/=MPI_SUCCESS.or.gamma_defect>10d0*tolerance)then
      message='periodic-phase sectors violate Gamma conjugate pairing';return
    endif
    if(self_conjugate)then
      sewing=matmul(conjg(transpose(aligned_rows)),spread(metric_weights,2,m)*generated)
      call MPI_Allreduce(MPI_IN_PLACE,sewing,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);if(ierr/=MPI_SUCCESS)return
      candidates=(0d0,0d0)
      do j=1,m
        candidates(j,j)=1d0;candidates(:,j)=candidates(:,j)+sewing(:,j)
        candidates(j,m+j)=cmplx(0d0,1d0,real64)
        candidates(:,m+j)=candidates(:,m+j)-cmplx(0d0,1d0,real64)*sewing(:,j)
      enddo
      gauge=(0d0,0d0);nfixed=0
      do j=1,2*m
        remote_row=candidates(:,j)
        do k=1,nfixed;remote_row=remote_row-gauge(:,k)*dot_product(gauge(:,k),remote_row);enddo
        norm_value=sqrt(max(0d0,real(dot_product(remote_row,remote_row),real64)))
        if(norm_value<=10d0*tolerance)cycle
        nfixed=nfixed+1;gauge(:,nfixed)=remote_row/norm_value
        if(nfixed==m)exit
      enddo
      if(nfixed/=m)then;message='self-conjugate Gamma fixed space is rank deficient';return;endif
      aligned_conjugate_rows=matmul(aligned_rows,gauge)
      generated=matmul(generated,conjg(gauge))
      gamma_defect=maxval(abs(aligned_conjugate_rows-generated))
      call MPI_Allreduce(gamma_defect,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr);gamma_defect=maximum
      if(ierr/=MPI_SUCCESS.or.gamma_defect>10d0*tolerance)then
        message='self-conjugate sector cannot be fixed to a Gamma-real gauge';return
      endif
    else
      aligned_conjugate_rows=generated
    endif
    ok=.true.
  end subroutine sew_dg_w90_periodic_phase_conjugate_sector


  subroutine align_dg_w90_cross_character_sector_gauge(comm,row_ids,reference_sector_rows,&
      target_sector_rows,w90_reference_rows,localization_weights,global_row_count,retained_state_count,&
      w90_frame_fingerprint,w90_unitarity_defect,tolerance,&
      aligned_target_rows,singular_values,&
      polar_defect,fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::reference_sector_rows(:,:),target_sector_rows(:,:),w90_reference_rows(:,:)
    real(real64),intent(in)::localization_weights(:),w90_unitarity_defect,tolerance
    integer,intent(in)::global_row_count,retained_state_count
    integer(int64),intent(in)::w90_frame_fingerprint
    complex(real64),allocatable,intent(out)::aligned_target_rows(:,:)
    real(real64),allocatable,intent(out)::singular_values(:)
    real(real64),intent(out)::polar_defect
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::reference_projection(:,:),target_projection(:,:),link(:,:),&
      svd_left(:,:),svd_right(:,:),polar(:,:),gram(:,:),svd_work(:),remote_row(:),projector_row(:)
    real(real64),allocatable::svd_rwork(:)
    integer,allocatable::owner(:),position(:),ownership_count(:)
    integer::nlocal,n,m,nw,i,j,k,rank,ierr,local_bad,global_bad,allocation_status
    integer::minimum_n,maximum_n,minimum_m,maximum_m,minimum_nw,maximum_nw,svd_info,svd_lwork
    integer::retained_singular_count
    real(real64)::minimum_tolerance,maximum_tolerance,singular_scale,safe_weight,&
      minimum_w90_defect,maximum_w90_defect
    integer(int64)::weight_bits,minimum_weight_bits,maximum_weight_bits,&
      minimum_frame_fingerprint,maximum_frame_fingerprint
    complex(real64)::projector_value
    integer(int64)::complex_elements,real_elements,integer_elements,byte_term
    logical::receipt_valid
    interface
      subroutine zgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)
        character,intent(in)::jobu,jobvt
        integer,intent(in)::m,n,lda,ldu,ldvt,lwork
        complex(8),intent(inout)::a(lda,*),work(*)
        real(8),intent(out)::s(*),rwork(*)
        complex(8),intent(out)::u(ldu,*),vt(ldvt,*)
        integer,intent(out)::info
      end subroutine
    end interface
    ok=.false.;message='';polar_defect=huge(1d0);fingerprint=0_int64;workspace_peak_bytes=0_int64
    nlocal=size(row_ids);m=size(reference_sector_rows,2);nw=size(w90_reference_rows,2);n=global_row_count
    local_bad=merge(0,1,n>=1.and.m>=1.and.int(nw,int64)>=2_int64*int(m,int64).and.&
      nw==retained_state_count.and.w90_frame_fingerprint/=0_int64.and.size(localization_weights)==nw.and.&
      size(reference_sector_rows,1)==nlocal.and.&
      all(shape(target_sector_rows)==[nlocal,m]).and.size(w90_reference_rows,1)==nlocal.and.&
      tolerance>=1d-15.and.tolerance<=1d-2.and.ieee_is_finite(tolerance).and.&
      w90_unitarity_defect>=0d0.and.w90_unitarity_defect<=tolerance.and.ieee_is_finite(w90_unitarity_defect).and.&
      all(row_ids>=1_int64).and.all(row_ids<=int(n,int64)).and.&
      all(ieee_is_finite(localization_weights)).and.&
      all(ieee_is_finite(real(reference_sector_rows))).and.all(ieee_is_finite(aimag(reference_sector_rows))).and.&
      all(ieee_is_finite(real(target_sector_rows))).and.all(ieee_is_finite(aimag(target_sector_rows))).and.&
      all(ieee_is_finite(real(w90_reference_rows))).and.all(ieee_is_finite(aimag(w90_reference_rows))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='invalid Wannier90 cross-character alignment contract';return
    endif
    safe_weight=sqrt(huge(1d0))/max(1d0,16d0*real(n,real64)*real(nw,real64))
    local_bad=merge(0,1,maxval(abs(reference_sector_rows))<=1d0+tolerance.and.&
      maxval(abs(target_sector_rows))<=1d0+tolerance.and.maxval(abs(w90_reference_rows))<=1d0+tolerance.and.&
      maxval(abs(localization_weights))<=safe_weight)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='Wannier90 cross-character input magnitude is unsafe';return
    endif
    call MPI_Allreduce(n,minimum_n,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character metadata reduction failed';return;endif
    call MPI_Allreduce(n,maximum_n,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character metadata reduction failed';return;endif
    call MPI_Allreduce(m,minimum_m,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character metadata reduction failed';return;endif
    call MPI_Allreduce(m,maximum_m,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character metadata reduction failed';return;endif
    call MPI_Allreduce(nw,minimum_nw,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character metadata reduction failed';return;endif
    call MPI_Allreduce(nw,maximum_nw,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character metadata reduction failed';return;endif
    call MPI_Allreduce(tolerance,minimum_tolerance,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character metadata reduction failed';return;endif
    call MPI_Allreduce(tolerance,maximum_tolerance,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character metadata reduction failed';return;endif
    call MPI_Allreduce(w90_unitarity_defect,minimum_w90_defect,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character metadata reduction failed';return;endif
    call MPI_Allreduce(w90_unitarity_defect,maximum_w90_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character metadata reduction failed';return;endif
    do i=1,nw
      weight_bits=transfer(localization_weights(i),weight_bits)
      call MPI_Allreduce(weight_bits,minimum_weight_bits,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='Wannier90 localization-weight agreement reduction failed';return;endif
      call MPI_Allreduce(weight_bits,maximum_weight_bits,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_weight_bits/=maximum_weight_bits)then
        message='Wannier90 localization weights disagree across ranks';return
      endif
    enddo
    if(minimum_n/=maximum_n.or.minimum_m/=maximum_m.or.minimum_nw/=maximum_nw.or.&
      minimum_tolerance/=maximum_tolerance)then
      message='Wannier90 cross-character metadata disagree across ranks';return
    endif
    if(minimum_w90_defect/=maximum_w90_defect)then
      message='Wannier90 localization metadata disagree across ranks';return
    endif
    call MPI_Allreduce(w90_frame_fingerprint,minimum_frame_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 frame provenance reduction failed';return;endif
    call MPI_Allreduce(w90_frame_fingerprint,maximum_frame_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_frame_fingerprint/=maximum_frame_fingerprint)then
      message='Wannier90 frame provenance disagrees across ranks';return
    endif
    if(m>huge(0)/m.or.m>huge(0)/5.or.m>huge(0)/nw.or.&
      (nlocal>0.and.m>huge(0)/nlocal))then
      message='Wannier90 cross-character extent overflows';return
    endif
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character rank query failed';return;endif
    allocate(owner(n),position(n),ownership_count(n),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character ownership allocation failed';return;endif
    owner=0;position=0;ownership_count=0
    do i=1,nlocal
      if(row_ids(i)<=int(n,int64))then
        owner(int(row_ids(i)))=rank+1;position(int(row_ids(i)))=i;ownership_count(int(row_ids(i)))=1
      endif
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,owner,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character ownership reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,position,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character ownership reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership_count/=1))then
      message='Wannier90 cross-character rows are not uniquely owned';return
    endif
    allocate(reference_projection(m,nw),target_projection(m,nw),link(m,m),svd_left(m,m),&
      svd_right(m,m),polar(m,m),gram(m,m),singular_values(m),aligned_target_rows(nlocal,m),&
      remote_row(m),projector_row(n),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character workspace allocation failed';return;endif
    gram=matmul(conjg(transpose(reference_sector_rows)),reference_sector_rows)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 reference-sector Gram reduction failed';return;endif
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    if(maxval(abs(gram))>10d0*tolerance)then;message='Wannier90 reference sector is not orthonormal';return;endif
    gram=matmul(conjg(transpose(target_sector_rows)),target_sector_rows)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 target-sector Gram reduction failed';return;endif
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    if(maxval(abs(gram))>10d0*tolerance)then;message='Wannier90 target sector is not orthonormal';return;endif
    reference_projection=matmul(conjg(transpose(reference_sector_rows)),w90_reference_rows)
    target_projection=matmul(conjg(transpose(target_sector_rows)),w90_reference_rows)
    call MPI_Allreduce(MPI_IN_PLACE,reference_projection,m*nw,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 reference projection reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,target_projection,m*nw,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 target projection reduction failed';return;endif
    ! The left target projection makes the polar transform contragredient under
    ! a target-frame rotation: target*R, link -> R^H*link.
    do j=1,nw;target_projection(:,j)=target_projection(:,j)*localization_weights(j);enddo
    link=matmul(target_projection,conjg(transpose(reference_projection)))
    allocate(svd_rwork(max(1,5*m)),svd_work(1),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character SVD allocation failed';return;endif
    svd_left=link;svd_lwork=-1
    call zgesvd('A','A',m,m,svd_left,m,singular_values,polar,m,svd_right,m,svd_work,svd_lwork,svd_rwork,svd_info)
    local_bad=merge(0,1,svd_info==0.and.ieee_is_finite(real(svd_work(1))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='Wannier90 cross-character SVD workspace query failed';return
    endif
    if(real(svd_work(1),real64)>real(huge(0),real64))then
      message='Wannier90 cross-character SVD workspace overflows';return
    endif
    svd_lwork=max(1,ceiling(real(svd_work(1),real64)));deallocate(svd_work)
    allocate(svd_work(svd_lwork),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character SVD allocation failed';return;endif
    svd_left=link
    call zgesvd('A','A',m,m,svd_left,m,singular_values,polar,m,svd_right,m,svd_work,svd_lwork,svd_rwork,svd_info)
    local_bad=merge(0,1,svd_info==0.and.all(ieee_is_finite(singular_values)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='Wannier90 cross-character localization-link SVD failed';return
    endif
    singular_scale=max(1d0,maxval(singular_values))
    retained_singular_count=count(singular_values>tolerance*singular_scale)
    if(retained_singular_count>0.and.retained_singular_count<m)then
      if(abs(singular_values(retained_singular_count)-singular_values(retained_singular_count+1))<=&
        10d0*tolerance*singular_scale)then
        message='Wannier90 cross-character rank threshold splits a degenerate cluster';return
      endif
    endif
    local_bad=merge(0,1,minval(singular_values)>tolerance*singular_scale)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='singular Wannier90 cross-character localization link';return
    endif
    polar=matmul(polar,svd_right)
    gram=matmul(conjg(transpose(polar)),polar)
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    polar_defect=maxval(abs(gram))
    local_bad=merge(0,1,polar_defect<=10d0*tolerance)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character polar is not unitary';return;endif
    aligned_target_rows=matmul(target_sector_rows,polar)
    fingerprint=ieor(int(z'3C6EF372FE94F82B',int64),w90_frame_fingerprint)
    do k=1,n
      remote_row=(0d0,0d0)
      if(rank==owner(k)-1)remote_row=aligned_target_rows(position(k),:)
      call MPI_Bcast(remote_row,m,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character projector stream failed';return;endif
      projector_row=(0d0,0d0)
      do i=1,nlocal
        projector_row(int(row_ids(i)))=sum(remote_row*conjg(aligned_target_rows(i,:)))
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,projector_row,n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character fingerprint reduction failed';return;endif
      if(rank==0)then
        do j=1,n
          projector_value=projector_row(j)
          fingerprint=ieor(ishftc(fingerprint,13),nint(real(projector_value)/(100d0*tolerance),int64))
          fingerprint=ieor(ishftc(fingerprint,13),nint(aimag(projector_value)/(100d0*tolerance),int64))
        enddo
      endif
    enddo
    call MPI_Bcast(fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 cross-character fingerprint broadcast failed';return;endif
    complex_elements=0_int64;real_elements=0_int64;integer_elements=0_int64;receipt_valid=.true.
    call checked_add(complex_elements,size(reference_sector_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(target_sector_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(w90_reference_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(reference_projection,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(target_projection,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(link,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(svd_left,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(svd_right,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(polar,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(gram,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(aligned_target_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(svd_work,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(remote_row,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(projector_row,kind=int64),receipt_valid)
    call checked_add(complex_elements,max(int(m,int64)*int(m,int64),int(nlocal,int64)*int(m,int64)),receipt_valid)
    call checked_add(real_elements,size(singular_values,kind=int64),receipt_valid)
    call checked_add(real_elements,size(svd_rwork,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(row_ids,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(owner,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(position,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(ownership_count,kind=int64),receipt_valid)
    if(receipt_valid)call checked_product([complex_elements,16_int64],workspace_peak_bytes,receipt_valid)
    if(receipt_valid)call checked_product([real_elements,8_int64],byte_term,receipt_valid)
    if(receipt_valid)call checked_add(workspace_peak_bytes,byte_term,receipt_valid)
    if(receipt_valid)call checked_product([integer_elements,8_int64],byte_term,receipt_valid)
    if(receipt_valid)call checked_add(workspace_peak_bytes,byte_term,receipt_valid)
    if(.not.receipt_valid.or.workspace_peak_bytes<=0_int64)then
      workspace_peak_bytes=0_int64;message='Wannier90 cross-character workspace receipt overflows';return
    endif
    ok=.true.
  end subroutine align_dg_w90_cross_character_sector_gauge

  subroutine align_dg_w90_character_sector_gauge(comm,row_ids,sector_rows,reference_rows,&
      gamma_rows,conjugate_rows,gamma_sewing_defect,tolerance,&
      aligned_rows,aligned_conjugate_rows,singular_values,&
      canonical_channel_keys,polar_defect,gamma_defect,fingerprint,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::sector_rows(:,:),reference_rows(:,:),gamma_rows(:,:),conjugate_rows(:,:)
    real(real64),intent(in)::gamma_sewing_defect,tolerance
    complex(real64),allocatable,intent(out)::aligned_rows(:,:),aligned_conjugate_rows(:,:)
    real(real64),allocatable,intent(out)::singular_values(:)
    integer(int64),allocatable,intent(out)::canonical_channel_keys(:)
    real(real64),intent(out)::polar_defect,gamma_defect
    integer(int64),intent(out)::fingerprint,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::link(:,:),polar(:,:),svd_left(:,:),svd_right(:,:),gram(:,:),&
      remote_row(:),generated(:,:),svd_work(:)
    complex(real64),allocatable::projector_row(:)
    complex(real64),allocatable::ordered_reference(:,:)
    real(real64),allocatable::svd_rwork(:)
    integer,allocatable::owner(:),position(:),ownership_count(:)
    integer,allocatable::channel_order(:)
    integer::nlocal,n,m,i,j,k,rank,ierr,svd_info,svd_lwork,local_bad,global_bad,allocation_status
    integer::minimum_n,maximum_n,minimum_m,maximum_m,retained_singular_count
    integer(int64)::complex_elements,real_elements,integer_elements,byte_term
    real(real64)::minimum_tolerance,maximum_tolerance,minimum_gamma_receipt,maximum_gamma_receipt,&
      singular_scale,pivot_magnitude,candidate_magnitude
    complex(real64)::pivot_value,stream_value,phase_factor
    logical::receipt_valid
    complex(real64)::projector_value
    interface
      subroutine zgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)
        character,intent(in)::jobu,jobvt
        integer,intent(in)::m,n,lda,ldu,ldvt,lwork
        complex(8),intent(inout)::a(lda,*),work(*)
        real(8),intent(out)::s(*),rwork(*)
        complex(8),intent(out)::u(ldu,*),vt(ldvt,*)
        integer,intent(out)::info
      end subroutine
    end interface
    ok=.false.;message='';polar_defect=huge(1d0);gamma_defect=huge(1d0)
    fingerprint=0_int64;workspace_peak_bytes=0_int64
    nlocal=size(row_ids);m=size(sector_rows,2);n=size(gamma_rows,2)
    local_bad=merge(0,1,n>=1.and.m>=1.and.nlocal>=0.and.size(sector_rows,1)==nlocal.and.&
        all(shape(reference_rows)==[nlocal,m]).and.all(shape(conjugate_rows)==[nlocal,m]).and.&
        gamma_sewing_defect>=0d0.and.gamma_sewing_defect<=tolerance.and.&
        ieee_is_finite(gamma_sewing_defect).and.&
        size(gamma_rows,1)==nlocal.and.tolerance>=1d-15.and.tolerance<=1d-2.and.&
        ieee_is_finite(tolerance).and.all(row_ids>=1_int64).and.all(row_ids<=int(n,int64)).and.&
        all(ieee_is_finite(real(sector_rows))).and.all(ieee_is_finite(aimag(sector_rows))).and.&
        all(ieee_is_finite(real(reference_rows))).and.all(ieee_is_finite(aimag(reference_rows))).and.&
        all(ieee_is_finite(real(gamma_rows))).and.all(ieee_is_finite(aimag(gamma_rows))).and.&
        all(ieee_is_finite(real(conjugate_rows))).and.all(ieee_is_finite(aimag(conjugate_rows))))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then
      message='invalid Wannier90 character-sector alignment contract';return
    endif
    call MPI_Allreduce(n,minimum_n,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 metadata agreement reduction failed';return;endif
    call MPI_Allreduce(n,maximum_n,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 metadata agreement reduction failed';return;endif
    call MPI_Allreduce(m,minimum_m,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 metadata agreement reduction failed';return;endif
    call MPI_Allreduce(m,maximum_m,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_n/=maximum_n.or.minimum_m/=maximum_m)then
      message='Wannier90 alignment metadata disagree across ranks';return
    endif
    call MPI_Allreduce(tolerance,minimum_tolerance,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 metadata agreement reduction failed';return;endif
    call MPI_Allreduce(tolerance,maximum_tolerance,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 metadata agreement reduction failed';return;endif
    call MPI_Allreduce(gamma_sewing_defect,minimum_gamma_receipt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 metadata agreement reduction failed';return;endif
    call MPI_Allreduce(gamma_sewing_defect,maximum_gamma_receipt,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_tolerance/=maximum_tolerance.or.&
        minimum_gamma_receipt/=maximum_gamma_receipt)then
      message='Wannier90 alignment metadata disagree across ranks';return
    endif
    if(m>0.and.m>huge(0)/m)then;message='Wannier90 sector-link MPI count overflows';return;endif
    if(m>huge(0)/5)then;message='Wannier90 SVD workspace extent overflows';return;endif
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 communicator rank query failed';return;endif
    allocate(owner(n),position(n),ownership_count(n),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='Wannier90 ownership allocation failed';return;endif
    owner=0;position=0;ownership_count=0
    do i=1,nlocal
      owner(int(row_ids(i)))=rank+1;position(int(row_ids(i)))=i;ownership_count(int(row_ids(i)))=1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,owner,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 ownership rank reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,position,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 ownership position reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,n,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership_count/=1))then
      message='Wannier90 sector rows are not uniquely owned';return
    endif
    allocate(link(m,m),polar(m,m),svd_left(m,m),svd_right(m,m),gram(m,m),singular_values(m),&
      ordered_reference(nlocal,m),channel_order(m),canonical_channel_keys(m),&
      aligned_rows(nlocal,m),aligned_conjugate_rows(nlocal,m),generated(nlocal,m),remote_row(m),&
      projector_row(n),&
      stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='Wannier90 sector workspace allocation failed';return;endif
    ordered_reference=reference_rows
    do j=1,m
      candidate_magnitude=maxval(abs(ordered_reference(:,j)))
      call MPI_Allreduce(candidate_magnitude,pivot_magnitude,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='Wannier90 reference pivot reduction failed';return;endif
      pivot_value=(0d0,0d0)
      do k=1,n
        stream_value=(0d0,0d0)
        if(rank==owner(k)-1)stream_value=ordered_reference(position(k),j)
        call MPI_Bcast(stream_value,1,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;message='Wannier90 reference phase stream failed';return;endif
        if(abs(stream_value)>=pivot_magnitude-10d0*tolerance*max(1d0,pivot_magnitude))then
          pivot_value=stream_value;exit
        endif
      enddo
      pivot_magnitude=abs(pivot_value)
      if(pivot_magnitude<=tolerance)then;message='singular Wannier90 reference channel';return;endif
      phase_factor=conjg(pivot_value)/pivot_magnitude
      ordered_reference(:,j)=ordered_reference(:,j)*phase_factor
      canonical_channel_keys(j)=int(z'BB67AE8584CAA73B',int64)
      do k=1,n
        stream_value=(0d0,0d0)
        if(rank==owner(k)-1)stream_value=ordered_reference(position(k),j)
        call MPI_Bcast(stream_value,1,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;message='Wannier90 reference key stream failed';return;endif
        canonical_channel_keys(j)=ieor(ishftc(canonical_channel_keys(j),11),&
          nint(real(stream_value,real64)/(1000d0*tolerance),int64))
        canonical_channel_keys(j)=ieor(ishftc(canonical_channel_keys(j),11),&
          nint(aimag(stream_value)/(1000d0*tolerance),int64))
      enddo
    enddo
    channel_order=[(i,i=1,m)]
    do i=2,m
      k=channel_order(i);j=i-1
      do while(j>=1)
        if(canonical_channel_keys(channel_order(j))<=canonical_channel_keys(k))exit
        channel_order(j+1)=channel_order(j);j=j-1
      enddo
      channel_order(j+1)=k
    enddo
    do i=2,m
      if(canonical_channel_keys(channel_order(i))==canonical_channel_keys(channel_order(i-1)))then
        message='Wannier90 canonical reference-channel fingerprints collide';return
      endif
    enddo
    ordered_reference=ordered_reference(:,channel_order)
    canonical_channel_keys=canonical_channel_keys(channel_order)
    gram=matmul(conjg(transpose(sector_rows)),sector_rows)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    if(ierr/=MPI_SUCCESS.or.maxval(abs(gram))>10d0*tolerance)then
      message='Wannier90 input sector frame is not orthonormal';return
    endif
    gram=matmul(conjg(transpose(ordered_reference)),ordered_reference)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    if(ierr/=MPI_SUCCESS.or.maxval(abs(gram))>10d0*tolerance)then
      message='Wannier90 reference frame is not orthonormal';return
    endif
    gram=matmul(conjg(transpose(conjugate_rows)),conjugate_rows)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    if(ierr/=MPI_SUCCESS.or.maxval(abs(gram))>10d0*tolerance)then
      message='Wannier90 conjugate sector frame is not orthonormal';return
    endif
    link=matmul(conjg(transpose(sector_rows)),ordered_reference)
    call MPI_Allreduce(MPI_IN_PLACE,link,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 localization-link reduction failed';return;endif
    allocate(svd_rwork(max(1,5*m)),svd_work(1),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='Wannier90 SVD query allocation failed';return;endif
    svd_left=link
    svd_lwork=-1
    call zgesvd('A','A',m,m,svd_left,m,singular_values,polar,m,svd_right,m,svd_work,&
      svd_lwork,svd_rwork,svd_info)
    if(svd_info/=0.or..not.ieee_is_finite(real(svd_work(1))))then
      message='Wannier90 localization-link SVD workspace query failed';return
    endif
    if(real(svd_work(1),real64)>real(huge(0),real64))then
      message='Wannier90 SVD workspace extent overflows';return
    endif
    svd_lwork=max(1,ceiling(real(svd_work(1),real64)));deallocate(svd_work)
    allocate(svd_work(svd_lwork),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='Wannier90 SVD workspace allocation failed';return;endif
    svd_left=link
    call zgesvd('A','A',m,m,svd_left,m,singular_values,polar,m,svd_right,m,svd_work,&
      svd_lwork,svd_rwork,svd_info)
    if(svd_info/=0.or..not.all(ieee_is_finite(singular_values)))then
      message='Wannier90 localization-link SVD failed';return
    endif
    singular_scale=max(1d0,maxval(singular_values))
    retained_singular_count=count(singular_values>tolerance*singular_scale)
    if(retained_singular_count>0.and.retained_singular_count<m)then
      if(abs(singular_values(retained_singular_count)-singular_values(retained_singular_count+1))<=&
          10d0*tolerance*singular_scale)then
        message='Wannier90 rank threshold splits a degenerate singular-value cluster';return
      endif
    endif
    if(minval(singular_values)<=tolerance*singular_scale)then
      message='singular Wannier90 localization link';return
    endif
    polar=matmul(polar,svd_right)
    gram=matmul(conjg(transpose(polar)),polar)
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    polar_defect=maxval(abs(gram));if(polar_defect>10d0*tolerance)then
      message='Wannier90 localization-link polar factor is not unitary';return
    endif
    aligned_rows=matmul(sector_rows,polar)
    generated=(0d0,0d0)
    do k=1,n
      remote_row=(0d0,0d0)
      if(rank==owner(k)-1)remote_row=conjg(aligned_rows(position(k),:))
      call MPI_Bcast(remote_row,m,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='Wannier90 Gamma row stream failed';return;endif
      do i=1,nlocal;generated(i,:)=generated(i,:)+gamma_rows(i,k)*remote_row;enddo
    enddo
    aligned_conjugate_rows=generated
    gram=matmul(conjg(transpose(generated)),generated)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    if(ierr/=MPI_SUCCESS.or.maxval(abs(gram))>10d0*tolerance)then
      message='Wannier90 Gamma image leaks outside a unitary conjugate-sector frame';return
    endif
    gram=matmul(conjg(transpose(conjugate_rows)),generated)
    call MPI_Allreduce(MPI_IN_PLACE,gram,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 Gamma overlap reduction failed';return;endif
    gram=matmul(conjg(transpose(gram)),gram)
    do i=1,m;gram(i,i)=gram(i,i)-1d0;enddo
    gamma_defect=maxval(abs(gram));if(gamma_defect>10d0*tolerance)then
      message='Wannier90 aligned sectors violate Gamma conjugate pairing';return
    endif
    fingerprint=int(z'6A09E667F3BCC909',int64)
    do k=1,n
      remote_row=(0d0,0d0)
      if(rank==owner(k)-1)remote_row=aligned_rows(position(k),:)
      call MPI_Bcast(remote_row,m,MPI_DOUBLE_COMPLEX,owner(k)-1,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='Wannier90 aligned-projector row stream failed';return;endif
      projector_row=(0d0,0d0)
      do i=1,nlocal
        projector_row(int(row_ids(i)))=sum(remote_row*conjg(aligned_rows(i,:)))
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,projector_row,n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='Wannier90 aligned-projector fingerprint reduction failed';return;endif
      if(rank==0)then
        do j=1,n
          projector_value=projector_row(j)
          fingerprint=ieor(ishftc(fingerprint,13),&
            nint(real(projector_value,real64)/(100d0*tolerance),int64))
          fingerprint=ieor(ishftc(fingerprint,13),nint(aimag(projector_value)/(100d0*tolerance),int64))
        enddo
      endif
    enddo
    call MPI_Bcast(fingerprint,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Wannier90 alignment fingerprint broadcast failed';return;endif
    complex_elements=0_int64;real_elements=0_int64;integer_elements=0_int64;receipt_valid=.true.
    call checked_add(complex_elements,size(sector_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(reference_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(gamma_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(conjugate_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(link,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(polar,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(svd_left,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(svd_right,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(svd_work,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(gram,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(ordered_reference,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(aligned_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(aligned_conjugate_rows,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(generated,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(remote_row,kind=int64),receipt_valid)
    call checked_add(complex_elements,size(projector_row,kind=int64),receipt_valid)
    ! Conservative allowance for the largest MATMUL/LAPACK temporary owned by this routine.
    call checked_add(complex_elements,max(int(m,int64)*int(m,int64),int(nlocal,int64)*int(m,int64)),receipt_valid)
    call checked_add(real_elements,size(singular_values,kind=int64),receipt_valid)
    call checked_add(real_elements,size(svd_rwork,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(row_ids,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(canonical_channel_keys,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(owner,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(position,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(ownership_count,kind=int64),receipt_valid)
    call checked_add(integer_elements,size(channel_order,kind=int64),receipt_valid)
    if(receipt_valid)call checked_product([complex_elements,16_int64],workspace_peak_bytes,receipt_valid)
    if(receipt_valid)call checked_product([real_elements,8_int64],byte_term,receipt_valid)
    if(receipt_valid)call checked_add(workspace_peak_bytes,byte_term,receipt_valid)
    if(receipt_valid)call checked_product([integer_elements,8_int64],byte_term,receipt_valid)
    if(receipt_valid)call checked_add(workspace_peak_bytes,byte_term,receipt_valid)
    if(.not.receipt_valid.or.workspace_peak_bytes<=0_int64)then
      workspace_peak_bytes=0_int64;message='Wannier90 alignment workspace receipt overflows';return
    endif
    ok=.true.;message=''
  end subroutine align_dg_w90_character_sector_gauge
#endif

  subroutine validate_dg_w90_convergence_log(path,maximum_iterations,iterations,ok,message)
    character(*),intent(in)::path
    integer,intent(in)::maximum_iterations
    integer,intent(out)::iterations
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::unit,io,parsed_iteration
    character(1024)::line
    logical::exists,have_final,have_convergence
    iterations=-1;ok=.false.;message='';have_final=.false.;have_convergence=.false.
    if(len_trim(path)==0.or.maximum_iterations<1)then
      message='invalid Wannier90 convergence-log contract';return
    endif
    inquire(file=trim(path),exist=exists)
    if(.not.exists)then;message='Wannier90 convergence log is missing';return;endif
    open(newunit=unit,file=trim(path),status='old',action='read',iostat=io)
    if(io/=0)then;message='Wannier90 convergence log cannot be opened';return;endif
    do
      read(unit,'(a)',iostat=io)line
      if(io/=0)exit
      if(index(line,'<-- CONV')>0)then
        read(line,*,iostat=io)parsed_iteration
        if(io==0)iterations=max(iterations,parsed_iteration)
        io=0
      endif
      if(index(line,'Wannierisation convergence criteria satisfied')>0)have_convergence=.true.
      if(index(adjustl(line),'Final State')==1)have_final=.true.
    enddo
    close(unit)
    if(.not.have_final)then;message='Wannier90 convergence log has no final state';return;endif
    if(.not.have_convergence)then;message='Wannier90 did not report convergence';return;endif
    if(iterations<0.or.iterations>=maximum_iterations)then
      message='Wannier90 exhausted its iteration limit';return
    endif
    ok=.true.;message=''
  end subroutine validate_dg_w90_convergence_log

  subroutine inherit_dg_w90_affine_receipts(transform,affine_subspace_defect,tolerance,&
      identity_defect,unitarity_defect,closure_defect,workspace_peak_bytes,ok,message)
    complex(real64),intent(in)::transform(:,:)
    real(real64),intent(in)::affine_subspace_defect,tolerance
    real(real64),intent(out)::identity_defect,unitarity_defect,closure_defect
    integer(int64),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::gram(:,:)
    integer::nstate,i
    ok=.false.;message='';identity_defect=huge(1d0);unitarity_defect=huge(1d0)
    closure_defect=huge(1d0);workspace_peak_bytes=0_int64;nstate=size(transform,1)
    if(nstate<1.or.size(transform,2)/=nstate.or.tolerance<=0d0.or.&
       affine_subspace_defect<0d0.or..not.ieee_is_finite(tolerance).or.&
       .not.ieee_is_finite(affine_subspace_defect).or.&
       .not.all(ieee_is_finite(real(transform))).or.&
       .not.all(ieee_is_finite(aimag(transform))))then
      message='invalid MLWF affine-receipt inheritance contract';return
    endif
    allocate(gram(nstate,nstate));gram=matmul(conjg(transpose(transform)),transform)
    do i=1,nstate;gram(i,i)=gram(i,i)-1d0;enddo
    unitarity_defect=maxval(abs(gram))
    identity_defect=affine_subspace_defect;closure_defect=affine_subspace_defect
    workspace_peak_bytes=int(storage_size((0d0,0d0))/8,int64)*int(size(gram),int64)
    ok=max(identity_defect,max(unitarity_defect,closure_defect))<=tolerance
    if(.not.ok)message='MLWF gauge cannot inherit the accepted affine proof'
  end subroutine inherit_dg_w90_affine_receipts

  subroutine apply_dg_w90_gamma_transform(comm,physical_ids,values,gradients,transform,centers,&
      tolerance,ok,message,spreads)
    integer,intent(in)::comm
    integer(int64),intent(in)::physical_ids(:)
    complex(real64),intent(inout)::values(:,:),transform(:,:)
    complex(real64),intent(inout),optional::gradients(:,:,:)
    real(real64),intent(inout)::centers(:,:)
    real(real64),intent(inout),optional::spreads(:)
    real(real64),intent(in)::tolerance
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nstate,npoint,i,j,k,axis,ierr,status
    integer,allocatable::order(:)
    logical,allocatable::used(:)
    complex(real64),allocatable::ordered_transform(:,:),point_values(:),point_gradients(:,:),gram(:,:)
    real(real64),allocatable::ordered_centers(:,:),ordered_spreads(:),local_maximum(:),global_maximum(:)
    integer(int64),allocatable::local_id(:),global_id(:)
    complex(real64),allocatable::local_pivot(:),global_pivot(:)
    real(real64)::scale
    logical::precedes
    ok=.false.;message='';status=0;nstate=size(values,1);npoint=size(values,2)
    if(nstate<=0.or.size(values,2)/=size(physical_ids).or.any(shape(transform)/=[nstate,nstate]).or.&
        any(shape(centers)/=[3,nstate]).or.tolerance<=0d0.or..not.ieee_is_finite(tolerance).or.&
        .not.all(ieee_is_finite(real(values))).or..not.all(ieee_is_finite(aimag(values))).or.&
        .not.all(ieee_is_finite(real(transform))).or..not.all(ieee_is_finite(aimag(transform))).or.&
        .not.all(ieee_is_finite(centers)).or.any(physical_ids<=0_int64))status=1
    if(present(gradients))then
      if(any(shape(gradients)/=[3,nstate,npoint]))status=1
      if(status==0)then
        if(.not.all(ieee_is_finite(real(gradients))).or.&
            .not.all(ieee_is_finite(aimag(gradients))))status=1
      endif
    endif
    if(present(spreads))then
      if(size(spreads)/=nstate.or..not.all(ieee_is_finite(spreads)))status=1
    endif
    if(maxval(abs(aimag(transform)))>tolerance*max(1d0,maxval(abs(transform))))status=1
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then;message='invalid Gamma MLWF transform contract';return;endif
    allocate(gram(nstate,nstate));gram=matmul(conjg(transpose(transform)),transform)
    do i=1,nstate;gram(i,i)=gram(i,i)-1d0;enddo
    if(maxval(abs(gram))>tolerance*max(1d0,real(nstate,real64)))status=2
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then;message='Gamma MLWF transform is not unitary';return;endif
    allocate(order(nstate),used(nstate));used=.false.
    do i=1,nstate
      order(i)=0
      do j=1,nstate
        if(used(j))cycle
        if(order(i)==0)then
          order(i)=j
        else
          precedes=.false.
          do axis=1,3
            if(modulo(centers(axis,j),1d0)<modulo(centers(axis,order(i)),1d0)-tolerance)then
              precedes=.true.;exit
            else if(modulo(centers(axis,j),1d0)>modulo(centers(axis,order(i)),1d0)+tolerance)then
              exit
            endif
          enddo
          if(.not.precedes.and.all(abs(modulo(centers(:,j),1d0)-&
              modulo(centers(:,order(i)),1d0))<=tolerance))then
            do k=1,nstate
              if(abs(transform(k,j))>abs(transform(k,order(i)))+tolerance)then
                precedes=.true.;exit
              else if(abs(transform(k,j))<abs(transform(k,order(i)))-tolerance)then
                exit
              endif
            enddo
          endif
          if(precedes)order(i)=j
        endif
      enddo
      used(order(i))=.true.
    enddo
    allocate(ordered_transform(nstate,nstate),ordered_centers(3,nstate))
    ordered_transform=transform(:,order);ordered_centers=centers(:,order)
    if(present(spreads))then
      allocate(ordered_spreads(nstate));ordered_spreads=spreads(order)
    endif
    allocate(point_values(nstate))
    if(present(gradients))allocate(point_gradients(3,nstate))
    do j=1,npoint
      point_values=matmul(transpose(ordered_transform),values(:,j))
      if(present(gradients))then
        do axis=1,3
          point_gradients(axis,:)=matmul(transpose(ordered_transform),gradients(axis,:,j))
        enddo
        gradients(:,:,j)=point_gradients
      endif
      values(:,j)=point_values
    enddo
    allocate(local_maximum(nstate),global_maximum(nstate),local_id(nstate),global_id(nstate),&
      local_pivot(nstate),global_pivot(nstate))
    do i=1,nstate
      if(npoint>0)then
        j=maxloc(abs(values(i,:)),dim=1);local_maximum(i)=abs(values(i,j))
      else
        local_maximum(i)=-1d0
      endif
    enddo
    call MPI_Allreduce(local_maximum,global_maximum,nstate,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    local_id=huge(0_int64)
    do i=1,nstate;do j=1,npoint
      scale=max(1d0,global_maximum(i))
      if(abs(abs(values(i,j))-global_maximum(i))<=tolerance*scale)&
        local_id(i)=min(local_id(i),physical_ids(j))
    enddo;enddo
    call MPI_Allreduce(local_id,global_id,nstate,MPI_INTEGER8,MPI_MIN,comm,ierr)
    local_pivot=(0d0,0d0)
    do i=1,nstate;do j=1,npoint
      if(physical_ids(j)==global_id(i))local_pivot(i)=values(i,j)
    enddo;enddo
    call MPI_Allreduce(local_pivot,global_pivot,nstate,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(global_maximum<=tolerance).or.&
        maxval(abs(aimag(global_pivot)))>tolerance*max(1d0,maxval(abs(global_pivot))))then
      message='cannot determine canonical Gamma MLWF signs';return
    endif
    do i=1,nstate
      if(real(global_pivot(i),real64)<0d0)then
        ordered_transform(:,i)=-ordered_transform(:,i);values(i,:)=-values(i,:)
        if(present(gradients))gradients(:,i,:)=-gradients(:,i,:)
      endif
    enddo
    transform=ordered_transform;centers=ordered_centers
    if(present(spreads))spreads=ordered_spreads
    ok=.true.
#else
    ok=.false.;message='Gamma MLWF transform application requires MPI'
#endif
  end subroutine apply_dg_w90_gamma_transform

  subroutine assemble_dg_w90_gamma_matrices(comm,values,anchors,weights,fractional,nncell,&
      coordinator_byte_limit,m_matrix,a_matrix,coordinator_bytes,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,nncell(:,:)
    complex(real64),intent(in)::values(:,:),anchors(:,:)
    real(real64),intent(in)::weights(:),fractional(:,:)
    integer(int64),intent(in)::coordinator_byte_limit
    complex(real64),allocatable,intent(out)::m_matrix(:,:,:),a_matrix(:,:)
    integer(int64),intent(out)::coordinator_bytes,workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,parameter::tile_size=32
    integer::rank,ierr,status,nband,nwann,npoint,nntot,m0,m1,n0,n1,b,p,m,n,count,allocation_status
    integer::local_dimensions(3),minimum_dimensions(3),maximum_dimensions(3)
    integer(int64)::output_elements,output_bytes,tile_bytes,complex_bytes,peak
    integer(int64)::minimum_limit,maximum_limit
    real(real64)::angle
    complex(real64)::phase
    complex(real64),allocatable::local_tile(:,:),reduced_tile(:,:)
    logical::arithmetic_ok
    ok=.false.;message='';coordinator_bytes=0_int64;workspace_peak_bytes=0_int64;status=0
    call MPI_Comm_rank(comm,rank,ierr)
    nband=size(values,1);npoint=size(values,2);nwann=size(anchors,1);nntot=size(nncell,2)
    if(ierr/=MPI_SUCCESS.or.nband<=0.or.nwann/=nband.or.npoint<0.or.&
        size(anchors,2)/=npoint.or.size(weights)/=npoint.or.&
        any(shape(fractional)/=[3,npoint]).or.size(nncell,1)/=3.or.nntot<=0.or.&
        coordinator_byte_limit<0_int64.or..not.all(ieee_is_finite(real(values))).or.&
        .not.all(ieee_is_finite(aimag(values))).or..not.all(ieee_is_finite(real(anchors))).or.&
        .not.all(ieee_is_finite(aimag(anchors))).or..not.all(ieee_is_finite(weights)).or.&
        .not.all(ieee_is_finite(fractional)).or.any(weights<0d0))status=1
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then
      allocate(m_matrix(0,0,0),a_matrix(0,0));message='invalid distributed Wannier90 matrix contract';return
    endif
    local_dimensions=[nband,nwann,nntot]
    call MPI_Allreduce(local_dimensions,minimum_dimensions,3,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(local_dimensions,maximum_dimensions,3,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(coordinator_byte_limit,minimum_limit,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(coordinator_byte_limit,maximum_limit,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(minimum_dimensions/=maximum_dimensions).or.minimum_limit/=maximum_limit)status=1
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then
      allocate(m_matrix(0,0,0),a_matrix(0,0));message='rank-inconsistent Wannier90 matrix contract';return
    endif
    call estimate_dg_w90_coordinator_bytes(nband,nwann,nntot,1,coordinator_bytes,&
      arithmetic_ok,message)
    if(.not.arithmetic_ok)status=2
    if(arithmetic_ok.and.coordinator_bytes>coordinator_byte_limit)status=3
    complex_bytes=int(storage_size((0d0,0d0))/8,int64)
    call checked_product([int(nband,int64),int(nband,int64),int(nntot,int64)],&
      output_elements,arithmetic_ok)
    if(arithmetic_ok)call checked_add(output_elements,int(nband,int64)*int(nwann,int64),arithmetic_ok)
    if(arithmetic_ok)call checked_product([output_elements,complex_bytes],output_bytes,arithmetic_ok)
    if(.not.arithmetic_ok)status=2
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then
      allocate(m_matrix(0,0,0),a_matrix(0,0));workspace_peak_bytes=0_int64
      if(status==3)message='Wannier90 coordinator byte limit exceeded'
      if(status==2)message='Wannier90 matrix byte estimate overflow'
      return
    endif
    allocation_status=0
    if(rank==0)then
      allocate(m_matrix(nband,nband,nntot),a_matrix(nband,nwann),stat=allocation_status)
      if(allocation_status==0)then;m_matrix=(0d0,0d0);a_matrix=(0d0,0d0);endif
    else
      allocate(m_matrix(0,0,0),a_matrix(0,0),stat=allocation_status)
    endif
    call MPI_Allreduce(MPI_IN_PLACE,allocation_status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(allocation_status/=0.or.ierr/=MPI_SUCCESS)then
      if(.not.allocated(m_matrix))allocate(m_matrix(0,0,0))
      if(.not.allocated(a_matrix))allocate(a_matrix(0,0))
      message='cannot allocate Wannier90 coordinator matrices';return
    endif
    peak=merge(output_bytes,0_int64,rank==0)
    do b=1,nntot
      do n0=1,nband,tile_size
        n1=min(n0+tile_size-1,nband)
        do m0=1,nband,tile_size
          m1=min(m0+tile_size-1,nband);allocate(local_tile(m1-m0+1,n1-n0+1));local_tile=(0d0,0d0)
          do p=1,npoint
            angle=-2d0*acos(-1d0)*dot_product(real(nncell(:,b),real64),fractional(:,p))
            phase=cmplx(cos(angle),sin(angle),real64)
            do n=n0,n1;do m=m0,m1
              local_tile(m-m0+1,n-n0+1)=local_tile(m-m0+1,n-n0+1)+&
                weights(p)*conjg(values(m,p))*phase*values(n,p)
            enddo;enddo
          enddo
          count=size(local_tile);allocate(reduced_tile(size(local_tile,1),size(local_tile,2)))
          call MPI_Reduce(local_tile,reduced_tile,count,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
          call checked_product([2_int64,int(count,int64),complex_bytes],tile_bytes,arithmetic_ok)
          if(arithmetic_ok)peak=max(peak,merge(output_bytes,0_int64,rank==0)+tile_bytes)
          if(rank==0.and.ierr==MPI_SUCCESS)m_matrix(m0:m1,n0:n1,b)=reduced_tile
          deallocate(local_tile,reduced_tile);if(ierr/=MPI_SUCCESS)status=4
        enddo
      enddo
    enddo
    do n0=1,nwann,tile_size
      n1=min(n0+tile_size-1,nwann)
      do m0=1,nband,tile_size
        m1=min(m0+tile_size-1,nband);allocate(local_tile(m1-m0+1,n1-n0+1));local_tile=(0d0,0d0)
        do p=1,npoint;do n=n0,n1;do m=m0,m1
          local_tile(m-m0+1,n-n0+1)=local_tile(m-m0+1,n-n0+1)+&
            weights(p)*conjg(values(m,p))*anchors(n,p)
        enddo;enddo;enddo
        count=size(local_tile);allocate(reduced_tile(size(local_tile,1),size(local_tile,2)))
        call MPI_Reduce(local_tile,reduced_tile,count,MPI_DOUBLE_COMPLEX,MPI_SUM,0,comm,ierr)
        call checked_product([2_int64,int(count,int64),complex_bytes],tile_bytes,arithmetic_ok)
        if(arithmetic_ok)peak=max(peak,merge(output_bytes,0_int64,rank==0)+tile_bytes)
        if(rank==0.and.ierr==MPI_SUCCESS)a_matrix(m0:m1,n0:n1)=reduced_tile
        deallocate(local_tile,reduced_tile);if(ierr/=MPI_SUCCESS)status=4
      enddo
    enddo
    workspace_peak_bytes=peak
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    ok=status==0.and.ierr==MPI_SUCCESS
    if(.not.ok)message='Wannier90 distributed matrix reduction failed'
#else
    ok=.false.;message='Wannier90 matrix assembly requires MPI';coordinator_bytes=0_int64
    workspace_peak_bytes=0_int64;allocate(m_matrix(0,0,0),a_matrix(0,0))
#endif
  end subroutine assemble_dg_w90_gamma_matrices

  subroutine setup_dg_w90_gamma_library(comm,seed,real_lattice,reciprocal_lattice,atom_symbols,&
      atoms_cart,nband,nwann,nntot,nncell,ok,message)
    integer,intent(in)::comm,nband,nwann
    character(*),intent(in)::seed
    real(real64),intent(in)::real_lattice(3,3),reciprocal_lattice(3,3),atoms_cart(:,:)
    character(*),intent(in)::atom_symbols(:)
    integer,intent(out)::nntot
    integer,allocatable,intent(out)::nncell(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#if defined(USE_MPI) && defined(USE_WANNIER90)
    integer,parameter::num_nnmax=12
    integer::rank,ierr,unit,io,axis,atom,num_bands_out,num_wann_out,status
    logical::dmn_exists
    integer::mp_grid(3),nnlist(1,num_nnmax),nncell_max(3,1,num_nnmax),exclude_bands(max(1,nband))
    integer::proj_l(max(1,nband)),proj_m(max(1,nband)),proj_radial(max(1,nband))
    integer::proj_s(max(1,nband))
    real(real64)::kpoint(3,1),proj_site(3,max(1,nband)),proj_z(3,max(1,nband)),&
      proj_x(3,max(1,nband)),proj_zona(max(1,nband)),proj_s_qaxis(3,max(1,nband))
    real(real64)::real_lattice_w90(3,3),reciprocal_lattice_w90(3,3),&
      atoms_cart_w90(3,size(atom_symbols))
    logical::geometry_ok
    character(len(message))::geometry_message
    interface
      subroutine wannier_setup(seed_name,mp_grid_loc,num_kpts_loc,real_lattice_loc,&
          recip_lattice_loc,kpt_latt_loc,num_bands_tot,num_atoms_loc,atom_symbols_loc,&
          atoms_cart_loc,gamma_only_loc,spinors_loc,nntot_loc,nnlist_loc,nncell_loc,&
          num_bands_loc,num_wann_loc,proj_site_loc,proj_l_loc,proj_m_loc,proj_radial_loc,&
          proj_z_loc,proj_x_loc,proj_zona_loc,exclude_bands_loc,proj_s_loc,proj_s_qaxis_loc)
        import real64,num_nnmax
        character(*),intent(in)::seed_name
        integer,intent(in)::mp_grid_loc(3),num_kpts_loc,num_bands_tot,num_atoms_loc
        real(real64),intent(in)::real_lattice_loc(3,3),recip_lattice_loc(3,3),kpt_latt_loc(3,num_kpts_loc)
        character(*),intent(in)::atom_symbols_loc(num_atoms_loc)
        real(real64),intent(in)::atoms_cart_loc(3,num_atoms_loc)
        logical,intent(in)::gamma_only_loc,spinors_loc
        integer,intent(out)::nntot_loc,nnlist_loc(num_kpts_loc,num_nnmax),&
          nncell_loc(3,num_kpts_loc,num_nnmax),num_bands_loc,num_wann_loc
        real(real64),intent(out)::proj_site_loc(3,num_bands_tot),proj_z_loc(3,num_bands_tot),&
          proj_x_loc(3,num_bands_tot),proj_zona_loc(num_bands_tot)
        integer,intent(out)::proj_l_loc(num_bands_tot),proj_m_loc(num_bands_tot),&
          proj_radial_loc(num_bands_tot),exclude_bands_loc(num_bands_tot),proj_s_loc(num_bands_tot)
        real(real64),intent(out)::proj_s_qaxis_loc(3,num_bands_tot)
      end subroutine wannier_setup
    end interface
    ok=.false.;message='';nntot=0;status=0
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS.or.nband<=0.or.nwann/=nband.or.size(atom_symbols)<=0.or.&
        any(shape(atoms_cart)/=[3,size(atom_symbols)]).or..not.all(ieee_is_finite(real_lattice)).or.&
        .not.all(ieee_is_finite(reciprocal_lattice)).or..not.all(ieee_is_finite(atoms_cart)))status=1
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then;message='invalid Gamma Wannier90 setup contract';return;endif
    mp_grid=[1,1,1];kpoint=0d0
    if(rank==0)then
      inquire(file=trim(seed)//'.dmn',exist=dmn_exists)
      if(.not.dmn_exists)status=4
    endif
    if(rank==0.and.status==0)then
      call convert_dg_w90_library_geometry(real_lattice,reciprocal_lattice,atoms_cart,&
        real_lattice_w90,reciprocal_lattice_w90,atoms_cart_w90,geometry_ok,geometry_message)
      if(.not.geometry_ok)status=5
    endif
    if(rank==0.and.status==0)then
      open(newunit=unit,file=trim(seed)//'.win',status='replace',action='write',iostat=io)
      if(io/=0)then
        status=2
      else
        write(unit,'(a,i0)')'num_bands = ',nband
        write(unit,'(a,i0)')'num_wann = ',nwann
        write(unit,'(a)')'num_iter = 400'
        write(unit,'(a)')'conv_tol = 1.d-10'
        write(unit,'(a)')'conv_window = 5'
        write(unit,'(a)')'gamma_only = true'
        write(unit,'(a)')'site_symmetry = .true.'
        write(unit,'(a)')'symmetrize_eps = 1.d-10'
        write(unit,'(a)')'begin unit_cell_cart';write(unit,'(a)')'bohr'
        do axis=1,3;write(unit,'(3(es24.16,1x))')real_lattice(:,axis);enddo
        write(unit,'(a)')'end unit_cell_cart'
        write(unit,'(a)')'begin atoms_cart';write(unit,'(a)')'bohr'
        do atom=1,size(atom_symbols)
          write(unit,'(a,1x,3(es24.16,1x))')trim(atom_symbols(atom)),atoms_cart(:,atom)
        enddo
        write(unit,'(a)')'end atoms_cart'
        ! Library mode receives the deterministic spectral trial overlap as A_matrix_loc.
        ! Omitting the projections block prevents setup from generating an unrelated
        ! random trial gauge; Wannier90 permits this and initializes num_proj=num_wann.
        write(unit,'(a)')'mp_grid = 1 1 1'
        write(unit,'(a)')'begin kpoints';write(unit,'(a)')'0.0 0.0 0.0'
        write(unit,'(a)')'end kpoints';close(unit)
        call wannier_setup(trim(seed),mp_grid,1,real_lattice_w90,reciprocal_lattice_w90,kpoint,nband,&
          size(atom_symbols),atom_symbols,atoms_cart_w90,.true.,.false.,nntot,nnlist,nncell_max,&
          num_bands_out,num_wann_out,proj_site,proj_l,proj_m,proj_radial,proj_z,proj_x,&
          proj_zona,exclude_bands,proj_s,proj_s_qaxis)
        if(nntot<1.or.nntot>num_nnmax.or.num_bands_out/=nband.or.num_wann_out/=nwann)status=3
      endif
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(nntot,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(nncell_max,3*num_nnmax,MPI_INTEGER,0,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then;message='Wannier90 Gamma library setup failed';return;endif
    allocate(nncell(3,nntot));nncell=nncell_max(:,1,1:nntot);ok=.true.
#else
    ok=.false.;message='Wannier90 Gamma setup requires MPI and USE_WANNIER90';nntot=0
#endif
  end subroutine setup_dg_w90_gamma_library

  subroutine run_dg_w90_gamma_library(comm,seed,real_lattice,reciprocal_lattice,atom_symbols,&
      atoms_cart,m_matrix,a_matrix,eigenvalues,initial_gauge_spread,tolerance,transform,centers,&
      spreads,spread,ok,message,convergence_iterations_out)
    integer,intent(in)::comm
    character(*),intent(in)::seed
    real(real64),intent(in)::real_lattice(3,3),reciprocal_lattice(3,3),atoms_cart(:,:),&
      eigenvalues(:),initial_gauge_spread,tolerance
    character(*),intent(in)::atom_symbols(:)
    complex(real64),intent(in)::m_matrix(:,:,:),a_matrix(:,:)
    complex(real64),allocatable,intent(out)::transform(:,:)
    real(real64),allocatable,intent(out)::centers(:,:),spreads(:)
    real(real64),intent(out)::spread(3)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,intent(out),optional::convergence_iterations_out
#if defined(USE_MPI) && defined(USE_WANNIER90)
    integer::rank,ierr,nband,nwann,nntot,status,mp_grid(3),matrix_dimensions(3),convergence_iterations
    real(real64)::kpoint(3,1)
    real(real64),parameter::bohr_to_angstrom=0.52917721067_real64
    real(real64)::real_lattice_w90(3,3),reciprocal_lattice_w90(3,3),&
      atoms_cart_w90(3,size(atom_symbols))
    complex(real64),allocatable::u(:,:,:),uopt(:,:,:),m4(:,:,:,:),a3(:,:,:)
    real(real64),allocatable::e2(:,:)
    logical,allocatable::lwindow(:,:)
    logical::geometry_ok
    character(len(message))::geometry_message
    interface
      subroutine wannier_run(seed_name,mp_grid_loc,num_kpts_loc,real_lattice_loc,&
          recip_lattice_loc,kpt_latt_loc,num_bands_loc,num_wann_loc,nntot_loc,num_atoms_loc,&
          atom_symbols_loc,atoms_cart_loc,gamma_only_loc,m_matrix_loc,a_matrix_loc,&
          eigenvalues_loc,u_matrix_loc,u_matrix_opt_loc,lwindow_loc,wann_centres_loc,&
          wann_spreads_loc,spread_loc)
        import real64
        character(*),intent(in)::seed_name
        integer,intent(in)::mp_grid_loc(3),num_kpts_loc,num_bands_loc,num_wann_loc,nntot_loc,num_atoms_loc
        real(real64),intent(in)::real_lattice_loc(3,3),recip_lattice_loc(3,3),kpt_latt_loc(3,num_kpts_loc)
        character(*),intent(in)::atom_symbols_loc(num_atoms_loc)
        real(real64),intent(in)::atoms_cart_loc(3,num_atoms_loc)
        logical,intent(in)::gamma_only_loc
        complex(real64),intent(in)::m_matrix_loc(num_bands_loc,num_bands_loc,nntot_loc,num_kpts_loc),&
          a_matrix_loc(num_bands_loc,num_wann_loc,num_kpts_loc)
        real(real64),intent(in)::eigenvalues_loc(num_bands_loc,num_kpts_loc)
        complex(real64),intent(out)::u_matrix_loc(num_wann_loc,num_wann_loc,num_kpts_loc),&
          u_matrix_opt_loc(num_bands_loc,num_wann_loc,num_kpts_loc)
        logical,intent(out)::lwindow_loc(num_bands_loc,num_kpts_loc)
        real(real64),intent(out)::wann_centres_loc(3,num_wann_loc),wann_spreads_loc(num_wann_loc),spread_loc(3)
      end subroutine wannier_run
    end interface
    ok=.false.;message='';spread=0d0;status=0;convergence_iterations=-1
    if(present(convergence_iterations_out))convergence_iterations_out=-1
    call MPI_Comm_rank(comm,rank,ierr)
    matrix_dimensions=0
    if(rank==0)matrix_dimensions=[size(m_matrix,1),size(a_matrix,2),size(m_matrix,3)]
    call MPI_Bcast(matrix_dimensions,3,MPI_INTEGER,0,comm,ierr)
    nband=matrix_dimensions(1);nwann=matrix_dimensions(2);nntot=matrix_dimensions(3)
    if(ierr/=MPI_SUCCESS.or.nband<=0.or.nwann/=nband.or.size(eigenvalues)/=nband.or.nntot<=0.or.&
        any(shape(atoms_cart)/=[3,size(atom_symbols)]).or.&
        .not.all(ieee_is_finite(eigenvalues)).or..not.all(ieee_is_finite(real_lattice)).or.&
        .not.all(ieee_is_finite(reciprocal_lattice)).or..not.all(ieee_is_finite(atoms_cart)))status=1
    if(rank==0)then
      if(size(m_matrix,2)/=nband.or.size(a_matrix,1)/=nband.or.&
          .not.all(ieee_is_finite(real(m_matrix))).or..not.all(ieee_is_finite(aimag(m_matrix))).or.&
          .not.all(ieee_is_finite(real(a_matrix))).or..not.all(ieee_is_finite(aimag(a_matrix))))status=1
    else if(size(m_matrix)/=0.or.size(a_matrix)/=0)then
      status=1
    endif
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then;message='invalid Gamma Wannier90 run contract';return;endif
    allocate(transform(nwann,nwann),centers(3,nwann),spreads(nwann));transform=(0d0,0d0)
    centers=0d0;spreads=0d0;mp_grid=[1,1,1];kpoint=0d0
    if(rank==0)then
      call convert_dg_w90_library_geometry(real_lattice,reciprocal_lattice,atoms_cart,&
        real_lattice_w90,reciprocal_lattice_w90,atoms_cart_w90,geometry_ok,geometry_message)
      if(.not.geometry_ok)then
        message=geometry_message;status=3
      endif
    endif
    if(rank==0.and.status==0)then
      allocate(u(nwann,nwann,1),uopt(nband,nwann,1),lwindow(nband,1),&
        m4(nband,nband,nntot,1),a3(nband,nwann,1),e2(nband,1))
      m4(:,:,:,1)=m_matrix;a3(:,:,1)=a_matrix;e2(:,1)=eigenvalues
      call wannier_run(trim(seed),mp_grid,1,real_lattice_w90,reciprocal_lattice_w90,kpoint,nband,nwann,&
        nntot,size(atom_symbols),atom_symbols,atoms_cart_w90,.true.,m4,a3,e2,u,uopt,lwindow,&
        centers,spreads,spread)
      centers=centers/bohr_to_angstrom
      spreads=spreads/(bohr_to_angstrom*bohr_to_angstrom)
      spread=spread/(bohr_to_angstrom*bohr_to_angstrom)
      transform=matmul(uopt(:,:,1),u(:,:,1))
      call validate_dg_w90_convergence_log(trim(seed)//'.wout',400,convergence_iterations,ok,message)
      if(ok)call validate_dg_w90_result(transform,centers,spreads,spread,initial_gauge_spread,&
        tolerance,ok,message)
      status=merge(0,2,ok)
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(convergence_iterations,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(message,len(message),MPI_CHARACTER,0,comm,ierr)
    call MPI_Bcast(transform,size(transform),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    call MPI_Bcast(centers,size(centers),MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(spreads,size(spreads),MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(spread,3,MPI_DOUBLE_PRECISION,0,comm,ierr)
    ok=status==0.and.ierr==MPI_SUCCESS
    if(present(convergence_iterations_out))convergence_iterations_out=convergence_iterations
    if(ok)message=''
#else
    ok=.false.;message='Wannier90 Gamma run requires MPI and USE_WANNIER90';spread=0d0
    if(present(convergence_iterations_out))convergence_iterations_out=-1
#endif
  end subroutine run_dg_w90_gamma_library

  subroutine estimate_dg_w90_coordinator_bytes(nband,nwann,nntot,nkpoint,nbytes,ok,message)
    integer,intent(in)::nband,nwann,nntot,nkpoint
    integer(int64),intent(out)::nbytes
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64)::nb,nw,nn,nk,complex_elements,real_elements,term
    ok=.false.;message='';nbytes=0_int64
    if(nband<=0.or.nwann<=0.or.nwann>nband.or.nntot<=0.or.nkpoint<=0)then
      message='invalid Wannier90 coordinator dimensions';return
    endif
    nb=int(nband,int64);nw=int(nwann,int64);nn=int(nntot,int64);nk=int(nkpoint,int64)
    complex_elements=0_int64;real_elements=0_int64
    ! Input plus the library-owned original/projected overlap copies.
    call checked_product([2_int64,nb,nb,nn,nk],term,ok)
    if(.not.ok)then;message='Wannier90 M-matrix byte estimate overflow';return;endif
    call checked_add(complex_elements,term,ok);if(.not.ok)goto 900
    ! A input, optimized subspace, returned U/Uopt, and SALMON scatter copy.
    call checked_product([3_int64,nb,nw,nk],term,ok)
    if(.not.ok)then;message='Wannier90 A/Uopt byte estimate overflow';return;endif
    call checked_add(complex_elements,term,ok);if(.not.ok)goto 900
    call checked_product([2_int64,nw,nw,nk],term,ok)
    if(.not.ok)then;message='Wannier90 U-matrix byte estimate overflow';return;endif
    call checked_add(complex_elements,term,ok);if(.not.ok)goto 900
    call checked_product([nb,nk],term,ok);if(.not.ok)goto 900
    call checked_add(real_elements,term,ok);if(.not.ok)goto 900
    call checked_product([4_int64,nw],term,ok);if(.not.ok)goto 900
    call checked_add(real_elements,term,ok);if(.not.ok)goto 900
    call checked_add(real_elements,3_int64,ok);if(.not.ok)goto 900
    call checked_product([complex_elements,int(storage_size((0d0,0d0))/8,int64)],term,ok)
    if(.not.ok)goto 900
    nbytes=term
    call checked_product([real_elements,int(storage_size(0d0)/8,int64)],term,ok)
    if(.not.ok)goto 900
    call checked_add(nbytes,term,ok);if(.not.ok)goto 900
    ok=.true.;return
900 message='Wannier90 coordinator byte estimate overflow';nbytes=0_int64
  end subroutine estimate_dg_w90_coordinator_bytes

  subroutine validate_dg_w90_result(transform,centers,spreads,spread,initial_gauge_spread,&
      tolerance,ok,message)
    complex(real64),intent(in)::transform(:,:)
    real(real64),intent(in)::centers(:,:),spreads(:),spread(:),initial_gauge_spread,tolerance
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::gram(:,:)
    real(real64)::scale,defect,imaginary_defect
    integer::i,nwann
    ok=.false.;message='';nwann=size(transform,1)
    if(nwann<=0.or.size(transform,2)/=nwann.or.any(shape(centers)/=[3,nwann]).or.&
        size(spreads)/=nwann.or.size(spread)/=3.or.tolerance<=0d0.or.&
        .not.ieee_is_finite(tolerance).or..not.ieee_is_finite(initial_gauge_spread).or.&
        initial_gauge_spread<0d0)then
      message='invalid Wannier90 result dimensions or tolerance';return
    endif
    if(.not.all(ieee_is_finite(real(transform))).or..not.all(ieee_is_finite(aimag(transform))).or.&
        .not.all(ieee_is_finite(centers)).or..not.all(ieee_is_finite(spreads)).or.&
        .not.all(ieee_is_finite(spread)))then
      message='Wannier90 result is nonfinite';return
    endif
    scale=max(1d0,max(maxval(abs(spreads)),maxval(abs(spread))))
    if(any(spreads < -tolerance*scale).or.any(spread < -tolerance*scale))then
      message='Wannier90 result has a physically negative spread';return
    endif
    scale=max(1d0,maxval(abs(transform)))
    imaginary_defect=maxval(abs(aimag(transform)))
    if(imaginary_defect>tolerance*scale)then
      message='Wannier90 transform violates the Gamma-real gauge';return
    endif
    allocate(gram(nwann,nwann));gram=matmul(conjg(transpose(transform)),transform)
    do i=1,nwann;gram(i,i)=gram(i,i)-1d0;enddo
    defect=maxval(abs(gram))
    if(defect>tolerance*max(1d0,real(nwann,real64)))then
      message='Wannier90 transform is not unitary';return
    endif
    if(spread(3)>initial_gauge_spread+tolerance*max(1d0,initial_gauge_spread))then
      message='Wannier90 increased the gauge-dependent spread';return
    endif
    ok=.true.
  end subroutine validate_dg_w90_result

  subroutine checked_product(factors,value,ok)
    integer(int64),intent(in)::factors(:)
    integer(int64),intent(out)::value
    logical,intent(out)::ok
    integer::i
    value=1_int64;ok=all(factors>=0_int64)
    if(.not.ok)return
    do i=1,size(factors)
      if(factors(i)==0_int64)then;value=0_int64;return;endif
      if(value>huge(value)/factors(i))then;value=0_int64;ok=.false.;return;endif
      value=value*factors(i)
    enddo
  end subroutine checked_product

  subroutine checked_add(value,increment,ok)
    integer(int64),intent(inout)::value
    integer(int64),intent(in)::increment
    logical,intent(inout)::ok
    if(.not.ok)return
    ok=value>=0_int64.and.increment>=0_int64.and.value<=huge(value)-increment
    if(ok)value=value+increment
  end subroutine checked_add
end module dg_overlapping_wannier_w90
