program test_dg_overlapping_wannier_eigenexa_mpi
  use mpi
  use structures,only:s_parallel_info
  use eigen_libs_mod
  use eigen_eigenexa,only:eigen_pdsyevd_ex_distributed_blocks
  use dg_overlapping_wannier_solver,only:solve_dg_overlapping_wannier_generalized_eigenexa
  use dg_overlapping_wannier_construction,only:build_dg_group_averaged_occupied_candidates_eigenexa,&
    build_dg_cocycle_averaged_occupied_candidates_eigenexa,measure_dg_rank_fixed_symmetry_residuals,&
    split_dg_translation_character_sector_eigenexa,validate_dg_translation_sector_cluster,&
    diagonalize_dg_spectral_basin_operator,select_dg_spectral_basin_channel_ranks
  use dg_overlapping_wannier_construction,only:propagate_dg_spectral_basin_orbit_channels
  use dg_overlapping_wannier_construction,only:build_dg_spectral_channel_generator_actions
  use dg_overlapping_wannier_construction,only:compose_dg_occupied_complement_trial_rows
  implicit none
  type(s_parallel_info)::info
  integer::comm,rank,nproc,ierr,i,p,nlocal
  integer(8),allocatable::row_ids(:)
  complex(8),allocatable::hrows(:,:),srows(:,:),coeff(:,:)
  real(8)::metric_diagonal(4),target_eigenvalues(4),residual,orthogonality,condition,gamma_defect
  real(8),allocatable::eigenvalues(:)
  integer(8)::workspace,signature
  logical::ok
  character(256)::message
  character(32)::case_name
  real(8)::solve_tolerance

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  case_name='normal';if(command_argument_count()>=1)call get_command_argument(1,case_name)
  if(index(trim(case_name),'sector')==1)then;call run_sector_case();call MPI_Finalize(ierr);stop;endif
  if(index(trim(case_name),'average')==1)then;call run_average_case();call MPI_Finalize(ierr);stop;endif
  if(trim(case_name)=='cocycle')then;call run_cocycle_case();call MPI_Finalize(ierr);stop;endif
  if(index(trim(case_name),'spectral_basin')==1)then
    call run_spectral_basin_case();call MPI_Finalize(ierr);stop
  endif
  call eigen_init(comm);call eigen_get_procs(p,info%nprow,info%npcol)
  call eigen_get_id(p,info%myrow,info%mycol);call eigen_get_matdims(4,info%nrow_local,info%ncol_local)
  info%flag_eigenexa_init=.true.
  solve_tolerance=1d-10
  nlocal=count([(mod(i-1,nproc)==rank,i=1,4)])
  allocate(row_ids(nlocal),hrows(nlocal,4),srows(nlocal,4),coeff(4,3),eigenvalues(3))
  metric_diagonal=[2d0,1.5d0,1.2d0,0.9d0];target_eigenvalues=[0.2d0,0.6d0,1.1d0,2d0]
  hrows=(0d0,0d0);srows=(0d0,0d0);p=0
  do i=1,4
    if(mod(i-1,nproc)/=rank)cycle
    p=p+1;row_ids(p)=i;srows(p,i)=metric_diagonal(i)
    hrows(p,i)=metric_diagonal(i)*target_eigenvalues(i)
  enddo
  select case(trim(case_name))
  case('degenerate')
    target_eigenvalues=[0d0,1d0,2d0,2d0]
    solve_tolerance=1d-6
    do p=1,nlocal;hrows(p,:)=0d0;hrows(p,int(row_ids(p)))=&
      metric_diagonal(int(row_ids(p)))*target_eigenvalues(int(row_ids(p)));enddo
  case('nonreal')
    if(nlocal>0.and.row_ids(1)==1_8)hrows(1,1)=hrows(1,1)+cmplx(0d0,1d-4,8)
  case('illmetric')
    do p=1,nlocal;if(row_ids(p)==4_8)srows(p,4)=1d-16;enddo
  case('residual')
    solve_tolerance=1d-30
    do p=1,nlocal
      if(row_ids(p)==1_8)hrows(p,2)=0.123456789d0
      if(row_ids(p)==2_8)hrows(p,1)=0.123456789d0
    enddo
  end select
  call solve_dg_overlapping_wannier_generalized_eigenexa(info,comm,row_ids,hrows,srows,3,&
    solve_tolerance,1d-12,1d-12,coeff,eigenvalues,residual,orthogonality,condition,gamma_defect,&
    workspace,ok,message)
  if(trim(case_name)/='normal')then
    call require(.not.ok,'negative generalized EigenExa case must reject')
    if(trim(case_name)=='degenerate')call require(index(message,'gap=')>0.and.&
      index(message,'threshold=')>0,'degenerate boundary rejection reports spectrum and threshold')
    if(rank==0)write(*,'(3a,i0)')'REJECT ',trim(case_name),' ranks=',nproc
    call eigen_free();call MPI_Finalize(ierr);stop
  endif
  call require(ok,trim(message));call require(maxval(abs(eigenvalues-target_eigenvalues(1:3)))<1d-10,&
    'known generalized EigenExa spectrum')
  call require(residual<1d-10.and.orthogonality<1d-10,'generalized EigenExa quality receipts')
  call require(workspace>0_8.and.gamma_defect==0d0,'measured workspace and Gamma-real receipts')
  signature=nint(sum(eigenvalues*[1d0,3d0,7d0])*1d12,8)
  if(rank==0)write(*,'(a,i0,a,i0)')'EIGENEXA ranks=',nproc,' signature=',signature
  call eigen_free();call MPI_Finalize(ierr)
contains
  subroutine run_spectral_basin_case()
    complex(8)::operator(4,4),rotation(4,4),rotated_operator(4,4)
    real(8),allocatable::spectrum(:),rotated_spectrum(:)
    integer,allocatable::block_offsets(:),rotated_offsets(:)
    integer::orbit_map(2,1),selected_ranks(2),payload_collectives
    integer::ii,jj,local_row
    real(8)::catalog_spectra(4,2)
    logical::block_ends(4,2)
    integer(8),allocatable::propagation_row_ids(:)
    complex(8),allocatable::propagation_generators(:,:,:),representative_vectors(:,:),trial_rows(:,:)
    complex(8),allocatable::target_action_rows(:,:,:),full_trial_rows(:,:),complement_rows(:,:)
    integer(8),allocatable::full_row_ids(:),complement_row_ids(:)
    real(8)::trial_gram_defect,trial_frame_defect,target_action_unitarity,target_action_block_defect
    real(8)::residual,rotated_residual,angle
    integer(8)::fingerprint,rotated_fingerprint,workspace
    logical::basin_ok
    character(256)::basin_message

    operator=(0d0,0d0);operator(1,1)=1d0;operator(2,2)=0.6d0
    operator(3,3)=0.6d0;operator(4,4)=0d0
    angle=0.37d0;rotation=(0d0,0d0)
    rotation(1,1)=cos(angle);rotation(1,2)=sin(angle)
    rotation(2,1)=-sin(angle);rotation(2,2)=cos(angle)
    rotation(3,3)=1d0;rotation(4,4)=1d0
    rotated_operator=matmul(conjg(transpose(rotation)),matmul(operator,rotation))
    if(trim(case_name)=='spectral_basin_nonhermitian')then
      operator(1,2)=cmplx(0.1d0,0.2d0,8)
      call diagonalize_dg_spectral_basin_operator(comm,operator,9901_8,1d-10,spectrum,&
        block_offsets,residual,fingerprint,workspace,basin_ok,basin_message)
      call require(.not.basin_ok,'non-Hermitian spectral basin operator is rejected')
      if(rank==0)write(*,'(a,i0)')'REJECT spectral_basin_nonhermitian ranks=',nproc
      return
    endif
    if(trim(case_name)=='spectral_basin_complement')then
      nlocal=count([(mod(ii-1,nproc)==rank,ii=1,6)])
      allocate(full_row_ids(nlocal));local_row=0
      do ii=1,6
        if(mod(ii-1,nproc)/=rank)cycle
        local_row=local_row+1;full_row_ids(local_row)=ii
      enddo
      allocate(complement_row_ids(count(full_row_ids>2_8)),&
        complement_rows(count(full_row_ids>2_8),4));local_row=0
      do ii=1,size(full_row_ids)
        if(full_row_ids(ii)<=2_8)cycle
        local_row=local_row+1;complement_row_ids(local_row)=full_row_ids(ii)-2_8
        complement_rows(local_row,:)=0d0
        complement_rows(local_row,int(complement_row_ids(local_row)))=1d0
      enddo
      call compose_dg_occupied_complement_trial_rows(comm,full_row_ids,2,complement_row_ids,&
        complement_rows,1d-10,full_trial_rows,trial_gram_defect,fingerprint,workspace,basin_ok,basin_message)
      call require(basin_ok.and.trial_gram_defect<1d-12.and.size(full_trial_rows,2)==6,&
        'occupied identity and localized complement compose one complete trial frame')
      trial_frame_defect=0d0
      do local_row=1,size(full_row_ids)
        ii=int(full_row_ids(local_row))
        do jj=1,6
          trial_frame_defect=max(trial_frame_defect,abs(full_trial_rows(local_row,jj)-&
            merge((1d0,0d0),(0d0,0d0),ii==jj)))
        enddo
      enddo
      call require(trial_frame_defect<1d-12,'occupied/complement composition preserves canonical rows')
      allocate(propagation_generators(nlocal,6,1));propagation_generators=0d0
      do local_row=1,nlocal
        ii=int(full_row_ids(local_row))
        select case(ii)
        case(1,2);propagation_generators(local_row,ii,1)=1d0
        case(3,4);propagation_generators(local_row,ii+2,1)=1d0
        case(5,6);propagation_generators(local_row,ii-2,1)=1d0
        end select
      enddo
      orbit_map(:,1)=[2,1];selected_ranks=[2,2]
      call build_dg_spectral_channel_generator_actions(comm,full_row_ids,propagation_generators,&
        full_trial_rows,orbit_map,selected_ranks,8801_8,fingerprint,1d-10,target_action_rows,&
        target_action_unitarity,target_action_block_defect,rotated_fingerprint,workspace,basin_ok,&
        basin_message,preserved_prefix=2)
      call require(basin_ok.and.target_action_unitarity<1d-12.and.target_action_block_defect<1d-12,&
        'target action preserves occupied prefix and permutes only complement basin blocks')
      if(rank==0)write(*,'(a,i0,a,i0)')'SPECTRAL_COMPLEMENT ranks=',nproc,' signature=',fingerprint
      return
    endif
    call diagonalize_dg_spectral_basin_operator(comm,operator,9901_8,1d-10,spectrum,&
      block_offsets,residual,fingerprint,workspace,basin_ok,basin_message)
    call require(basin_ok,trim(basin_message))
    call require(maxval(abs(spectrum-[1d0,0.6d0,0.6d0,0d0]))<1d-12.and.&
      all(block_offsets==[1,2,4,5]).and.residual<1d-12,&
      'spectral basin eigensystem preserves complete unresolved blocks')
    call diagonalize_dg_spectral_basin_operator(comm,rotated_operator,9901_8,1d-10,rotated_spectrum,&
      rotated_offsets,rotated_residual,rotated_fingerprint,workspace,basin_ok,basin_message)
    call require(basin_ok.and.maxval(abs(rotated_spectrum-spectrum))<1d-12.and.&
      all(rotated_offsets==block_offsets).and.rotated_fingerprint==fingerprint,&
      'spectral basin blocks are retained-frame gauge covariant')
    catalog_spectra(:,1)=[1d0,0.6d0,0d0,0d0];catalog_spectra(:,2)=catalog_spectra(:,1)
    block_ends=.false.;block_ends(1,:)=.true.;block_ends(2,:)=.true.;block_ends(4,:)=.true.
    orbit_map(:,1)=[2,1]
    call select_dg_spectral_basin_channel_ranks(comm,catalog_spectra,block_ends,orbit_map,4,1d-10,&
      selected_ranks,fingerprint,workspace,basin_ok,basin_message,payload_collectives)
    call require(basin_ok.and.all(selected_ranks==[2,2]).and.payload_collectives==6,&
      'spectral basin catalog spans retained rank with equal orbit ranks')
    nlocal=count([(mod(ii-1,nproc)==rank,ii=1,4)])
    allocate(propagation_row_ids(nlocal),propagation_generators(nlocal,4,1),representative_vectors(4,2))
    propagation_generators=(0d0,0d0);representative_vectors=(0d0,0d0)
    representative_vectors(1,1)=1d0;representative_vectors(2,2)=1d0;local_row=0
    do ii=1,4
      if(mod(ii-1,nproc)/=rank)cycle
      local_row=local_row+1;propagation_row_ids(local_row)=ii
      select case(ii)
      case(1);propagation_generators(local_row,3,1)=1d0
      case(2);propagation_generators(local_row,4,1)=1d0
      case(3);propagation_generators(local_row,1,1)=1d0
      case(4);propagation_generators(local_row,2,1)=1d0
      end select
    enddo
    call propagate_dg_spectral_basin_orbit_channels(comm,propagation_row_ids,propagation_generators,&
      orbit_map,selected_ranks,representative_vectors,8801_8,fingerprint,1d-10,trial_rows,&
      trial_gram_defect,rotated_fingerprint,workspace,basin_ok,basin_message)
    call require(basin_ok.and.trial_gram_defect<1d-12,&
      'one representative eigenspace propagates to a complete row-owned trial frame')
    trial_frame_defect=0d0
    do local_row=1,nlocal
      ii=int(propagation_row_ids(local_row))
      trial_frame_defect=max(trial_frame_defect,&
        maxval(abs(trial_rows(local_row,:)-[(merge((1d0,0d0),(0d0,0d0),jj==ii),jj=1,4)])))
    enddo
    call require(trial_frame_defect<1d-12,&
      'propagated spectral basin trial frame has canonical basin-column order')
    call build_dg_spectral_channel_generator_actions(comm,propagation_row_ids,propagation_generators,&
      trial_rows,orbit_map,selected_ranks,8801_8,rotated_fingerprint,1d-10,target_action_rows,&
      target_action_unitarity,target_action_block_defect,fingerprint,workspace,basin_ok,basin_message)
    call require(basin_ok.and.target_action_unitarity<1d-12.and.target_action_block_defect<1d-12.and.&
      maxval(abs(target_action_rows-propagation_generators))<1d-12,&
      'streamed target action has the known basin permutation and no off-block leakage')
    if(trim(case_name)=='spectral_basin_split')then
      catalog_spectra(:,1)=[1d0,0.6d0,0.6d0,0d0];catalog_spectra(:,2)=catalog_spectra(:,1)
      block_ends=.false.;block_ends(1,:)=.true.;block_ends(3,:)=.true.;block_ends(4,:)=.true.
      call select_dg_spectral_basin_channel_ranks(comm,catalog_spectra,block_ends,orbit_map,4,1d-10,&
        selected_ranks,fingerprint,workspace,basin_ok,basin_message)
      call require(.not.basin_ok,'spectral basin catalog never splits an unresolved local block')
      if(rank==0)write(*,'(a,i0)')'REJECT spectral_basin_split ranks=',nproc
      return
    endif
    fingerprint=ieor(fingerprint,int(size(block_offsets),8))
    if(rank==0)write(*,'(a,i0,a,i0)')'SPECTRAL_BASIN_EIGEN ranks=',nproc,' signature=',fingerprint
  end subroutine run_spectral_basin_case

  subroutine run_sector_case()
    use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
    integer,parameter::n=8,ncharacter=4,ngenerator=2,multiplicity=2
    integer::ii,jj,kk,local_row,sector,sector_rank,generator_orders(ngenerator)
    integer(8),allocatable::sector_row_ids(:)
    complex(8),allocatable::generator_rows(:,:,:),gamma_rows(:,:),sector_rows(:,:),full_sector(:,:),&
      gauge(:,:),dense_generator(:,:,:),gamma_sewing(:,:),trivial_generators(:,:,:),&
      trivial_gamma(:,:),trivial_characters(:,:)
    complex(8)::characters(ncharacter,ngenerator),phase
    real(8)::identity_defect,unitarity_defect,commutator_defect,order_defect,gamma_pairing_defect
    real(8)::projector_defect,local_projector_defect
    real(8)::cluster_fixture(6)
    real(8),allocatable::cluster_matrix(:,:),cluster_vectors(:,:)
    integer::character_conjugates(ncharacter),element_words(ncharacter,ngenerator)
    integer,allocatable::trivial_orders(:),trivial_words(:,:),trivial_conjugates(:)
    integer(8)::sector_fingerprint,sector_workspace,sector_signature
    logical::sector_ok
    character(256)::sector_message

    nlocal=count([(mod(ii-1,nproc)==rank,ii=1,n)])
    allocate(sector_row_ids(nlocal),generator_rows(nlocal,n,ngenerator),gamma_rows(nlocal,n),&
      gauge(n,n),dense_generator(n,n,ngenerator),gamma_sewing(n,n))
    characters(:,1)=[(1d0,0d0),(1d0,0d0),(-1d0,0d0),(-1d0,0d0)]
    characters(:,2)=[(1d0,0d0),(-1d0,0d0),(1d0,0d0),(-1d0,0d0)]
    generator_orders=2;character_conjugates=[1,2,3,4];gauge=(0d0,0d0)
    element_words=reshape([0,0,1,1,0,1,0,1],[ncharacter,ngenerator])
    if(trim(case_name)=='sector_split_cluster')then
      cluster_fixture=[0d0,0d0,0d0,0d0,5d-11,1d0]
      call eigen_init(comm);call eigen_get_procs(p,info%nprow,info%npcol)
      call eigen_get_id(p,info%myrow,info%mycol);call eigen_get_matdims(6,info%nrow_local,info%ncol_local)
      info%flag_eigenexa_init=.true.
      allocate(cluster_matrix(info%nrow_local,info%ncol_local),&
        cluster_vectors(info%nrow_local,info%ncol_local));cluster_matrix=0d0
      do ii=1,6
        if(eigen_owner_node(ii,info%nprow,info%myrow)==info%myrow.and.&
            eigen_owner_node(ii,info%npcol,info%mycol)==info%mycol)&
          cluster_matrix(eigen_translate_g2l(ii,info%nprow,info%myrow),&
            eigen_translate_g2l(ii,info%npcol,info%mycol))=cluster_fixture(ii)
      enddo
      call eigen_pdsyevd_ex_distributed_blocks(info,6,cluster_matrix,cluster_fixture,cluster_vectors,&
        sector_ok,sector_message)
      call require(sector_ok,'split-cluster distributed EigenExa solve')
      call validate_dg_translation_sector_cluster(cluster_fixture,4,1d-10,sector_ok,sector_message)
      call require(.not.sector_ok.and.index(trim(sector_message),'boundary splits')>0,&
        'split-cluster fixture reaches the spectral boundary gate')
      if(rank==0)write(*,'(a,i0)')'REJECT sector_split_cluster ranks=',nproc
      call eigen_free()
      return
    endif
    do jj=1,n
      do ii=1,n
        phase=exp(cmplx(0d0,2d0*acos(-1d0)*real((ii-1)*(jj-1),8)/real(n,8),8))
        gauge(ii,jj)=phase/sqrt(real(n,8))*exp(cmplx(0d0,0.071d0*real(ii,8),8))
      enddo
    enddo
    gamma_sewing=matmul(gauge,transpose(gauge))
    dense_generator=(0d0,0d0)
    do kk=1,ngenerator
      do jj=1,n
        sector=(jj-1)/multiplicity+1
        dense_generator(:,:,kk)=dense_generator(:,:,kk)+characters(sector,kk)*&
          spread(gauge(:,jj),2,n)*spread(conjg(gauge(:,jj)),1,n)
      enddo
    enddo
    local_row=0
    do ii=1,n
      if(mod(ii-1,nproc)/=rank)cycle
      local_row=local_row+1;sector_row_ids(local_row)=ii
      generator_rows(local_row,:,:)=dense_generator(ii,:,:)
      gamma_rows(local_row,:)=gamma_sewing(ii,:)
    enddo
    select case(trim(case_name))
    case('sector_noncommuting')
      do local_row=1,nlocal
        ii=int(sector_row_ids(local_row))
        generator_rows(local_row,:,2)=generator_rows(local_row,:,2)+&
          (cos(0.4d0)-1d0)*gauge(ii,1)*conjg(gauge(:,1))+&
          (1d0-cos(0.4d0))*gauge(ii,7)*conjg(gauge(:,7))-&
          sin(0.4d0)*(gauge(ii,1)*conjg(gauge(:,7))+gauge(ii,7)*conjg(gauge(:,1)))
      enddo
    case('sector_nonunitary')
      generator_rows(:,:,1)=1.01d0*generator_rows(:,:,1)
    case('sector_wrong_order')
      generator_rows(:,:,1)=exp(cmplx(0d0,0.2d0,8))*generator_rows(:,:,1)
    case('sector_rank_losing')
      generator_rows(:,:,2)=generator_rows(:,:,1)
    case('sector_nonfinite')
      if(rank==0)generator_rows(1,1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
    case('sector_gamma_nonunitary')
      gamma_rows=1.01d0*gamma_rows
    case('sector_gamma_noninvolutory')
      gamma_rows=(0d0,0d0)
      do local_row=1,nlocal
        ii=int(sector_row_ids(local_row))
        if(ii==1)then;gamma_rows(local_row,1)=cos(0.2d0);gamma_rows(local_row,2)=-sin(0.2d0)
        elseif(ii==2)then;gamma_rows(local_row,1)=sin(0.2d0);gamma_rows(local_row,2)=cos(0.2d0)
        else;gamma_rows(local_row,ii)=1d0
        endif
      enddo
    case('sector_gamma_covariance')
      gamma_rows=(0d0,0d0)
      do local_row=1,nlocal;gamma_rows(local_row,int(sector_row_ids(local_row)))=1d0;enddo
    end select
    call eigen_init(comm);call eigen_get_procs(p,info%nprow,info%npcol)
    call eigen_get_id(p,info%myrow,info%mycol);call eigen_get_matdims(n,info%nrow_local,info%ncol_local)
    info%flag_eigenexa_init=.true.
    if(trim(case_name)=='sector_trivial')then
      allocate(trivial_generators(nlocal,n,0),trivial_gamma(nlocal,n),trivial_characters(1,0),&
        trivial_orders(0),trivial_words(1,0),trivial_conjugates(1))
      trivial_gamma=(0d0,0d0);trivial_conjugates=1
      do local_row=1,nlocal;trivial_gamma(local_row,int(sector_row_ids(local_row)))=1d0;enddo
      call split_dg_translation_character_sector_eigenexa(info,comm,sector_row_ids,trivial_generators,&
        trivial_gamma,trivial_characters,trivial_orders,trivial_words,trivial_conjugates,1,1d-10,81231_8,&
        sector_rows,sector_rank,identity_defect,unitarity_defect,commutator_defect,order_defect,&
        gamma_pairing_defect,sector_fingerprint,sector_workspace,sector_ok,sector_message)
      call require(sector_ok.and.sector_rank==n.and.all(shape(sector_rows)==[nlocal,n]).and.&
        gamma_pairing_defect<1d-12,'trivial translation group returns the full row-owned sector')
      if(rank==0)write(*,'(a,i0,a,i0)')'SECTOR_TRIVIAL ranks=',nproc,' signature=',sector_fingerprint
      call eigen_free();return
    endif
    call split_dg_translation_character_sector_eigenexa(info,comm,sector_row_ids,generator_rows,gamma_rows,&
      characters,generator_orders,element_words,character_conjugates,4,1d-10,77123_8,sector_rows,sector_rank,&
      identity_defect,unitarity_defect,&
      commutator_defect,order_defect,gamma_pairing_defect,sector_fingerprint,sector_workspace,&
      sector_ok,sector_message)
    if(trim(case_name)/='sector')then
      call require(.not.sector_ok,'adverse translation-sector generator case must reject')
      select case(trim(case_name))
      case('sector_noncommuting')
        call require(commutator_defect>1d-10.and.unitarity_defect<1d-10.and.order_defect<1d-10,&
          'noncommuting fixture reaches only the commutator gate')
      case('sector_nonunitary')
        call require(unitarity_defect>1d-10,'nonunitary fixture reaches the unitarity gate')
      case('sector_wrong_order')
        call require(order_defect>1d-10,'wrong-order fixture reaches the finite-order gate')
      case('sector_rank_losing')
        call require(index(trim(sector_message),'multiplicity')>0,&
          'rank-losing and split-cluster fixtures reach multiplicity gates')
      case('sector_nonfinite')
        call require(index(trim(sector_message),'invalid')>0,'nonfinite fixture reaches contract gate')
      case('sector_gamma_nonunitary','sector_gamma_noninvolutory','sector_gamma_covariance')
        call require(index(trim(sector_message),'Gamma')>0,'invalid Gamma sewing reaches a Gamma gate')
      end select
      if(rank==0)write(*,'(3a,i0)')'REJECT ',trim(case_name),' ranks=',nproc
      call eigen_free();return
    endif
    call require(sector_ok,trim(sector_message))
    call require(sector_rank==multiplicity.and.size(sector_rows,1)==nlocal.and.&
      size(sector_rows,2)==multiplicity,'known equal translation-character multiplicity')
    allocate(full_sector(n,multiplicity));full_sector=(0d0,0d0)
    do local_row=1,nlocal;full_sector(int(sector_row_ids(local_row)),:)=sector_rows(local_row,:);enddo
    call MPI_Allreduce(MPI_IN_PLACE,full_sector,n*multiplicity,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    local_projector_defect=0d0
    do ii=1,n
      do jj=1,n
        projector_defect=abs(sum(full_sector(ii,:)*conjg(full_sector(jj,:)))-&
          sum(gauge(ii,7:8)*conjg(gauge(jj,7:8))))
        local_projector_defect=max(local_projector_defect,projector_defect)
      enddo
    enddo
    call MPI_Allreduce(local_projector_defect,projector_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call require(projector_defect<1d-8.and.identity_defect<1d-10.and.unitarity_defect<1d-10.and.&
      commutator_defect<1d-10.and.order_defect<1d-10.and.gamma_pairing_defect<1d-10.and.&
      sector_workspace>0_8,'translation-sector projector and quality receipts')
    sector_signature=sector_fingerprint+int(sector_rank,8)
    if(rank==0)write(*,'(a,i0,a,i0)')'SECTOR ranks=',nproc,' signature=',sector_signature
    call eigen_free()
  end subroutine run_sector_case

  subroutine run_average_case()
    complex(8),allocatable::average_occupied(:,:),average_candidates(:,:)
    real(8),allocatable::average_spectrum(:)
    real(8),allocatable::average_weights(:)
    integer(8),allocatable::average_maps(:,:)
    integer::average_product(2,2),ii,global_point,global_count,average_rank,average_requested
    real(8)::average_trace,average_closure,average_gamma,average_selected_edge,&
      average_rejected_edge,average_cluster_gap
    integer(8)::average_workspace,average_signature
    logical::average_ok
    character(256)::average_message
    if(index(trim(case_name),'average_hamiltonian')==1)then
      call run_average_hamiltonian_case()
      return
    endif
    allocate(average_occupied(1,2),average_weights(2),average_maps(2,2))
    average_occupied=(0d0,0d0);average_weights=1d0;global_count=2*nproc
    do ii=1,2
      global_point=2*rank+ii
      average_maps(ii,1)=global_point
      average_maps(ii,2)=modulo(global_point-1+global_count/2,global_count)+1
      if(global_point==1)average_occupied(1,ii)=1d0
    enddo
    if(trim(case_name)=='average_unique')average_occupied=1d0/sqrt(real(global_count,8))
    if(trim(case_name)=='average_nonorthogonal')average_occupied=2d0/sqrt(real(global_count,8))
    average_product=reshape([1,2,2,1],[2,2])
    average_requested=merge(2,1,trim(case_name)=='average')
    call eigen_init(comm);call eigen_get_procs(p,info%nprow,info%npcol)
    call eigen_get_id(p,info%myrow,info%mycol);call eigen_get_matdims(2,info%nrow_local,info%ncol_local)
    info%flag_eigenexa_init=.true.
    call build_dg_group_averaged_occupied_candidates_eigenexa(info,comm,average_occupied,&
      average_weights,average_maps,average_product,1,average_requested,1d-12,average_candidates,average_spectrum,&
      average_rank,average_trace,average_closure,average_gamma,average_workspace,average_ok,average_message,&
      average_selected_edge,average_rejected_edge,average_cluster_gap)
    if(trim(case_name)=='average_split')then
      call require(.not.average_ok,'group-average selection must reject a split degenerate block')
      if(rank==0)write(*,'(a,i0)')'REJECT average_split ranks=',nproc
      call eigen_free();return
    endif
    if(trim(case_name)=='average_nonorthogonal')then
      call require(.not.average_ok,'group-average input occupied space must be metric orthonormal')
      if(rank==0)write(*,'(a,i0)')'REJECT average_nonorthogonal ranks=',nproc
      call eigen_free();return
    endif
    if(trim(case_name)=='average_unique')then
      call require(average_ok.and.average_rank==1.and.abs(average_spectrum(1)-1d0)<1d-12.and.&
        abs(average_trace-1d0)<1d-12.and.average_closure<1d-12.and.&
        abs(average_selected_edge-1d0)<1d-12.and.abs(average_rejected_edge)<1d-12.and.&
        abs(average_cluster_gap-1d0)<1d-12,&
        'unique invariant group-average rank is accepted')
      average_signature=nint(average_spectrum(1)*1d12,8)
      if(rank==0)write(*,'(a,i0,a,i0)')'AVERAGE_UNIQUE ranks=',nproc,' signature=',average_signature
      call eigen_free();return
    endif
    call require(average_ok,trim(average_message))
    call require(average_rank==2.and.maxval(abs(average_spectrum-[0.5d0,0.5d0]))<1d-12,&
      'distributed group-average spectrum')
    call require(abs(average_trace-1d0)<1d-12.and.average_closure<1d-12.and.&
      average_gamma==0d0.and.average_workspace>0_8,'distributed group-average receipts')
    average_signature=nint(sum(average_spectrum*[1d0,3d0])*1d12,8)
    if(rank==0)write(*,'(a,i0,a,i0)')'AVERAGE ranks=',nproc,' signature=',average_signature
    if(rank==0)average_maps(:,2)=1_8
    call build_dg_group_averaged_occupied_candidates_eigenexa(info,comm,average_occupied,&
      average_weights,average_maps,average_product,1,2,1d-12,average_candidates,average_spectrum,&
      average_rank,average_trace,average_closure,average_gamma,average_workspace,average_ok,average_message)
    call require(.not.average_ok,'distributed group-average rejects a non-group point action')
    call eigen_free()
  end subroutine run_average_case

  subroutine run_average_hamiltonian_case()
      use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
      complex(8)::occupied_hamiltonian(2,2)
      complex(8),allocatable::hamiltonian_occupied(:,:),hamiltonian_candidates(:,:)
      real(8),allocatable::hamiltonian_weights(:),hamiltonian_spectrum(:)
      integer(8),allocatable::hamiltonian_maps(:,:)
      real(8)::hamiltonian_trace,hamiltonian_closure,hamiltonian_gamma,&
        primary_selected,primary_rejected,primary_gap,secondary_selected,&
        secondary_rejected,secondary_gap,secondary_residual,local_defect,global_defect
      integer::hamiltonian_product(2,2),hamiltonian_rank,boundary_dimension,local_index,global_point
      integer(8)::hamiltonian_workspace,hamiltonian_signature,hamiltonian_fingerprint
      logical::average_ok
      character(256)::average_message

      allocate(hamiltonian_occupied(2,2),hamiltonian_weights(2),hamiltonian_maps(2,2))
      hamiltonian_occupied=(0d0,0d0);hamiltonian_weights=1d0
      do local_index=1,2
        global_point=2*rank+local_index
        hamiltonian_maps(local_index,:)=global_point
        if(global_point==1)hamiltonian_occupied(1,local_index)=1d0
        if(global_point==2)hamiltonian_occupied(2,local_index)=1d0
      enddo
      occupied_hamiltonian=(0d0,0d0)
      occupied_hamiltonian(1,1)=0d0;occupied_hamiltonian(2,2)=2d0
      select case(trim(case_name))
      case('average_hamiltonian_nonhermitian')
        occupied_hamiltonian(1,2)=cmplx(0.2d0,0.1d0,8)
      case('average_hamiltonian_nonfinite')
        occupied_hamiltonian(2,2)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,8)
      case('average_hamiltonian_disagree')
        if(rank==0)occupied_hamiltonian(2,2)=3d0
      case('average_hamiltonian_degenerate')
        occupied_hamiltonian=0d0
      end select
      hamiltonian_product=reshape([1,2,2,1],[2,2])
      call eigen_init(comm);call eigen_get_procs(p,info%nprow,info%npcol)
      call eigen_get_id(p,info%myrow,info%mycol);call eigen_get_matdims(4,info%nrow_local,info%ncol_local)
      info%flag_eigenexa_init=.true.
      call build_dg_group_averaged_occupied_candidates_eigenexa(info,comm,hamiltonian_occupied,&
        hamiltonian_weights,hamiltonian_maps,hamiltonian_product,1,1,1d-12,&
        hamiltonian_candidates,hamiltonian_spectrum,hamiltonian_rank,hamiltonian_trace,&
        hamiltonian_closure,hamiltonian_gamma,hamiltonian_workspace,average_ok,average_message,&
        primary_selected,primary_rejected,primary_gap,occupied_hamiltonian=occupied_hamiltonian,&
        secondary_selected_edge=secondary_selected,secondary_rejected_edge=secondary_rejected,&
        secondary_cluster_gap=secondary_gap,primary_boundary_dimension=boundary_dimension,&
        hamiltonian_fingerprint=hamiltonian_fingerprint,&
        secondary_eigensystem_residual=secondary_residual)
      if(trim(case_name)/='average_hamiltonian')then
        call require(.not.average_ok,'invalid occupied Hamiltonian tiebreak must reject collectively')
        if(rank==0)write(*,'(3a,i0)')'REJECT ',trim(case_name),' ranks=',nproc
        call eigen_free();return
      endif
      call require(average_ok,trim(average_message))
      call require(hamiltonian_rank==1.and.boundary_dimension==2.and.&
        abs(primary_selected-1d0)<1d-12.and.abs(primary_rejected-1d0)<1d-12.and.&
        abs(secondary_selected)<1d-12.and.abs(secondary_rejected-2d0)<1d-12.and.&
        abs(secondary_gap-2d0)<1d-12.and.secondary_residual<1d-12.and.&
        hamiltonian_fingerprint/=0_8.and.hamiltonian_workspace>0_8,&
        'occupied Hamiltonian resolves only the primary boundary-degenerate block')
      local_defect=0d0
      do local_index=1,2
        global_point=2*rank+local_index
        local_defect=max(local_defect,abs(hamiltonian_candidates(1,local_index)-&
          merge((1d0,0d0),(0d0,0d0),global_point==1)))
      enddo
      call MPI_Allreduce(local_defect,global_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      call require(ierr==MPI_SUCCESS.and.global_defect<1d-12,&
        'occupied Hamiltonian selects the known low-energy candidate')
      hamiltonian_signature=nint(1d12*sum(abs(hamiltonian_candidates)),8)
      call MPI_Allreduce(MPI_IN_PLACE,hamiltonian_signature,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
      hamiltonian_signature=hamiltonian_signature+int(boundary_dimension,8)
      hamiltonian_signature=ieor(hamiltonian_signature,hamiltonian_fingerprint)
      if(rank==0)write(*,'(a,i0,a,i0)')'AVERAGE_HAMILTONIAN ranks=',nproc,&
        ' signature=',hamiltonian_signature
      call eigen_free()
  end subroutine run_average_hamiltonian_case

  subroutine run_cocycle_case()
    complex(8),allocatable::occupied(:,:),candidates(:,:)
    real(8),allocatable::weights(:),candidate_spectrum(:)
    integer(8),allocatable::translation_maps(:,:),representative_maps(:,:)
    integer(8),allocatable::full_maps(:,:)
    real(8),allocatable::full_total(:),full_boundary(:),full_interior(:)
    logical,allocatable::no_boundary(:)
    integer::point_product(2,2),cocycle(2,2),ii,global_point,global_count,candidate_rank
    real(8)::projector_trace,closure,gamma
    integer(8)::workspace,cocycle_signature
    logical::cocycle_ok
    character(256)::cocycle_message

    global_count=4*nproc
    allocate(occupied(1,4),weights(4),translation_maps(4,2),representative_maps(4,2),&
      full_maps(4,4),full_total(4),full_boundary(4),full_interior(4),no_boundary(4))
    occupied=(0d0,0d0);weights=1d0
    do ii=1,4
      global_point=4*rank+ii
      translation_maps(ii,1)=global_point
      translation_maps(ii,2)=modulo(global_point-1+global_count/2,global_count)+1
      representative_maps(ii,1)=global_point
      representative_maps(ii,2)=modulo(global_point-1+global_count/4,global_count)+1
      do p=1,4
        full_maps(ii,p)=modulo(global_point-1+(p-1)*global_count/4,global_count)+1
      enddo
      if(global_point==1.or.global_point==1+global_count/2)occupied(1,ii)=1d0/sqrt(2d0)
    enddo
    point_product=reshape([1,2,2,1],[2,2])
    cocycle=reshape([1,1,1,2],[2,2])
    call eigen_init(comm);call eigen_get_procs(p,info%nprow,info%npcol)
    call eigen_get_id(p,info%myrow,info%mycol);call eigen_get_matdims(2,info%nrow_local,info%ncol_local)
    info%flag_eigenexa_init=.true.
    call build_dg_cocycle_averaged_occupied_candidates_eigenexa(info,comm,occupied,weights,&
      translation_maps,representative_maps,point_product,cocycle,1,2,1d-12,candidates,&
      candidate_spectrum,candidate_rank,projector_trace,closure,gamma,workspace,cocycle_ok,cocycle_message)
    call require(cocycle_ok,trim(cocycle_message))
    call require(candidate_rank==2.and.maxval(abs(candidate_spectrum-[0.5d0,0.5d0]))<1d-12,&
      'cocycle representative average matches the explicit affine orbit spectrum')
    call require(abs(projector_trace-1d0)<1d-12.and.closure<1d-12.and.gamma==0d0.and.workspace>0_8,&
      'cocycle representative average receipts')
    no_boundary=.false.
    call measure_dg_rank_fixed_symmetry_residuals(comm,candidates,weights,full_maps,no_boundary,&
      total_residual=full_total,boundary_residual=full_boundary,interior_residual=full_interior,&
      ok=cocycle_ok,message=cocycle_message)
    call require(cocycle_ok.and.maxval(full_total)<1d-12,&
      'cocycle representative average equals the explicit full-affine projector')
    cocycle_signature=nint(sum(candidate_spectrum*[1d0,3d0])*1d12,8)
    if(rank==0)write(*,'(a,i0,a,i0)')'COCYCLE ranks=',nproc,' signature=',cocycle_signature
    cocycle(2,2)=1
    call build_dg_cocycle_averaged_occupied_candidates_eigenexa(info,comm,occupied,weights,&
      translation_maps,representative_maps,point_product,cocycle,1,2,1d-12,candidates,&
      candidate_spectrum,candidate_rank,projector_trace,closure,gamma,workspace,cocycle_ok,cocycle_message)
    call require(.not.cocycle_ok,'corrupt representative cocycle must reject')
    call eigen_free()
  end subroutine

  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_failure,global_failure,error
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,comm,error)
    if(global_failure/=0)error stop label
  end subroutine
end program
