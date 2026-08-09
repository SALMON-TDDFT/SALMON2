module dg_overlapping_wannier_w90
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::estimate_dg_w90_coordinator_bytes,validate_dg_w90_result
  public::setup_dg_w90_gamma_library,run_dg_w90_gamma_library
  public::assemble_dg_w90_gamma_matrices
contains
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
    integer::mp_grid(3),nnlist(1,num_nnmax),nncell_max(3,1,num_nnmax),exclude_bands(max(1,nband))
    integer::proj_l(max(1,nband)),proj_m(max(1,nband)),proj_radial(max(1,nband))
    integer::proj_s(max(1,nband))
    real(real64)::kpoint(3,1),proj_site(3,max(1,nband)),proj_z(3,max(1,nband)),&
      proj_x(3,max(1,nband)),proj_zona(max(1,nband)),proj_s_qaxis(3,max(1,nband))
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
      open(newunit=unit,file=trim(seed)//'.win',status='replace',action='write',iostat=io)
      if(io/=0)then
        status=2
      else
        write(unit,'(a,i0)')'num_bands = ',nband
        write(unit,'(a,i0)')'num_wann = ',nwann
        write(unit,'(a)')'num_iter = 200'
        write(unit,'(a)')'conv_tol = 1.d-12'
        write(unit,'(a)')'conv_window = 5'
        write(unit,'(a)')'gamma_only = true'
        write(unit,'(a)')'begin unit_cell_cart';write(unit,'(a)')'bohr'
        do axis=1,3;write(unit,'(3(es24.16,1x))')real_lattice(:,axis);enddo
        write(unit,'(a)')'end unit_cell_cart'
        write(unit,'(a)')'begin atoms_cart';write(unit,'(a)')'bohr'
        do atom=1,size(atom_symbols)
          write(unit,'(a,1x,3(es24.16,1x))')trim(atom_symbols(atom)),atoms_cart(:,atom)
        enddo
        write(unit,'(a)')'end atoms_cart'
        write(unit,'(a)')'begin projections';write(unit,'(a)')'random'
        write(unit,'(a)')'end projections'
        write(unit,'(a)')'mp_grid = 1 1 1'
        write(unit,'(a)')'begin kpoints';write(unit,'(a)')'0.0 0.0 0.0'
        write(unit,'(a)')'end kpoints';close(unit)
        call wannier_setup(trim(seed),mp_grid,1,real_lattice,reciprocal_lattice,kpoint,nband,&
          size(atom_symbols),atom_symbols,atoms_cart,.true.,.false.,nntot,nnlist,nncell_max,&
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
      spreads,spread,ok,message)
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
#if defined(USE_MPI) && defined(USE_WANNIER90)
    integer::rank,ierr,nband,nwann,nntot,status,mp_grid(3)
    real(real64)::kpoint(3,1)
    complex(real64),allocatable::u(:,:,:),uopt(:,:,:),m4(:,:,:,:),a3(:,:,:)
    real(real64),allocatable::e2(:,:)
    logical,allocatable::lwindow(:,:)
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
    ok=.false.;message='';spread=0d0;status=0
    call MPI_Comm_rank(comm,rank,ierr)
    nband=size(m_matrix,1);nwann=size(a_matrix,2);nntot=size(m_matrix,3)
    if(ierr/=MPI_SUCCESS.or.nband<=0.or.size(m_matrix,2)/=nband.or.nwann/=nband.or.&
        size(a_matrix,1)/=nband.or.size(eigenvalues)/=nband.or.nntot<=0.or.&
        any(shape(atoms_cart)/=[3,size(atom_symbols)]).or.&
        .not.all(ieee_is_finite(real(m_matrix))).or..not.all(ieee_is_finite(aimag(m_matrix))).or.&
        .not.all(ieee_is_finite(real(a_matrix))).or..not.all(ieee_is_finite(aimag(a_matrix))).or.&
        .not.all(ieee_is_finite(eigenvalues)).or..not.all(ieee_is_finite(real_lattice)).or.&
        .not.all(ieee_is_finite(reciprocal_lattice)).or..not.all(ieee_is_finite(atoms_cart)))status=1
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(status/=0.or.ierr/=MPI_SUCCESS)then;message='invalid Gamma Wannier90 run contract';return;endif
    allocate(transform(nwann,nwann),centers(3,nwann),spreads(nwann));transform=(0d0,0d0)
    centers=0d0;spreads=0d0;mp_grid=[1,1,1];kpoint=0d0
    if(rank==0)then
      allocate(u(nwann,nwann,1),uopt(nband,nwann,1),lwindow(nband,1),&
        m4(nband,nband,nntot,1),a3(nband,nwann,1),e2(nband,1))
      m4(:,:,:,1)=m_matrix;a3(:,:,1)=a_matrix;e2(:,1)=eigenvalues
      call wannier_run(trim(seed),mp_grid,1,real_lattice,reciprocal_lattice,kpoint,nband,nwann,&
        nntot,size(atom_symbols),atom_symbols,atoms_cart,.true.,m4,a3,e2,u,uopt,lwindow,&
        centers,spreads,spread)
      transform=u(:,:,1)
      call validate_dg_w90_result(transform,centers,spreads,spread,initial_gauge_spread,&
        tolerance,ok,message)
      status=merge(0,2,ok)
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(transform,size(transform),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    call MPI_Bcast(centers,size(centers),MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(spreads,size(spreads),MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(spread,3,MPI_DOUBLE_PRECISION,0,comm,ierr)
    ok=status==0.and.ierr==MPI_SUCCESS
    if(ok)then;message='';else;message='Wannier90 Gamma library run failed validation';endif
#else
    ok=.false.;message='Wannier90 Gamma run requires MPI and USE_WANNIER90';spread=0d0
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
    logical,intent(out)::ok
    ok=value>=0_int64.and.increment>=0_int64.and.value<=huge(value)-increment
    if(ok)value=value+increment
  end subroutine checked_add
end module dg_overlapping_wannier_w90
