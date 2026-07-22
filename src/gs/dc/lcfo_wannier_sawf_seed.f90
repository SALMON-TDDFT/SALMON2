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
  public :: build_sawf_local_band_representation
  public :: canonicalize_sawf_wannier_center
  public :: transform_sawf_wannier_occupation
  public :: build_sawf_projected_wannier
  public :: build_sawf_projected_wannier_from_overlap
  public :: apply_sawf_projected_wannier_transform
  public :: apply_sawf_projected_wannier_gradient_transform
  public :: build_sawf_projected_buffer_gradients
  public :: canonicalize_sawf_bond_identity
  public :: build_sawf_wannier_density
  public :: qualify_sawf_wannier_density_projection
  public :: diagnose_sawf_discrete_wannier_spread
  public :: assemble_sawf_diagonal_periodic_links
  public :: normalize_sawf_projected_wannier_columns

contains

  subroutine normalize_sawf_projected_wannier_columns(norm,polar_transform,core_values,p_values,info)
    real(8),intent(in)::norm(:)
    complex(8),intent(inout)::polar_transform(:,:),core_values(:,:),p_values(:,:)
    integer,intent(out)::info
    real(8)::scale(size(norm))
    integer::iw

    info=1
    if(size(norm)<=0.or.any(shape(polar_transform)/=[size(norm),size(norm)]).or.&
        size(core_values,1)<=0.or.size(p_values,1)<=0.or.&
        size(core_values,2)/=size(norm).or.size(p_values,2)/=size(norm).or.&
        any(norm<=0d0).or..not.all(ieee_is_finite(norm)).or.&
        .not.all(ieee_is_finite(real(polar_transform))).or.&
        .not.all(ieee_is_finite(aimag(polar_transform))).or.&
        .not.all(ieee_is_finite(real(core_values))).or.&
        .not.all(ieee_is_finite(aimag(core_values))).or.&
        .not.all(ieee_is_finite(real(p_values))).or.&
        .not.all(ieee_is_finite(aimag(p_values))))return
    scale=1d0/sqrt(norm)
    if(.not.all(ieee_is_finite(scale)))return
    do iw=1,size(norm)
      polar_transform(:,iw)=polar_transform(:,iw)*scale(iw)
      core_values(:,iw)=core_values(:,iw)*scale(iw)
      p_values(:,iw)=p_values(:,iw)*scale(iw)
    enddo
    info=0
  end subroutine normalize_sawf_projected_wannier_columns

  subroutine assemble_sawf_diagonal_periodic_links(values,coordinates,bvec,quadrature_weight,&
      norm,diagonal_link,info)
    complex(8),intent(in)::values(:,:)
    real(8),intent(in)::coordinates(:,:),bvec(:,:),quadrature_weight
    real(8),intent(out)::norm(:)
    complex(8),intent(out)::diagonal_link(:,:)
    integer,intent(out)::info
    integer::point,neighbor
    complex(8)::phase

    info=1;norm=0d0;diagonal_link=(0d0,0d0)
    if(size(values,1)<=0.or.size(values,2)<=0.or.size(coordinates,1)/=3.or.&
        size(coordinates,2)/=size(values,1).or.size(bvec,1)/=3.or.size(bvec,2)<=0.or.&
        size(norm)/=size(values,2).or.&
        any(shape(diagonal_link)/=[size(values,2),size(bvec,2)]).or.&
        quadrature_weight<=0d0.or..not.ieee_is_finite(quadrature_weight).or.&
        .not.all(ieee_is_finite(real(values))).or..not.all(ieee_is_finite(aimag(values))).or.&
        .not.all(ieee_is_finite(coordinates)).or..not.all(ieee_is_finite(bvec)))return
    do point=1,size(values,1)
      norm=norm+abs(values(point,:))**2*quadrature_weight
      do neighbor=1,size(bvec,2)
        phase=exp(cmplx(0d0,-dot_product(bvec(:,neighbor),coordinates(:,point)),8))
        diagonal_link(:,neighbor)=diagonal_link(:,neighbor)+&
          abs(values(point,:))**2*phase*quadrature_weight
      enddo
    enddo
    if(any(norm<=0d0).or..not.all(ieee_is_finite(norm)).or.&
        .not.all(ieee_is_finite(real(diagonal_link))).or.&
        .not.all(ieee_is_finite(aimag(diagonal_link))))then
      norm=0d0;diagonal_link=(0d0,0d0);return
    endif
    info=0
  end subroutine assemble_sawf_diagonal_periodic_links

  subroutine diagnose_sawf_discrete_wannier_spread(diagonal_link,norm,bvec,weight,&
      center,omega,center_valid,info,require_unit_norm)
    complex(8),intent(in)::diagonal_link(:,:)
    real(8),intent(in)::norm(:),bvec(:,:),weight(:)
    real(8),intent(out)::center(:,:),omega(:)
    logical,intent(out)::center_valid(:)
    integer,intent(out)::info
    logical,intent(in),optional::require_unit_norm
    real(8)::metric(3,3),identity(3,3),phase(size(weight)),scale,pair_tolerance
    complex(8)::link(size(weight))
    integer::iw,ib,jb,axis,pair_count

    info=1;center=0d0;omega=huge(1d0);center_valid=.false.
    if(size(diagonal_link,1)<=0.or.size(diagonal_link,2)<=0.or.&
        size(norm)/=size(diagonal_link,1).or.size(bvec,1)/=3.or.&
        size(bvec,2)/=size(diagonal_link,2).or.size(weight)/=size(diagonal_link,2).or.&
        any(shape(center)/=[3,size(diagonal_link,1)]).or.&
        size(omega)/=size(diagonal_link,1).or.size(center_valid)/=size(diagonal_link,1).or.&
        any(norm<=0d0).or.any(weight<=0d0).or.&
        .not.all(ieee_is_finite(real(diagonal_link))).or.&
        .not.all(ieee_is_finite(aimag(diagonal_link))).or.&
        .not.all(ieee_is_finite(norm)).or..not.all(ieee_is_finite(bvec)).or.&
        .not.all(ieee_is_finite(weight)))return
    if(present(require_unit_norm))then
      if(require_unit_norm.and.any(abs(norm-1d0)>1d-8))return
    endif
    metric=0d0
    do ib=1,size(weight)
      do axis=1,3
        metric(axis,:)=metric(axis,:)+weight(ib)*bvec(axis,ib)*bvec(:,ib)
      enddo
    enddo
    identity=0d0
    do axis=1,3;identity(axis,axis)=1d0;enddo
    if(maxval(abs(metric-identity))>1000d0*epsilon(1d0))return
    do iw=1,size(diagonal_link,1)
      link=diagonal_link(iw,:)/norm(iw)
      scale=max(1d0,maxval(abs(link)))
      pair_tolerance=1000d0*epsilon(1d0)*scale
      do ib=1,size(weight)
        pair_count=0
        do jb=1,size(weight)
          if(maxval(abs(bvec(:,jb)+bvec(:,ib)))<=1000d0*epsilon(1d0)*&
              max(1d0,maxval(abs(bvec))))then
            pair_count=pair_count+1
            if(abs(weight(jb)-weight(ib))>1000d0*epsilon(1d0)*max(1d0,weight(ib)).or.&
                abs(link(jb)-conjg(link(ib)))>pair_tolerance)return
          endif
        enddo
        if(pair_count/=1)return
      enddo
      if(any(abs(link)<=sqrt(epsilon(1d0))))cycle
      phase=atan2(aimag(link),real(link,8))
      center(:,iw)=0d0
      do ib=1,size(weight)
        center(:,iw)=center(:,iw)-weight(ib)*bvec(:,ib)*phase(ib)
      enddo
      omega(iw)=0d0
      do ib=1,size(weight)
        omega(iw)=omega(iw)+weight(ib)*(1d0-abs(link(ib))**2+&
          (phase(ib)+dot_product(bvec(:,ib),center(:,iw)))**2)
      enddo
      if(omega(iw)<-1000d0*epsilon(1d0)*max(1d0,maxval(weight)))return
      omega(iw)=max(0d0,omega(iw));center_valid(iw)=.true.
    enddo
    info=0
  end subroutine diagnose_sawf_discrete_wannier_spread

  subroutine canonicalize_sawf_wannier_center(center,cell_length,core_lower,core_upper,tolerance,&
      wrapped,image,owned,info)
    real(8),intent(in)::center(3),cell_length(3),core_lower(3),core_upper(3),tolerance
    real(8),intent(out)::wrapped(3)
    integer,intent(out)::image(3),info
    logical,intent(out)::owned
    integer::axis

    info=1;wrapped=0d0;image=0;owned=.false.
    if(any(cell_length<=0d0).or.any(core_lower<0d0).or.any(core_upper>cell_length).or.&
        any(core_upper<=core_lower).or.tolerance<=0d0.or.&
        .not.all(ieee_is_finite(center)).or..not.all(ieee_is_finite(cell_length)).or.&
        .not.all(ieee_is_finite(core_lower)).or..not.all(ieee_is_finite(core_upper)).or.&
        .not.ieee_is_finite(tolerance))return
    do axis=1,3
      image(axis)=floor(center(axis)/cell_length(axis))
      wrapped(axis)=center(axis)-dble(image(axis))*cell_length(axis)
      if(wrapped(axis)<0d0)then
        wrapped(axis)=wrapped(axis)+cell_length(axis);image(axis)=image(axis)-1
      endif
      if(abs(wrapped(axis)-cell_length(axis))<=tolerance)then
        wrapped(axis)=0d0;image(axis)=image(axis)+1
      elseif(abs(wrapped(axis))<=tolerance)then
        wrapped(axis)=0d0
      endif
      if(abs(wrapped(axis)-core_lower(axis))<=tolerance)wrapped(axis)=core_lower(axis)
      if(abs(wrapped(axis)-core_upper(axis))<=tolerance)wrapped(axis)=core_upper(axis)
    enddo
    owned=all(wrapped>=core_lower.and.wrapped<core_upper)
    info=0
  end subroutine canonicalize_sawf_wannier_center

  subroutine canonicalize_sawf_bond_identity(atom_a,atom_b,directed_image,atoms,canonical_image,info)
    integer,intent(in)::atom_a,atom_b,directed_image(3)
    integer,intent(out)::atoms(2),canonical_image(3),info
    info=1;atoms=0;canonical_image=0
    if(atom_a<=0.or.atom_b<=0.or.atom_a==atom_b)return
    if(atom_a<atom_b)then
      atoms=[atom_a,atom_b];canonical_image=directed_image
    else
      atoms=[atom_b,atom_a];canonical_image=-directed_image
    endif
    info=0
  end subroutine canonicalize_sawf_bond_identity

  subroutine build_sawf_wannier_density(values,occupation,weight,density,charge,info)
    complex(8),intent(in)::values(:,:),occupation(:,:)
    real(8),intent(in)::weight(:)
    real(8),intent(out)::density(:),charge
    integer,intent(out)::info
    integer::point,nsource
    real(8)::scale
    complex(8)::rho

    info=1;density=0d0;charge=0d0;nsource=size(values,2)
    if(size(values,1)<=0.or.nsource<=0.or.any(shape(occupation)/=[nsource,nsource]).or.&
        size(weight)/=size(values,1).or.size(density)/=size(values,1).or.any(weight<=0d0).or.&
        .not.all(ieee_is_finite(real(values))).or..not.all(ieee_is_finite(aimag(values))).or.&
        .not.all(ieee_is_finite(real(occupation))).or..not.all(ieee_is_finite(aimag(occupation))).or.&
        .not.all(ieee_is_finite(weight)))return
    scale=max(1d0,maxval(abs(occupation)))
    if(maxval(abs(occupation-conjg(transpose(occupation))))>100d0*epsilon(1d0)*scale)return
    do point=1,size(values,1)
      rho=sum(values(point,:)*matmul(occupation,conjg(values(point,:))))
      if(abs(aimag(rho))>100d0*epsilon(1d0)*max(1d0,abs(real(rho,8))))return
      density(point)=real(rho,8)
      if(density(point)<-100d0*epsilon(1d0)*scale)return
      density(point)=max(0d0,density(point))
    enddo
    charge=sum(weight*density)
    if(.not.ieee_is_finite(charge))return
    info=0
  end subroutine build_sawf_wannier_density

  subroutine qualify_sawf_wannier_density_projection(source,projected,normalized,weight,expected_charge,&
      captured_norm,tolerance,projection_residual,normalization_residual,charge_error,info)
    real(8),intent(in)::source(:),projected(:),normalized(:),weight(:)
    real(8),intent(in)::expected_charge,captured_norm,tolerance
    real(8),intent(out)::projection_residual,normalization_residual,charge_error
    integer,intent(out)::info
    real(8)::source_norm,projected_norm

    info=1;projection_residual=huge(1d0);normalization_residual=huge(1d0);charge_error=huge(1d0)
    if(size(source)<=0.or.size(projected)/=size(source).or.size(normalized)/=size(source).or.&
        size(weight)/=size(source).or.any(weight<=0d0).or.any(source<0d0).or.any(projected<0d0).or.&
        any(normalized<0d0).or.expected_charge<=0d0.or.tolerance<=0d0.or.&
        .not.all(ieee_is_finite(source)).or..not.all(ieee_is_finite(projected)).or.&
        .not.all(ieee_is_finite(normalized)).or..not.all(ieee_is_finite(weight)).or.&
        .not.all(ieee_is_finite([expected_charge,captured_norm,tolerance])))return
    source_norm=sqrt(sum(weight*source**2));projected_norm=sqrt(sum(weight*projected**2))
    if(source_norm<=0d0.or.projected_norm<=0d0)return
    projection_residual=sqrt(sum(weight*(projected-source)**2))/source_norm
    normalization_residual=sqrt(sum(weight*(normalized-projected)**2))/projected_norm
    charge_error=abs(sum(weight*projected)-expected_charge)
    if(max(projection_residual,max(normalization_residual,max(charge_error,&
      abs(1d0-captured_norm))))>tolerance)return
    info=0
  end subroutine qualify_sawf_wannier_density_projection

  subroutine transform_sawf_wannier_occupation(transform,f_source,f_normalized,info)
    complex(8),intent(in)::transform(:,:),f_source(:,:)
    complex(8),intent(out)::f_normalized(:,:)
    integer,intent(out)::info
    complex(8),allocatable::inverse(:,:),work(:)
    integer,allocatable::pivot(:)
    integer::n,lwork
    external::zgetrf,zgetri

    info=1;f_normalized=(0d0,0d0);n=size(transform,1)
    if(n<=0.or.size(transform,2)/=n.or.any(shape(f_source)/=[n,n]).or.&
        any(shape(f_normalized)/=[n,n]).or..not.all(ieee_is_finite(real(transform))).or.&
        .not.all(ieee_is_finite(aimag(transform))).or..not.all(ieee_is_finite(real(f_source))).or.&
        .not.all(ieee_is_finite(aimag(f_source))))return
    allocate(inverse(n,n),pivot(n),work(1));inverse=transform
    call zgetrf(n,n,inverse,n,pivot,info);if(info/=0)return
    lwork=-1;call zgetri(n,inverse,n,pivot,work,lwork,info)
    if(info/=0.or..not.ieee_is_finite(real(work(1))).or..not.ieee_is_finite(aimag(work(1))))return
    lwork=max(1,int(real(work(1))));deallocate(work);allocate(work(lwork))
    call zgetri(n,inverse,n,pivot,work,lwork,info);if(info/=0)return
    f_normalized=matmul(inverse,matmul(f_source,conjg(transpose(inverse))))
    f_normalized=0.5d0*(f_normalized+conjg(transpose(f_normalized)))
    if(.not.all(ieee_is_finite(real(f_normalized))).or.&
        .not.all(ieee_is_finite(aimag(f_normalized))))then;info=1;return;endif
    info=0
  end subroutine transform_sawf_wannier_occupation

  subroutine build_sawf_projected_wannier(occupied,trial,weight,rank_cutoff,wannier,condition,info)
    complex(8),intent(in)::occupied(:,:),trial(:,:)
    real(8),intent(in)::weight,rank_cutoff
    complex(8),intent(out)::wannier(:,:)
    real(8),intent(out)::condition
    integer,intent(out)::info
    complex(8),allocatable::projected(:,:),gram(:,:),inverse_sqrt(:,:),work(:),identity(:,:)
    real(8),allocatable::eigenvalue(:),rwork(:)
    integer::npoint,noccupied,ntrial,i,j,lwork
    external::zheev

    info=1;wannier=(0d0,0d0);condition=huge(1d0)
    npoint=size(occupied,1);noccupied=size(occupied,2);ntrial=size(trial,2)
    if(npoint<=0.or.noccupied<ntrial.or.ntrial<=0.or.size(trial,1)/=npoint.or.&
        any(shape(wannier)/=[npoint,ntrial]).or.weight<=0d0.or.rank_cutoff<=0d0.or.&
        .not.ieee_is_finite(weight).or..not.ieee_is_finite(rank_cutoff).or.&
        .not.all(ieee_is_finite(real(occupied))).or..not.all(ieee_is_finite(aimag(occupied))).or.&
        .not.all(ieee_is_finite(real(trial))).or..not.all(ieee_is_finite(aimag(trial))))return
    allocate(identity(noccupied,noccupied));identity=0
    do i=1,noccupied;identity(i,i)=1;enddo
    if(maxval(abs(weight*matmul(conjg(transpose(occupied)),occupied)-identity))>&
        100d0*rank_cutoff)return
    allocate(projected(npoint,ntrial),gram(ntrial,ntrial),inverse_sqrt(ntrial,ntrial),&
      eigenvalue(ntrial),rwork(max(1,3*ntrial-2)),work(1))
    projected=matmul(occupied,weight*matmul(conjg(transpose(occupied)),trial))
    gram=weight*matmul(conjg(transpose(projected)),projected)
    lwork=-1;call zheev('V','U',ntrial,gram,ntrial,eigenvalue,work,lwork,rwork,info)
    if(info/=0.or..not.ieee_is_finite(real(work(1))))return
    lwork=max(1,int(real(work(1))));deallocate(work);allocate(work(lwork))
    call zheev('V','U',ntrial,gram,ntrial,eigenvalue,work,lwork,rwork,info)
    if(info/=0.or..not.all(ieee_is_finite(eigenvalue)).or.&
        eigenvalue(1)<=rank_cutoff*max(1d0,eigenvalue(ntrial)))then;info=1;return;endif
    condition=eigenvalue(ntrial)/eigenvalue(1);inverse_sqrt=0
    do j=1,ntrial;do i=1,ntrial
      inverse_sqrt(i,j)=sum(gram(i,:)*conjg(gram(j,:))/sqrt(eigenvalue(:)))
    enddo;enddo
    wannier=matmul(projected,inverse_sqrt)
    if(.not.all(ieee_is_finite(real(wannier))).or..not.all(ieee_is_finite(aimag(wannier))))then
      info=1;return
    endif
    info=0
  end subroutine build_sawf_projected_wannier

  subroutine build_sawf_projected_wannier_from_overlap(occupied_values,projection_overlap,weight,&
      rank_cutoff,wannier_values,polar_transform,condition,info)
    complex(8),intent(in)::occupied_values(:,:),projection_overlap(:,:)
    real(8),intent(in)::weight,rank_cutoff
    complex(8),intent(out)::wannier_values(:,:),polar_transform(:,:)
    real(8),intent(out)::condition
    integer,intent(out)::info
    complex(8),allocatable::gram(:,:),work(:)
    real(8),allocatable::eigenvalue(:),rwork(:)
    integer::noccupied,ntrial,i,j,lwork
    external::zheev

    info=1;wannier_values=(0d0,0d0);polar_transform=(0d0,0d0);condition=huge(1d0)
    noccupied=size(occupied_values,2);ntrial=size(projection_overlap,2)
    if(size(occupied_values,1)<=0.or.noccupied<=0.or.ntrial<=0.or.&
        size(projection_overlap,1)/=noccupied.or.&
        any(shape(wannier_values)/=[size(occupied_values,1),ntrial]).or.&
        any(shape(polar_transform)/=[ntrial,ntrial]).or.weight<=0d0.or.rank_cutoff<=0d0.or.&
        .not.ieee_is_finite(weight).or..not.ieee_is_finite(rank_cutoff).or.&
        .not.all(ieee_is_finite(real(occupied_values))).or.&
        .not.all(ieee_is_finite(aimag(occupied_values))).or.&
        .not.all(ieee_is_finite(real(projection_overlap))).or.&
        .not.all(ieee_is_finite(aimag(projection_overlap))))return
    allocate(gram(ntrial,ntrial),eigenvalue(ntrial),rwork(max(1,3*ntrial-2)),work(1))
    gram=matmul(conjg(transpose(projection_overlap)),projection_overlap)
    lwork=-1;call zheev('V','U',ntrial,gram,ntrial,eigenvalue,work,lwork,rwork,info)
    if(info/=0.or..not.ieee_is_finite(real(work(1))))return
    lwork=max(1,int(real(work(1))));deallocate(work);allocate(work(lwork))
    call zheev('V','U',ntrial,gram,ntrial,eigenvalue,work,lwork,rwork,info)
    if(info/=0.or..not.all(ieee_is_finite(eigenvalue)).or.&
        eigenvalue(1)<=rank_cutoff*max(1d0,eigenvalue(ntrial)))then;info=1;return;endif
    condition=eigenvalue(ntrial)/eigenvalue(1);polar_transform=(0d0,0d0)
    do j=1,ntrial;do i=1,ntrial
      polar_transform(i,j)=sum(gram(i,:)*conjg(gram(j,:))/sqrt(eigenvalue(:)))
    enddo;enddo
    wannier_values=matmul(matmul(occupied_values,projection_overlap),polar_transform)
    if(.not.all(ieee_is_finite(real(wannier_values))).or.&
        .not.all(ieee_is_finite(aimag(wannier_values))))then;info=1;return;endif
    info=0
  end subroutine build_sawf_projected_wannier_from_overlap

  subroutine apply_sawf_projected_wannier_transform(occupied_samples,projection_overlap,&
      polar_transform,wannier_samples,info)
    complex(8),intent(in)::occupied_samples(:,:),projection_overlap(:,:),polar_transform(:,:)
    complex(8),intent(out)::wannier_samples(:,:)
    integer,intent(out)::info
    integer::noccupied,ntrial

    info=1;wannier_samples=(0d0,0d0)
    noccupied=size(occupied_samples,2);ntrial=size(projection_overlap,2)
    if(size(occupied_samples,1)<=0.or.noccupied<=0.or.ntrial<=0.or.&
        size(projection_overlap,1)/=noccupied.or.&
        any(shape(polar_transform)/=[ntrial,ntrial]).or.&
        any(shape(wannier_samples)/=[size(occupied_samples,1),ntrial]).or.&
        .not.all(ieee_is_finite(real(occupied_samples))).or.&
        .not.all(ieee_is_finite(aimag(occupied_samples))).or.&
        .not.all(ieee_is_finite(real(projection_overlap))).or.&
        .not.all(ieee_is_finite(aimag(projection_overlap))).or.&
        .not.all(ieee_is_finite(real(polar_transform))).or.&
        .not.all(ieee_is_finite(aimag(polar_transform))))return
    wannier_samples=matmul(matmul(occupied_samples,projection_overlap),polar_transform)
    if(.not.all(ieee_is_finite(real(wannier_samples))).or.&
        .not.all(ieee_is_finite(aimag(wannier_samples))))then
      wannier_samples=(0d0,0d0);return
    endif
    info=0
  end subroutine apply_sawf_projected_wannier_transform

  subroutine apply_sawf_projected_wannier_gradient_transform(occupied_gradients,&
      projection_overlap,polar_transform,wannier_gradients,info)
    complex(8),intent(in)::occupied_gradients(:,:,:),projection_overlap(:,:),polar_transform(:,:)
    complex(8),intent(out)::wannier_gradients(:,:,:)
    integer,intent(out)::info
    integer::axis,local_info,npoint,noccupied,ntrial

    info=1;wannier_gradients=(0d0,0d0)
    npoint=size(occupied_gradients,2);noccupied=size(occupied_gradients,3)
    ntrial=size(projection_overlap,2)
    if(size(occupied_gradients,1)/=3.or.npoint<=0.or.noccupied<=0.or.ntrial<=0.or.&
        size(projection_overlap,1)/=noccupied.or.any(shape(polar_transform)/=[ntrial,ntrial]).or.&
        any(shape(wannier_gradients)/=[3,npoint,ntrial]).or.&
        .not.all(ieee_is_finite(real(occupied_gradients))).or.&
        .not.all(ieee_is_finite(aimag(occupied_gradients))))return
    do axis=1,3
      call apply_sawf_projected_wannier_transform(occupied_gradients(axis,:,:),projection_overlap,&
        polar_transform,wannier_gradients(axis,:,:),local_info)
      if(local_info/=0)then
        wannier_gradients=(0d0,0d0);return
      endif
    enddo
    info=0
  end subroutine apply_sawf_projected_wannier_gradient_transform

  subroutine build_sawf_projected_buffer_gradients(occupied_stencil,gradient_coefficients,&
      projection_overlap,polar_transform,wannier_gradients,info)
    complex(8),intent(in)::occupied_stencil(:,:,:,:),projection_overlap(:,:),polar_transform(:,:)
    real(8),intent(in)::gradient_coefficients(:,:)
    complex(8),intent(out)::wannier_gradients(:,:,:)
    integer,intent(out)::info
    complex(8),allocatable::occupied_gradient(:,:)
    integer::radius,extent(3),noccupied,ntrial,ix,iy,iz,point,axis,dist,local_info

    info=1;wannier_gradients=(0d0,0d0);radius=size(gradient_coefficients,1)
    extent=[size(occupied_stencil,1),size(occupied_stencil,2),size(occupied_stencil,3)]-2*radius
    noccupied=size(occupied_stencil,4);ntrial=size(projection_overlap,2)
    if(radius<=0.or.size(gradient_coefficients,2)/=3.or.any(extent<=0).or.noccupied<=0.or.&
        ntrial<=0.or.size(projection_overlap,1)/=noccupied.or.&
        any(shape(polar_transform)/=[ntrial,ntrial]).or.&
        any(shape(wannier_gradients)/=[3,product(extent),ntrial]).or.&
        .not.all(ieee_is_finite(gradient_coefficients)).or.&
        .not.all(ieee_is_finite(real(occupied_stencil))).or.&
        .not.all(ieee_is_finite(aimag(occupied_stencil))))return
    allocate(occupied_gradient(product(extent),noccupied))
    do axis=1,3
      occupied_gradient=(0d0,0d0);point=0
      do iz=1,extent(3);do iy=1,extent(2);do ix=1,extent(1)
        point=point+1
        do dist=1,radius
          select case(axis)
          case(1)
            occupied_gradient(point,:)=occupied_gradient(point,:)+gradient_coefficients(dist,axis)*&
              (occupied_stencil(ix+radius+dist,iy+radius,iz+radius,:)-&
               occupied_stencil(ix+radius-dist,iy+radius,iz+radius,:))
          case(2)
            occupied_gradient(point,:)=occupied_gradient(point,:)+gradient_coefficients(dist,axis)*&
              (occupied_stencil(ix+radius,iy+radius+dist,iz+radius,:)-&
               occupied_stencil(ix+radius,iy+radius-dist,iz+radius,:))
          case(3)
            occupied_gradient(point,:)=occupied_gradient(point,:)+gradient_coefficients(dist,axis)*&
              (occupied_stencil(ix+radius,iy+radius,iz+radius+dist,:)-&
               occupied_stencil(ix+radius,iy+radius,iz+radius-dist,:))
          end select
        enddo
      enddo;enddo;enddo
      call apply_sawf_projected_wannier_transform(occupied_gradient,projection_overlap,&
        polar_transform,wannier_gradients(axis,:,:),local_info)
      if(local_info/=0)return
    enddo
    info=0
  end subroutine build_sawf_projected_buffer_gradients

  subroutine build_sawf_local_band_representation(states,point_map,weight,tolerance,representation,ok,message)
    complex(8),intent(in)::states(:,:)
    integer,intent(in)::point_map(:,:)
    real(8),intent(in)::weight,tolerance
    complex(8),intent(out)::representation(:,:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::transformed(:,:),identity(:,:)
    logical,allocatable::seen(:)
    real(8)::unitarity
    integer::operation,source,target,state,npoint,nstate

    ok=.false.;message='';representation=(0d0,0d0);npoint=size(states,1);nstate=size(states,2)
    if(npoint<=0.or.nstate<=0.or.size(point_map,1)/=npoint.or.size(point_map,2)<=0.or. &
        size(representation,1)/=nstate.or.size(representation,2)/=nstate.or. &
        size(representation,3)/=size(point_map,2).or..not.ieee_is_finite(weight).or.weight<=0d0.or. &
        .not.ieee_is_finite(tolerance).or.tolerance<=0d0.or. &
        .not.all(ieee_is_finite(real(states))).or..not.all(ieee_is_finite(aimag(states))))then
      message='SAWF local D_band inputs are invalid';return
    end if
    allocate(transformed(npoint,nstate),seen(npoint),identity(nstate,nstate))
    identity=(0d0,0d0);do state=1,nstate;identity(state,state)=(1d0,0d0);end do
    do operation=1,size(point_map,2)
      transformed=(0d0,0d0);seen=.false.
      do source=1,npoint
        target=point_map(source,operation)
        if(target<1.or.target>npoint.or.seen(target))then
          message='SAWF local stabilizer grid map is not a permutation';return
        end if
        seen(target)=.true.;transformed(target,:)=states(source,:)
      end do
      if(.not.all(seen))then;message='SAWF local stabilizer grid map is incomplete';return;end if
      representation(:,:,operation)=weight*matmul(conjg(transpose(states)),transformed)
      unitarity=maxval(abs(matmul(conjg(transpose(representation(:,:,operation))), &
        representation(:,:,operation))-identity))
      if(unitarity>tolerance)then
        message='SAWF local D_band is not unitary on the retained eigenspace';return
      end if
    end do
    ok=.true.
  end subroutine build_sawf_local_band_representation

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
      states,energies,ok,message,coefficients)
    real(8),intent(in)::buffer_basis(:,:),weight,h_basis(:,:),rank_tolerance
    complex(8),allocatable,intent(out)::states(:,:)
    real(8),allocatable,intent(out)::energies(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8),allocatable,intent(out),optional::coefficients(:,:)
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
    if(present(coefficients))then
      allocate(coefficients(size(coeff,1),size(coeff,2)));coefficients=coeff
    end if
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
