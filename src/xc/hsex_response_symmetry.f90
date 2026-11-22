! Conservative target-k little group on an exact fractional permutation grid.
module hsex_response_symmetry
  implicit none
  private
  public :: s_response_stars,response_q_stars
  type :: s_response_stars
    integer :: nsym=0,nrep=0,h_probes=0
    integer,allocatable :: representative(:),operation(:),column(:,:)
    real(8) :: max_h_covariance=0d0,max_projector_error=0d0,max_energy_error=0d0
  end type
  abstract interface
    subroutine h_action(x,y,status)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  subroutine response_q_stars(n,mesh,qmesh,cell,rec,atoms,species,frac,levels,occupied, &
      destination,q,g,counts,sources,apply,stars,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    integer,intent(in) :: n(3),mesh(3),qmesh(3),species(:),destination,g(:,:,:),counts(:),sources(:)
    real(8),intent(in) :: cell(3,3),rec(3,3),atoms(:,:),frac(:,:),levels(:,:),q(:,:)
    complex(8),intent(in) :: occupied(:,:,:)
    procedure(h_action) :: apply
    type(s_response_stars),intent(out) :: stars
    integer,intent(out) :: status
    integer,parameter :: permutations(3,6)=reshape([1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1],[3,6])
    real(8),parameter :: projector_tol=1d-6,energy_tol=1d-7,covariance_tol=1d-7
    integer :: rotations(3,3,6),rotation(3,3),ns,ng,nv,nk,nq,op,i,j,a,b,site,coord(3),old(3),delta(3)
    integer :: iq,jq,source,dest_source,ig,jg,trial,info,found,used_atoms(size(species)),mapped(3),capacity
    integer,allocatable :: grid(:,:),columns(:)
    complex(8),allocatable :: phase(:,:),transformed(:,:),overlap(:,:),probe(:),action(:),moved(:),moved_action(:)
    real(8),allocatable :: positions(:,:)
    real(8) :: metric(3,3),difference(3),xyz(3),h_error,p_error,e_error,q_delta(3),scale
    logical :: valid
    status=1;ng=product(n);nv=size(occupied,2);nk=size(frac,2);nq=size(q,2)
    if(any(n<1).or.any(mesh<1).or.any(qmesh<1).or.nv<1.or.nq<1)return
    if(destination<1.or.destination>nk.or.size(occupied,1)/=ng.or.size(occupied,3)/=nk)return
    if(size(frac,1)/=3.or.size(q,1)/=3.or.size(counts)/=nq.or.size(sources)/=nq)return
    if(any(counts<1).or.any(counts>size(g,2)).or.size(g,1)/=3.or.size(g,3)/=nq)return
    if(size(levels,1)<nv.or.size(levels,2)/=nk.or.any(sources<1).or.any(sources>nk))return
    if(size(atoms,1)/=3.or.size(atoms,2)/=size(species))return
    if(.not.screen_finite(cell).or..not.screen_finite(rec).or..not.screen_finite(atoms))return
    if(.not.screen_finite(frac).or..not.screen_finite(q).or..not.screen_finite(levels))return
    if(.not.screen_finite(occupied))return
    capacity=size(g,2)
    allocate(grid(ng,6),phase(ng,6),positions(3,size(species)),transformed(ng,nv),overlap(nv,nv), &
      probe(ng),action(ng),moved(ng),moved_action(ng),columns(capacity), &
      stars%representative(nq),stars%operation(nq),stars%column(capacity,nq),stat=info)
    if(info/=0)return
    positions=matmul(transpose(rec),atoms)/(2d0*acos(-1d0))
    metric=matmul(transpose(cell),cell);ns=0
    do op=1,6
      rotation=0
      do i=1,3
        rotation(i,permutations(i,op))=1
      enddo
      if(any(n/=n(permutations(:,op))).or.any(mesh/=mesh(permutations(:,op))).or. &
          any(qmesh/=qmesh(permutations(:,op))))cycle
      if(maxval(abs(matmul(transpose(rotation),matmul(metric,rotation))-metric))> &
          1d-10*max(1d0,maxval(abs(metric))))cycle
      difference=matmul(rotation,frac(:,destination))-frac(:,destination)
      if(maxval(abs(difference-anint(difference)))>1d-10)cycle
      valid=.true.;used_atoms=0
      do a=1,size(species)
        xyz=matmul(rotation,positions(:,a));found=0
        do b=1,size(species)
          if(species(a)/=species(b).or.used_atoms(b)/=0)cycle
          difference=xyz-positions(:,b);difference=difference-anint(difference)
          if(maxval(abs(difference))<1d-8)then
            found=b
            exit
          endif
        enddo
        if(found==0)then
          valid=.false.
          exit
        endif
        used_atoms(found)=1
      enddo
      if(.not.valid)cycle
      do site=0,ng-1
        coord=[mod(site,n(1)),mod(site/n(1),n(2)),site/(n(1)*n(2))]
        old=matmul(transpose(rotation),coord)
        grid(site+1,ns+1)=1+old(1)+n(1)*(old(2)+n(2)*old(3))
        xyz=dble(coord)/n
        phase(site+1,ns+1)=exp(cmplx(0d0,2d0*acos(-1d0)* &
          dot_product(matmul(rotation,frac(:,destination))-frac(:,destination),xyz),8))
      enddo
      call check_projector(destination,destination,rotation,ns+1,p_error,e_error)
      if(p_error>projector_tol.or.e_error>energy_tol)cycle
      h_error=0d0
      if(op/=1)then
        do trial=1,2
          do site=1,ng
            probe(site)=cmplx(sin(dble(site*trial)*sqrt(2d0)),cos(dble(site*(trial+1))*sqrt(3d0)),8)
          enddo
          probe=probe/norm2(abs(probe))
          moved=phase(:,ns+1)*probe(grid(:,ns+1))
          call apply(probe,action,info);stars%h_probes=stars%h_probes+1
          if(info/=0.or..not.screen_finite(action))return
          call apply(moved,moved_action,info);stars%h_probes=stars%h_probes+1
          if(info/=0.or..not.screen_finite(moved_action))return
          scale=max(norm2(abs(action)),tiny(1d0))
          if(.not.screen_finite([scale]))return
          h_error=max(h_error,norm2(abs(moved_action-phase(:,ns+1)*action(grid(:,ns+1))))/scale)
          if(.not.screen_finite([h_error]))return
        enddo
      endif
      write(*,'(a,i3,3es15.6)')'HSEX_REDUCED_SYMMETRY_CANDIDATE permutation/H/projector/energy ',op,h_error,p_error,e_error
      if(h_error>covariance_tol)cycle
      ns=ns+1;rotations(:,:,ns)=rotation
      stars%max_h_covariance=max(stars%max_h_covariance,h_error)
      stars%max_projector_error=max(stars%max_projector_error,p_error)
      stars%max_energy_error=max(stars%max_energy_error,e_error)
    enddo
    if(ns<1)return
    stars%nsym=ns;stars%representative=0;stars%operation=0;stars%column=0
    do iq=1,nq
      stars%representative(iq)=iq;stars%operation(iq)=1
      do ig=1,counts(iq)
        stars%column(ig,iq)=ig
      enddo
      candidate_representatives: do jq=1,iq-1
        if(stars%representative(jq)/=jq)cycle
        if(counts(iq)/=counts(jq))cycle
        source=sources(jq);dest_source=sources(iq)
        do op=1,ns
          rotation=rotations(:,:,op)
          q_delta=matmul(rotation,q(:,jq))-q(:,iq)
          if(maxval(abs(q_delta-anint(q_delta)))>1d-10)cycle
          difference=matmul(rotation,frac(:,source))-frac(:,dest_source)
          if(maxval(abs(difference-anint(difference)))>1d-10)cycle
          call check_projector(source,dest_source,rotation,op,p_error,e_error)
          if(p_error>projector_tol.or.e_error>energy_tol)cycle
          valid=.true.;columns=0;delta=nint(q_delta)
          do ig=1,counts(jq)
            mapped=matmul(rotation,g(:,ig,jq))+delta;found=0
            do jg=1,counts(iq)
              if(all(g(:,jg,iq)==mapped))then
                found=jg
                exit
              endif
            enddo
            if(found==0)then
              valid=.false.
              exit
            endif
            columns(found)=ig
          enddo
          if(.not.valid.or.any(columns(:counts(iq))==0))cycle
          stars%representative(iq)=jq;stars%operation(iq)=op;stars%column(:,iq)=columns
          stars%max_projector_error=max(stars%max_projector_error,p_error)
          stars%max_energy_error=max(stars%max_energy_error,e_error)
          exit candidate_representatives
        enddo
      enddo candidate_representatives
      if(stars%representative(iq)==iq)stars%nrep=stars%nrep+1
      write(*,'(a,3i6)')'HSEX_REDUCED_SYMMETRY_STAR q/representative/operation ',iq, &
        stars%representative(iq),stars%operation(iq)
    enddo
    status=0
  contains
    subroutine check_projector(from,to,rot,operation_slot,projector_error,energy_error)
      implicit none
      integer,intent(in) :: from,to,rot(3,3),operation_slot
      real(8),intent(out) :: projector_error,energy_error
      integer :: point,band,left
      complex(8) :: factor
      real(8) :: position(3),fold(3)
      fold=matmul(rot,frac(:,from))-frac(:,to)
      do point=0,ng-1
        position=dble([mod(point,n(1)),mod(point/n(1),n(2)),point/(n(1)*n(2))])/n
        factor=exp(cmplx(0d0,2d0*acos(-1d0)*dot_product(fold,position),8))
        transformed(point+1,:)=factor*occupied(grid(point+1,operation_slot),:,from)
      enddo
      overlap=matmul(conjg(transpose(occupied(:,:,to))),transformed)
      transformed=transformed-matmul(occupied(:,:,to),overlap)
      projector_error=norm2(reshape(abs(transformed),[ng*nv]))/sqrt(dble(nv))
      do band=1,nv
        do left=1,nv
          overlap(left,band)=overlap(left,band)*(levels(left,to)-levels(band,from))
        enddo
      enddo
      energy_error=norm2(reshape(abs(overlap),[nv*nv]))/sqrt(dble(nv))
      if(.not.screen_finite([projector_error,energy_error]))then
        projector_error=huge(1d0);energy_error=huge(1d0)
      endif
    end subroutine
  end subroutine
end module
