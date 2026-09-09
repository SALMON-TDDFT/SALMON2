module dg_canonical_pp_fingerprint
  ! Canonical provenance for the supported periodic, scalar-relativistic DG Hybrid route.
  ! Hash only final production tables and their declared meaning domains; large reader
  ! workspaces and inactive padding are deliberately outside this schema.
  use iso_fortran_env,only:int64
  use structures,only:s_pp_info
  use dg_portable_sha256,only:s_dg_sha256_context,dg_sha256_schema,dg_sha256_init,&
    dg_sha256_update_int64,dg_sha256_update_integer,dg_sha256_update_real64,&
    dg_sha256_update_logical,dg_sha256_update_character,dg_sha256_final
  implicit none
  private
  integer(int64),parameter::canonical_pp_schema=1_int64
  public::canonical_pp_fingerprint,canonical_pp_digest,canonical_pp_valence_sum
contains
  pure integer(int64) function canonical_pp_fingerprint(pp) result(fingerprint)
    type(s_pp_info),intent(in)::pp
    integer::element,radial,channel,angular,nprojector,radial_extent,nlcc_end
    fingerprint=0_int64
    if(.not.valid_canonical_pp(pp))return
    fingerprint=int(z'6A09E667F3BCC909',int64)
    call mix_character(fingerprint,'SALMON-DG-CANONICAL-PP')
    call mix_int64(fingerprint,canonical_pp_schema)
    call mix_integer(fingerprint,size(pp%zps));call mix_integer(fingerprint,pp%lmax)
    call mix_integer(fingerprint,pp%nrmax);call mix_logical(fingerprint,pp%flag_nlcc)
    do element=1,size(pp%zps)
      call mix_character(fingerprint,pp%atom_symbol(element))
      call mix_integer(fingerprint,pp%zps(element));call mix_integer(fingerprint,pp%mr(element))
      call mix_integer(fingerprint,pp%lref(element));call mix_integer(fingerprint,pp%mlps(element))
      call mix_integer(fingerprint,pp%nrps(element));call mix_integer(fingerprint,pp%nrloc(element))
      call mix_real(fingerprint,pp%rps(element));call mix_real(fingerprint,pp%rloc(element))
      do angular=0,pp%mlps(element)
        call mix_integer(fingerprint,pp%nproj(angular,element))
      enddo
      nprojector=sum(pp%nproj(0:pp%mlps(element),element))
      do channel=0,nprojector-1
        call mix_real(fingerprint,pp%anorm(channel,element))
        call mix_integer(fingerprint,pp%inorm(channel,element))
        call mix_integer(fingerprint,pp%inorm_so(channel,element))
        if(pp%inorm_so(channel,element)/=0)&
          call mix_real(fingerprint,pp%anorm_so(channel,element))
      enddo
      nlcc_end=0;if(pp%flag_nlcc)nlcc_end=nlcc_meaning_end(pp,element)
      radial_extent=max(pp%mr(element)+1,pp%nrloc(element),pp%nrps(element),&
        pp%nrps_ao(element),nlcc_end)
      do radial=1,radial_extent
        call mix_real(fingerprint,pp%rad(radial,element))
      enddo
      do radial=1,pp%nrloc(element)
        call mix_real(fingerprint,pp%vloctbl(radial,element))
        call mix_real(fingerprint,pp%dvloctbl(radial,element))
      enddo
      do radial=1,pp%nrps(element)
        call mix_real(fingerprint,pp%radnl(radial,element))
        do channel=0,nprojector-1
          call mix_real(fingerprint,pp%udvtbl(radial,channel,element))
          call mix_real(fingerprint,pp%dudvtbl(radial,channel,element))
          if(pp%inorm_so(channel,element)/=0)then
            call mix_real(fingerprint,pp%udvtbl_so(radial,channel,element))
            call mix_real(fingerprint,pp%dudvtbl_so(radial,channel,element))
          endif
        enddo
      enddo
      call mix_integer(fingerprint,pp%nrps_ao(element))
      call mix_real(fingerprint,pp%rps_ao(element))
      do radial=1,pp%nrps_ao(element)
        do channel=0,nprojector-1
          call mix_real(fingerprint,pp%upptbl_ao(radial,channel,element))
        enddo
      enddo
      do radial=1,pp%mr(element)
        call mix_real(fingerprint,pp%rho_pp_tbl(radial,element))
      enddo
      if(pp%flag_nlcc)then
        do radial=1,nlcc_end
          call mix_real(fingerprint,pp%rho_nlcc_tbl(radial,element))
          call mix_real(fingerprint,pp%tau_nlcc_tbl(radial,element))
        enddo
      endif
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
  end function canonical_pp_fingerprint

  pure function canonical_pp_digest(pp) result(digest)
    type(s_pp_info),intent(in)::pp
    integer(int64)::digest(4)
    type(s_dg_sha256_context)::hash
    integer::element,radial,channel,angular,nprojector,radial_extent,nlcc_end
    digest=0_int64;if(.not.valid_canonical_pp(pp))return
    call dg_sha256_init(hash);call dg_sha256_update_character(hash,'SALMON-DG-CANONICAL-PP-AUTH')
    call dg_sha256_update_int64(hash,dg_sha256_schema)
    call dg_sha256_update_integer(hash,size(pp%zps));call dg_sha256_update_integer(hash,pp%lmax)
    call dg_sha256_update_integer(hash,pp%nrmax);call dg_sha256_update_logical(hash,pp%flag_nlcc)
    do element=1,size(pp%zps)
      call dg_sha256_update_character(hash,pp%atom_symbol(element))
      call dg_sha256_update_integer(hash,pp%zps(element));call dg_sha256_update_integer(hash,pp%mr(element))
      call dg_sha256_update_integer(hash,pp%lref(element));call dg_sha256_update_integer(hash,pp%mlps(element))
      call dg_sha256_update_integer(hash,pp%nrps(element));call dg_sha256_update_integer(hash,pp%nrloc(element))
      call dg_sha256_update_real64(hash,pp%rps(element));call dg_sha256_update_real64(hash,pp%rloc(element))
      do angular=0,pp%mlps(element);call dg_sha256_update_integer(hash,pp%nproj(angular,element));enddo
      nprojector=sum(pp%nproj(0:pp%mlps(element),element))
      do channel=0,nprojector-1
        call dg_sha256_update_real64(hash,pp%anorm(channel,element))
        call dg_sha256_update_integer(hash,pp%inorm(channel,element))
        call dg_sha256_update_integer(hash,pp%inorm_so(channel,element))
        if(pp%inorm_so(channel,element)/=0)call dg_sha256_update_real64(hash,pp%anorm_so(channel,element))
      enddo
      nlcc_end=0;if(pp%flag_nlcc)nlcc_end=nlcc_meaning_end(pp,element)
      radial_extent=max(pp%mr(element)+1,pp%nrloc(element),pp%nrps(element),pp%nrps_ao(element),nlcc_end)
      do radial=1,radial_extent;call dg_sha256_update_real64(hash,pp%rad(radial,element));enddo
      do radial=1,pp%nrloc(element)
        call dg_sha256_update_real64(hash,pp%vloctbl(radial,element))
        call dg_sha256_update_real64(hash,pp%dvloctbl(radial,element))
      enddo
      do radial=1,pp%nrps(element)
        call dg_sha256_update_real64(hash,pp%radnl(radial,element))
        do channel=0,nprojector-1
          call dg_sha256_update_real64(hash,pp%udvtbl(radial,channel,element))
          call dg_sha256_update_real64(hash,pp%dudvtbl(radial,channel,element))
          if(pp%inorm_so(channel,element)/=0)then
            call dg_sha256_update_real64(hash,pp%udvtbl_so(radial,channel,element))
            call dg_sha256_update_real64(hash,pp%dudvtbl_so(radial,channel,element))
          endif
        enddo
      enddo
      call dg_sha256_update_integer(hash,pp%nrps_ao(element));call dg_sha256_update_real64(hash,pp%rps_ao(element))
      do radial=1,pp%nrps_ao(element);do channel=0,nprojector-1
        call dg_sha256_update_real64(hash,pp%upptbl_ao(radial,channel,element))
      enddo;enddo
      do radial=1,pp%mr(element);call dg_sha256_update_real64(hash,pp%rho_pp_tbl(radial,element));enddo
      if(pp%flag_nlcc)then;do radial=1,nlcc_end
        call dg_sha256_update_real64(hash,pp%rho_nlcc_tbl(radial,element))
        call dg_sha256_update_real64(hash,pp%tau_nlcc_tbl(radial,element))
      enddo;endif
    enddo
    call dg_sha256_final(hash,digest)
  end function canonical_pp_digest

  pure real(8) function canonical_pp_valence_sum(pp) result(valence)
    type(s_pp_info),intent(in)::pp
    valence=0d0
    if(allocated(pp%zps))valence=sum(real(pp%zps,8))
  end function canonical_pp_valence_sum

  pure logical function valid_canonical_pp(pp) result(valid)
    type(s_pp_info),intent(in)::pp
    integer::element,nelement,nprojector,max_projector,max_radial,max_mr,nlcc_end
    valid=.false.
    if(pp%lmax<0.or.pp%nrmax<1)return
    if(.not.allocated(pp%atom_symbol).or..not.allocated(pp%mr).or..not.allocated(pp%lref).or.&
      .not.allocated(pp%nrps).or..not.allocated(pp%mlps).or..not.allocated(pp%nproj).or.&
      .not.allocated(pp%zps).or..not.allocated(pp%nrloc).or..not.allocated(pp%rloc).or.&
      .not.allocated(pp%rps).or..not.allocated(pp%anorm).or..not.allocated(pp%inorm).or.&
      .not.allocated(pp%anorm_so).or..not.allocated(pp%inorm_so).or.&
      .not.allocated(pp%rad).or..not.allocated(pp%radnl).or..not.allocated(pp%vloctbl).or.&
      .not.allocated(pp%dvloctbl).or..not.allocated(pp%udvtbl).or..not.allocated(pp%dudvtbl).or.&
      .not.allocated(pp%udvtbl_so).or..not.allocated(pp%dudvtbl_so).or.&
      .not.allocated(pp%rho_pp_tbl).or..not.allocated(pp%rho_nlcc_tbl).or.&
      .not.allocated(pp%tau_nlcc_tbl).or..not.allocated(pp%nrps_ao).or.&
      .not.allocated(pp%rps_ao).or..not.allocated(pp%upptbl_ao))return
    nelement=size(pp%zps);if(nelement<1)return
    if(size(pp%atom_symbol)<nelement.or.size(pp%mr)<nelement.or.size(pp%lref)<nelement.or.&
      size(pp%nrps)<nelement.or.size(pp%mlps)<nelement.or.size(pp%nrloc)<nelement.or.&
      size(pp%rloc)<nelement.or.size(pp%rps)<nelement.or.size(pp%nrps_ao)<nelement.or.&
      size(pp%rps_ao)<nelement)return
    if(lbound(pp%nproj,1)>0.or.ubound(pp%nproj,1)<pp%lmax.or.size(pp%nproj,2)<nelement)return
    if(pp%flag_nlcc)then
      if(lbound(pp%rho_nlcc_tbl,1)>1.or.ubound(pp%rho_nlcc_tbl,1)<pp%nrmax.or.&
        size(pp%rho_nlcc_tbl,2)<nelement)return
      if(lbound(pp%tau_nlcc_tbl,1)>1.or.ubound(pp%tau_nlcc_tbl,1)<pp%nrmax.or.&
        size(pp%tau_nlcc_tbl,2)<nelement)return
    endif
    max_projector=0;max_radial=0;max_mr=0
    do element=1,nelement
      if(pp%mlps(element)<0.or.pp%mlps(element)>pp%lmax.or.pp%mr(element)<1.or.&
        pp%nrps(element)<1.or.pp%nrloc(element)<1.or.pp%nrps_ao(element)<1)return
      if(pp%mr(element)+1>pp%nrmax.or.pp%nrps(element)>pp%nrmax.or.pp%nrloc(element)>pp%nrmax)return
      if(pp%nrps_ao(element)>pp%nrmax)return
      if(any(pp%nproj(0:pp%mlps(element),element)<0))return
      nprojector=sum(pp%nproj(0:pp%mlps(element),element));if(nprojector<1)return
      nlcc_end=0
      if(pp%flag_nlcc)then
        nlcc_end=nlcc_meaning_end(pp,element);if(nlcc_end<1)return
      endif
      max_projector=max(max_projector,nprojector)
      max_radial=max(max_radial,pp%mr(element)+1,pp%nrps(element),pp%nrloc(element),&
        pp%nrps_ao(element),nlcc_end)
      max_mr=max(max_mr,pp%mr(element))
    enddo
    if(lbound(pp%anorm,1)>0.or.ubound(pp%anorm,1)<max_projector-1.or.size(pp%anorm,2)<nelement)return
    if(lbound(pp%inorm,1)>0.or.ubound(pp%inorm,1)<max_projector-1.or.size(pp%inorm,2)<nelement)return
    if(lbound(pp%anorm_so,1)>0.or.ubound(pp%anorm_so,1)<max_projector-1.or.&
      size(pp%anorm_so,2)<nelement)return
    if(lbound(pp%inorm_so,1)>0.or.ubound(pp%inorm_so,1)<max_projector-1.or.&
      size(pp%inorm_so,2)<nelement)return
    if(lbound(pp%rad,1)>1.or.ubound(pp%rad,1)<max_radial.or.size(pp%rad,2)<nelement)return
    if(lbound(pp%radnl,1)>1.or.ubound(pp%radnl,1)<maxval(pp%nrps(1:nelement)).or.&
      size(pp%radnl,2)<nelement)return
    if(lbound(pp%vloctbl,1)>1.or.ubound(pp%vloctbl,1)<maxval(pp%nrloc(1:nelement)).or.&
      size(pp%vloctbl,2)<nelement)return
    if(lbound(pp%dvloctbl,1)>1.or.ubound(pp%dvloctbl,1)<maxval(pp%nrloc(1:nelement)).or.&
      size(pp%dvloctbl,2)<nelement)return
    if(lbound(pp%udvtbl,1)>1.or.ubound(pp%udvtbl,1)<maxval(pp%nrps(1:nelement)).or.&
      lbound(pp%udvtbl,2)>0.or.ubound(pp%udvtbl,2)<max_projector-1.or.size(pp%udvtbl,3)<nelement)return
    if(lbound(pp%dudvtbl,1)>1.or.ubound(pp%dudvtbl,1)<maxval(pp%nrps(1:nelement)).or.&
      lbound(pp%dudvtbl,2)>0.or.ubound(pp%dudvtbl,2)<max_projector-1.or.size(pp%dudvtbl,3)<nelement)return
    if(lbound(pp%udvtbl_so,1)>1.or.ubound(pp%udvtbl_so,1)<maxval(pp%nrps(1:nelement)).or.&
      lbound(pp%udvtbl_so,2)>0.or.ubound(pp%udvtbl_so,2)<max_projector-1.or.&
      size(pp%udvtbl_so,3)<nelement)return
    if(lbound(pp%dudvtbl_so,1)>1.or.ubound(pp%dudvtbl_so,1)<maxval(pp%nrps(1:nelement)).or.&
      lbound(pp%dudvtbl_so,2)>0.or.ubound(pp%dudvtbl_so,2)<max_projector-1.or.&
      size(pp%dudvtbl_so,3)<nelement)return
    if(lbound(pp%upptbl_ao,1)>1.or.ubound(pp%upptbl_ao,1)<maxval(pp%nrps_ao(1:nelement)).or.&
      lbound(pp%upptbl_ao,2)>0.or.ubound(pp%upptbl_ao,2)<max_projector-1.or.&
      size(pp%upptbl_ao,3)<nelement)return
    if(lbound(pp%rho_pp_tbl,1)>1.or.ubound(pp%rho_pp_tbl,1)<max_mr.or.size(pp%rho_pp_tbl,2)<nelement)return
    if(lbound(pp%rho_nlcc_tbl,1)>1.or.ubound(pp%rho_nlcc_tbl,1)<max_mr.or.&
      size(pp%rho_nlcc_tbl,2)<nelement)return
    if(lbound(pp%tau_nlcc_tbl,1)>1.or.ubound(pp%tau_nlcc_tbl,1)<max_mr.or.&
      size(pp%tau_nlcc_tbl,2)<nelement)return
    valid=.true.
  end function valid_canonical_pp

  pure integer function nlcc_meaning_end(pp,element) result(last)
    type(s_pp_info),intent(in)::pp
    integer,intent(in)::element
    integer::radial
    last=0
    do radial=1,pp%nrmax
      if(pp%rho_nlcc_tbl(radial,element)+pp%tau_nlcc_tbl(radial,element)<1d-6)then
        last=min(pp%nrmax,radial+1)
        return
      endif
    enddo
  end function nlcc_meaning_end

  pure subroutine mix_integer(hash,value)
    integer(int64),intent(inout)::hash
    integer,intent(in)::value
    call mix_int64(hash,int(value,int64))
  end subroutine mix_integer

  pure subroutine mix_int64(hash,value)
    integer(int64),intent(inout)::hash
    integer(int64),intent(in)::value
    integer::byte
    ! This is the bit-exact schema-1 transform.  ISHFTC and IEOR are defined
    ! bit operations; unlike signed integer multiplication they cannot overflow.
    do byte=0,7
      hash=ieor(ishftc(hash,7),int(ibits(value,8*byte,8),int64))
    enddo
  end subroutine mix_int64

  pure subroutine mix_real(hash,value)
    integer(int64),intent(inout)::hash
    real(8),intent(in)::value
    integer(int64)::bits
    bits=transfer(value,bits);call mix_int64(hash,bits)
  end subroutine mix_real

  pure subroutine mix_logical(hash,value)
    integer(int64),intent(inout)::hash
    logical,intent(in)::value
    call mix_integer(hash,merge(1,0,value))
  end subroutine mix_logical

  pure subroutine mix_character(hash,value)
    integer(int64),intent(inout)::hash
    character(*),intent(in)::value
    integer::position
    call mix_integer(hash,len(value))
    do position=1,len(value)
      hash=ieor(ishftc(hash,7),int(iachar(value(position:position)),int64))
    enddo
  end subroutine mix_character
end module dg_canonical_pp_fingerprint
