module dg_portable_sha256
  ! Incremental SHA-256 over canonical bytes.  32-bit words are represented by
  ! int64 values and every arithmetic result is masked, avoiding signed overflow.
  use,intrinsic::iso_fortran_env,only:int64,real64
  implicit none
  private
  integer(int64),parameter::mask32=int(z'FFFFFFFF',int64)
  integer(int64),parameter,public::dg_sha256_schema=4_int64
  integer(int64),parameter::initial_state(8)=[int(z'6A09E667',int64),int(z'BB67AE85',int64),&
    int(z'3C6EF372',int64),int(z'A54FF53A',int64),int(z'510E527F',int64),&
    int(z'9B05688C',int64),int(z'1F83D9AB',int64),int(z'5BE0CD19',int64)]
  integer(int64),parameter::round_constant(64)=[&
    int(z'428A2F98',int64),int(z'71374491',int64),int(z'B5C0FBCF',int64),int(z'E9B5DBA5',int64),&
    int(z'3956C25B',int64),int(z'59F111F1',int64),int(z'923F82A4',int64),int(z'AB1C5ED5',int64),&
    int(z'D807AA98',int64),int(z'12835B01',int64),int(z'243185BE',int64),int(z'550C7DC3',int64),&
    int(z'72BE5D74',int64),int(z'80DEB1FE',int64),int(z'9BDC06A7',int64),int(z'C19BF174',int64),&
    int(z'E49B69C1',int64),int(z'EFBE4786',int64),int(z'0FC19DC6',int64),int(z'240CA1CC',int64),&
    int(z'2DE92C6F',int64),int(z'4A7484AA',int64),int(z'5CB0A9DC',int64),int(z'76F988DA',int64),&
    int(z'983E5152',int64),int(z'A831C66D',int64),int(z'B00327C8',int64),int(z'BF597FC7',int64),&
    int(z'C6E00BF3',int64),int(z'D5A79147',int64),int(z'06CA6351',int64),int(z'14292967',int64),&
    int(z'27B70A85',int64),int(z'2E1B2138',int64),int(z'4D2C6DFC',int64),int(z'53380D13',int64),&
    int(z'650A7354',int64),int(z'766A0ABB',int64),int(z'81C2C92E',int64),int(z'92722C85',int64),&
    int(z'A2BFE8A1',int64),int(z'A81A664B',int64),int(z'C24B8B70',int64),int(z'C76C51A3',int64),&
    int(z'D192E819',int64),int(z'D6990624',int64),int(z'F40E3585',int64),int(z'106AA070',int64),&
    int(z'19A4C116',int64),int(z'1E376C08',int64),int(z'2748774C',int64),int(z'34B0BCB5',int64),&
    int(z'391C0CB3',int64),int(z'4ED8AA4A',int64),int(z'5B9CCA4F',int64),int(z'682E6FF3',int64),&
    int(z'748F82EE',int64),int(z'78A5636F',int64),int(z'84C87814',int64),int(z'8CC70208',int64),&
    int(z'90BEFFFA',int64),int(z'A4506CEB',int64),int(z'BEF9A3F7',int64),int(z'C67178F2',int64)]
  type,public::s_dg_sha256_context
    integer(int64)::state(8)=initial_state
    integer(int64)::buffer(64)=0_int64,total_bytes=0_int64
    integer::buffered=0
  end type s_dg_sha256_context
  public::dg_sha256_init,dg_sha256_update_bytes,dg_sha256_update_int64,&
    dg_sha256_update_integer,dg_sha256_update_real64,dg_sha256_update_logical,&
    dg_sha256_update_character,dg_sha256_final
contains
  pure subroutine dg_sha256_init(context)
    type(s_dg_sha256_context),intent(out)::context
    context%state=initial_state;context%buffer=0_int64
    context%total_bytes=0_int64;context%buffered=0
  end subroutine
  pure subroutine dg_sha256_update_bytes(context,bytes)
    type(s_dg_sha256_context),intent(inout)::context
    integer(int64),intent(in)::bytes(:)
    integer::i
    do i=1,size(bytes)
      context%buffered=context%buffered+1;context%buffer(context%buffered)=iand(bytes(i),255_int64)
      context%total_bytes=context%total_bytes+1_int64
      if(context%buffered==64)then;call compress(context%state,context%buffer);context%buffered=0;endif
    enddo
  end subroutine
  pure subroutine dg_sha256_update_int64(context,value)
    type(s_dg_sha256_context),intent(inout)::context;integer(int64),intent(in)::value
    integer(int64)::bytes(8);integer::i
    do i=1,8;bytes(i)=ibits(value,8*(i-1),8);enddo
    call dg_sha256_update_bytes(context,bytes)
  end subroutine
  pure subroutine dg_sha256_update_integer(context,value)
    type(s_dg_sha256_context),intent(inout)::context;integer,intent(in)::value
    call dg_sha256_update_int64(context,int(value,int64))
  end subroutine
  pure subroutine dg_sha256_update_real64(context,value)
    type(s_dg_sha256_context),intent(inout)::context;real(real64),intent(in)::value
    integer(int64)::bits;bits=transfer(value,bits);call dg_sha256_update_int64(context,bits)
  end subroutine
  pure subroutine dg_sha256_update_logical(context,value)
    type(s_dg_sha256_context),intent(inout)::context;logical,intent(in)::value
    call dg_sha256_update_int64(context,merge(1_int64,0_int64,value))
  end subroutine
  pure subroutine dg_sha256_update_character(context,value)
    type(s_dg_sha256_context),intent(inout)::context;character(*),intent(in)::value
    integer(int64),allocatable::bytes(:);integer::i
    call dg_sha256_update_integer(context,len(value));allocate(bytes(len(value)))
    do i=1,len(value);bytes(i)=iachar(value(i:i));enddo
    call dg_sha256_update_bytes(context,bytes)
  end subroutine
  pure subroutine dg_sha256_final(context,digest)
    type(s_dg_sha256_context),intent(in)::context;integer(int64),intent(out)::digest(4)
    type(s_dg_sha256_context)::work;integer(int64)::length_bits;integer::i
    work=context;length_bits=shiftl(context%total_bytes,3);call append(work,128_int64)
    do while(work%buffered/=56);call append(work,0_int64);enddo
    do i=7,0,-1;call append(work,ibits(length_bits,8*i,8));enddo
    do i=1,4;digest(i)=ior(shiftl(work%state(2*i-1),32),work%state(2*i));enddo
  contains
    pure subroutine append(target,value)
      type(s_dg_sha256_context),intent(inout)::target;integer(int64),intent(in)::value
      target%buffered=target%buffered+1;target%buffer(target%buffered)=value
      if(target%buffered==64)then;call compress(target%state,target%buffer);target%buffered=0;endif
    end subroutine
  end subroutine
  pure subroutine compress(state,block)
    integer(int64),intent(inout)::state(8);integer(int64),intent(in)::block(64)
    integer(int64)::w(64),a,b,c,d,e,f,g,h,t1,t2;integer::i,j
    do i=1,16
      j=4*(i-1);w(i)=ior(shiftl(block(j+1),24),ior(shiftl(block(j+2),16),&
        ior(shiftl(block(j+3),8),block(j+4))))
    enddo
    do i=17,64;w(i)=low32(sigma1(w(i-2))+w(i-7)+sigma0(w(i-15))+w(i-16));enddo
    a=state(1);b=state(2);c=state(3);d=state(4);e=state(5);f=state(6);g=state(7);h=state(8)
    do i=1,64
      t1=low32(h+capsigma1(e)+choose(e,f,g)+round_constant(i)+w(i));t2=low32(capsigma0(a)+majority(a,b,c))
      h=g;g=f;f=e;e=low32(d+t1);d=c;c=b;b=a;a=low32(t1+t2)
    enddo
    state=[low32(state(1)+a),low32(state(2)+b),low32(state(3)+c),low32(state(4)+d),&
      low32(state(5)+e),low32(state(6)+f),low32(state(7)+g),low32(state(8)+h)]
  end subroutine
  pure integer(int64) function low32(x);integer(int64),intent(in)::x;low32=iand(x,mask32);end function
  pure integer(int64) function rotr(x,n);integer(int64),intent(in)::x;integer,intent(in)::n
    rotr=iand(ior(shiftr(iand(x,mask32),n),shiftl(iand(x,mask32),32-n)),mask32);end function
  pure integer(int64) function choose(x,y,z);integer(int64),intent(in)::x,y,z
    choose=iand(ieor(iand(x,y),iand(not(x),z)),mask32);end function
  pure integer(int64) function majority(x,y,z);integer(int64),intent(in)::x,y,z
    majority=iand(ieor(ieor(iand(x,y),iand(x,z)),iand(y,z)),mask32);end function
  pure integer(int64) function capsigma0(x);integer(int64),intent(in)::x
    capsigma0=ieor(ieor(rotr(x,2),rotr(x,13)),rotr(x,22));end function
  pure integer(int64) function capsigma1(x);integer(int64),intent(in)::x
    capsigma1=ieor(ieor(rotr(x,6),rotr(x,11)),rotr(x,25));end function
  pure integer(int64) function sigma0(x);integer(int64),intent(in)::x
    sigma0=ieor(ieor(rotr(x,7),rotr(x,18)),shiftr(iand(x,mask32),3));end function
  pure integer(int64) function sigma1(x);integer(int64),intent(in)::x
    sigma1=ieor(ieor(rotr(x,17),rotr(x,19)),shiftr(iand(x,mask32),10));end function
end module dg_portable_sha256
