module dg_portable_sha256
  ! Portable SHA-256 compression for 16-byte chained receipts.  Every 32-bit
  ! word is held in int64 and explicitly masked, so no signed overflow occurs.
  use,intrinsic::iso_fortran_env,only:int64,real64
  implicit none
  private
  integer(int64),parameter::mask32=int(z'FFFFFFFF',int64)
  integer(int64),parameter,public::dg_sha256_schema=3_int64
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
  public::dg_sha256_mix_int64,dg_sha256_mix_integer,dg_sha256_mix_real64,&
    dg_sha256_mix_logical,dg_sha256_mix_character
contains
  pure subroutine dg_sha256_mix_int64(receipt,value)
    integer(int64),intent(inout)::receipt
    integer(int64),intent(in)::value
    integer(int64)::w(64),h(8),a,b,c,d,e,f,g,hh,t1,t2,bytes(16)
    integer::i,j
    h=[int(z'6A09E667',int64),int(z'BB67AE85',int64),int(z'3C6EF372',int64),&
      int(z'A54FF53A',int64),int(z'510E527F',int64),int(z'9B05688C',int64),&
      int(z'1F83D9AB',int64),int(z'5BE0CD19',int64)]
    do i=1,8
      bytes(i)=int(ibits(receipt,8*(i-1),8),int64)
      bytes(8+i)=int(ibits(value,8*(i-1),8),int64)
    enddo
    w=0_int64
    do i=1,4
      j=4*(i-1)
      w(i)=ior(shiftl(bytes(j+1),24),ior(shiftl(bytes(j+2),16),&
        ior(shiftl(bytes(j+3),8),bytes(j+4))))
    enddo
    w(5)=int(z'80000000',int64);w(16)=128_int64
    do i=17,64
      w(i)=low32(small_sigma1(w(i-2))+w(i-7)+small_sigma0(w(i-15))+w(i-16))
    enddo
    a=h(1);b=h(2);c=h(3);d=h(4);e=h(5);f=h(6);g=h(7);hh=h(8)
    do i=1,64
      t1=low32(hh+big_sigma1(e)+choose(e,f,g)+round_constant(i)+w(i))
      t2=low32(big_sigma0(a)+majority(a,b,c))
      hh=g;g=f;f=e;e=low32(d+t1);d=c;c=b;b=a;a=low32(t1+t2)
    enddo
    h=[low32(h(1)+a),low32(h(2)+b),low32(h(3)+c),low32(h(4)+d),&
      low32(h(5)+e),low32(h(6)+f),low32(h(7)+g),low32(h(8)+hh)]
    receipt=ior(shiftl(h(1),32),h(2))
  contains
    pure integer(int64) function low32(x)
      integer(int64),intent(in)::x
      low32=iand(x,mask32)
    end function low32
    pure integer(int64) function rotate_right(x,n)
      integer(int64),intent(in)::x
      integer,intent(in)::n
      rotate_right=iand(ior(shiftr(iand(x,mask32),n),shiftl(iand(x,mask32),32-n)),mask32)
    end function rotate_right
    pure integer(int64) function choose(x,y,z)
      integer(int64),intent(in)::x,y,z
      choose=iand(ieor(iand(x,y),iand(not(x),z)),mask32)
    end function choose
    pure integer(int64) function majority(x,y,z)
      integer(int64),intent(in)::x,y,z
      majority=iand(ieor(ieor(iand(x,y),iand(x,z)),iand(y,z)),mask32)
    end function majority
    pure integer(int64) function big_sigma0(x)
      integer(int64),intent(in)::x
      big_sigma0=ieor(ieor(rotate_right(x,2),rotate_right(x,13)),rotate_right(x,22))
    end function big_sigma0
    pure integer(int64) function big_sigma1(x)
      integer(int64),intent(in)::x
      big_sigma1=ieor(ieor(rotate_right(x,6),rotate_right(x,11)),rotate_right(x,25))
    end function big_sigma1
    pure integer(int64) function small_sigma0(x)
      integer(int64),intent(in)::x
      small_sigma0=ieor(ieor(rotate_right(x,7),rotate_right(x,18)),shiftr(iand(x,mask32),3))
    end function small_sigma0
    pure integer(int64) function small_sigma1(x)
      integer(int64),intent(in)::x
      small_sigma1=ieor(ieor(rotate_right(x,17),rotate_right(x,19)),shiftr(iand(x,mask32),10))
    end function small_sigma1
  end subroutine dg_sha256_mix_int64

  pure subroutine dg_sha256_mix_integer(receipt,value)
    integer(int64),intent(inout)::receipt
    integer,intent(in)::value
    call dg_sha256_mix_int64(receipt,int(value,int64))
  end subroutine dg_sha256_mix_integer
  pure subroutine dg_sha256_mix_real64(receipt,value)
    integer(int64),intent(inout)::receipt
    real(real64),intent(in)::value
    integer(int64)::bits
    bits=transfer(value,bits);call dg_sha256_mix_int64(receipt,bits)
  end subroutine dg_sha256_mix_real64
  pure subroutine dg_sha256_mix_logical(receipt,value)
    integer(int64),intent(inout)::receipt
    logical,intent(in)::value
    call dg_sha256_mix_int64(receipt,merge(1_int64,0_int64,value))
  end subroutine dg_sha256_mix_logical
  pure subroutine dg_sha256_mix_character(receipt,value)
    integer(int64),intent(inout)::receipt
    character(*),intent(in)::value
    integer::i
    call dg_sha256_mix_integer(receipt,len(value))
    do i=1,len(value);call dg_sha256_mix_integer(receipt,iachar(value(i:i)));enddo
  end subroutine dg_sha256_mix_character
end module dg_portable_sha256
