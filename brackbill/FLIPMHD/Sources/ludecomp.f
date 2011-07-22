      subroutine ludecomp(ncells,
     &     a11,a12,a13,a14,
     &     a21,a22,a23,a24,
     &     a31,a32,a33,a34,
     &     a41,a42,a43,a44)
c
      real
     &     a11(*),a12(*),a13(*),a14(*),
     &     a21(*),a22(*),a23(*),a24(*),
     &     a31(*),a32(*),a33(*),a34(*),
     &     a41(*),a42(*),a43(*),a44(*)
c
c
      do n=1,ncells
c
c
      a11(n)=1./(a11(n)+1.e-30)
c
      a21(n)=a21(n)*a11(n)
      a22(n)=a22(n)-a21(n)*a12(n)
      a23(n)=a23(n)-a21(n)*a13(n)
      a24(n)=a24(n)-a21(n)*a14(n)
c
      a31(n)=a31(n)*a11(n)
      b32=a32(n)-a31(n)*a12(n)
      b33=a33(n)-a31(n)*a13(n)
      b34=a34(n)-a31(n)*a14(n)
c
      a41(n)=a41(n)*a11(n)
      b42=a42(n)-a41(n)*a12(n)
      b43=a43(n)-a41(n)*a13(n)
      b44=a44(n)-a41(n)*a14(n)
c
      a22(n)=1./(a22(n)+1.e-30)
c
      a32(n)=b32*a22(n)
      a33(n)=b33-a32(n)*a23(n)
      a34(n)=b34-a32(n)*a24(n)
c
      a42(n)=b42*a22(n)
      c43=b43-a42(n)*a23(n)
      c44=b44-a42(n)*a24(n)
c
      a33(n)=1./(a33(n)+1.e-30)
c
      a43(n)=c43*a33(n)
      a44(n)=c44-a43(n)*a34(n)
c
      a44(n)=1./(a44(n)+1.e-30)
c
c
      enddo
c
c
      return
      end
