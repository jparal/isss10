      module blcom_com_M
      use vast_kind_param, only:  double
      use corgan_com_M
!     TYPE DECLARATIONS
!
      character(len=12), dimension(nspecies) :: SpeciesName
!
      integer, dimension(nreg) :: npcelx,npcely,npcelz,icoi
      integer, dimension(0:nreg) :: numtot
!
      real(double) ::      &
             fibar , fjbar , fkbar , rfibar, rfjbar, rfkbar, rbjbkb

!
      real(double) ::&
          eps,resistivity,&
          modex(nreg), modey(nreg),     &
          uvi(nreg),vvi(nreg),wvi(nreg),utordrft(nreg),&
          xcenter(nreg),ycenter(nreg),zcenter(nreg),&
          rex(nreg),rey(nreg),rez(nreg),rix(nreg),riy(nreg),riz(nreg),&
          siep(nreg), rhr(nreg),qom(nreg),&
          bxr(nreg),byr(nreg),bzr(nreg),&
          xvi(8,nreg),yvi(8,nreg),zvi(8,nreg)
!
!
      real(double) ::&
          wate(itdim,29),wght(8),number(itdim),&
          sie(itdim),rho(itdim),pl(itdim),&
          numberv(itdim),color(itdim)
!
!
      real(double) ::&
          u(itdim),v(itdim),w(itdim),&
          ul(itdim),vl(itdim),wl(itdim),&
          xptilde(itdim),yptilde(itdim),zptilde(itdim),&
          axv1(itdim),ayv1(itdim),azv1(itdim),&
          mc(itdim),sie1p(itdim),&
          ax(itdim),ay(itdim),az(itdim),&
          mv(itdim),umom(itdim),vmom(itdim),wmom(itdim),&
          work(itdim),pdv(itdim),csq(itdim),&
          vasq(itdim),&
          divu(itdim),gradx(itdim),grady(itdim),gradz(itdim),&
          exx(itdim),exy(itdim),exz(itdim),&
          eyy(itdim),eyz(itdim),ezz(itdim)
      integer :: imp
!
      real(double) ::&
          dudx(itdim),dudy(itdim),dudz(itdim),&
          dvdx(itdim),dvdy(itdim),dvdz(itdim),&
          dwdx(itdim),dwdy(itdim),dwdz(itdim)
!
      real(double) ::&
          pixx(itdim),pixy(itdim),pixz(itdim),&
          piyy(itdim),piyz(itdim),pizz(itdim),&
          divpix(itdim),divpiy(itdim),divpiz(itdim)
!
      real(double) ::&
          delta,&
          courant,lagrange,&
          bxn(itdim),byn(itdim),bzn(itdim),&
          bxl(itdim),byl(itdim),bzl(itdim),&
          bmagx(itdim),bmagy(itdim),bmagz(itdim),&
          divuphix(itdim),divuphiy(itdim),divuphiz(itdim),  &
          dbvx(itdim),dbvy(itdim),dbvz(itdim),              &
          jx(itdim),jy(itdim),jz(itdim),&
          ex(itdim),ey(itdim),ez(itdim),&
          EdotJ(itdim),&
          CurlEx(itdim),CurlEy(itdim),CurlEz(itdim),&
          Ohmic_heating(itdim),&
          bxv(itdim),byv(itdim),bzv(itdim),&
          p(itdim)
!
      integer ::&
          itmax,itmag,numit,iter_pois
!
      real(double) ::&
          error,&
          a11(itdim),a12(itdim),a13(itdim),a14(itdim),&
          a21(itdim),a22(itdim),a23(itdim),a24(itdim),&
          a31(itdim),a32(itdim),a33(itdim),a34(itdim),&
          a41(itdim),a42(itdim),a43(itdim),a44(itdim)

      real(double) ::&
          efnrg(nhst),ebnrg(nhst),eknrg(nhst),&
          charge(nhst),thistry(nhst),reconnected_flux(0:nhst)

!
      integer ::&
          iphead(itdim),ijkcell(itdim),ijkvtx(itdim),&
          ijktmp2(itdim),ijktmp3(itdim),ijktmp4(itdim),&
          ijkctmp(itdim),istep(idzg)
!
      logical :: adaptg,explicit,fields,cartesian,&
          periodic_x,periodic_y,periodic_z,resist
!
      real(double) ::&
          lama,lamb,lamc,&
          xl,xr,yb,yt,ze,zf,&
          x(itdim),y(itdim),z(itdim),&
          vol(itdim),vvol(itdim),&
          tsix(itdim),tsiy(itdim),tsiz(itdim),&
          etax(itdim),etay(itdim),etaz(itdim),&
          nux(itdim),nuy(itdim),nuz(itdim),&
          area_x(itdim),area_y(itdim),area_z(itdim)
!
!      real(double) ::
!           c1x(itdim),c2x(itdim),c3x(itdim),c4x(itdim),
!            c5x(itdim),c6x(itdim),c7x(itdim),c8x(itdim),
!            c1y(itdim),c2y(itdim),c3y(itdim),c4y(itdim),
!            c5y(itdim),c6y(itdim),c7y(itdim),c8y(itdim),
!            c1z(itdim),c2z(itdim),c3z(itdim),c4z(itdim),
!            c5z(itdim),c6z(itdim),c7z(itdim),c8z(itdim)

      real(double) ::&
          wgrid(itdim)
!
      integer ::&
          nptotl,npsampl,&
          npsav(nlist),&
          ico(0:npart),&
          link(0:npart),&
          iphd2(itdim)

      integer, dimension(0:npart) :: ijkParticle
!
      real(double) ::&
          px(0:npart),py(0:npart),pz(0:npart),&
          pxi(0:npart),peta(0:npart),pzta(0:npart),&
          up(0:npart),vp(0:npart),wp(0:npart),&
          ep(0:npart),&
          mass(0:npart),&
          mupx(0:npart),mupy(0:npart),mupz(0:npart),&
          uptilde(itdim),vptilde(itdim),wptilde(itdim)
!
      end module blcom_com_M
