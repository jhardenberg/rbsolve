#include "config.h"
c ______________________________________________________________________

      subroutine FOURIER1_TMPT (u,uo)

#include "param.h"
      integer*4 indice,nw2,ik1,ik2
      real*8 uo(Nx,Nz,0:Nyl+1),a(Nx+2,Nz+2)
      complex*16 u(nwax+1,2*nwaz+1,0:Nyl+1)

      do j=Nylmin0,Nylmax0

        do k=1,Nz+2
	  do i=1,Nx+2
	    a(i,k)=0.d0
          end do
        end do

c - coeff. di Fourier C_k  (k=0, k=1,nwaz) -
        do k=1,nwaz+1
          indice=0
	  do i=1,nwax+1
            indice=indice+2
            a(indice-1,k)=dreal(u(i,k,j))
            a(indice,k)=dimag(u(i,k,j))
          end do
        end do

c - coeff. di Fourier C_(Nz-k)  (k=1,nwaz) -
	nw2=2*nwaz+2
        do k=1,nwaz
          indice=0
          ik1=Nz-(k-1)
          ik2=nw2-k
	  do i=1,nwax+1
            indice=indice+2
            a(indice-1,ik1)=dreal(u(i,ik2,j))
            a(indice,ik1)=dimag(u(i,ik2,j))
          end do
        enddo

        call rffi(a,Nx,Nz)
     
        do k=1,Nz
          do i=1,Nx
            uo(i,k,j)=a(i,k)
          end do
        end do

      end do

      return
      

      end      


      subroutine INITFOUR_TMPT (iest)
#include "param.h"
      integer*4 ifax1(13),ifax2(13),ierror
      integer iest
      real*8 ex1(2*(Nx+2)),ex2(2*(Nz+2)),work(2*(Nx+2)*(Nz+2)) 
      common /work/work
      common /ifax1/ifax1
      common /ifax2/ifax2
      common /ex1/ex1
      common /ex2/ex2
 
      call rffp(ex1,ex2,ifax1,ifax2,Nx,Nz)
 
      if((ifax1(1).eq.-99).or.(ifax1(2).eq.-99)) goto 991

      return

 991  ierror=1
      write(6,9910) ierror
 9910 format('IERROR=',I2,'Nx or Nz CONTAIN FACTORS OTHER THAN 2,3,5')

      end



      subroutine rffp(ex1,ex2,ifax1,ifax2,n1,n2)
******---------------***************************
      implicit double precision(a-h,o-z)
      call rfftst(ex1,ifax1,n1)
      call cfftst(ex2,ifax2,n2)
      return
      end
      subroutine rffd(u,n1,n2)
******---------------***************************
#include "param.h"
      integer*4 ifax1(13),ifax2(13)
      real*8 ex1(2*(Nx+2)),ex2(2*(Nz+2)),work(2*(Nx+2)*(Nz+2)) 
      real*8 u(*)
      integer*4 n1,n2
      common /work/work
      common /ifax1/ifax1
      common /ifax2/ifax2
      common /ex1/ex1
      common /ex2/ex2
      call rfftsb(u,work,ex1,ifax1,1,n1+2,n1,n2,-1)
      call cfftsb(u,work,ex2,ifax2,(n1+2)/2,1,n2,(n1+2)/2,-1)
      return
      end
      subroutine rffi(u,n1,n2)
******---------------***************************
#include "param.h"
      integer*4 ifax1(13),ifax2(13)
      real*8 ex1(2*(Nx+2)),ex2(2*(Nz+2)),work(2*(Nx+2)*(Nz+2)) 
      real*8 u(*)
      integer*4 n1,n2
      common /work/work
      common /ifax1/ifax1
      common /ifax2/ifax2
      common /ex1/ex1
      common /ex2/ex2

      call cfftsb(u,work,ex2,ifax2,(n1+2)/2,1,n2,(n1+2)/2,1)
      call rfftsb(u,work,ex1,ifax1,1,n1+2,n1,n2,1)
      return
      end
      
************---------------------------------------********************
*     SBFOR : Any machine using a Fortran Version of FFT991 and CFFT99
******-----------------------------------------------------------******
*       Standard Basic subroutines version 2.01  18-01-90
*
*  1. Content :
*
*            rfftst : initialisations for rfftsb
*            rfftsb : real Fast Fourier transform
*
*            cfftst : initialisations for cfftsb
*            cfftsb : complex fast Fourier Transform
*
*            rcosst : initialisations for rcossb
*            rcossb : real fast cosine transform
*
*            ccosst : initialisations for ccossb
*            ccossb : complex fast cosine transform
*
*            rsinst : initialisations for rsinsb
*            rsinsb : real  fast sine transform
*
*            rsymst : initialisations for rsymsb
*            rsymsb : real symetric fast Fourier transform
*
*            csymst : initialisations for csymsb
*            csymsb : complex symetric fast Fourier transform
*
*            rantst : initialisations for rantsb
*            rantsb : real antisymetric fast Fourier transform
*
*  2. Internal subroutine used :
*       the package Temperton FFT991 and CFFT99 for real  and complex ffts
*
*
*
*
*  3. NB : in this version priviledging  real ftt fft991 it is choosen:
*       rfftsb uses the real Temperton fft991
*       cfftsb uses the complex Temperton cfft99
*       rcossb uses rsymsb. (Thual)
*       ccossb uses rcossb.
*       rsinsb uses rantsb. (Thual)
*       rsymsb uses rfftsb. (Brachet)
*       rantb uses  rfftsb. (Brachet)
*
*
***********************************************************************
      subroutine rfftst(ex,ifax,n)
******-----------------***********
      implicit double precision(a-h,o-z)
      call fftfax(n,ifax,ex)
      return
      end
      subroutine rfftsb(u,work,ex,ifax,inc,jump,n,nft,isign)
******-----------------*************************************
      implicit double precision(a-h,o-z)
      real*8 u(*)
      if (isign.eq.-1)  then
      do 10 i=1,nft
      i1=(i-1)*jump+n*inc+inc
      i2=i1+inc
      u(i1)=0.d0
      u(i2)=0.d0
 10   continue
      endif
      call fft991(u,work,ex,ifax,inc,jump,n,nft,isign)
      return
      end
      subroutine cfftst(ex,ifax,n)
******-----------------***********
      implicit double precision(a-h,o-z)
      call cftfax(n,ifax,ex)
      return
      end

      subroutine cfftsb(u,work,ex,ifax,inc,jump,n,nft,isign)
******-----------------***********
      implicit double precision(a-h,o-z)
      call cfft99(u,work,ex,ifax,inc,jump,n,nft,isign)
      if(isign.eq.-1) then
      call ftnorm(u,n,inc,jump,nft,isign)
      endif
      return
      end
      subroutine rcosst (ex,ifax,ssini,ssinq,npts)
******-----------------******************************
      implicit double precision(a-h,o-z)
c initialisation pour rsymsb
      dimension ssini(*),ssinq(*)
      call fftfax(npts,ifax,ex)
      pi=  4.*atan(1.)
      pih= pi/2.
      do 10 i=2,npts
      ssini(i)= .25/ sin(pi*float(i-1)/float(npts))
      ssinq(i)=      sin(pi*float(i-1)/float(npts))
  10  continue
      ssinq(1)=0.
      ssinq(npts+1)=0.
      return
      end
      subroutine rcossb(a,work,ex,ifax,ssini,
     1                 ssinq,inc,jump,npts,nft,isign)
******-----------------******************************
      implicit double precision(a-h,o-z)
      real*8 a(*),work(*),ssini(*),ssinq(*)
c physique spectral isign=1 :  ui=sum an cos nxi vers an
c spectral physique isign=-1:  an vers ui
*  physique spectral isign=-1
      if (isign.eq.-1) then
      call rsymsb(a,work,ex,ifax,ssini,ssinq,inc,jump,npts,nft,isign)
      do 100 i=2,npts
CDIR$ IVDEP
      do 11 iv=1,nft
      ivec=(iv-1)*jump
      ia=(i-1)*inc+1+ivec
      a(ia)=2.*a(ia)
  11  continue
 100  continue
      return
      endif
* inverse spectral physique isign=1
      if(isign.eq.1) then
      do 200 i=2,npts
CDIR$ IVDEP
      do 21 iv=1,nft
      ivec=(iv-1)*jump
      ia=(i-1)*inc+1+ivec
      a(ia)=.5*a(ia)
  21  continue
 200  continue
      call rsymsb(a,work,ex,ifax,ssini,ssinq,inc,jump,npts,nft,isign)
      return
      endif
      end
      subroutine ccosst (ex,ifax,ssini,ssinq,npts)
******-----------------******************************
      implicit double precision(a-h,o-z)
      dimension ssini(*),ssinq(*)
      call fftfax(npts,ifax,ex)
      pi=  4.*atan(1.)
      pih= pi/2.
      do 10 i=2,npts
      ssini(i)= .25/ sin(pi*float(i-1)/float(npts))
      ssinq(i)=      sin(pi*float(i-1)/float(npts))
  10  continue
      ssinq(1)=0.
      ssinq(npts+1)=0.
      return
      end
      subroutine ccossb(a,work,ex,ifax,ssini,
     1                 ssinq,inc,jump,npts,nft,isign)
******-----------------******************************
      implicit double precision(a-h,o-z)
      real*8 a(*),work(*),ssini(*)
      real*8 ssinq(*)
      call rcossb(a(1),work,ex,ifax,ssini,ssinq,
     1                2*inc,2*jump,npts,nft,isign)
      call rcossb(a(2),work,ex,ifax,ssini,ssinq,
     1                2*inc,2*jump,npts,nft,isign)
      return
      end
      subroutine rsinst (ex,ifax,ssini,ssinq,npts)
******-----------------******************************
      implicit double precision(a-h,o-z)
c initialisation pour rsymsb
      dimension ssini(*),ssinq(*)
      call fftfax(npts,ifax,ex)
      pi=  4.*atan(1.)
      pih= pi/2.
      do 10 i=2,npts
      ssini(i)= .25/ sin(pi*float(i-1)/float(npts))
      ssinq(i)=      sin(pi*float(i-1)/float(npts))
  10  continue
      ssinq(1)=0.
      ssinq(npts+1)=0.
      return
      end
      subroutine rsinsb(a,work,ex,ifax,ssini,
     1                 ssinq,inc,jump,npts,nft,isign)
******-----------------******************************
      implicit double precision(a-h,o-z)
      real*8 a(*),work(*),ssini(*),ssinq(*)
c physique spectral isign=1 :  ui=sum an sin nxi vers an
c spectral physique isign=-1:  an vers ui
*  physique spectral isign=-1
      if (isign.eq.-1) then
      call rantsb(a,work,ex,ifax,ssini,ssinq,inc,jump,npts,nft,isign)
      do 100 i=2,npts-1
CDIR$ IVDEP
      do 11 iv=1,nft
      ivec=(iv-1)*jump
      ia=(i-1)*inc+1+ivec
      a(ia)=2.*a(ia)
  11  continue
 100  continue
      return
      endif
* inverse spectral physique isign=1
      if(isign.eq.1) then
      do 200 i=2,npts-1
CDIR$ IVDEP
      do 21 iv=1,nft
      ivec=(iv-1)*jump
      ia=(i-1)*inc+1+ivec
      a(ia)=.5*a(ia)
  21  continue
 200  continue
      call rantsb(a,work,ex,ifax,ssini,ssinq,inc,jump,npts,nft,isign)
      return
      endif
      end
      subroutine rsymst (ex,ifax,ssini,ssinq,npts)
******-----------------******************************
      implicit double precision(a-h,o-z)
      dimension ssini(*),ssinq(*)
      call fftfax(npts,ifax,ex)
      pi=  4.*atan(1.)
      pih= pi/2.
      do 10 i=2,npts
      ssini(i)= .25/ sin(pi*float(i-1)/float(npts))
      ssinq(i)=      sin(pi*float(i-1)/float(npts))
  10  continue
      ssinq(1)=0.
      ssinq(npts+1)=0.
      return
      end
      subroutine rsymsb(a,work,ex,ifax,ssini,
     1                 ssinq,inc,jump,npts,nft,isign)
******-----------------******************************
      implicit double precision(a-h,o-z)
      dimension a(*),work(*),ssini(*),ssinq(*)
c transformation cosinus reel basee sur temperton reel complex
c attention! ne marche que si la dimension du champ >=npts+3
      if (isign.eq.-1) goto 1000
c cas spectral --> physique
CDIR$ IVDEP
      do 10 iv=1,nft
      iv2=iv+nft
      ivec=(iv-1)*jump
      ia=1+inc+ivec
      work(iv)=a(ia)
      work(iv2)=a(ia)
      a(ia)=0.
  10  continue
      do 100 i=4,npts,2
CDIR$ IVDEP
      do 11 iv=1,nft
      iv2=iv+nft
      iv3=iv2+nft
      ivec=(iv-1)*jump
      ia=(i-1)*inc+1+ivec
      work(iv)=work(iv)+a(ia)
      work(iv3)=a(ia)
      a(ia)=a(ia)-work(iv2)
      work(iv2)=work(iv3)
  11  continue
 100  continue
CDIR$ IVDEP
      do 12 iv=1,nft
      ivec=(iv-1)*jump
      a((npts+1)*inc+ivec+1)=0.
      a((npts+2)*inc+ivec+1)=2.*work(iv)
  12  continue
c
      call fft991(a,work,ex,ifax,inc,jump,npts,nft,isign)
c
CDIR$ IVDEP
      do 13 iv=1,nft
      ivec=(iv-1)*jump
      a((npts)*inc+ivec+1)=a(1+ivec)-a((npts+2)*inc+ivec+1)
      a(1+ivec)=a(1+ivec)+a((npts+2)*inc+ivec+1)
  13  continue
      do 110 i=2,npts/2
      it =(i-1)*inc+1
      isym=(npts+1-i)*inc+1
CDIR$ IVDEP
      do 14 iv=1,nft
      iv2=iv+nft
      ivec=(iv-1)*jump
      work(iv)=.5*(a(it+ivec)+a(isym+ivec))
      work(iv2)=ssini(i)*(a(it+ivec)-a(isym+ivec))
      a(it+ivec)=work(iv)+work(iv2)
      a(isym+ivec)=work(iv)-work(iv2)
  14  continue
 110  continue
CDIR$ IVDEP
      do 15 iv=1,nft
      ivec=(iv-1)*jump
      a((npts+1)*inc+ivec+1)=0.
      a((npts+2)*inc+ivec+1)=0.
  15  continue
      return
 1000 continue
c cas physique --> spectral
CDIR$ IVDEP
      do 23 iv=1,nft
      ivec=(iv-1)*jump
      a((npts+2)*inc+ivec+1)=.25*(a(1+ivec)-a(npts*inc+ivec+1))
      a(1+ivec)=.5*(a(1+ivec)+a(npts*inc+ivec+1))
   23 continue
      do 210 i=2,npts/2
      it=(i-1)*inc+1
      isym=(npts+1-i)*inc+1
CDIR$ IVDEP
      do 24 iv=1,nft
      iv2=iv+nft
      ivec=(iv-1)*jump
      work(iv)=.5*(a(it+ivec)+a(isym+ivec))
      work(iv2)=ssinq(i)*(a(it+ivec)-a(isym+ivec))
      a(it+ivec)=work(iv)+work(iv2)
      a(isym+ivec)=work(iv)-work(iv2)
  24  continue
 210  continue
c
      call fft991(a,work,ex,ifax,inc,jump,npts,nft,isign)
c
CDIR$ IVDEP
      do 20 iv=1,nft
      iv2=iv+nft
      ivec=(iv-1)*jump
      ia=(npts-1)*inc+ivec+1
      work(iv)=-a(ia)
      work(iv2)=-a(ia)
  20  continue
      do 200 i=(npts-2),4,-2
CDIR$ IVDEP
      do 21 iv=1,nft
      ivec=(iv-1)*jump
      iv2=iv+nft
      ia=(i-1)*inc+1+ivec
      work(iv)=work(iv)-a(ia)
      work(iv2)=work(iv2)+work(iv)
  21  continue
 200  continue
CDIR$ IVDEP
      do 22 iv=1,nft
      iv2=iv+nft
      ivec=(iv-1)*jump
      work(iv)=(a(inc*(npts+2)+ivec+1)-work(iv2))/float(npts/2)
  22  continue
      do 205 i=npts,4,-2
CDIR$ IVDEP
      do 25 iv=1,nft
      iv2=iv+nft
      ivec=(iv-1)*jump
      ia=(i-1)*inc+ivec+1
      work(iv2)=a(ia)
      a(ia)=work(iv)
      work(iv)=work(iv)-work(iv2)
  25  continue
 205  continue
CDIR$ IVDEP
      do 26,iv=1,nft
      ivec=(iv-1)*jump
      ia=inc+ivec+1
      a(ia)=work(iv)
      a((npts+1)*inc+ivec+1)=0.
      a((npts+2)*inc+ivec+1)=0.
  26  continue
      return
      end
      subroutine csymst(ex,ifax,ssini,ssinq,npts)
******-----------------******************************
      implicit double precision(a-h,o-z)
c initialisation pour rsymsb
      dimension ssini(*),ssinq(*)
      call fftfax(npts,ifax,ex)
      pi=  4.*atan(1.)
      pih= pi/2.
      do 10 i=2,npts
      ssini(i)= .25/ sin(pi*float(i-1)/float(npts))
      ssinq(i)=      sin(pi*float(i-1)/float(npts))
  10  continue
      ssinq(1)=0.
      ssinq(npts+1)=0.
      return
      end
      subroutine csymsb(a,work,ex,ifax,ssini,
     1                 ssinq,inc,jump,m,nft,isign)
******-----------------***************************
      implicit double precision(a-h,o-z)
      real*8 a(*),work(*),ssini(*),ssinq(*)
      call rsymsb(a(1),work,ex,ifax,ssini,
     1           ssinq,2*inc,2*jump,npts,nft,isign)
      call rsymsb(a(2),work,ex,ifax,ssini,
     1           ssinq,2*inc,2*jump,npts,nft,isign)
      return
      end
      subroutine rantst (ex,ifax,ssini,ssinq,npts)
******-----------------******************************
      implicit double precision(a-h,o-z)
c initialisation pour rsymsb
      dimension ssini(*),ssinq(*)
      call fftfax(npts,ifax,ex)
      pi=  4.*atan(1.)
      pih= pi/2.
      do 10 i=2,npts
      ssini(i)= .25/ sin(pi*float(i-1)/float(npts))
      ssinq(i)=      sin(pi*float(i-1)/float(npts))
  10  continue
      ssinq(1)=0.
      ssinq(npts+1)=0.
      return
      end
      subroutine rantsb(a,work,ex,ifax,ssini,
     1                 ssinq,inc,jump,npts,nft,isign)
******-----------------******************************
      implicit double precision(a-h,o-z)
      dimension a(*),work(*),ssini(*),ssinq(*)
c transformation sinus reel basee sur temperton reel complex
c attention! ne marche que si  la dimension du champ>=npts+2
c (comme pour temperton reel complexe)
      if (isign.eq.-1) goto 1000
c cas spectral --> physique
CDIR$ IVDEP
      do 10 iv=1,nft
      iv2=iv+nft
      ivec=(iv-1)*jump
      ia=1+inc+ivec
      work(iv2)=a(ia)
      a(ia)=0.
      a(1+ivec)=2.*work(iv2)
  10  continue
      do 100 i=3,npts-1,2
CDIR$ IVDEP
      do 11 iv=1,nft
      iv2=iv+nft
      ivec=(iv-1)*jump
      ia=(i-1)*inc+1+ivec
      work(iv)=a(ia)
      a(ia)=a(ia+inc)-work(iv2)
      work(iv2)=a(ia+inc)
      a(ia+inc)=-work(iv)
  11  continue
 100  continue
CDIR$ IVDEP
      do 12 iv=1,nft
      iv2=iv+nft
      ivec=(iv-1)*jump
      a(npts*inc+ivec+1)=-2.*work(iv2)
      a((npts+1)*inc+ivec+1)=0.
  12  continue
c
      call fft991(a,work,ex,ifax,inc,jump,npts,nft,isign)
c
CDIR$ IVDEP
      do 13 iv=1,nft
      ivec=(iv-1)*jump
      a(npts*inc+ivec+1)=0.
      a(1+ivec)=0.
  13  continue
      do 110 i=2,npts/2
      it =(i-1)*inc+1
      isym=(npts+1-i)*inc+1
CDIR$ IVDEP
      do 14 iv=1,nft
      iv2=iv+nft
      ivec=(iv-1)*jump
      work(iv)=.5*(a(it+ivec)-a(isym+ivec))
      work(iv2)=ssini(i)*(a(it+ivec)+a(isym+ivec))
      a(it+ivec)=work(iv)+work(iv2)
      a(isym+ivec)=-work(iv)+work(iv2)
  14  continue
 110  continue
CDIR$ IVDEP
      do 15 iv=1,nft
      ivec=(iv-1)*jump
      a((npts+1)*inc+ivec+1)=0.
      a(inc*npts/2+ivec+1)=.5*a(inc*npts/2+ivec+1)
  15  continue
      return
 1000 continue
c cas physique --> spectral
      do 23 iv=1,nft
      ivec=(iv-1)*jump
      a(1+ivec)=0.
   23 continue
      do 210 i=2,npts/2
      it=(i-1)*inc+1
      isym=(npts+1-i)*inc+1
CDIR$ IVDEP
      do 24 iv=1,nft
      iv2=iv+nft
      ivec=(iv-1)*jump
      work(iv)=.5*(a(it+ivec)-a(isym+ivec))
      work(iv2)=ssinq(i)*(a(it+ivec)+a(isym+ivec))
      a(it+ivec)=work(iv)+work(iv2)
      a(isym+ivec)=-work(iv)+work(iv2)
  24  continue
 210  continue
      do 27 iv=1,nft
      ivec=(iv-1)*jump
      ia=inc*npts/2+ivec+1
      a(ia)=2.*a(ia)
  27  continue
c
      call fft991(a,work,ex,ifax,inc,jump,npts,nft,isign)
c
CDIR$ IVDEP
      do 20 iv=1,nft
      iv2=iv+nft
      ivec=(iv-1)*jump
      ia=1+inc+ivec
      a(ia)=.5*a(1+ivec)
      a(1+ivec)=0.
  20  continue
      do 205 i=3,npts-1,2
CDIR$ IVDEP
      do 25 iv=1,nft
      iv2=iv+nft
      ivec=(iv-1)*jump
      ia=(i-1)*inc+1+ivec
      work(iv)=a(ia+inc)
      a(ia+inc)=a(ia)+a(ia-inc)
      a(ia)=-work(iv)
  25  continue
 205  continue
CDIR$ IVDEP
      do 26 iv=1,nft
      ivec=(iv-1)*jump
      a(npts*inc+ivec+1)=0.
      a((npts+1)*inc+ivec+1)=0.
  26  continue
      return
      end
      subroutine ftnorm(a,n,inc,jump,nft,isign)
******-----------------***********
      implicit double precision(a-h,o-z)
      common/normal/normal
c ?      real a(2,0:n)
      dimension a(2,0:n)
      if(isign.eq.-1)then
      f=1.d0/dfloat(n)
      do 1 j=0,jump*(nft-1),jump
      do 1 k=j,j+(n-1)*inc,inc
      a(1,k)=a(1,k)*f
1     a(2,k)=a(2,k)*f
      endif
      return
      end
******************************************************
* now is the fortran source of  cray FFT991 and CFFT99
* got from Ncar xmp 18-01-90
******************************************************
      SUBROUTINE FFT99(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
      implicit double precision(a-h,o-z)
C
C PURPOSE      PERFORMS MULTIPLE FAST FOURIER TRANSFORMS.  THIS PACKAGE
C              WILL PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX
C              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
C              TRANSFORMS, I.E.  GIVEN A SET OF REAL DATA VECTORS, THE
C              PACKAGE RETURNS A SET OF 'HALF-COMPLEX' FOURIER
C              COEFFICIENT VECTORS, OR VICE VERSA.  THE LENGTH OF THE
C              TRANSFORMS MUST BE AN EVEN NUMBER GREATER THAN 4 THAT HAS
C              NO OTHER FACTORS EXCEPT POSSIBLY POWERS OF 2, 3, AND 5.
C              THIS IS AN ALL FORTRAN VERSION OF THE CRAYLIB PACKAGE
C              THAT IS MOSTLY WRITTEN IN CAL.
C
C              THE PACKAGE FFT99F CONTAINS SEVERAL USER-LEVEL ROUTINES:
C
C            SUBROUTINE FFTFAX
C                AN INITIALIZATION ROUTINE THAT MUST BE CALLED ONCE
C                BEFORE A SEQUENCE OF CALLS TO THE FFT ROUTINES
C                (PROVIDED THAT N IS NOT CHANGED).
C
C            SUBROUTINES FFT99 AND FFT991
C                TWO FFT ROUTINES THAT RETURN SLIGHTLY DIFFERENT
C                ARRANGEMENTS OF THE DATA IN GRIDPOINT SPACE.
C
C
C ACCESS       THIS FORTRAN VERSION MAY BE ACCESSED WITH
C
C                   *FORTRAN,P=XLIB,SN=FFT99F
C
C              TO ACCESS THE CRAY OBJECT CODE, CALLING THE USER ENTRY
C              POINTS FROM A CRAY PROGRAM IS SUFFICIENT.  THE SOURCE
C              FORTRAN AND CAL CODE FOR THE CRAYLIB VERSION MAY BE
C              ACCESSED USING
C
C                   FETCH P=CRAYLIB,SN=FFT99
C                   FETCH P=CRAYLIB,SN=CAL99
C
C USAGE        LET N BE OF THE FORM 2**P * 3**Q * 5**R, WHERE P .GE. 1,
C              Q .GE. 0, AND R .GE. 0.  THEN A TYPICAL SEQUENCE OF
C              CALLS TO TRANSFORM A GIVEN SET OF REAL VECTORS OF LENGTH
C              N TO A SET OF 'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS
C              OF LENGTH N IS
C
C                   DIMENSION IFAX(13),TRIGS(3*N/2+1),A(M*(N+2)),
C                  +          WORK(M*(N+1))
C
C                   CALL FFTFAX (N, IFAX, TRIGS)
C                   CALL FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
C
C              SEE THE INDIVIDUAL WRITE-UPS FOR FFTFAX, FFT99, AND
C              FFT991 BELOW, FOR A DETAILED DESCRIPTION OF THE
C              ARGUMENTS.
C
C HISTORY      THE PACKAGE WAS WRITTEN BY CLIVE TEMPERTON AT ECMWF IN
C              NOVEMBER, 1978.  IT WAS MODIFIED, DOCUMENTED, AND TESTED
C              FOR NCAR BY RUSS REW IN SEPTEMBER, 1980.
C
C-----------------------------------------------------------------------
C
C SUBROUTINE FFTFAX (N,IFAX,TRIGS)
C
C PURPOSE      A SET-UP ROUTINE FOR FFT99 AND FFT991.  IT NEED ONLY BE
C              CALLED ONCE BEFORE A SEQUENCE OF CALLS TO THE FFT
C              ROUTINES (PROVIDED THAT N IS NOT CHANGED).
C
C ARGUMENT     IFAX(13),TRIGS(3*N/2+1)
C DIMENSIONS
C
C ARGUMENTS
C
C ON INPUT     N
C               AN EVEN NUMBER GREATER THAN 4 THAT HAS NO PRIME FACTOR
C               GREATER THAN 5.  N IS THE LENGTH OF THE TRANSFORMS (SEE
C               THE DOCUMENTATION FOR FFT99 AND FFT991 FOR THE
C               DEFINITIONS OF THE TRANSFORMS).
C
C              IFAX
C               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED
C               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING
C               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN A MILLION.
C
C              TRIGS
C               A FLOATING POINT ARRAY OF DIMENSION 3*N/2 IF N/2 IS
C               EVEN, OR 3*N/2+1 IF N/2 IS ODD.
C
C ON OUTPUT    IFAX
C               CONTAINS THE FACTORIZATION OF N/2.  IFAX(1) IS THE
C               NUMBER OF FACTORS, AND THE FACTORS THEMSELVES ARE STORED
C               IN IFAX(2),IFAX(3),...  IF FFTFAX IS CALLED WITH N ODD,
C               OR IF N HAS ANY PRIME FACTORS GREATER THAN 5, IFAX(1)
C               IS SET TO -99.
C
C              TRIGS
C               AN ARRAY OF TRIGNOMENTRIC FUNCTION VALUES SUBSEQUENTLY
C               USED BY THE FFT ROUTINES.
C
C-----------------------------------------------------------------------
C
C SUBROUTINE FFT991 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
C                       AND
C SUBROUTINE FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
C
C PURPOSE      PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX
C              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
C              TRANSFORMS, USING ORDINARY SPATIAL ORDER OF GRIDPOINT
C              VALUES (FFT991) OR EXPLICIT CYCLIC CONTINUITY IN THE
C              GRIDPOINT VALUES (FFT99).  GIVEN A SET
C              OF REAL DATA VECTORS, THE PACKAGE RETURNS A SET OF
C              'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS, OR VICE
C              VERSA.  THE LENGTH OF THE TRANSFORMS MUST BE AN EVEN
C              NUMBER THAT HAS NO OTHER FACTORS EXCEPT POSSIBLY POWERS
C              OF 2, 3, AND 5.  THESE VERSION OF FFT991 AND FFT99 ARE
C              OPTIMIZED FOR USE ON THE CRAY-1.
C
C ARGUMENT     A(M*(N+2)), WORK(M*(N+1)), TRIGS(3*N/2+1), IFAX(13)
C DIMENSIONS
C
C ARGUMENTS
C
C ON INPUT     A
C               AN ARRAY OF LENGTH M*(N+2) CONTAINING THE INPUT DATA
C               OR COEFFICIENT VECTORS.  THIS ARRAY IS OVERWRITTEN BY
C               THE RESULTS.
C
C              WORK
C               A WORK ARRAY OF DIMENSION M*(N+1)
C
C              TRIGS
C               AN ARRAY SET UP BY FFTFAX, WHICH MUST BE CALLED FIRST.
C
C              IFAX
C               AN ARRAY SET UP BY FFTFAX, WHICH MUST BE CALLED FIRST.
C
C              INC
C               THE INCREMENT (IN WORDS) BETWEEN SUCCESSIVE ELEMENTS OF
C               EACH DATA OR COEFFICIENT VECTOR (E.G.  INC=1 FOR
C               CONSECUTIVELY STORED DATA).
C
C              JUMP
C               THE INCREMENT (IN WORDS) BETWEEN THE FIRST ELEMENTS OF
C               SUCCESSIVE DATA OR COEFFICIENT VECTORS.  ON THE CRAY-1,
C               TRY TO ARRANGE DATA SO THAT JUMP IS NOT A MULTIPLE OF 8
C               (TO AVOID MEMORY BANK CONFLICTS).  FOR CLARIFICATION OF
C               INC AND JUMP, SEE THE EXAMPLES BELOW.
C
C              N
C               THE LENGTH OF EACH TRANSFORM (SEE DEFINITION OF
C               TRANSFORMS, BELOW).
C
C              M
C               THE NUMBER OF TRANSFORMS TO BE DONE SIMULTANEOUSLY.
C
C              ISIGN
C               = +1 FOR A TRANSFORM FROM FOURIER COEFFICIENTS TO
C                    GRIDPOINT VALUES.
C               = -1 FOR A TRANSFORM FROM GRIDPOINT VALUES TO FOURIER
C                    COEFFICIENTS.
C
C ON OUTPUT    A
C               IF ISIGN = +1, AND M COEFFICIENT VECTORS ARE SUPPLIED
C               EACH CONTAINING THE SEQUENCE:
C
C               A(0),B(0),A(1),B(1),...,A(N/2),B(N/2)  (N+2 VALUES)
C
C               THEN THE RESULT CONSISTS OF M DATA VECTORS EACH
C               CONTAINING THE CORRESPONDING N+2 GRIDPOINT VALUES:
C
C               FOR FFT991, X(0), X(1), X(2),...,X(N-1),0,0.
C               FOR FFT99, X(N-1),X(0),X(1),X(2),...,X(N-1),X(0).
C                   (EXPLICIT CYCLIC CONTINUITY)
C
C               WHEN ISIGN = +1, THE TRANSFORM IS DEFINED BY:
C                 X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
C                 WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
C                 AND I=SQRT (-1)
C
C               IF ISIGN = -1, AND M DATA VECTORS ARE SUPPLIED EACH
C               CONTAINING A SEQUENCE OF GRIDPOINT VALUES X(J) AS
C               DEFINED ABOVE, THEN THE RESULT CONSISTS OF M VECTORS
C               EACH CONTAINING THE CORRESPONDING FOURIER COFFICIENTS
C               A(K), B(K), 0 .LE. K .LE N/2.
C
C               WHEN ISIGN = -1, THE INVERSE TRANSFORM IS DEFINED BY:
C                 C(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*EXP(-2*I*J*K*PI/N))
C                 WHERE C(K)=A(K)+I*B(K) AND I=SQRT(-1)
C
C               A CALL WITH ISIGN=+1 FOLLOWED BY A CALL WITH ISIGN=-1
C               (OR VICE VERSA) RETURNS THE ORIGINAL DATA.
C
C               NOTE: THE FACT THAT THE GRIDPOINT VALUES X(J) ARE REAL
C               IMPLIES THAT B(0)=B(N/2)=0.  FOR A CALL WITH ISIGN=+1,
C               IT IS NOT ACTUALLY NECESSARY TO SUPPLY THESE ZEROS.
C
C EXAMPLES      GIVEN 19 DATA VECTORS EACH OF LENGTH 64 (+2 FOR EXPLICIT
C               CYCLIC CONTINUITY), COMPUTE THE CORRESPONDING VECTORS OF
C               FOURIER COEFFICIENTS.  THE DATA MAY, FOR EXAMPLE, BE
C               ARRANGED LIKE THIS:
C
C FIRST DATA   A(1)=    . . .                A(66)=             A(70)
C VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS)
C
C SECOND DATA  A(71)=   . . .                                  A(140)
C VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS)
C
C               AND SO ON.  HERE INC=1, JUMP=70, N=64, M=19, ISIGN=-1,
C               AND FFT99 SHOULD BE USED (BECAUSE OF THE EXPLICIT CYCLIC
C               CONTINUITY).
C
C               ALTERNATIVELY THE DATA MAY BE ARRANGED LIKE THIS:
C
C                FIRST         SECOND                          LAST
C                DATA          DATA                            DATA
C                VECTOR        VECTOR                          VECTOR
C
C                 A(1)=         A(2)=                           A(19)=
C
C                 X(63)         X(63)       . . .               X(63)
C        A(20)=   X(0)          X(0)        . . .               X(0)
C        A(39)=   X(1)          X(1)        . . .               X(1)
C                  .             .                               .
C                  .             .                               .
C                  .             .                               .
C
C               IN WHICH CASE WE HAVE INC=19, JUMP=1, AND THE REMAINING
C               PARAMETERS ARE THE SAME AS BEFORE.  IN EITHER CASE, EACH
C               COEFFICIENT VECTOR OVERWRITES THE CORRESPONDING INPUT
C               DATA VECTOR.
C
C-----------------------------------------------------------------------
      DIMENSION A(N),WORK(N),TRIGS(N),IFAX(1)
C
C     SUBROUTINE "FFT99" - MULTIPLE FAST REAL PERIODIC TRANSFORM
C     CORRESPONDING TO OLD SCALAR ROUTINE FFT9
C     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
C     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
C     (1970), 315-337)
C
C     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
C     WORK IS AN AREA OF SIZE (N+1)*LOT
C     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
C     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
C     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
C         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
C     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
C     N IS THE LENGTH OF THE DATA VECTORS
C     LOT IS THE NUMBER OF DATA VECTORS
C     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
C           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
C
C     ORDERING OF COEFFICIENTS:
C         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
C         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
C
C     ORDERING OF DATA:
C         X(N-1),X(0),X(1),X(2),...,X(N),X(0)
C         I.E. EXPLICIT CYCLIC CONTINUITY; (N+2) LOCATIONS REQUIRED
C
C     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
C     PARALLEL
C
C     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
C
C     DEFINITION OF TRANSFORMS:
C     -------------------------
C
C     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
C         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
C
C     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
C               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
C
C
C
C
      NFAX=IFAX(1)
      NX=N+1
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30
C
C     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=INC+1
      JBASE=1
      DO 20 L=1,LOT
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 10 M=1,N
      WORK(J)=A(I)
      I=I+INC
      J=J+1
   10 CONTINUE
      IBASE=IBASE+JUMP
      JBASE=JBASE+NX
   20 CONTINUE
C
      IGO=60
      GO TO 40
C
C     PREPROCESSING (ISIGN=+1)
C     ------------------------
C
   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60
C
C     COMPLEX TRANSFORM
C     -----------------
C
   40 CONTINUE
      IA=INC+1
      LA=1
      DO 80 K=1,NFAX
      IF (IGO.EQ.60) GO TO 60
   50 CONTINUE
      CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS,
     *   INK,2,JUMP,NX,LOT,NH,IFAX(K+1),LA)
      IGO=60
      GO TO 70
   60 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS,
     *    2,INK,NX,JUMP,LOT,NH,IFAX(K+1),LA)
      IGO=50
   70 CONTINUE
      LA=LA*IFAX(K+1)
   80 CONTINUE
C
      IF (ISIGN.EQ.-1) GO TO 130
C
C     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=IA
      DO 100 L=1,LOT
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 90 M=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
   90 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  100 CONTINUE
C
C     FILL IN CYCLIC BOUNDARY POINTS
  110 CONTINUE
      IA=1
      IB=N*INC+1
CDIR$ IVDEP
      DO 120 L=1,LOT
      A(IA)=A(IB)
      A(IB+INC)=A(IA+INC)
      IA=IA+JUMP
      IB=IB+JUMP
  120 CONTINUE
      GO TO 140
C
C     POSTPROCESSING (ISIGN=-1):
C     --------------------------
C
  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
C
  140 CONTINUE
      RETURN
      END
      SUBROUTINE FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      implicit double precision(a-h,o-z)
      DIMENSION A(N),WORK(N),TRIGS(N)
C
C     SUBROUTINE FFT99A - PREPROCESSING STEP FOR FFT99, ISIGN=+1
C     (SPECTRAL TO GRIDPOINT TRANSFORM)
C
      NH=N/2
      NX=N+1
      INK=INC+INC
C
C     A(0) AND A(N/2)
      IA=1
      IB=N*INC+1
      JA=1
      JB=2
CDIR$ IVDEP
      DO 10 L=1,LOT
      WORK(JA)=A(IA)+A(IB)
      WORK(JB)=A(IA)-A(IB)
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   10 CONTINUE
C
C     REMAINING WAVENUMBERS
      IABASE=2*INC+1
      IBBASE=(N-2)*INC+1
      JABASE=3
      JBBASE=N-1
C
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
CDIR$ IVDEP
      DO 20 L=1,LOT
      WORK(JA)=(A(IA)+A(IB))-
     *    (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JB)=(A(IA)+A(IB))+
     *    (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JA+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))+
     *    (A(IA+INC)-A(IB+INC))
      WORK(JB+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))-
     *    (A(IA+INC)-A(IB+INC))
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   20 CONTINUE
      IABASE=IABASE+INK
      IBBASE=IBBASE-INK
      JABASE=JABASE+2
      JBBASE=JBBASE-2
   30 CONTINUE
C
      IF (IABASE.NE.IBBASE) GO TO 50
C     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
CDIR$ IVDEP
      DO 40 L=1,LOT
      WORK(JA)=2.d0*A(IA)
      WORK(JA+1)=-2.d0*A(IA+INC)
      IA=IA+JUMP
      JA=JA+NX
   40 CONTINUE
C
   50 CONTINUE
      RETURN
      END
      SUBROUTINE FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
      implicit double precision(a-h,o-z)
      DIMENSION WORK(N),A(N),TRIGS(N)
C
C     SUBROUTINE FFT99B - POSTPROCESSING STEP FOR FFT99, ISIGN=-1
C     (GRIDPOINT TO SPECTRAL TRANSFORM)
C
      NH=N/2
      NX=N+1
      INK=INC+INC
C
C     A(0) AND A(N/2)
      SCALE=1.d0/DFLOAT(N)
      IA=1
      IB=2
      JA=1
      JB=N*INC+1
CDIR$ IVDEP
      DO 10 L=1,LOT
      A(JA)=SCALE*(WORK(IA)+WORK(IB))
      A(JB)=SCALE*(WORK(IA)-WORK(IB))
      A(JA+INC)=0.d0
      A(JB+INC)=0.d0
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   10 CONTINUE
C
C     REMAINING WAVENUMBERS
      SCALE=0.5d0*SCALE
      IABASE=3
      IBBASE=N-1
      JABASE=2*INC+1
      JBBASE=(N-2)*INC+1
C
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
CDIR$ IVDEP
      DO 20 L=1,LOT
      A(JA)=SCALE*((WORK(IA)+WORK(IB))
     *   +(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JB)=SCALE*((WORK(IA)+WORK(IB))
     *   -(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JA+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))
     *    +(WORK(IB+1)-WORK(IA+1)))
      A(JB+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))
     *    -(WORK(IB+1)-WORK(IA+1)))
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   20 CONTINUE
      IABASE=IABASE+2
      IBBASE=IBBASE-2
      JABASE=JABASE+INK
      JBBASE=JBBASE-INK
   30 CONTINUE
C
      IF (IABASE.NE.IBBASE) GO TO 50
C     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
      SCALE=2.d0*SCALE
CDIR$ IVDEP
      DO 40 L=1,LOT
      A(JA)=SCALE*WORK(IA)
      A(JA+INC)=-SCALE*WORK(IA+1)
      IA=IA+NX
      JA=JA+JUMP
   40 CONTINUE
C
   50 CONTINUE
      RETURN
      END
      SUBROUTINE FFT991(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
      implicit double precision(a-h,o-z)
      DIMENSION A(N),WORK(N),TRIGS(N),IFAX(1)
C
C     SUBROUTINE "FFT991" - MULTIPLE REAL/HALF-COMPLEX PERIODIC
C     FAST FOURIER TRANSFORM
C
C     SAME AS FFT99 EXCEPT THAT ORDERING OF DATA CORRESPONDS TO
C     THAT IN MRFFT2
C
C     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
C     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
C     (1970), 315-337)
C
C     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
C     WORK IS AN AREA OF SIZE (N+1)*LOT
C     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
C     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
C     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
C         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
C     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
C     N IS THE LENGTH OF THE DATA VECTORS
C     LOT IS THE NUMBER OF DATA VECTORS
C     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
C           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
C
C     ORDERING OF COEFFICIENTS:
C         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
C         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
C
C     ORDERING OF DATA:
C         X(0),X(1),X(2),...,X(N-1)
C
C     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
C     PARALLEL
C
C     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
C
C     DEFINITION OF TRANSFORMS:
C     -------------------------
C
C     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
C         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
C
C     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
C               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
C
C
C
      NFAX=IFAX(1)
      NX=N+1
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30
C
C     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=1
      JBASE=1
      DO 20 L=1,LOT
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 10 M=1,N
      WORK(J)=A(I)
      I=I+INC
      J=J+1
   10 CONTINUE
      IBASE=IBASE+JUMP
      JBASE=JBASE+NX
   20 CONTINUE
C
      IGO=60
      GO TO 40
C
C     PREPROCESSING (ISIGN=+1)
C     ------------------------
C
   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60
C
C     COMPLEX TRANSFORM
C     -----------------
C
   40 CONTINUE
      IA=1
      LA=1
      DO 80 K=1,NFAX
      IF (IGO.EQ.60) GO TO 60
   50 CONTINUE
      CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS,
     *   INK,2,JUMP,NX,LOT,NH,IFAX(K+1),LA)
      IGO=60
      GO TO 70
   60 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS,
     *    2,INK,NX,JUMP,LOT,NH,IFAX(K+1),LA)
      IGO=50
   70 CONTINUE
      LA=LA*IFAX(K+1)
   80 CONTINUE
C
      IF (ISIGN.EQ.-1) GO TO 130
C
C     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=1
      DO 100 L=1,LOT
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 90 M=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
   90 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  100 CONTINUE
C
C     FILL IN ZEROS AT END
  110 CONTINUE
      IB=N*INC+1
CDIR$ IVDEP
      DO 120 L=1,LOT
      A(IB)=0.d0
      A(IB+INC)=0.d0
      IB=IB+JUMP
  120 CONTINUE
      GO TO 140
C
C     POSTPROCESSING (ISIGN=-1):
C     --------------------------
C
  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
C
  140 CONTINUE
      RETURN
      END
      SUBROUTINE FFTFAX(N,IFAX,TRIGS)
      implicit double precision(a-h,o-z)
      DIMENSION IFAX(13),TRIGS(1)
C
C MODE 3 IS USED FOR REAL/HALF-COMPLEX TRANSFORMS.  IT IS POSSIBLE
C TO DO COMPLEX/COMPLEX TRANSFORMS WITH OTHER VALUES OF MODE, BUT
C DOCUMENTATION OF THE DETAILS WERE NOT AVAILABLE WHEN THIS ROUTINE
C WAS WRITTEN.
C
      DATA MODE /3/
      CALL FAX (IFAX, N, MODE)
      I = IFAX(1)
      IF (IFAX(I+1) .GT. 5 .OR. N .LE. 4) IFAX(1) = -99
      IF (IFAX(1) .LE. 0 )CALL ULIBER
      CALL FFTRIG (TRIGS, N, MODE)
      RETURN
      END
      SUBROUTINE FAX(IFAX,N,MODE)
      implicit double precision(a-h,o-z)
      DIMENSION IFAX(10)
      NN=N
      IF (IABS(MODE).EQ.1) GO TO 10
      IF (IABS(MODE).EQ.8) GO TO 10
      NN=N/2
      IF ((NN+NN).EQ.N) GO TO 10
      IFAX(1)=-99
      RETURN
   10 K=1
C     TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
      K=K+1
      IFAX(K)=4
      NN=NN/4
      IF (NN.EQ.1) GO TO 80
      GO TO 20
C     TEST FOR EXTRA FACTOR OF 2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
C     TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40
C     NOW FIND REMAINING FACTORS
   50 L=5
      INC=2
C     INC ALTERNATELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L
      IF (NN.EQ.1) GO TO 80
      GO TO 60
   70 L=L+INC
      INC=6-INC
      GO TO 60
   80 IFAX(1)=K-1
C     IFAX(1) CONTAINS NUMBER OF FACTORS
      NFAX=IFAX(1)
C     SORT FACTORS INTO ASCENDING ORDER
      IF (NFAX.EQ.1) GO TO 110
      DO 100 II=2,NFAX
      ISTOP=NFAX+2-II
      DO 90 I=2,ISTOP
      IF (IFAX(I+1).GE.IFAX(I)) GO TO 90
      ITEM=IFAX(I)
      IFAX(I)=IFAX(I+1)
      IFAX(I+1)=ITEM
   90 CONTINUE
  100 CONTINUE
  110 CONTINUE
      RETURN
      END
      SUBROUTINE FFTRIG(TRIGS,N,MODE)
      implicit double precision(a-h,o-z)
      DIMENSION TRIGS(1)
      PI=2.d0*ASIN(1.d0)
      IMODE=IABS(MODE)
      NN=N
      IF (IMODE.GT.1.AND.IMODE.LT.6) NN=N/2
      DEL=(PI+PI)/FLOAT(NN)
      L=NN+NN
      DO 10 I=1,L,2
      ANGLE=0.5d0*DFLOAT(I-1)*DEL
      TRIGS(I)=DCOS(ANGLE)
      TRIGS(I+1)=DSIN(ANGLE)
   10 CONTINUE
      IF (IMODE.EQ.1) RETURN
      IF (IMODE.EQ.8) RETURN
      DEL=0.5d0*DEL
      NH=(NN+1)/2
      L=NH+NH
      LA=NN+NN
      DO 20 I=1,L,2
      ANGLE=0.5d0*DFLOAT(I-1)*DEL
      TRIGS(LA+I)=DCOS(ANGLE)
      TRIGS(LA+I+1)=DSIN(ANGLE)
   20 CONTINUE
      IF (IMODE.LE.3) RETURN
      DEL=0.5d0*DEL
      LA=LA+NN
      IF (MODE.EQ.5) GO TO 40
      DO 30 I=2,NN
      ANGLE=DFLOAT(I-1)*DEL
      TRIGS(LA+I)=2.d0*DSIN(ANGLE)
   30 CONTINUE
      RETURN
   40 CONTINUE
      DEL=0.5d0*DEL
      DO 50 I=2,N
      ANGLE=DFLOAT(I-1)*DEL
      TRIGS(LA+I)=DSIN(ANGLE)
   50 CONTINUE
      RETURN
      END
      SUBROUTINE CFFT99(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
      implicit double precision(a-h,o-z)
C
C PURPOSE      PERFORMS MULTIPLE FAST FOURIER TRANSFORMS.  THIS PACKAGE
C              WILL PERFORM A NUMBER OF SIMULTANEOUS COMPLEX PERIODIC
C              FOURIER TRANSFORMS OR CORRESPONDING INVERSE TRANSFORMS.
C              THAT IS, GIVEN A SET OF COMPLEX GRIDPOINT VECTORS, THE
C              PACKAGE RETURNS A SET OF COMPLEX FOURIER
C              COEFFICIENT VECTORS, OR VICE VERSA.  THE LENGTH OF THE
C              TRANSFORMS MUST BE A NUMBER GREATER THAN 1 THAT HAS
C              NO PRIME FACTORS OTHER THAN 2, 3, AND 5.
C
C              THE PACKAGE CFFT99 CONTAINS SEVERAL USER-LEVEL ROUTINES:
C
C            SUBROUTINE CFTFAX
C                AN INITIALIZATION ROUTINE THAT MUST BE CALLED ONCE
C                BEFORE A SEQUENCE OF CALLS TO CFFT99
C                (PROVIDED THAT N IS NOT CHANGED).
C
C            SUBROUTINE CFFT99
C                THE ACTUAL TRANSFORM ROUTINE ROUTINE, CABABLE OF
C                PERFORMING BOTH THE TRANSFORM AND ITS INVERSE.
C                HOWEVER, AS THE TRANSFORMS ARE NOT NORMALIZED,
C                THE APPLICATION OF A TRANSFORM FOLLOWED BY ITS
C                INVERSE WILL YIELD THE ORIGINAL VALUES MULTIPLIED
C                BY N.
C
C
C ACCESS       *FORTRAN,P=XLIB,SN=CFFT99
C
C
C USAGE        LET N BE OF THE FORM 2**P * 3**Q * 5**R, WHERE P .GE. 0,
C              Q .GE. 0, AND R .GE. 0.  THEN A TYPICAL SEQUENCE OF
C              CALLS TO TRANSFORM A GIVEN SET OF COMPLEX VECTORS OF
C              LENGTH N TO A SET OF (UNSCALED) COMPLEX FOURIER
C              COEFFICIENT VECTORS OF LENGTH N IS
C
C                   DIMENSION IFAX(13),TRIGS(2*N)
C                   COMPLEX A(...), WORK(...)
C
C                   CALL CFTFAX (N, IFAX, TRIGS)
C                   CALL CFFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
C
C              THE OUTPUT VECTORS OVERWRITE THE INPUT VECTORS, AND
C              THESE ARE STORED IN A.  WITH APPROPRIATE CHOICES FOR
C              THE OTHER ARGUMENTS, THESE VECTORS MAY BE CONSIDERED
C              EITHER THE ROWS OR THE COLUMNS OF THE ARRAY A.
C              SEE THE INDIVIDUAL WRITE-UPS FOR CFTFAX AND
C              CFFT99 BELOW, FOR A DETAILED DESCRIPTION OF THE
C              ARGUMENTS.
C
C HISTORY      THE PACKAGE WAS WRITTEN BY CLIVE TEMPERTON AT ECMWF IN
C              NOVEMBER, 1978.  IT WAS MODIFIED, DOCUMENTED, AND TESTED
C              FOR NCAR BY RUSS REW IN SEPTEMBER, 1980.  IT WAS
C              FURTHER MODIFIED FOR THE FULLY COMPLEX CASE BY DAVE
C              FULKER IN NOVEMBER, 1980.
C
C-----------------------------------------------------------------------
C
C SUBROUTINE CFTFAX (N,IFAX,TRIGS)
C
C PURPOSE      A SET-UP ROUTINE FOR CFFT99.  IT NEED ONLY BE
C              CALLED ONCE BEFORE A SEQUENCE OF CALLS TO CFFT99,
C              PROVIDED THAT N IS NOT CHANGED.
C
C ARGUMENT     IFAX(13),TRIGS(2*N)
C DIMENSIONS
C
C ARGUMENTS
C
C ON INPUT     N
C               AN EVEN NUMBER GREATER THAN 1 THAT HAS NO PRIME FACTOR
C               GREATER THAN 5.  N IS THE LENGTH OF THE TRANSFORMS (SEE
C               THE DOCUMENTATION FOR CFFT99 FOR THE DEFINITION OF
C               THE TRANSFORMS).
C
C              IFAX
C               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED
C               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING
C               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN 1 MILLION.
C
C              TRIGS
C               A REAL ARRAY OF DIMENSION 2*N
C
C ON OUTPUT    IFAX
C               CONTAINS THE FACTORIZATION OF N.  IFAX(1) IS THE
C               NUMBER OF FACTORS, AND THE FACTORS THEMSELVES ARE STORED
C               IN IFAX(2),IFAX(3),...  IF N HAS ANY PRIME FACTORS
C               GREATER THAN 5, IFAX(1) IS SET TO -99.
C
C              TRIGS
C               AN ARRAY OF TRIGONOMETRIC FUNCTION VALUES SUBSEQUENTLY
C               USED BY THE CFT ROUTINES.
C
C-----------------------------------------------------------------------
C
C SUBROUTINE CFFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
C
C PURPOSE      PERFORM A NUMBER OF SIMULTANEOUS (UNNORMALIZED) COMPLEX
C              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
C              TRANSFORMS.  GIVEN A SET OF COMPLEX GRIDPOINT
C              VECTORS, THE PACKAGE RETURNS A SET OF
C              COMPLEX FOURIER COEFFICIENT VECTORS, OR VICE
C              VERSA.  THE LENGTH OF THE TRANSFORMS MUST BE A
C              NUMBER HAVING NO PRIME FACTORS OTHER THAN
C              2, 3, AND 5.  THIS ROUTINE IS
C              OPTIMIZED FOR USE ON THE CRAY-1.
C
C ARGUMENT     COMPLEX A(N*INC+(LOT-1)*JUMP), WORK(N*LOT)
C DIMENSIONS   REAL TRIGS(2*N), INTEGER IFAX(13)
C
C ARGUMENTS
C
C ON INPUT     A
C               A COMPLEX ARRAY OF LENGTH N*INC+(LOT-1)*JUMP CONTAINING
C               THE INPUT GRIDPOINT OR COEFFICIENT VECTORS.  THIS ARRAY
C               OVERWRITTEN BY THE RESULTS.
C
C               N.B. ALTHOUGH THE ARRAY A IS USUALLY CONSIDERED TO BE OF
C               TYPE COMPLEX IN THE CALLING PROGRAM, IT IS TREATED AS
C               REAL WITHIN THE TRANSFORM PACKAGE.  THIS REQUIRES THAT
C               SUCH TYPE CONFLICTS ARE PERMITTED IN THE USER"S
C               ENVIRONMENT, AND THAT THE STORAGE OF COMPLEX NUMBERS
C               MATCHES THE ASSUMPTIONS OF THIS ROUTINE.  THIS ROUTINE
C               ASSUMES THAT THE REAL AND IMAGINARY PORTIONS OF A
C               COMPLEX NUMBER OCCUPY ADJACENT ELEMENTS OF MEMORY.  IF
C               THESE CONDITIONS ARE NOT MET, THE USER MUST TREAT THE
C               ARRAY A AS REAL (AND OF TWICE THE ABOVE LENGTH), AND
C               WRITE THE CALLING PROGRAM TO TREAT THE REAL AND
C               IMAGINARY PORTIONS EXPLICITLY.
C
C              WORK
C               A COMPLEX WORK ARRAY OF LENGTH N*LOT OR A REAL ARRAY
C               OF LENGTH 2*N*LOT.  SEE N.B. ABOVE.
C
C              TRIGS
C               AN ARRAY SET UP BY CFTFAX, WHICH MUST BE CALLED FIRST.
C
C              IFAX
C               AN ARRAY SET UP BY CFTFAX, WHICH MUST BE CALLED FIRST.
C
C
C               N.B. IN THE FOLLOWING ARGUMENTS, INCREMENTS ARE MEASURED
C               IN WORD PAIRS, BECAUSE EACH COMPLEX ELEMENT IS ASSUMED
C               TO OCCUPY AN ADJACENT PAIR OF WORDS IN MEMORY.
C
C              INC
C               THE INCREMENT (IN WORD PAIRS) BETWEEN SUCCESSIVE ELEMENT
C               OF EACH (COMPLEX) GRIDPOINT OR COEFFICIENT VECTOR
C               (E.G.  INC=1 FOR CONSECUTIVELY STORED DATA).
C
C              JUMP
C               THE INCREMENT (IN WORD PAIRS) BETWEEN THE FIRST ELEMENTS
C               OF SUCCESSIVE DATA OR COEFFICIENT VECTORS.  ON THE CRAY-
C               TRY TO ARRANGE DATA SO THAT JUMP IS NOT A MULTIPLE OF 8
C               (TO AVOID MEMORY BANK CONFLICTS).  FOR CLARIFICATION OF
C               INC AND JUMP, SEE THE EXAMPLES BELOW.
C
C              N
C               THE LENGTH OF EACH TRANSFORM (SEE DEFINITION OF
C               TRANSFORMS, BELOW).
C
C              LOT
C               THE NUMBER OF TRANSFORMS TO BE DONE SIMULTANEOUSLY.
C
C              ISIGN
C               = -1 FOR A TRANSFORM FROM GRIDPOINT VALUES TO FOURIER
C                    COEFFICIENTS.
C               = +1 FOR A TRANSFORM FROM FOURIER COEFFICIENTS TO
C                    GRIDPOINT VALUES.
C
C ON OUTPUT    A
C               IF ISIGN = -1, AND LOT GRIDPOINT VECTORS ARE SUPPLIED,
C               EACH CONTAINING THE COMPLEX SEQUENCE:
C
C               G(0),G(1), ... ,G(N-1)  (N COMPLEX VALUES)
C
C               THEN THE RESULT CONSISTS OF LOT COMPLEX VECTORS EACH
C               CONTAINING THE CORRESPONDING N COEFFICIENT VALUES:
C
C               C(0),C(1), ... ,C(N-1)  (N COMPLEX VALUES)
C
C               DEFINED BY:
C                 C(K) = SUM(J=0,...,N-1)( G(J)*EXP(-2*I*J*K*PI/N) )
C                 WHERE I = SQRT(-1)
C
C
C               IF ISIGN = +1, AND LOT COEFFICIENT VECTORS ARE SUPPLIED,
C               EACH CONTAINING THE COMPLEX SEQUENCE:
C
C               C(0),C(1), ... ,C(N-1)  (N COMPLEX VALUES)
C
C               THEN THE RESULT CONSISTS OF LOT COMPLEX VECTORS EACH
C               CONTAINING THE CORRESPONDING N GRIDPOINT VALUES:
C
C               G(0),G(1), ... ,G(N-1)  (N COMPLEX VALUES)
C
C               DEFINED BY:
C                 G(J) = SUM(K=0,...,N-1)( G(K)*EXP(+2*I*J*K*PI/N) )
C                 WHERE I = SQRT(-1)
C
C
C               A CALL WITH ISIGN=-1 FOLLOWED BY A CALL WITH ISIGN=+1
C               (OR VICE VERSA) RETURNS THE ORIGINAL DATA, MULTIPLIED
C               BY THE FACTOR N.
C
C
C EXAMPLE       GIVEN A 64 BY 9 GRID OF COMPLEX VALUES, STORED IN
C               A 66 BY 9 COMPLEX ARRAY, A, COMPUTE THE TWO DIMENSIONAL
C               FOURIER TRANSFORM OF THE GRID.  FROM TRANSFORM THEORY,
C               IT IS KNOWN THAT A TWO DIMENSIONAL TRANSFORM CAN BE
C               OBTAINED BY FIRST TRANSFORMING THE GRID ALONG ONE
C               DIRECTION, THEN TRANSFORMING THESE RESULTS ALONG THE
C               ORTHOGONAL DIRECTION.
C
C               COMPLEX A(66,9), WORK(64,9)
C               REAL TRIGS1(128), TRIGS2(18)
C               INTEGER IFAX1(13), IFAX2(13)
C
C               SET UP THE IFAX AND TRIGS ARRAYS FOR EACH DIRECTION:
C
C               CALL CFTFAX(64, IFAX1, TRIGS1)
C               CALL CFTFAX( 9, IFAX2, TRIGS2)
C
C               IN THIS CASE, THE COMPLEX VALUES OF THE GRID ARE
C               STORED IN MEMORY AS FOLLOWS (USING U AND V TO
C               DENOTE THE REAL AND IMAGINARY COMPONENTS, AND
C               ASSUMING CONVENTIONAL FORTRAN STORAGE):
C
C   U(1,1), V(1,1), U(2,1), V(2,1),  ...  U(64,1), V(64,1), 4 NULLS,
C
C   U(1,2), V(1,2), U(2,2), V(2,2),  ...  U(64,2), V(64,2), 4 NULLS,
C
C   .       .       .       .         .   .        .        .
C   .       .       .       .         .   .        .        .
C   .       .       .       .         .   .        .        .
C
C   U(1,9), V(1,9), U(2,9), V(2,9),  ...  U(64,9), V(64,9), 4 NULLS.
C
C               WE CHOOSE (ARBITRARILY) TO TRANSORM FIRST ALONG THE
C               DIRECTION OF THE FIRST SUBSCRIPT.  THUS WE DEFINE
C               THE LENGTH OF THE TRANSFORMS, N, TO BE 64, THE
C               NUMBER OF TRANSFORMS, LOT, TO BE 9, THE INCREMENT
C               BETWEEN ELEMENTS OF EACH TRANSFORM, INC, TO BE 1,
C               AND THE INCREMENT BETWEEN THE STARTING POINTS
C               FOR EACH TRANSFORM, JUMP, TO BE 66 (THE FIRST
C               DIMENSION OF A).
C
C               CALL CFFT99( A, WORK, TRIGS1, IFAX1, 1, 66, 64, 9, -1)
C
C               TO TRANSFORM ALONG THE DIRECTION OF THE SECOND SUBSCRIPT
C               THE ROLES OF THE INCREMENTS ARE REVERSED.  THUS WE DEFIN
C               THE LENGTH OF THE TRANSFORMS, N, TO BE 9, THE
C               NUMBER OF TRANSFORMS, LOT, TO BE 64, THE INCREMENT
C               BETWEEN ELEMENTS OF EACH TRANSFORM, INC, TO BE 66,
C               AND THE INCREMENT BETWEEN THE STARTING POINTS
C               FOR EACH TRANSFORM, JUMP, TO BE 1
C
C               CALL CFFT99( A, WORK, TRIGS2, IFAX2, 66, 1, 9, 64, -1)
C
C               THESE TWO SEQUENTIAL STEPS RESULTS IN THE TWO-DIMENSIONA
C               FOURIER COEFFICIENT ARRAY OVERWRITING THE INPUT
C               GRIDPOINT ARRAY, A.  THE SAME TWO STEPS APPLIED AGAIN
C               WITH ISIGN = +1 WOULD RESULT IN THE RECONSTRUCTION OF
C               THE GRIDPOINT ARRAY (MULTIPLIED BY A FACTOR OF 64*9).
C
C
C-----------------------------------------------------------------------
      DIMENSION A(*),WORK(*),TRIGS(*),IFAX(*)
C
C     SUBROUTINE "CFFT99" - MULTIPLE FAST COMPLEX FOURIER TRANSFORM
C
C     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
C     WORK IS AN AREA OF SIZE N*LOT
C     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
C     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N
C     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
C         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
C     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
C     N IS THE LENGTH OF THE DATA VECTORS
C     LOT IS THE NUMBER OF DATA VECTORS
C     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
C           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
C
C
C     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
C     PARALLEL.
C
C
C
C
C
      NN = N+N
      INK=INC+INC
      JUM = JUMP+JUMP
      NFAX=IFAX(1)
      JNK = 2
      JST = 2
      IF (ISIGN.GE.0) GO TO 30
C
C     THE INNERMOST TEMPERTON ROUTINES HAVE NO FACILITY FOR THE
C     FORWARD (ISIGN = -1) TRANSFORM.  THEREFORE, THE INPUT MUST BE
C     REARRANGED AS FOLLOWS:
C
C     THE ORDER OF EACH INPUT VECTOR,
C
C     G(0), G(1), G(2), ... , G(N-2), G(N-1)
C
C     IS REVERSED (EXCLUDING G(0)) TO YIELD
C
C     G(0), G(N-1), G(N-2), ... , G(2), G(1).
C
C     WITHIN THE TRANSFORM, THE CORRESPONDING EXPONENTIAL MULTIPLIER
C     IS THEN PRECISELY THE CONJUGATE OF THAT FOR THE NORMAL
C     ORDERING.  THUS THE FORWARD (ISIGN = -1) TRANSFORM IS
C     ACCOMPLISHED
C
C     FOR NFAX ODD, THE INPUT MUST BE TRANSFERRED TO THE WORK ARRAY,
C     AND THE REARRANGEMENT CAN BE DONE DURING THE MOVE.
C
      JNK = -2
      JST = NN-2
      IF (MOD(NFAX,2).EQ.1) GOTO 40
C
C     FOR NFAX EVEN, THE REARRANGEMENT MUST BE APPLIED DIRECTLY TO
C     THE INPUT ARRAY.  THIS CAN BE DONE BY SWAPPING ELEMENTS.
C
      IBASE = 1
      ILAST = (N-1)*INK
      NH = N/2
      DO 20 L=1,LOT
      I1 = IBASE+INK
      I2 = IBASE+ILAST
CDIR$ IVDEP
      DO 10 M=1,NH
C     SWAP REAL AND IMAGINARY PORTIONS
      HREAL = A(I1)
      HIMAG = A(I1+1)
      A(I1) = A(I2)
      A(I1+1) = A(I2+1)
      A(I2) = HREAL
      A(I2+1) = HIMAG
      I1 = I1+INK
      I2 = I2-INK
   10 CONTINUE
      IBASE = IBASE+JUM
   20 CONTINUE
      GOTO 100
C
   30 CONTINUE
      IF (MOD(NFAX,2).EQ.0) GOTO 100
C
   40 CONTINUE
C
C     DURING THE TRANSFORM PROCESS, NFAX STEPS ARE TAKEN, AND THE
C     RESULTS ARE STORED ALTERNATELY IN WORK AND IN A.  IF NFAX IS
C     ODD, THE INPUT DATA ARE FIRST MOVED TO WORK SO THAT THE FINAL
C     RESULT (AFTER NFAX STEPS) IS STORED IN ARRAY A.
C
      IBASE=1
      JBASE=1
      DO 60 L=1,LOT
C     MOVE REAL AND IMAGINARY PORTIONS OF ELEMENT ZERO
      WORK(JBASE) = A(IBASE)
      WORK(JBASE+1) = A(IBASE+1)
      I=IBASE+INK
      J=JBASE+JST
CDIR$ IVDEP
      DO 50 M=2,N
C     MOVE REAL AND IMAGINARY PORTIONS OF OTHER ELEMENTS (POSSIBLY IN
C     REVERSE ORDER, DEPENDING ON JST AND JNK)
      WORK(J) = A(I)
      WORK(J+1) = A(I+1)
      I=I+INK
      J=J+JNK
   50 CONTINUE
      IBASE=IBASE+JUM
      JBASE=JBASE+NN
   60 CONTINUE
C
  100 CONTINUE
C
C     PERFORM THE TRANSFORM PASSES, ONE PASS FOR EACH FACTOR.  DURING
C     EACH PASS THE DATA ARE MOVED FROM A TO WORK OR FROM WORK TO A.
C
C     FOR NFAX EVEN, THE FIRST PASS MOVES FROM A TO WORK
      IGO = 110
C     FOR NFAX ODD, THE FIRST PASS MOVES FROM WORK TO A
      IF (MOD(NFAX,2).EQ.1) IGO = 120
      LA=1
      DO 140 K=1,NFAX
      IF (IGO.EQ.120) GO TO 120
  110 CONTINUE
      CALL VPASSM(A(1),A(2),WORK(1),WORK(2),TRIGS,
     *   INK,2,JUM,NN,LOT,N,IFAX(K+1),LA)
      IGO=120
      GO TO 130
  120 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(1),A(2),TRIGS,
     *    2,INK,NN,JUM,LOT,N,IFAX(K+1),LA)
      IGO=110
  130 CONTINUE
      LA=LA*IFAX(K+1)
  140 CONTINUE
C
C     AT THIS POINT THE FINAL TRANSFORM RESULT IS STORED IN A.
C
      RETURN
      END
      SUBROUTINE CFTFAX(N,IFAX,TRIGS)
      implicit double precision(a-h,o-z)
      DIMENSION IFAX(13),TRIGS(1)
C
C     THIS ROUTINE WAS MODIFIED FROM TEMPERTON"S ORIGINAL
C     BY DAVE FULKER.  IT NO LONGER PRODUCES FACTORS IN ASCENDING
C     ORDER, AND THERE ARE NONE OF THE ORIGINAL 'MODE' OPTIONS.
C
C ON INPUT     N
C               THE LENGTH OF EACH COMPLEX TRANSFORM TO BE PERFORMED
C
C               N MUST BE GREATER THAN 1 AND CONTAIN NO PRIME
C               FACTORS GREATER THAN 5.
C
C ON OUTPUT    IFAX
C               IFAX(1)
C                 THE NUMBER OF FACTORS CHOSEN OR -99 IN CASE OF ERROR
C               IFAX(2) THRU IFAX( IFAX(1)+1 )
C                 THE FACTORS OF N IN THE FOLLOWIN ORDER:  APPEARING
C                 FIRST ARE AS MANY FACTORS OF 4 AS CAN BE OBTAINED.
C                 SUBSEQUENT FACTORS ARE PRIMES, AND APPEAR IN
C                 ASCENDING ORDER, EXCEPT FOR MULTIPLE FACTORS.
C
C              TRIGS
C               2N SIN AND COS VALUES FOR USE BY THE TRANSFORM ROUTINE
C
      CALL FACT(N,IFAX)
      K = IFAX(1)
      IF (K .LT. 1 .OR. IFAX(K+1) .GT. 5) IFAX(1) = -99
      IF (IFAX(1) .LE. 0 )CALL ULIBER
      CALL CFTRIG (N, TRIGS)
      RETURN
      END
      SUBROUTINE FACT(N,IFAX)
      implicit double precision(a-h,o-z)
C     FACTORIZATION ROUTINE THAT FIRST EXTRACTS ALL FACTORS OF 4
      DIMENSION IFAX(13)
      IF (N.GT.1) GO TO 10
      IFAX(1) = 0
      IF (N.LT.1) IFAX(1) = -99
      RETURN
   10 NN=N
      K=1
C     TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
      K=K+1
      IFAX(K)=4
      NN=NN/4
      IF (NN.EQ.1) GO TO 80
      GO TO 20
C     TEST FOR EXTRA FACTOR OF 2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
C     TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40
C     NOW FIND REMAINING FACTORS
   50 L=5
      MAX = SQRT(FLOAT(NN))
      INC=2
C     INC ALTERNATELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L
      IF (NN.EQ.1) GO TO 80
      GO TO 60
   70 IF (L.GT.MAX) GO TO 75
      L=L+INC
      INC=6-INC
      GO TO 60
   75 K = K+1
      IFAX(K) = NN
   80 IFAX(1)=K-1
C     IFAX(1) NOW CONTAINS NUMBER OF FACTORS
      RETURN
      END
      SUBROUTINE CFTRIG(N,TRIGS)
      implicit double precision(a-h,o-z)
      DIMENSION TRIGS(1)
      PI=2.d0*ASIN(1.d0)
      DEL=(PI+PI)/DFLOAT(N)
      L=N+N
      DO 10 I=1,L,2
      ANGLE=0.5d0*DFLOAT(I-1)*DEL
      TRIGS(I)=DCOS(ANGLE)
      TRIGS(I+1)=DSIN(ANGLE)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE VPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA)
      implicit double precision(a-h,o-z)
      DIMENSION A(N),B(N),C(N),D(N),TRIGS(N)
*      do i=1,n
*       write(*,*) a(i),b(i),c(i),d(i)
*      end do
C
C     SUBROUTINE "VPASSM" - MULTIPLE VERSION OF "VPASSA"
C     PERFORMS ONE PASS THROUGH DATA
C     AS PART OF MULTIPLE COMPLEX (INVERSE) FFT ROUTINE
C     A IS FIRST REAL INPUT VECTOR
C     B IS FIRST IMAGINARY INPUT VECTOR
C     C IS FIRST REAL OUTPUT VECTOR
C     D IS FIRST IMAGINARY OUTPUT VECTOR
C     TRIGS IS PRECALCULATED TABLE OF SINES & COSINES
C     INC1 IS ADDRESSING INCREMENT FOR A AND B
C     INC2 IS ADDRESSING INCREMENT FOR C AND D
C     INC3 IS ADDRESSING INCREMENT BETWEEN A"S & B"S
C     INC4 IS ADDRESSING INCREMENT BETWEEN C"S & D"S
C     LOT IS THE NUMBER OF VECTORS
C     N IS LENGTH OF VECTORS
C     IFAC IS CURRENT FACTOR OF N
C     LA IS PRODUCT OF PREVIOUS FACTORS
C
      DATA SIN36/0.587785252292473/,COS36/0.809016994374947/,
     *     SIN72/0.951056516295154/,COS72/0.309016994374947/,
     *     SIN60/0.866025403784437/
C
      M=N/IFAC
      IINK=M*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.GT.4) RETURN
      GO TO (10,50,90,130),IGO
C
C     CODING FOR FACTOR 2
C
   10 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      DO 20 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 15 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=A(IA+I)-A(IB+I)
      D(JB+J)=B(IA+I)-B(IB+I)
      I=I+INC3
      J=J+INC4
   15 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   20 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 40 K=LA1,M,LA
      KB=K+K-2
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      DO 30 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 25 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)-B(IB+I))
      D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)-B(IB+I))
      I=I+INC3
      J=J+INC4
   25 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   30 CONTINUE
      JBASE=JBASE+JUMP
   40 CONTINUE
      RETURN
C
C     CODING FOR FACTOR 3
C
C-----------------------------------------------------------------------

   50 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      DO 60 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 55 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=(A(IA+I)-0.5d0*(A(IB+I)+A(IC+I)))-
     & (SIN60*(B(IB+I)-B(IC+I)))
      C(JC+J)=(A(IA+I)-0.5d0*(A(IB+I)+A(IC+I)))+
     & (SIN60*(B(IB+I)-B(IC+I)))
      D(JB+J)=(B(IA+I)-0.5d0*(B(IB+I)+B(IC+I)))+
     & (SIN60*(A(IB+I)-A(IC+I)))
      D(JC+J)=(B(IA+I)-0.5d0*(B(IB+I)+B(IC+I)))-
     & (SIN60*(A(IB+I)-A(IC+I)))
      I=I+INC3
      J=J+INC4
   55 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   60 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 80 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      DO 70 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 65 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=
     *    C1*((A(IA+I)-0.5d0*(A(IB+I)+A(IC+I)))-
     *   (SIN60*(B(IB+I)-B(IC+I))))
     *   -S1*((B(IA+I)-0.5d0*(B(IB+I)+B(IC+I)))+
     *   (SIN60*(A(IB+I)-A(IC+I))))
      D(JB+J)=
     *    S1*((A(IA+I)-0.5d0*(A(IB+I)+A(IC+I)))-
     *    (SIN60*(B(IB+I)-B(IC+I))))
     *    +C1*((B(IA+I)-0.5d0*(B(IB+I)+B(IC+I)))+
     *    (SIN60*(A(IB+I)-A(IC+I))))
      C(JC+J)=
     *    C2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+
     *    (SIN60*(B(IB+I)-B(IC+I))))
     *    -S2*((B(IA+I)-0.5d0*(B(IB+I)+B(IC+I)))-
     *    (SIN60*(A(IB+I)-A(IC+I))))
      D(JC+J)=
     *    S2*((A(IA+I)-0.5d0*(A(IB+I)+A(IC+I)))+
     *    (SIN60*(B(IB+I)-B(IC+I))))
     *    +C2*((B(IA+I)-0.5d0*(B(IB+I)+B(IC+I)))-
     *    (SIN60*(A(IB+I)-A(IC+I))))
      I=I+INC3
      J=J+INC4
   65 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   70 CONTINUE
      JBASE=JBASE+JUMP
   80 CONTINUE
      RETURN
C
C     CODING FOR FACTOR 4
C
   90 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      DO 100 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 95 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      D(JC+J)=(B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I))
      C(JB+J)=(A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))
      C(JD+J)=(A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))
      D(JB+J)=(B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I))
      D(JD+J)=(B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I))
      I=I+INC3
      J=J+INC4
   95 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  100 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 120 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      DO 110 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 105 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      C(JC+J)=
     *    C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     *   -S2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      D(JC+J)=
     *    S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     *   +C2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      C(JB+J)=
     *    C1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))
     *   -S1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      D(JB+J)=
     *    S1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))
     *   +C1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      C(JD+J)=
     *    C3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))
     *   -S3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      D(JD+J)=
     *    S3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))
     *   +C3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  105 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  110 CONTINUE
      JBASE=JBASE+JUMP
  120 CONTINUE
      RETURN
C
C     CODING FOR FACTOR 5
C
  130 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      IE=ID+IINK
      JE=JD+JINK
      DO 140 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 135 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *  -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      C(JE+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *  +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      D(JB+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *  +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      D(JE+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *  -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      C(JC+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *  -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      C(JD+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *  +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      D(JC+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *  +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      D(JD+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *  -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  135 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  140 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 160 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      DO 150 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 145 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=
     *    C1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   -S1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JB+J)=
     *    S1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   +C1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JE+J)=
     *    C4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   -S4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JE+J)=
     *    S4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   +C4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JC+J)=
     *    C2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   -S2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JC+J)=
     *    S2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   +C2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      C(JD+J)=
     *    C3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   -S3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JD+J)=
     *    S3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   +C3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      I=I+INC3
      J=J+INC4
  145 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  150 CONTINUE
      JBASE=JBASE+JUMP
  160 CONTINUE
      RETURN
      END
******************************
* you need theses subroutines
* not coming from xmp
******************************
* missing subroutines
      subroutine uliber
      implicit double precision(a-h,o-z)
*     print * , ' uliber'
      end





