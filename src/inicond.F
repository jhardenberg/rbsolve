c	Creates an initial linear conductive condition
c	perturbed with white noise on all wavenumbers

#include "config.h"

      program main

#include "param.h"

      complex*16 t(nwax+1,2*nwaz+1,0:Nylmem+1)


      real*8 yp(0:Ny+1),yv(0:Ny+1)
      real*8 qp(Ny-1),rp(Ny-1),sp(Ny-1),qv(Ny-2),rv(Ny-2),sv(Ny-2)

      character*12 file0
      character*45 file

      character cc*1,str*30

      real*8 ptan,qtan,ran2
      real*8 DD,Lx,Lz,Re,invRe,dt,dx,dz
      real*8 kxmin,kzmin,Ra,Kscalar,Pr

      integer*4 nsave,ttot
      real*8 alfa(3),beta(3),gamma(3),zeta(3)

      common/ALLPARA/ptan,qtan,DD,Re,Lx,Lz,invRe,dt,dx,dz,
     +   kxmin,kzmin,nsave,ttot,alfa,beta,gamma,zeta,Ra,Kscalar,Pr

      real*8 A(nwax+1,2*nwaz+1)
      real*8 kxm(nwax+1),kzm(2*nwaz+1)
      complex*16 mm(nwax+1),nn(2*nwaz+1)
      real*8 K2(nwax+1,2*nwaz+1)

      integer ik,idum,l
      real*8 b,c,K0,sigmaL,mediaL
      real*8 T_L((Ny+1)/2),B_L((Ny+1)/2),rr,ReT,ImT,ampl
      real*8 pi2
      parameter (pi2=6.28318530717959d0)


#ifndef NOMPI 
      print *,"Please compile with the NOMPI option!"
      stop
#endif      
 
      iam=1

      b=18
      c=6
      K0=40
      sigmaL=2
      mediaL=10
      ampl=0.001

      idum=-32433

      rr=ran2(idum)

      call slice()
      call read_para(0)
      call Crea_wave(kxm,kzm,K2,mm,nn,kxmin,kzmin)
      call punti (DD,ptan,qtan,yp,yv,qp,rp,sp,qv,rv,sv)

      do j=1,Ny
	do i=1,nwax+1
           do k=1,2*nwaz+1
		t(i,k,j)=(0.d0,0.d0)
           enddo
        enddo
      enddo

      call savefield(t,'u',0,0.d0)
      call savefield(t,'v',0,0.d0)
      call savefield(t,'w',0,0.d0)

c     Crea lo spettro orrizontale
      do i=1,nwax+1
       do k=1,2*nwaz+1
        A(i,k)=dsqrt((2.d0*(dsqrt(K2(i,k))/K0)**c)
     &    /(1+(dsqrt(K2(i,k))/K0)**b))
       enddo
      enddo
c     Crea lo spettro verticale (una gaussiana)





      do l=1,(Ny+1)/2

         B_L(l)=(1/(sigmaL*dsqrt(pi2)))*
     &   dexp(-0.5d0*((l-mediaL)/sigmaL)**2.d0)

         print *,l
      do j=1,Ny

        T_L(l)=dsin(l*pi2*yv(j))

        do i=1,nwax+1
          do k=1,2*nwaz+1
            rr=ran2(idum)
            ReT=A(i,k)*dcos(rr*pi2)
            ImT=A(i,k)*dsin(rr*pi2)
            t(i,k,j)=t(i,k,j)+ampl*B_L(l)*dcmplx(ReT,ImT)*T_L(l)    
         enddo
        enddo
      enddo
      enddo
      do j=1,Ny
          t(1,1,j)=(yv(j)/DD+DD/2)*
     & (TEMPERATURE_TOP-TEMPERATURE_BOTTOM)+TEMPERATURE_BOTTOM      
      enddo
      call savefield(t,'t',0,0.d0)
  
      end



