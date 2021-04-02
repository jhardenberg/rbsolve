#include "config.h"
#ifdef PRINF
c ______________________________________________________________________
c Pressure Helmholtz eq. for pressure
c ______________________________________________________________________

      subroutine prinf_helmp (t,D,E,F,yv,yp,qp,sp,
     +     kxm,kzm,p,xlh,xuh,reda,redb,redc,redr,Ra,pbot,ptop)
 
#include "param.h"
      real*8 qp(Ny),sp(Ny),yv(0:Ny+1),yp(0:Ny+1)
      real*8 D(nwax+1,2*nwaz+1,Nyl),F(nwax+1,2*nwaz+1,Nyl)
      real*8 E(nwax+1,2*nwaz+1,Nyl),pbot,ptop
      complex*16 t(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 v(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 p(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 TN(nwax+1,2*nwaz+1,Nylmem)

      real*8 kxm(nwax+1),kzm(2*nwaz+1),Ra
      real*8 rp(Ny),cff

#define NX ((nwax+1)*(2*nwaz+1))

      real*8 xlh(NX,Nyl)
      real*8 xuh(NX,Nyl)

#define NNR (2*NPROC-2)

      real*8 reda(NX,0:NNR+1)
      real*8 redb(NX,0:NNR+1)
      real*8 redc(NX,0:NNR+1)
      complex*16 redr(NX,0:NNR+1)

      p(:,:,:)=-9e25
      call TN_prinf_helmp (t,yv,yp,qp,sp,kxm,kzm,TN,Ra,pbot,ptop)
      call tridmpi(D,E,F,TN,reda,redb,redc,redr,
     &     xlh,xuh,p(1,1,1),
     &     Nylmax1,NX,NPROC,iam)
      p(1,1,:)=p(1,1,:)*(1.d0-pbot)*(1.d0-ptop)

      return
      end
#endif
