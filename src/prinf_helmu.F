#include "config.h"
c ______________________________________________________________________
c Pressure Helmholtz eq. for velocity
c ______________________________________________________________________

      subroutine prinf_helmu (p,t,D,E,F,yp,
     +          kxm,kzm,u,v,w,l,
     +     xlh,xuh,reda,redb,redc,redr,Ra )
 
#include "param.h"
      integer*4 l
      real*8 D(nwax+1,2*nwaz+1,Nyl),F(nwax+1,2*nwaz+1,Nyl),yp(0:Ny+1)
      real*8 E(nwax+1,2*nwaz+1,Nyl)
      complex*16 u(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 v(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 w(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 p(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 t(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 TNU(nwax+1,2*nwaz+1,Nylmem)
      complex*16 TNV(nwax+1,2*nwaz+1,Nylmem)
      complex*16 TNW(nwax+1,2*nwaz+1,Nylmem)

      real*8 kxm(nwax+1),kzm(2*nwaz+1),Ra

#define NX ((nwax+1)*(2*nwaz+1))

      real*8 xlh(NX,Nyl)
      real*8 xuh(NX,Nyl)

#define NNR (2*NPROC-2)

      real*8 reda(NX,0:NNR+1)
      real*8 redb(NX,0:NNR+1)
      real*8 redc(NX,0:NNR+1)
      complex*16 redr(NX,0:NNR+1)

      call TN_prinf_helmu (p,t,yp,kxm,kzm,TNU,TNV,TNW,l,Ra)
      call tridmpi(D,E,F,TNU,reda,redb,redc,redr,
     &     xlh,xuh,u(1,1,1),
     &     Nylmax1,NX,NPROC,iam)
      call tridmpi(D,E,F,TNW,reda,redb,redc,redr,
     &     xlh,xuh,w(1,1,1),
     &     Nylmax1,NX,NPROC,iam)
      call tridmpi(D,E,F,TNV,reda,redb,redc,redr,
     &     xlh,xuh,v(1,1,1),
     &     Nylmax1,NX,NPROC,iam)
         v(1,1,:)=0.0

      return
      end
