#include "config.h"
c ______________________________________________________________________

      subroutine Poiss_imp (ustar,vstar,wstar,D,E,F,yv,
     +          kxm,kzm,p,dt,alfa,beta,l,
     +     xlh,xuh,reda,redb,redc,redr )
      
#include "param.h"
      integer*4 l
      real*8 D(nwax+1,2*nwaz+1,Nyl),F(nwax+1,2*nwaz+1,Nyl),yv(0:Ny+1)
      real*8 E(nwax+1,2*nwaz+1,Nyl)
      complex*16 ustar(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 vstar(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 wstar(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 p(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 TN(nwax+1,2*nwaz+1,Nylmem)

      real*8 kxm(nwax+1),kzm(2*nwaz+1),dt,alfa(3),beta(3)

#define NX ((nwax+1)*(2*nwaz+1))

      real*8 xlh(NX,Nyl)
      real*8 xuh(NX,Nyl)

#define NNR (2*NPROC-2)

      real*8 reda(NX,0:NNR+1)
      real*8 redb(NX,0:NNR+1)
      real*8 redc(NX,0:NNR+1)
      complex*16 redr(NX,0:NNR+1)

      call TN_Poiss (ustar,vstar,wstar,yv,kxm,kzm,TN,
     +                dt,alfa,beta,l)
      call tridmpi(D,E,F,TN,reda,redb,redc,redr,
     &     xlh,xuh,p(1,1,1),
     &     Nylmax1,NX,NPROC,iam)

cxxxxx  these are most probably useless!!! jh & gp 07/2002
#ifndef FREE_SLIP_BOTTOM
      if(iam.eq.1) then 
         p(:,:,0)=p(:,:,1)
      endif
#else
      if(iam.eq.1) then
         p(:,:,0)=p(:,:,1)
      endif

#endif
#ifndef FREE_SLIP_TOP
      if(iam.eq.NPROC) then
         p(:,:,Nyl)=p(:,:,Nyl-1)
      endif
#else
      if(iam.eq.NPROC) then
         p(:,:,Nyl)=p(:,:,Nyl-1)
      endif
#endif

      return
      end
