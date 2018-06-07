#include "config.h"
c ______________________________________________________________________
c Pressure Helmholtz eq. for pressure
c ______________________________________________________________________

      subroutine prinf_helmp (t,D,E,F,yv,
     +          kxm,kzm,p,l,
     +     xlh,xuh,reda,redb,redc,redr,Ra )
 
#include "param.h"
      integer*4 l
      real*8 D(nwax+1,2*nwaz+1,Nyl),F(nwax+1,2*nwaz+1,Nyl),yv(0:Ny+1)
      real*8 E(nwax+1,2*nwaz+1,Nyl)
      complex*16 t(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 p(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 TN(nwax+1,2*nwaz+1,Nylmem)

      real*8 kxm(nwax+1),kzm(2*nwaz+1),Ra

#define NX ((nwax+1)*(2*nwaz+1))

      real*8 xlh(NX,Nyl)
      real*8 xuh(NX,Nyl)

#define NNR (2*NPROC-2)

      real*8 reda(NX,0:NNR+1)
      real*8 redb(NX,0:NNR+1)
      real*8 redc(NX,0:NNR+1)
      complex*16 redr(NX,0:NNR+1)

      call TN_prinf_helmp (t,yv,kxm,kzm,TN,l,Ra)
      call tridmpi(D,E,F,TN,reda,redb,redc,redr,
     &     xlh,xuh,p(1,1,1),
     &     Nylmax1,NX,NPROC,iam)

cxxxxx  these are most probably useless!!! jh & gp 07/2002
#ifndef FREE_SLIP_BOTTOM
      if(iam.eq.1) then 
       do k=1,2*nwaz+1
        do i=1,nwax+1
         p(i,k,0)=p(i,k,1)
        enddo
       enddo
      endif
#else
      if(iam.eq.1) then
       do k=1,2*nwaz+1
        do i=1,nwax+1
         p(i,k,0)=p(i,k,1)
        enddo
       enddo
      endif

#endif
#ifndef FREE_SLIP_TOP
      if(iam.eq.NPROC) then
       do k=1,2*nwaz+1
        do i=1,nwax+1
         p(i,k,Nyl)=p(i,k,Nyl-1)
        enddo
       enddo
      endif
#else
      if(iam.eq.NPROC) then
       do k=1,2*nwaz+1
        do i=1,nwax+1
         p(i,k,Nyl)=p(i,k,Nyl-1)
        enddo
       enddo
      endif
#endif


      return
      end
