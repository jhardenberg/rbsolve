#include "config.h"

c ______________________________________________________________________

      subroutine TN_Poiss (ustar,vstar,wstar,yv,kxm,kzm,TN,dt
     +                     ,alfa,beta,l)
							      
#include "param.h"
      integer*4 l,jj
      complex*16 ustar(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 vstar(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 wstar(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 TN(nwax+1,2*nwaz+1,Nyl),aa,bb
      real*8 kxm(nwax+1),kzm(2*nwaz+1),yv(0:Ny+1),dt,alfa(3),beta(3)
      real*8 cff1,aff
      
      aff=1.0/(dt*(alfa(l)+beta(l))) 
      do j=1,Nylmax
       jj=j+Nybase
       cff1=1.0/(yv(jj+1)-yv(jj))*aff
       do k=1,2*nwaz+1
        bb=(0.d0,1.d0)*kzm(k)*0.5*aff
        do i=1,nwax+1
          aa=(0.d0,1.d0)*kxm(i)*0.5*aff
	  TN(i,k,j)=aa*(ustar(i,k,j+1)+ustar(i,k,j))+
     &           (vstar(i,k,j+1)-vstar(i,k,j))*cff1+
     &           bb*(wstar(i,k,j+1)+wstar(i,k,j))
        enddo
       enddo
      enddo

	
#ifdef PRESSUREFIX
c     The famous pressure fix! GP& JH 07/2002
c     We impose p(1,1,Np)=0 to avoid singularity of the Helmholtz prooblem

c    ???? 11/07 JvH 
      if(iam.eq.NPROC) TN(1,1,Nylmax1)=0.d0
#endif

      return
      end


