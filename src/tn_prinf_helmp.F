#include "config.h"

c ______________________________________________________________________

      subroutine TN_prinf_helmp (t,yv,kxm,kzm,TN,l,Ra)
							      
#include "param.h"
      integer*4 l,jj
      complex*16 t(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 TN(nwax+1,2*nwaz+1,Nyl)
      real*8 kxm(nwax+1),kzm(2*nwaz+1),yv(0:Ny+1)
      real*8 cff1,Ra
      
      do j=1,Nylmax
       jj=j+Nybase
       cff1=1.0/(yv(jj+1)-yv(jj))*Ra
       do k=1,2*nwaz+1
        do i=1,nwax+1
	  TN(i,k,j)=(t(i,k,j+1)-t(i,k,j))*cff1
        enddo
       enddo
      enddo

	
#ifdef PRESSUREFIX
c     The famous pressure fix! GP& JH 07/2002
c     We impose p(1,1,Np)=0 to avoid singularity of the Helmholtz prooblem

c    ???? 11/07 JvH 
c      if(iam.eq.NPROC) TN(1,1,Nylmax1)=0.d0
#endif

      return
      end


