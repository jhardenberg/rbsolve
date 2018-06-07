#include "config.h"

c ______________________________________________________________________

      subroutine TN_prinf_helmu (p,t,yp,kxm,kzm,TNU,TNV,TNW,l,Ra)
							      
#include "param.h"
      integer*4 l,jj
      complex*16 TNU(nwax+1,2*nwaz+1,Nyl),aa,bb
      complex*16 TNW(nwax+1,2*nwaz+1,Nyl)
      complex*16 TNV(nwax+1,2*nwaz+1,Nyl)
      complex*16 p(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 t(nwax+1,2*nwaz+1,0:Nylmem+1)
      real*8 kxm(nwax+1),kzm(2*nwaz+1),yp(0:Ny+1)
      real*8 cff1,Ra
      
      do j=1,Nylmax
       jj=j+Nybase
       cff1=1.0/(yp(jj)-yp(jj-1))
       do k=1,2*nwaz+1
        bb=(0.d0,1.d0)*kzm(k)*0.5
        do i=1,nwax+1
          aa=(0.d0,1.d0)*kxm(i)*0.5
	  TNV(i,k,j)= -Ra*t(i,k,j)+
     &          (p(i,k,j)-p(i,k,j-1))*cff1
	  TNU(i,k,j)=aa*(p(i,k,j)+p(i,k,j-1))
	  TNW(i,k,j)=bb*(p(i,k,j)+p(i,k,j-1))
c CHECK +1 and -1!
        enddo
       enddo
      enddo

      return
      end


