
#include "config.h"

c---------------------------------------------------------------------
      subroutine diffus (u_n,invRe,qv,rv,sv,K2,Lu_n)
      
#include "param.h"
      
      complex*16 u_n(nwax+1,2*nwaz+1,0:Nyl+1)
      real*8 invRe,qv(Ny),rv(Ny),sv(Ny),rv1,sv1,qv1
      real*8 K2(nwax+1,2*nwaz+1)
      complex*16 Lu_n(nwax+1,2*nwaz+1,0:Nyl+1)
      integer jj

c Calcolo il termine diffusivo

      do j=Nylmin,Nylmax
       jj=j+Nybase
       qv1=invRe*2.d0*qv(jj)
       rv1=invRe*2.d0*rv(jj)
       sv1=invRe*2.d0*sv(jj)
       do k=1,2*nwaz+1   
        do i=1,nwax+1
       Lu_n(i,k,j)=qv1*u_n(i,k,j-1)-(rv1
     +   +invRe*K2(i,k))*u_n(i,k,j)+sv1*u_n(i,k,j+1) 

*       Lu_n(i,k,j)=invRe*(2.d0*qv(jj-1)*u_n(i,k,j-1)-(2.d0*rv(jj-1)
*     +             +K2(i,k))*u_n(i,k,j)+2.d0*sv(jj-1)*u_n(i,k,j+1))  
        enddo
       enddo
      enddo

      return
      end
      


      subroutine addrotv (u_n,v_n,ff,Lu_n,Lv_n)

#include "param.h"

      complex*16 u_n(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 v_n(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 Lu_n(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 Lv_n(nwax+1,2*nwaz+1,0:Nyl+1)
      real*8 ff
      integer jj

!      print *,Lu_n(128,128,Nyl/2),ff*v_n(128,128,Nyl/2),ff
!	stop
      do j=Nylmin,Nylmax
       do k=1,2*nwaz+1
        do i=1,nwax+1
           Lu_n(i,k,j)=Lu_n(i,k,j)-ff*v_n(i,k,j)
           Lv_n(i,k,j)=Lv_n(i,k,j)+ff*u_n(i,k,j)
        enddo
       enddo
      enddo 
    
      return
      end 


