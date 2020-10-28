
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
       Lu_n(:,:,j)=qv1*u_n(:,:,j-1)-(rv1
     +   +invRe*K2(:,:))*u_n(:,:,j)+sv1*u_n(:,:,j+1) 
      enddo

      return
      end

      subroutine addrotw (u_n,v_n,ff,Lu_n,Lv_n)

#include "param.h"

      complex*16 u_n(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 v_n(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 Lu_n(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 Lv_n(nwax+1,2*nwaz+1,0:Nyl+1)
      real*8 ff
      integer jj

      do j=Nylmin,Nylmax
           Lu_n(:,:,j)=Lu_n(:,:,j)-ff*v_n(:,:,j)
           Lv_n(:,:,j)=Lv_n(:,:,j)+ff*u_n(:,:,j)
      enddo 
    
      return
      end 

      subroutine addrot (u_n,v_n,w_n,om_y,om_z,Lu_n,Lv_n,Lw_n)
 
#include "param.h"
 
       complex*16 u_n(nwax+1,2*nwaz+1,0:Nyl+1)
       complex*16 v_n(nwax+1,2*nwaz+1,0:Nyl+1)
       complex*16 w_n(nwax+1,2*nwaz+1,0:Nyl+1)
       complex*16 Lu_n(nwax+1,2*nwaz+1,0:Nyl+1)
       complex*16 Lv_n(nwax+1,2*nwaz+1,0:Nyl+1)
       complex*16 Lw_n(nwax+1,2*nwaz+1,0:Nyl+1)
       real*8 om_y,om_z
       integer jj
 
       do j=Nylmin,Nylmax
            Lu_n(:,:,j)=Lu_n(:,:,j)+om_y*w_n(:,:,j)-om_z*v_n(:,:,j)
            Lv_n(:,:,j)=Lv_n(:,:,j)+om_z*u_n(:,:,j)
            Lw_n(:,:,j)=Lw_n(:,:,j)-om_y*u_n(:,:,j)
       enddo 
     
       return
       end 


#ifdef LAPSERATE
      subroutine add_lapse (Lt_n,v_1,lapse)

#include "param.h"

      complex*16 Lt_n(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 v_1(nwax+1,2*nwaz+1,0:Nyl+1)
      real*8 lapse

      do j=Nylmin,Nylmax
         Lt_n(:,:,j)= Lt_n(:,:,j)-lapse*v_1(:,:,j)
      enddo

      return
      end
#endif


#ifdef RADIATION
      subroutine add_radiation (Lt_n,qrad)

#include "param.h"

      complex*16 Lt_n(nwax+1,2*nwaz+1,0:Nyl+1)
      real*8 qrad

      do j=Nylmin,Nylmax
       Lt_n(1,1,j)= Lt_n(1,1,j)-qrad
      enddo

      return
      end
#endif



