#include "config.h"
c ______________________________________________________________________

      subroutine velocita (ustar,vstar,wstar,p,mm,nn,yp,yv
     +                     ,ut,vt,wt,dt,alfa,beta,l)

#include "param.h"
      integer*4 l,jj
      complex*16 ustar(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 vstar(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 wstar(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 ut(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 vt(nwax+1,2*nwaz+1,0:Nyl+1),wt(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 p(nwax+1,2*nwaz+1,0:Nyl+1),mm(nwax+1),nn(2*nwaz+1)
      real*8 yp(0:Ny+1),yv(0:Ny+1),dt,alfa(3),beta(3)
      real*8 cffa,cff1,cff2,cff0
      complex*16 pint
      
      cffa=dt*(alfa(l)+beta(l))

c 'pint' is p interpolated on yv(j)

      do j=Nylmin,Nylmax

      jj=j+Nybase
      
      cff1=1.0/(yp(jj)-yv(jj))*(yv(jj)-yp(jj-1))/(yp(jj)-yp(jj-1))
      cff2=1.0/(yv(jj)-yp(jj-1))*(yp(jj)-yv(jj))/(yp(jj)-yp(jj-1))
      cff0=1.0/(yp(jj)-yp(jj-1))*(yv(jj)-yp(jj-1))

      do k=1,2*nwaz+1
       do i=1,nwax+1
        pint=p(i,k,j-1)+(p(i,k,j)-p(i,k,j-1))*cff0
        ut(i,k,j)=ustar(i,k,j)-mm(i)*pint*cffa
        wt(i,k,j)=wstar(i,k,j)-nn(k)*pint*cffa
        vt(i,k,j)=vstar(i,k,j)-(((p(i,k,j)-pint)*cff1+ 
     +    (pint-p(i,k,j-1))*cff2))*cffa
        enddo
       enddo 

*       wt(1,1,j)=(0.d0,0.d0)
       wt(1,1,j)=wstar(1,1,j)

c	Does this solve our pressure problem ?
       vt(1,1,j)=(0.d0,0.d0)
cxxxx Should we impose vertical vel ave to 0 ?
       ut(1,1,j)=ustar(1,1,j)
      enddo 

      return
      end  


