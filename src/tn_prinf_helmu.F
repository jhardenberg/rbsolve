#include "config.h"

c ______________________________________________________________________

      subroutine TN_prinf_helmu (p,t,yp,yv,kxm,kzm,TNU,TNV,TNW,Ra)
							      
#include "param.h"
      integer*4 jj
      complex*16 TNU(nwax+1,2*nwaz+1,Nylmem)
      complex*16 TNW(nwax+1,2*nwaz+1,Nylmem)
      complex*16 TNV(nwax+1,2*nwaz+1,Nylmem)
      complex*16 p(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 t(nwax+1,2*nwaz+1,0:Nyl+1)
      real*8 kxm(nwax+1),kzm(2*nwaz+1),yp(0:Ny+1),yv(0:Ny+1)
      real*8 Ra,cff0,cff1,cff2,cff
      complex*16 mm,nn,pint

      do j=1,Nyl
       jj=j+Nybase
       cff0=1.0/(yp(jj)-yp(jj-1))*(yv(jj)-yp(jj-1))
       cff1=1.0/(yp(jj)-yv(jj))*(yv(jj)-yp(jj-1))/(yp(jj)-yp(jj-1))
       cff2=1.0/(yv(jj)-yp(jj-1))*(yp(jj)-yv(jj))/(yp(jj)-yp(jj-1))
       cff=1.d0/(yp(jj)-yp(jj-1))
       do k=1,2*nwaz+1
        nn=(0.d0,1.d0)*kzm(k)
        do i=1,nwax+1
          mm=(0.d0,1.d0)*kxm(i)
          pint=p(i,k,j-1)+(p(i,k,j)-p(i,k,j-1))*cff0
	  TNV(i,k,j)= -Ra*t(i,k,j)+
     &      (((p(i,k,j)-pint)*cff1+
     &        (pint-p(i,k,j-1))*cff2))
          TNU(i,k,j)=mm*pint
          TNW(i,k,j)=nn*pint

!!          TNV(i,k,j)= -Ra*t(i,k,j)+
!!    &       (p(i,k,j)-p(i,k,j-1))*cff
!!	  TNU(i,k,j)=mm*(p(i,k,j)+p(i,k,j-1))*0.5
!!          TNW(i,k,j)=nn*(p(i,k,j)+p(i,k,j-1))*0.5

c CHECK +1 and -1!
        enddo
       enddo
      enddo

      if(iam.eq.1) then ! No slip!
         TNU(:,:,1)=0.d0
         TNV(:,:,1)=0.d0
         TNW(:,:,1)=0.d0
      endif
      if(iam.eq.NPROC) then
         TNU(:,:,Nyl)=0.d0
         TNV(:,:,Nyl)=0.d0
         TNW(:,:,Nyl)=0.d0
      endif
      return
      end


      subroutine TN_prinf_helmv (p,t,yp,yv,kxm,kzm,TNU,TNV,TNW,Ra)
							      
#include "param.h"
      integer*4 jj
      complex*16 TNU(nwax+1,2*nwaz+1,Nylmem)
      complex*16 TNW(nwax+1,2*nwaz+1,Nylmem)
      complex*16 TNV(nwax+1,2*nwaz+1,Nylmem)
      complex*16 p(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 t(nwax+1,2*nwaz+1,0:Nyl+1)
      real*8 kxm(nwax+1),kzm(2*nwaz+1),yp(0:Ny+1),yv(0:Ny+1)
      real*8 Ra,cff0,cff1,cff2,a,b,cff
      complex*16 mm,nn,pint

      do j=1,Nyl
       jj=j+Nybase
       cff=1.d0/(yp(jj)-yp(jj-1))
       cff0=1.0/(yp(jj)-yp(jj-1))*(yv(jj)-yp(jj-1))
       cff1=1.0/(yp(jj)-yv(jj))*(yv(jj)-yp(jj-1))/(yp(jj)-yp(jj-1))
       cff2=1.0/(yv(jj)-yp(jj-1))*(yp(jj)-yv(jj))/(yp(jj)-yp(jj-1))
       do k=1,2*nwaz+1
        nn=(0.d0,1.d0)*kzm(k)
        do i=1,nwax+1
          pint=p(i,k,j-1)+(p(i,k,j)-p(i,k,j-1))*cff0
          mm=(0.d0,1.d0)*kxm(i)
!          TNV(i,k,j)= (p(i,k,j)-p(i,k,j-1))*cff
          TNU(i,k,j)=mm*pint
          TNW(i,k,j)=nn*pint
          TNV(i,k,j)=(p(i,k,j)-pint)*cff1+(pint-p(i,k,j-1))*cff2
        enddo
       enddo
      enddo
      if(iam.eq.1) then ! No slip!
         TNV(:,:,1)=0.d0
         TNU(:,:,1)=0.d0
         TNW(:,:,1)=0.d0
      endif
      if(iam.eq.NPROC) then
         TNV(:,:,Nyl)=0.d0
         TNU(:,:,Nyl)=0.d0
         TNW(:,:,Nyl)=0.d0
      endif
      return
      end


