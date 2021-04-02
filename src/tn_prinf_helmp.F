#include "config.h"

#ifdef PRINF
c ______________________________________________________________________
      subroutine TN_prinf_helmp (t,yv,yp,qp,sp,kxm,kzm,TN,Ra,pbot,ptop)
							      
#include "param.h"
      integer*4 jj
      real*8 qp(Ny),sp(Ny),pbot,ptop
      complex*16 t(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 TN(nwax+1,2*nwaz+1,Nyl)
      real*8 kxm(nwax+1),kzm(2*nwaz+1),yv(0:Ny+1),yp(0:Ny+1)
      real*8 cff1,Ra
     
      if ((pbot == 0.d0) .and. (ptop == 0.d0)) then
       do j=1,Nylmax1
         jj=j+Nybase
        cff1=1.0/(yv(jj+1)-yv(jj))*Ra
        do k=1,2*nwaz+1
         do i=1,nwax+1
  	  TN(i,k,j)=(t(i,k,j+1)-t(i,k,j))*cff1
         enddo
        enddo
       enddo
!  Impose vertical derivative of RaT for the mean
!  Impose 0 gradient at top and bottom for other modes
       if(iam.eq.1) then
!    TN(:,:,1)=TN(:,:,1)+2.d0*(yp(1)-yp(0))*qp(1)*pbot
!        cff1=1.0/(yv(2)-yv(1))*Ra
        TN(1,1,1)=TN(1,1,1)+
     &  2.d0*Ra*TEMPERATURE_BOTTOM*(yp(1)-yp(0))*qp(1)
       endif
       if (iam.eq.NPROC) then
!        TN(:,:,Nyl-1)=TN(:,:,Nyl-1)+2.d0*(yp(Ny)-yp(Ny-1))*qp(Ny-1)*ptop
!         cff1=1.0/(yv(Ny)-yv(Ny-1))*Ra
         TN(1,1,Nyl-1)=TN(1,1,Nyl-1)+
     &  2.d0*Ra*TEMPERATURE_BOTTOM*(yp(Ny)-yp(Ny-1))*sp(Ny-1)
       endif
      else 
       TN(:,:,:)=(0.d0,0.d0)
!  Impose 0 gradient at top and bottom for the mean
!  Impose pbot and ptop for other modes
       if(iam.eq.1) then
        TN(:,:,1)=TN(:,:,1)+2.d0*(yp(1)-yp(0))*qp(1)*pbot   
        TN(1,1,1)=0.d0
       endif
       if (iam.eq.NPROC) then
        TN(:,:,Nyl-1)=TN(:,:,Nyl-1)-2.d0*(yp(Ny)-yp(Ny-1))*qp(Ny-1)*ptop
        TN(1,1,Nyl-1)=0.d0
       endif
      endif

      return
      end
#endif
