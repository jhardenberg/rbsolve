#include "config.h"

c ______________________________________________________________________

      subroutine Crea_E_prinf (qp,rp,sp,yp,K2,E,Ra)
      
#include "param.h"
      real*8 K2(nwax+1,2*nwaz+1),E(nwax+1,2*nwaz+1,Nyl)
      real*8 qp(Ny),rp(Ny),sp(Ny),Ra
      real*8 yp(0:Ny+1)
      
      do k=1,2*nwaz+1
       do i=1,nwax+1
        do j=1,Nylmax1
	  E(i,k,j)=-(K2(i,k)+2.d0*rp(j+Nybase))
        enddo
       enddo
      enddo

c Impose zero gradient at top and bottom
      if(iam.eq.1) then
      do k=1,2*nwaz+1
       do i=1,nwax+1
        E(i,k,1)=E(i,k,1)+2.d0*qp(1)
       enddo
      enddo
c Impose zero gradient at top and Ra at bottom for the mean
      E(1,1,1)=E(1,1,1)-qp(1)*Ra*(yp(1)-yp(0))

      endif
      if (iam.eq.NPROC) then
       do k=1,2*nwaz+1
        do i=1,nwax+1
         E(i,k,Nyl-1)=E(i,k,Nyl-1)+2.d0*sp(Ny-1)
        enddo
       enddo
      endif
      E(1,1,Nyl-1)=E(1,1,Nyl-1)+sp(1)*Ra*(yp(1)-yp(0))
c     The famous pressure fix! GP& JH 07/2002
c     We impose p(1,1,Np)=0 to avoid singularity of the Helmholtz prooblem
c      if(iam.eq.NPROC) E(1,1,Nyl-1)=1.d0

      return 
      end


c ______________________________________________________________________

      subroutine Crea_E_prinf_u (qv,rv,sv,yv,K2,E,Ra)
      
#include "param.h"
      real*8 K2(nwax+1,2*nwaz+1),E(nwax+1,2*nwaz+1,Nyl)
      real*8 qv(Ny),rv(Ny),sv(Ny),Ra
      real*8 yv(0:Ny+1)
      
      do k=1,2*nwaz+1
       do i=1,nwax+1
        do j=1,Nylmax1
	  E(i,k,j)=-(K2(i,k)+2.d0*rv(j+Nybase))
        enddo
       enddo
      enddo

      return 
      end


