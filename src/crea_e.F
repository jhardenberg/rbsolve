#include "config.h"

c ______________________________________________________________________

      subroutine Crea_E (qp,rp,sp,K2,E)
      
#include "param.h"
      real*8 K2(nwax+1,2*nwaz+1),E(nwax+1,2*nwaz+1,Nyl)
      real*8 qp(Ny),rp(Ny),sp(Ny)
      
      do k=1,2*nwaz+1
       do i=1,nwax+1
        do j=1,Nylmax1
	  E(i,k,j)=-(K2(i,k)+2.d0*rp(j+Nybase))
        enddo
       enddo
      enddo

      if(iam.eq.1) then
      do k=1,2*nwaz+1
       do i=1,nwax+1
        E(i,k,1)=E(i,k,1)+2.d0*qp(1)
       enddo
      enddo
      endif
      if (iam.eq.NPROC) then
       do k=1,2*nwaz+1
        do i=1,nwax+1
         E(i,k,Nyl-1)=E(i,k,Nyl-1)+2.d0*sp(Ny-1)
        enddo
       enddo
      endif
c     The famous pressure fix! GP& JH 07/2002
c     We impose p(1,1,Np)=0 to avoid singularity of the Helmholtz prooblem
      if(iam.eq.NPROC) E(1,1,Nyl-1)=1.d0

      return 
      end
