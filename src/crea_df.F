
#include "config.h"

c ______________________________________________________________________

      subroutine Crea_DF (qp,sp,D,F)
      
#include "param.h"
      real*8 D(nwax+1,2*nwaz+1,Nyl),F(nwax+1,2*nwaz+1,Nyl)
      real*8 qp(Ny),sp(Ny)
      integer jj

      do i=1,nwax+1
      do k=1,2*nwaz+1
      do j=1,Nylmax1
	  jj=j+Nybase 
	  D(i,k,j)=2.d0*qp(jj)
	  F(i,k,j)=2.d0*sp(jj)
      enddo
	    
      if(iam.eq.1) D(i,k,1)=0.d0
      if(iam.eq.NPROC) F(i,k,Nyl-1)=0.d0
     
      enddo
      enddo
c     The famous pressure fix! GP& JH 07/2002
c     We impose p(1,1,Np)=0 to avoid singularity of the Helmholtz prooblem
      if(iam.eq.NPROC) then
         D(1,1,Nyl-1)=0.d0
cxxx  This assumes that > 1 slice is assigned to each CPU!
c       F(1,1,Nyl-2)=0.d0
c       JH11 commented this!
      endif
      return
      end
