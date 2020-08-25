
#include "config.h"

c ______________________________________________________________________

      subroutine Crea_DF_prinf (qp,sp,D,F)
      
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
       enddo
      enddo
c     The famous pressure fix! GP& JH 07/2002
c     We impose p(1,1,Ny-1)=0 to avoid singularity of the Helmholtz prooblem
c      if(iam.eq.NPROC) then
c         D(1,1,Nyl-1)=0.d0
c      endif
      return
      end

      subroutine Crea_DF_prinf_u (qv,sv,D,F)
#include "param.h"
      real*8 D(nwax+1,2*nwaz+1,Nyl),F(nwax+1,2*nwaz+1,Nyl)
      real*8 qv(Ny),sv(Ny)
      integer jj

      do i=1,nwax+1
       do k=1,2*nwaz+1
        do j=1,Nyl
          jj=j+Nybase
          D(i,k,j)=2.d0*qv(jj)
          F(i,k,j)=2.d0*sv(jj)
        enddo
        if(iam.eq.1) F(i,k,1)=0.d0      ! No slip !
        if(iam.eq.NPROC) D(i,k,Nyl)=0.d0
       enddo
      enddo
      return
      end

