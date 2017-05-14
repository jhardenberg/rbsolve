#include "config.h"
#ifdef HAVE_FFTW

* Some compilers (eg Intel Fortran Compiler) require these to be defined with an additional underscore
* for linking correctly (see the Makefile)

#ifdef FFTWUNDER
#define REAL_TO_COMPLEX dfftw_execute_
#else
#define REAL_TO_COMPLEX dfftw_execute
#endif
   
      subroutine FOURIER2_FFTW (uo,vo,uv,Nli,Nls)
      
#include "param.h"  
#include "fftw3.h"  
         
      integer*4  ik2 , nw2,ik1,Nli,Nls
      real*8 uo(Nx,Nz,0:Nyl+1),vo(Nx,Nz,0:Nyl+1)
      complex*16 uv(nwax+1,2*nwaz+1,0:Nyl+1)
      real*8 scale
      real*4 tdat0(2),tdat(2)
      real*4 ttot1,ttot0
      real*8 ain(Nx,Nz)
      complex*16 aout(Nx/2+1,Nz)

      integer mul
    
      real*4 gdat0(2),gdat1(2)
      real*4 tglob1,tglob0
      common/FFTW3/ain,aout	

      real*8 Kmask(nwax+1,nwaz*2+1)
      common/ALIASMASK/Kmask


c     tglob0=etime(gdat0)
c      do mul=1,100
*     The FFTW transforms are not normalized ! fix this
      scale=1./(Nx*Nz)

      do j=Nli,Nls

        do k=1,Nz
          do i=1,Nx
            ain(i,k)=uo(i,k,j)*vo(i,k,j)
          end do
        end do


c     call REAL_TO_COMPLEX(f_plan,ain,aout)


      call REAL_TO_COMPLEX(f_plan)

c -  C_k  (k=0, k=1,nwaz) -
        do k=1,nwaz+1 
          do i=1,nwax+1
            uv(i,k,j)=aout(i,k)*Kmask(i,k)
          end do
        end do
      
c -   C_(Nz-k)  (k=1,nwaz) -
        nw2=2*nwaz+2

        do k=1,nwaz
          ik2=Nz-(k-1)
	  ik1=nw2-k
	  do i=1,nwax+1
            uv(i,ik1,j)=aout(i,ik2)*Kmask(i,ik1)
          end do
        end do

      end do


c      enddo
c       tglob1=etime(gdat1)
c       print *,'TOTAL (',iam,'): t:',tglob1-tglob0,
c     &                         ' u:',gdat1(1)-gdat0(1),
c     &                         ' s:',gdat1(2)-gdat0(2)


            
      return
      end


#endif


c       Some stupid compilers (HPUX) need this!
        subroutine dummy2()
        end
