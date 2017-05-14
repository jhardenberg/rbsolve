#include "config.h"
#ifdef HAVE_FFTW

* Some compilers (eg Intel Fortran) require these to be defined with an additional underscore
* for linking correctly (see the Makefile)

#ifdef FFTWUNDER
#define COMPLEX_TO_REAL  dfftw_execute_
#define CREATE_PLAN_R2C      dfftw_plan_dft_r2c_2d_
#define CREATE_PLAN_C2R      dfftw_plan_dft_c2r_2d_
#else
#define COMPLEX_TO_REAL  dfftw_execute
#define CREATE_PLAN_R2C      dfftw_plan_dft_r2c_2d
#define CREATE_PLAN_C2R      dfftw_plan_dft_c2r_2d
#endif

      subroutine FOURIER1_FFTW (u,uo)

#include "param.h"

#include "fftw3.h"

      integer*4 indice,nw2,ik1,ik2,howmany,idist,odist
      real*8 uo(Nx,Nz,0:Nyl+1)
      complex*16 u(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 aout(Nx/2+1,Nz)
      real*8  ain(Nx,Nz)
 
      integer mul
      real*4 gdat0(2),gdat1(2)
      real*4 tglob1,tglob0
      common/FFTW3/ain,aout

c      tglob0=etime(gdat0)
c - coeff. di Fourier C_k  (k=0, k=1,nwaz) -
c      do mul=1,100

      do j=Nylmin0,Nylmax0

        do k=nwaz+2,Nz-nwaz
          do i=1,nwax+1
            aout(i,k)=(0.d0,0.d0)
          enddo
        enddo

        do k=1,Nz
          do i=nwax+2,Nx/2+1
            aout(i,k)=(0.d0,0.d0)
          enddo
        enddo

	do k=1,nwaz+1
          do i=1,nwax+1
            aout(i,k)=u(i,k,j)
          end do
        end do

c - coeff. di Fourier C_(Nz-k)  (k=1,nwaz) -
	nw2=2*nwaz+2
        do k=1,nwaz
          ik1=Nz-(k-1)
          ik2=nw2-k
	  do i=1,nwax+1
            aout(i,ik1)=u(i,ik2,j)
          end do
        enddo
       
        call COMPLEX_TO_REAL(i_plan)
       
c        call c2r(i_plan,aout,uo(1,1,j))
        do k=1,Nz
          do i=1,Nx
            uo(i,k,j)=ain(i,k)
          end do
        end do

      enddo


c      enddo
c       tglob1=etime(gdat1)
c       print *,'TOTAL1 (',iam,'): t:',tglob1-tglob0,
c     &                         ' u:',gdat1(1)-gdat0(1),
c     &                         ' s:',gdat1(2)-gdat0(2)

  
      return
     
      end      




      subroutine dealias (u)

#include "param.h"
#include "fftw_f77.h"

      complex*16 u(nwax+1,2*nwaz+1,0:Nyl+1)
      real*8 Kmask(nwax+1,nwaz*2+1)
      common/ALIASMASK/Kmask

      do j=0,Nyl+1
        do k=1,2*nwaz+1
          do i=1,nwax+1
            if(Kmask(i,k).eq.0) u(i,k,j)=(0.,0.)
          end do
        end do
      end do

      return
      end





      subroutine INITFOUR_FFTW (iest) 
#include "param.h"

#include "fftw3.h"
      integer iest
      real*8 ain(Nx,Nz)
      complex*16 aout(Nx/2+1,Nz)
      common/FFTW3/ain,aout


      real*8 Kmask(nwax+1,nwaz*2+1),scale
      integer*8 kk,kkx,kkz
      integer*8 cffx,cffz
      common/ALIASMASK/Kmask


      scale=1./(Nx*Nz)
      cffx=(nwax*nwax)
      cffz=(nwaz*nwaz)
      do k=1,nwaz+1
         kkz=k-1
         do i=1,nwax+1
            kkx=i-1
            Kmask(i,k)=scale
            if((kkx*kkx*cffz+kkz*kkz*cffx).gt.(cffx*cffz)) Kmask(i,k)=0.0
         enddo
      enddo

      do k=nwaz+2,nwaz*2+1
         kkz=k-(nwaz*2+1)-1
        print *,kkz
         do i=1,nwax+1
            kkx=i-1
            Kmask(i,k)=scale
            if((kkx*kkx*cffz+kkz*kkz*cffx).gt.(cffx*cffz)) Kmask(i,k)=0.0
         enddo
      enddo

      if(iest.eq.0) then
      call CREATE_PLAN_R2C(f_plan,Nx,Nz,ain,aout,
     +         FFTW_ESTIMATE)
      call CREATE_PLAN_C2R(i_plan,Nx,Nz,aout,ain,
     +         FFTW_ESTIMATE)
      elseif(iest.eq.1) then
c Faster variant for analysis tools
      call CREATE_PLAN_R2C(f_plan,Nx,Nz,ain,aout,
     +        FFTW_MEASURE)
      call CREATE_PLAN_C2R(i_plan,Nx,Nz,aout,ain,
     +        FFTW_MEASURE)
      elseif(iest.eq.2) then
c Slooow variant for long runs   
      call CREATE_PLAN_R2C(f_plan,Nx,Nz,ain,aout,
     +        FFTW_PATIENT)
      call CREATE_PLAN_C2R(i_plan,Nx,Nz,aout,ain,
     +        FFTW_PATIENT)
      elseif(iest.eq.3) then
c Even slowe variant for long runs   
      call CREATE_PLAN_R2C(f_plan,Nx,Nz,ain,aout,
     +        FFTW_EXHAUSTIVE)
      call CREATE_PLAN_C2R(i_plan,Nx,Nz,aout,ain,
     +        FFTW_EXHAUSTIVE)

      endif
      return

      end


#endif

c	Some stupid compilers (HPUX) need this!
	subroutine dummy1()
	end
