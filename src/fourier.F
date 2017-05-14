#define Nx 256
#define Ny 256
#define nwax 85
#define nway 85
#define NKx (nwax+1)
#define NKy (nway*2)

      subroutine fourier_c2r (u,uo)

#include "fftw3.h"

      real*8     uo(Nx,Ny)
      complex*16 u(NKx,NKy)
      complex*16 aout(Nx/2+1,Ny)
      real*8     ain(Nx,Ny)
 
      integer*8 f_plan,  i_plan
      common/FFTW3/ain,aout,f_plan,i_plan

      do k=nway+2,Ny-nway
         do i=1,nwax+1
            aout(i,k)=(0.d0,0.d0)
         enddo
      enddo

      do k=1,Ny
         do i=nwax+2,Nx/2+1
            aout(i,k)=(0.d0,0.d0)
         enddo
      enddo

      do k=1,nway+1
         do i=1,nwax+1
            aout(i,k)=u(i,k)
         end do
      end do

c - coeff.  Fourier C_(Ny-k)  (k=1,nway) -
      nw2=2*nway+1
      do k=1,nway-1
         ik1=Ny-(k-1)
         ik2=nw2-k
         do i=1,nwax+1
            aout(i,ik1)=u(i,ik2)
         end do
      enddo
       
      call dfftw_execute(i_plan)
       
      do k=1,Ny
         do i=1,Nx
            uo(i,k)=ain(i,k)
         end do
      end do

      return
     
      end      


      subroutine fourier_r2c (uo,u)
      
#include "fftw3.h"  
         
      complex*16 u(NKx,NKy)
      real*8     uo(Nx,Ny)
      real*8     ain(Nx,Ny)
      complex*16 aout(Nx/2+1,Ny)
      real*8 scale
      integer*8 f_plan,  i_plan
      common/FFTW3/ain,aout,f_plan,i_plan

      scale=1./(Nx*Ny)

      do k=1,Ny
         do i=1,Nx
            ain(i,k)=uo(i,k)
         end do
      end do

      call dfftw_execute(f_plan)

c -  C_k  (k=0, k=1,nway) -
      do k=1,nway+1 
         do i=1,nwax+1
            u(i,k)=aout(i,k)*scale
         end do
      end do
      
c -   C_(Ny-k)  (k=1,nway) -
      nw2=2*nway+1
      do k=1,nway-1
         ik2=Ny-(k-1)
         ik1=nw2-k
	 do i=1,nwax+1
            u(i,ik1)=aout(i,ik2)*scale
         end do
      end do

      return
      end


      subroutine initfour (iest) 

#include "fftw3.h"

      real*8     ain(Nx,Ny)
      complex*16 aout(Nx/2+1,Ny)

      integer*8 f_plan,  i_plan
      common/FFTW3/ain,aout,f_plan,i_plan

      if(iest.eq.0) then
      call dfftw_plan_dft_r2c_2d(f_plan,Nx,Ny,ain,aout,
     +         FFTW_MEASURE)
      call dfftw_plan_dft_c2r_2d(i_plan,Nx,Ny,aout,ain,
     +         FFTW_MEASURE)

      elseif(iest.eq.1) then
c Faster variant for analysis tools
      call dfftw_plan_dft_r2c_2d(f_plan,Nx,Ny,ain,aout,
     +        FFTW_ESTIMATE)
      call dfftw_plan_dft_c2r_2d(i_plan,Nx,Ny,aout,ain,
     +        FFTW_ESTIMATE)
      else
c Slooow variant for long runs   

      call dfftw_plan_dft_r2c_2d(f_plan,Nx,Ny,ain,aout,
     +        FFTW_PATIENT)
      call dfftw_plan_dft_c2r_2d(i_plan,Nx,Ny,aout,ain,
     +        FFTW_PATIENT)
      endif
      return

      end



