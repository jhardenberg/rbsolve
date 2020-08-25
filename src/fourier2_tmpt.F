
#include "config.h"
c ______________________________________________________________________
	       
      subroutine FOURIER2_TMPT (uo,vo,uv,Nli,Nls)
      
#include "param.h"  
         
      integer*4  indice , ik2 , nw2,ik1,Nli,Nls
      real*8 a(Nx+2,Nz+2)
      real*8 uo(Nx,Nz,0:Nyl+1),vo(Nx,Nz,0:Nyl+1)
      complex*16 uv(nwax+1,2*nwaz+1,0:Nyl+1)



      do k=1,Nz+2
          do i=1,Nx+2
             a(i,k)=0.d0
          end do
      end do
	
      do j=Nli,Nls
              
        do k=1,Nz
          do i=1,Nx
            a(i,k)=uo(i,k,j)*vo(i,k,j)
          end do
        end do

        call rffd(a,Nx,Nz)

c - coeff. di Fourier C_k  (k=0, k=1,nwaz) -
        do k=1,nwaz+1 
        indice=0
          do i=1,nwax+1
            indice=indice+2
            uv(i,k,j)=dcmplx(a(indice-1,k),a(indice,k))
          end do
        end do
      
c - coeff. di Fourier C_(Nz-k)  (k=1,nwaz) -
        nw2=2*nwaz+2
        do k=1,nwaz
          indice=0
          ik2=Nz-(k-1)
	  ik1=nw2-k
	  do i=1,nwax+1
            indice=indice+2
            uv(i,ik1,j)=
     1        dcmplx(a(indice-1,ik2),a(indice,ik2))
          end do
        end do

      end do

            
      return
      end

c le subroutine 'rffp', 'rffd', 'rrfi' e tutte le altre necessarie
c si trovano nel file 'fourier1.f'.

      subroutine FOURIER2_TMPT_SINGLE2D (uo,uv)
      
#include "param.h"  
         
      integer*4  indice , ik2 , nw2,ik1,Nli,Nls
      real*8 a(Nx+2,Nz+2)
      real*8 uo(Nx,Nz)
      complex*16 uv(nwax+1,2*nwaz+1)

      do k=1,Nz+2
          do i=1,Nx+2
             a(i,k)=0.d0
          end do
      end do
	
      do k=1,Nz
        do i=1,Nx
          a(i,k)=uo(i,k)
        end do
      end do

      call rffd(a,Nx,Nz)

c - coeff. di Fourier C_k  (k=0, k=1,nwaz) -
      do k=1,nwaz+1 
        indice=0
          do i=1,nwax+1
            indice=indice+2
            uv(i,k)=dcmplx(a(indice-1,k),a(indice,k))
        end do
      end do
      
c - coeff. di Fourier C_(Nz-k)  (k=1,nwaz) -
      nw2=2*nwaz+2
      do k=1,nwaz
        indice=0
        ik2=Nz-(k-1)
        ik1=nw2-k
        do i=1,nwax+1
          indice=indice+2
          uv(i,ik1)=dcmplx(a(indice-1,ik2),a(indice,ik2))
        end do
      end do

      return
      end

