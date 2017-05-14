#include "config.h"

c----------------------------------------------------------------------
      subroutine conv (uo,vo,wo,yv,mm,nn,Nu_n,Nv_n,Nw_n,
     &                  w1,w2,w3)

*     w1,w2,w3 are workspace arrays

#include "param.h"

      complex*16 Nu_n(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 Nv_n(nwax+1,2*nwaz+1,0:Nyl+1) 
      complex*16 Nw_n(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 w1(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 w2(nwax+1,2*nwaz+1,0:Nylmem+1) 
      complex*16 w3(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 mm(nwax+1),nn(2*nwaz+1) ,der1,nn1

      real*8 uo(Nx,Nz,0:Nylmem+1)
      real*8 vo(Nx,Nz,0:Nylmem+1)
      real*8 wo(Nx,Nz,0:Nylmem+1)

      real*8 yv(0:Ny+1)
      real*8 cff1,cff2
      integer jj

#define UU w1
#define WW w2
#define UW w3

c     Compute the double products
   
      call Fourier2 (uo,uo,UU,Nylmin,Nylmax)
      call Fourier2 (wo,wo,WW,Nylmin,Nylmax)
      call Fourier2 (uo,wo,UW,Nylmin,Nylmax)

c     Convective terms: sum of the double products

      do j=Nylmin,Nylmax
       do k=1,2*nwaz+1
         nn1=nn(k)
         do i=1,nwax+1
       Nu_n(i,k,j)=mm(i)*UU(i,k,j)+nn1*UW(i,k,j)
       Nw_n(i,k,j)=mm(i)*UW(i,k,j)+nn1*WW(i,k,j)
        end do
       end do
      end do

#define UV w1
#define VW w2
#define VV w3


*     Computes double products including v also on the ghost levels

      call Fourier2 (uo,vo,UV,Nylmin0,Nylmax0)
      call Fourier2 (vo,vo,VV,Nylmin0,Nylmax0)
      call Fourier2 (vo,wo,VW,Nylmin0,Nylmax0)

c     Convective terms: sum of the double products

      do j=Nylmin,Nylmax
       jj=j+Nybase
       cff1=(yv(jj)-yv(jj-1))/( (yv(jj+1)-yv(jj))*(yv(jj+1)-yv(jj-1)))
       cff2=(yv(jj+1)-yv(jj))/( (yv(jj)-yv(jj-1))*(yv(jj+1)-yv(jj-1)))

       do k=1,2*nwaz+1
         nn1=nn(k)
         do i=1,nwax+1
       Nu_n(i,k,j)=Nu_n(i,k,j)+
     +     (UV(i,k,j+1)-UV(i,k,j))*cff1+(UV(i,k,j)-UV(i,k,j-1))*cff2
       Nw_n(i,k,j)=Nw_n(i,k,j)+
     +     (VW(i,k,j+1)-VW(i,k,j))*cff1+(VW(i,k,j)-VW(i,k,j-1))*cff2     
       Nv_n(i,k,j)=mm(i)*UV(i,k,j)+nn1*VW(i,k,j)+
     +     (VV(i,k,j+1)-VV(i,k,j))*cff1+(VV(i,k,j)-VV(i,k,j-1))*cff2      
        end do
       end do
      end do

      
      return
      end

      

      subroutine conv_scalar (uo,vo,wo,to,yv,mm,nn,Nt_n,ut,vt,wt)
* Computes advective terms

* ut,vt,wt are temporary storage

#include "param.h"

      complex*16 Nt_n(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 ut(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 vt(nwax+1,2*nwaz+1,0:Nylmem+1) 
      complex*16 wt(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 mm(nwax+1),nn(2*nwaz+1) ,der1,nn1

      real*8 uo(Nx,Nz,0:Nylmem+1)
      real*8 vo(Nx,Nz,0:Nylmem+1)
      real*8 wo(Nx,Nz,0:Nylmem+1)
      real*8 to(Nx,Nz,0:Nylmem+1)

      real*8 yv(0:Ny+1)
      real*8 cff1,cff2
      integer jj

c     Compute the double products
   
      call Fourier2 (uo,to,ut,Nylmin,Nylmax)
      call Fourier2 (wo,to,wt,Nylmin,Nylmax)

*     Computes double products including v also on the ghost levels
      
      call Fourier2 (vo,to,vt,Nylmin0,Nylmax0)

c     Convective terms: sum of the double products

      do j=Nylmin,Nylmax
       jj=j+Nybase
       cff1=(yv(jj)-yv(jj-1))/( (yv(jj+1)-yv(jj))*(yv(jj+1)-yv(jj-1)))
       cff2=(yv(jj+1)-yv(jj))/( (yv(jj)-yv(jj-1))*(yv(jj+1)-yv(jj-1)))

       do k=1,2*nwaz+1
         nn1=nn(k)
         do i=1,nwax+1
       Nt_n(i,k,j)=mm(i)*ut(i,k,j)+nn1*wt(i,k,j)+
     +     (vt(i,k,j+1)-vt(i,k,j))*cff1+(vt(i,k,j)-vt(i,k,j-1))*cff2
        end do
       end do

      end do

      return
      end

