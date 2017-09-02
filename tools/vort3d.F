#include "config.h"


      program main

#include "param.h"

      complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 v(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 w(nwax+1,2*nwaz+1,0:Nylmem+1)
c      complex*16 wx(nwax+1,2*nwaz+1,0:Nylmem+1)
c      complex*16 uz(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 vv(nwax+1,2*nwaz+1,0:Nylmem+1)
      real*8 ss(Nx/2+1)
      integer*8 n(Nx/2+1)
      real*8 vvm
      real*8 ddm
      real*8 ttm
      real xmin,xmax
	real*8 sum,sumk
      integer istart,iend,istep,it,jmax,jmin,jj,iargc
      character*12 file0
      character*45 file

      character cc*1,str*30

      real*8 ptan,qtan
      real*8 DD,Lx,Lz,Re,invRe,dt,dx,dz
      real*8 kxmin,kzmin,Ra,Kscalar,Pr

      integer*4 nsave,ttot
      real*8 alfa(3),beta(3),gamma(3),zeta(3)

      common/ALLPARA/ptan,qtan,DD,Re,Lx,Lz,invRe,dt,dx,dz,
     +   kxmin,kzmin,nsave,ttot,alfa,beta,gamma,zeta,Ra,Kscalar,Pr

      real*8 kxm(nwax+1),kzm(2*nwaz+1),phi
      complex*16 mm(nwax+1),nn(2*nwaz+1)
      real*8 K2(nwax+1,2*nwaz+1)
      real*8 yp(0:Ny+1),yv(0:Ny+1)
      real*8 qp(Ny),rp(Ny),sp(Ny),qv(Ny),rv(Ny),sv(Ny)
      real*8 cff1,cff2

      integer ik
      real*8 pi2
      parameter (pi2=6.28318530717959d0)


#ifndef NOMPI 
      print *,"Please compile with the NOMPI option!"
      stop
#endif      

      iam=1

      if(iargc().lt.2) then
         print *, 'Usage: vort3d phi istart [iend istep]'
	 print *,'Computes vorticity component along vector 
     &   with angle phi in the y/z plane'
         stop
      endif
      call initfour(0)

      call slice()

      call getarg(1,str)
      read(str,*) phi

      call getarg(2,str)
      read(str,*) istart

      if(iargc().gt.2) then
       call getarg(3,str)
       read(str,*) iend
       call getarg(4,str)
       read(str,*) istep
      else
        iend=istart
        istep=1
      endif

      call read_para(0)
      call Crea_wave(kxm,kzm,K2,mm,nn,kxmin,kzmin)
      call punti (DD,ptan,qtan,yp,yv,qp,rp,sp,qv,rv,sv)

      do it=istart,iend,istep
         call readfield(u,'u',it)
         call readfield(v,'v',it)
         call readfield(w,'w',it)

	 do i=1,nwax+1
           do k=1,2*nwaz+1
             do j=0,Ny+1
c               uz(i,k,j)=nn(k)*u(i,k,j)
c               wx(i,k,j)=mm(i)*w(i,k,j)
c               vv(i,k,j)=wx(i,k,j)-uz(i,k,j)
         vv(i,k,j)=(mm(i)*w(i,k,j)-nn(k)*u(i,k,j))*sin(phi*pi2/360.)
             enddo
           enddo
         enddo
         do j=1,Ny
         cff1=(yv(j)-yv(j-1))/( (yv(j+1)-yv(j))*(yv(j+1)-yv(j-1)))
         cff2=(yv(j+1)-yv(j))/( (yv(j)-yv(j-1))*(yv(j+1)-yv(j-1)))
           print *, cff1,cff2,cos(phi*pi2/360.)
           do k=1,2*nwaz+1
             do i=1,nwax+1
c               vv(i,k,j)=uy(i,k,j)-vx(i,k,j)
c               vx(i,k,j)=mm(i)*v(i,k,j)

             vv(i,k,j)=vv(i,k,j)    + cos(phi*pi2/360.)*(
     &         - mm(i)*v(i,k,j)  
     &      + (u(i,k,j+1)-u(i,k,j))*cff1+(u(i,k,j)-u(i,k,j-1))*cff2 ) 
             end do
           end do
         end do

c	 vv(:,:,Ny+1)=vv(:,:,Ny)
c	 vv(:,:,0)=vv(:,:,1)
  	 call savefield(vv,'z',it,-1)
	  
      enddo
        print *,'DONE'

      end



