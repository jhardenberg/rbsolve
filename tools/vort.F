#include "config.h"


      program main

#include "param.h"

      complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 t(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 v(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 ux(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 vx(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 uy(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 vy(nwax+1,2*nwaz+1,0:Nylmem+1)
      real*8 ss(Nx/2+1)
      integer*8 n(Nx/2+1)
      real*8 uo(Nx,Nz,0:Nylmem+1)
      real*8 vo(Nx,Nz,0:Nylmem+1)
      real*8 uxo(Nx,Nz,0:Nylmem+1)
      real*8 vxo(Nx,Nz,0:Nylmem+1)
      real*8 uyo(Nx,Nz,0:Nylmem+1)
      real*8 vyo(Nx,Nz,0:Nylmem+1)
      real*8 to(Nx,Nz,0:Nylmem+1)
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

      real*8 kxm(nwax+1),kzm(2*nwaz+1)
      complex*16 mm(nwax+1),nn(2*nwaz+1)
      real*8 K2(nwax+1,2*nwaz+1)
	integer ik



#ifndef NOMPI 
      print *,"Please compile with the NOMPI option!"
      stop
#endif      

      iam=1

      if(iargc().lt.2) then
         print *, 'Usage: kmax y istart [iend istep]'
	 print *,'Computes the horizontal spectrum of variable yy'
	 print *, 'y<=0 means all levels'
         stop
      endif
      call initfour(0)

      call slice()

      call getarg(1,str)
      read(str,*) jj

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

      do it=istart,iend,istep
      call makenamebob(file0,'z',it,'.dat')
        print *,file0
            call readfield(u,'u',it)
            call readfield(v,'w',it)
            call readfield(t,'t',it)
	    do i=1,nwax+1
              do k=1,2*nwaz+1
                 do j=0,Ny+1
                  ux(i,k,j)=mm(i)*u(i,k,j)
                  uy(i,k,j)=nn(k)*u(i,k,j)
                  vx(i,k,j)=mm(i)*v(i,k,j)
                  vy(i,k,j)=nn(k)*v(i,k,j)
                 enddo
              enddo
           enddo
      open(unit=22,file=file0,form='formatted')
	  call Fourier1 (ux,uxo)
	  call Fourier1 (uy,uyo)
	  call Fourier1 (vx,vxo)
	  call Fourier1 (vy,vyo)
	  call Fourier1 (t,to)
        print *,'writing'    
            do k=1,Nz
               do i=1,Nx
                       vvm=vxo(i,k,jj)-uyo(i,k,jj)
                       ddm=vyo(i,k,jj)+uxo(i,k,jj)
                       ttm=to(i,k,jj)
	        write(22,*) vvm,ddm,ttm
	       enddo
	     enddo

	  close(22) 
      enddo

      end



