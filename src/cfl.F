#include "config.h"

        program main

#include "param.h"
#include "allpara.h"
        real*8 xmin,xmax,xnorm,xmino,xmaxo
        character cc*1
        character*30 str
        integer istart,iend,istep,it,jmax,jmin,jj,iargc
        complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1)
        complex*16 v(nwax+1,2*nwaz+1,0:Nylmem+1)
        complex*16 w(nwax+1,2*nwaz+1,0:Nylmem+1)
        real*8 uo(Nx,Nz,0:Nylmem+1),umax,umin
        real*8 vo(Nx,Nz,0:Nylmem+1)
        real*8 wo(Nx,Nz,0:Nylmem+1)
        character*45 file
        character*12 file0
        real*8 kxm(nwax+1),kzm(2*nwaz+1)
        complex*16 mm(nwax+1),nn(2*nwaz+1)
        real*8 K2(nwax+1,2*nwaz+1)

        real*8 yp(0:Ny+1),yv(0:Ny+1)
        real*8 qp(Ny),rp(Ny),sp(Ny),qv(Ny),rv(Ny),sv(Ny)
        real*8 tauw




#ifndef NOMPI
      print *,"Please compile with the NOMPI option!"
      stop
#endif
      iam=1

      if(iargc().lt.3) then
        print *, 'Usage: cfl istart iend istep'
        stop
      endif

       call initfour(0)
       call slice()

       call getarg(1,str)
       read(str,*) istart

       call getarg(2,str)
       read(str,*) iend

       call getarg(3,str)
       read(str,*) istep


      call read_para(1)
      call Crea_wave(kxm,kzm,K2,mm,nn,kxmin,kzmin)
      call punti (DD,ptan,qtan,yp,yv,qp,rp,sp,qv,rv,sv)



      do it=istart,iend,istep


      call readfield(u,'u',it)
      call readfield(v,'v',it)
      call readfield(w,'w',it)

      call dealias(u)
      call dealias(v)
      call dealias(w)
 
      call Fourier1(u,uo)
      call Fourier1(v,vo)
      call Fourier1(w,wo)

      call check_CFL(uo,vo,wo,it,dx,dz,yv,dt)


      enddo
   
      end


      subroutine check_CFL(uo,vo,wo,t,dx,dz,yv,dt)

#include "param.h"
      integer*4 t,jj
      real*8 uo(Nx,Nz,0:Nylmem+1)
      real*8 vo(Nx,Nz,0:Nylmem+1)
      real*8 wo(Nx,Nz,0:Nylmem+1)
      real*8 yv(0:Ny+1)
      real*8 dt,dx,dz,cflx,cflz,cfly,cflymax,umax,wmax,dy,cdy
      real*8 enx,eny,enz
      real*8 cfl(3,MAXCPU),ic(3),ene(3,MAXCPU)
      integer*4 icoord(3,MAXCPU)

*     Find local CFLs
*     this could also be done with an MPI reduce operation
*     but we want also a scalar NOMPI version

      wmax=0
      umax=0
      cflymax=0

      enx=0.
      eny=0.
      enz=0.
      do j=1,Nyl
         jj=j+Nybase
         dy=abs(yv(jj)-yv(jj-1))
         if(abs(yv(jj+1)-yv(jj)).lt.dy) dy=abs(yv(jj+1)-yv(jj))
         cdy=dt/dy
         do k=1,Nz
            do i=1,Nx  
               enx=enx+uo(i,k,j)*uo(i,k,j)
               eny=eny+vo(i,k,j)*vo(i,k,j)
               enz=enz+wo(i,k,j)*wo(i,k,j)
               if(abs(uo(i,k,j)).gt.umax) umax=abs(uo(i,k,j))
               if(abs(wo(i,k,j)).gt.wmax) wmax=abs(wo(i,k,j))
               cfly=abs(vo(i,k,j))*cdy
               if(cfly.gt.cflymax) cflymax=cfly
            enddo
         enddo
      enddo
      cfl(1,iam)=umax*dt/dx
      cfl(2,iam)=cflymax
      cfl(3,iam)=wmax*dt/dz

      ene(1,iam)=enx
      ene(2,iam)=eny
      ene(3,iam)=enz

      if(iam.eq.1) then
         cflx=0
         cfly=0
         cflz=0
         enx=0
         eny=0
         enz=0
         do i=1,NPROC
            if(cfl(1,i).gt.cflx) cflx=cfl(1,i)
            if(cfl(2,i).gt.cfly) cfly=cfl(2,i)
            if(cfl(3,i).gt.cflz) cflz=cfl(3,i)
            enx=enx+ene(1,i)
            eny=eny+ene(2,i)
            enz=enz+ene(3,i)
         enddo
         enx=enx/(2*Nx*Ny*Nz)
         eny=eny/(2*Nx*Ny*Nz)
         enz=enz/(2*Nx*Ny*Nz)

      write(*,'(i6,x,f11.6,x,a4,x,3f9.6,x,a4,x,3e13.5)')
     &         t,t*dt,'CFL:',cflx,cfly,cflz,'ENE:',enx,eny,enz



      endif
      end




