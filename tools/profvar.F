c Computes the variance of a profile in time (used for the rb_shear work)
c RB output files

#include "config.h"

      program main

#include "param.h"


      complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 v(nwax+1,2*nwaz+1,0:Nylmem+1)
      real*8 uo(Nx,Nz,0:Nylmem+1),umax,umin,yv(0:Ny+1)
      real*8 vo(Nx,Nz,0:Nylmem+1),dum,rt,rs,rr
      real*8 xmean(Ny),xvar,xskew,xkurt,xx,xerr,x1,xstd
      REAL*8 xtmean(Ny),xtvar(Ny),xtskew(Ny),xtkurt(Ny),xm,xv
      integer istart,iend,istep,it,jmax,jmin,nn,jj,iargc
      character cc*1,str*30,sa*2,sb*3
      character*12 file0

#ifndef NOMPI 
      print *,"Please compile with the NOMPI option!"
      stop
#endif      

      iam=1

      if(iargc().lt.2) then
         print *, 'Usage: prof cc istart [iend istep]'
         stop
      endif
      call initfour(0)

      open(unit=10,file='coord')
      do j=0,Ny+1
         read(10,*) sa,jj,sb,yv(j)
      enddo
      close(10)


      call slice()

      call getarg(1,cc)

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

      nn=Nx*Nz
      call makenamebob(file0,cc,istart,'.var')
      open(unit=20,file=file0)

	  do it=istart,iend,istep
             call readfield(u,cc,it)
             call Fourier1 (u,uo)
             xm=0.
             do j=1,Ny      
                xm=xm+u(1,1,j)
	     enddo
             xm=xm/Ny
             xv=0. 
             do j=1,Ny      
               xv=xv+(u(1,1,j)-xm)**2
	     enddo
	     xv=xv/Ny
	     write(20,'(i7,x,e16.7)') it,xv
      	   enddo

      close(20)
      end








