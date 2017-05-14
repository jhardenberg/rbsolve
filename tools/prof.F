c Computes mean,variance,skewness and kurtosis of 
c RB output files

#include "config.h"

      program main

#include "param.h"


      complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 v(nwax+1,2*nwaz+1,0:Nylmem+1)
      real*8 uo(Nx,Nz,0:Nylmem+1),umax,umin,yv(0:Ny+1)
      real*8 vo(Nx,Nz,0:Nylmem+1),dum,rt,rs,rr
      real*8 xmean(Ny),xvar,xskew,xkurt,xx,xerr,x1,xstd
      REAL*8 xtmean(Ny),xtvar(Ny),xtskew(Ny),xtkurt(Ny)
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

      if(cc.eq.'b') then
       open(unit=10,file='param0')
       do j=1,12
         read(10,*) dum
       enddo
       read(10,*) rt
       read(10,*) rs
       rr=rt/rs
       close(10)
       print *,rt,rs,rr
      endif

      do j=1,Ny
         xtmean(j)=0.
         xtvar(j)=0.
         xtskew(j)=0.
         xtkurt(j)=0.
      enddo

      nn=Nx*Nz

	  do it=istart,iend,istep

	     if(cc.ne.'b') then
                call readfield(u,cc,it)
                call Fourier1 (u,uo)
                do j=1,Ny      
                  xmean(j)=u(1,1,j)
	        enddo
	     else
                call readfield(u,'t',it)
                call Fourier1 (u,uo)
		call readfield(v,'s',it)
                call Fourier1 (v,vo)
                do j=1,Ny      
                 do k=1,Nz
                  do i=1,Nx
                    uo(i,k,j)=uo(i,k,j)*rr-vo(i,k,j)
	          enddo
	         enddo
                enddo
                do j=1,Ny      
                  xmean(j)=u(1,1,j)*rr-v(1,1,j)
	        enddo
	     endif
              
             do j=1,Ny      
               xvar=0.
               xskew=0.
               xkurt=0.
               xerr=0.
               do k=1,Nz
                 do i=1,Nx
                    x1=uo(i,k,j)-xmean(j)
                    xerr=xerr+x1
                    xx=x1*x1
                    xvar=xvar+xx
                    xx=xx*x1
                    xskew=xskew+xx
                    xx=xx*x1
                    xkurt=xkurt+xx
                 enddo
               enddo
               xvar=(xvar-xerr*xerr/nn)/(nn-1)
               xstd=sqrt(xvar)
               if(xvar.ne.0.) then
                  xskew=xskew/(nn*xvar*xstd) 
                  xkurt=xkurt/(nn*xvar*xvar)-3.0 
               else
                  xskew=0.
                  xkurt=0.
               endif
                 write(33,*) yv(j),xmean(j)
                 xtmean(j)=xtmean(j)+xmean(j)
                 xtvar(j)=xtvar(j)+xvar
                 xtskew(j)=xtskew(j)+xskew
                 xtkurt(j)=xtkurt(j)+xkurt
	     write(22,'(f10.2)') xvar
	  enddo
          write(34,*) xmean(1),xmean(Ny)

      enddo

      if(istart.ne.iend) then
         call makenamebob(file0,cc,istart,'.prm')
         nn=(iend-istart)/istep+1
      else
         call makenamebob(file0,cc,istart,'.prf')
         nn=1
      endif

      open(unit=20,file=file0)
      do j=1,Ny
         write(20,'(5e16.7)') yv(j),xtmean(j)/nn,
     +         xtvar(j)/nn,xtskew(j)/nn,xtkurt(j)/nn 
      enddo
	
      close(20)
      close(22)
      end








