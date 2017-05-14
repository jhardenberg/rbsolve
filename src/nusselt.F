c Converts NS/RB files to bob format
#include "config.h"

      program main

#include "param.h"
#include "allpara.h"

      complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1,2)
      complex*16 t(nwax+1,2*nwaz+1,0:Nylmem+1,2)
      real*8 uo(Nx,Nz,0:Nylmem+1),umax,umin
      real*8 to(Nx,Nz,0:Nylmem+1),ww
       character*30 str
      real xmin,xmax,ss,ss2,yc(0:Ny+1),sst,ss2t
      integer istart,iend,istep,it,jmax,jmin,iargc,jj,idum
      real*8 kxm(nwax+1),kzm(2*nwaz+1)
      complex*16 mm(nwax+1),nn(2*nwaz+1)
      real*8 K2(nwax+1,2*nwaz+1)

      real*8 yp(0:Ny+1),yv(0:Ny+1)
      real*8 qp(Ny),rp(Ny),sp(Ny),qv(Ny),rv(Ny),sv(Ny)


#ifndef NOMPI 
      print *,"Please compile with the NOMPI option!"
      stop
#endif      

      iam=1

      if(iargc().lt.1) then
         print *, 'Usage: nusselt istart [iend istep]'
         stop
      endif
      call initfour(0)
      call slice()

      call getarg(1,str)
      read(str,*) istart

      if(iargc().gt.1) then
       call getarg(2,str)
       read(str,*) iend
       call getarg(3,str)
       read(str,*) istep
      else
        iend=istart
        istep=1
      endif


      call read_para(1)
      call Crea_wave(kxm,kzm,K2,mm,nn,kxmin,kzmin)
      call punti (DD,ptan,qtan,yp,yv,qp,rp,sp,qv,rv,sv)

      sst=0
      ss2t=0
      do it=istart,iend,istep
            call readfield(u,'v',it)
            call readfield(t,'t',it)
            call Fourier1 (u,uo)
            call Fourier1 (t,to)
	    
	    ss2=0.
          ww=0.

	    do jj=2,Ny-1
		 ss=0.
             do k=1,Nz
               do i=1,Nx
		    ss=ss+uo(i,k,jj)*(to(i,k,jj)-real(t(1,1,jj,1)))
	       enddo
	      enddo
	      ss2=ss2+(ss/(Nx*Nz))*(yp(jj)-yp(jj-1))
	      ww=ww+(yp(jj)-yp(jj-1))
	    enddo
	    ww=ww+(yp(2)-yp(1))
	    write(20,*) it,ss2/ww+1,ss2/(ww-(yp(2)-yp(1)))+1
      enddo

      end

