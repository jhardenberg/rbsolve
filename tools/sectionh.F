c Converts NS/RB files to bob format
#include "config.h"

      program main

#include "param.h"


      complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1,2)
      complex*16 t(nwax+1,2*nwaz+1,0:Nylmem+1,2)
      real*8 uo(Nx,Nz,0:Nylmem+1),umax,umin
      real*8 to(Nx,Nz,0:Nylmem+1)
      real xmin,xmax
      integer istart,iend,istep,it,jmax,jmin,jj,iargc
      character*12 file0
	character*45 file
      character cc*1,str*30


#ifndef NOMPI 
      print *,"Please compile with the NOMPI option!"
      stop
#endif      

      iam=1

      if(iargc().lt.3) then
         print *, 'Usage: sectionh cc y istart [iend istep]'
         stop
      endif
      call initfour(0)
      call slice()


      call getarg(1,cc)

      call getarg(2,str)
      read(str,*) jj

      call getarg(3,str)
      read(str,*) istart

      if(iargc().gt.3) then
       call getarg(4,str)
       read(str,*) iend
       call getarg(5,str)
       read(str,*) istep
      else
        iend=istart
        istep=1
      endif

      do it=istart,iend,istep
             call makenamebob(file0,cc,it,'.dat')

        if(cc=='b') then
         call readfield(u,'s',it)
         call Fourier1 (u,uo)
         call readfield(t,'t',it)
         call Fourier1 (t,to)
         do j=1,Ny
               do k=1,Nz
                 do i=1,Nx
                   uo(i,k,j)=to(i,k,j)*1.025-uo(i,k,j)
                 enddo
               enddo
          enddo
        else
         call readfield(u,cc,it)
         call Fourier1 (u,uo)
        endif

	file=file0
	
	  open(unit=20,file=file,form='formatted')
               do k=1,Nz
                 do i=1,Nx
		write(20,'(e16.7)') uo(i,k,jj)
                 enddo
               enddo
	close(20)

      enddo

      end








