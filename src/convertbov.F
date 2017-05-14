c Converts NS/RB files to bob format
#include "config.h"

      program main

#include "param.h"


      complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 v(nwax+1,2*nwaz+1,0:Nylmem+1)
      real*8 uo(Nx,Nz,0:Nylmem+1),umax,umin
      real*8 vo(Nx,Nz,0:Nylmem+1),rt,rs,rr,dum
      real xmin,xmax
      integer istart,iend,istep,it,jmax,jmin,iargc
      character cc*1,str*30


#ifndef NOMPI 
      print *,"Please compile with the NOMPI option!"
      stop
#endif      

      iam=1

      if(iargc().lt.2) then
         print *, 'Usage: convert cc istart [iend istep]'
         stop
      endif
      call initfour(1)
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


      do it=istart,iend,istep

             if(cc.ne.'b') then
                call readfield(u,cc,it)
                call Fourier1 (u,uo)
             else
                call readfield(u,'t',it)
                do j=1,Ny
                   u(1,1,j)=0.
                enddo

                call Fourier1 (u,uo)
                call readfield(v,'s',it)
                do j=1,Ny
                   v(1,1,j)=0.
                enddo

                call Fourier1 (v,vo)
                do j=1,Ny
                 do k=1,Nz
                  do i=1,Nx
                    uo(i,k,j)=uo(i,k,j)*rr-vo(i,k,j)
                  enddo
                 enddo
                enddo
             endif

             call savebov(uo,cc,it,xmin,xmax)
      enddo

      end


      subroutine savebov (uo,ch,nsalv)

#include "param.h"

      character*1 ch
      real*8 uo(Nx,Nz,0:Nylmem+1)
      real*4 uo4(Nx,Nz,Ny)
      real center,cff
      integer*4 nsalv,ix
      character*12 file0
      integer ist

      call makenamebob(file0,ch,nsalv,'.bbv')
      print *,'saving ', file0

      do j=1,Ny
         do k=1,Nz
            do i=1,Nx
               uo4(i,k,j)=uo(i,k,j)
            enddo
         enddo
       enddo

      open(unit=12,file=file0,form='unformatted',status='unknown')
      write(12) uo4
      close(12)

      return
      end


