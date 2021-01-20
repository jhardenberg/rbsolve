c Converts NS/RB files to bob format
#include "config.h"

      program main

#include "param.h"


      complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1,2)
      complex*16 t(nwax+1,2*nwaz+1,0:Nylmem+1,2)
      real*8 uo(Nx,Nz,0:Nylmem+1),umax,umin
      real*8 uom(Nx,Nz), uov(Nx,Nz)
      real*8 to(Nx,Nz,0:Nylmem+1)
      real xmin,xmax
      integer istart,iend,istep,it,jmax,jmin,jj,iargc,nstep
      character*12 file0
      character*45 file
      character cc*1,str*30
      logical fmean


#ifndef NOMPI 
      print *,"Please compile with the NOMPI option!"
      stop
#endif      

      iam=1
      fmean=.false.
      uom(:,:)=0.d0

      if(iargc().lt.3) then
         print *,'Usage: sectionh cc y istart [iend istep]'
         print *,' Use a negative istep to produce mean and variance'
         print *,' instead of individual frames over time'
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
      if (istep.lt.0) then
         istep=-istep
         fmean=.true.
      endif
      nstep=0 
      do it=istart,iend,istep

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
        elseif(cc=='f') then
         call readfield(u,'v',it)
         call Fourier1 (u,uo)
         call readfield(t,'t',it)
         call Fourier1 (t,to)
         to(1,1,:)=0.0
         do j=1,Ny
               do k=1,Nz
                 do i=1,Nx
                   uo(i,k,j)=uo(i,k,j)*to(i,k,j)
                 enddo
               enddo
         enddo
        else
         call readfield(u,cc,it)
         call Fourier1 (u,uo)
        endif

        if (fmean) then
           uom(:,:)=uom(:,:)+uo(:,:,jj)
           uov(:,:)=uov(:,:)+uo(:,:,jj)*uo(:,:,jj)
        else
           call makenamebob(file0,cc,it,'.dat')
	   open(unit=20,file=file0,form='formatted')
           do k=1,Nz
              do i=1,Nx
	         write(20,'(e16.7)') uo(i,k,jj)
              enddo
           enddo
	   close(20)
        endif
        nstep=nstep+1
      enddo
      if (fmean) then
         uom(:,:)=uom(:,:)/nstep
         uov(:,:)=uov(:,:)/nstep
         uov(:,:)=uov(:,:)-uom(:,:)*uom(:,:)
         call makenamebob(file0,cc,it,'.dam')
         open(unit=20,file=file0,form='formatted')
         do k=1,Nz
            do i=1,Nx
               write(20,'(e16.7)') uom(i,k)
            enddo
         enddo
         close(20)
         call makenamebob(file0,cc,it,'.dav')
         open(unit=20,file=file0,form='formatted')
         do k=1,Nz
            do i=1,Nx
               write(20,'(e16.7)') uov(i,k)
            enddo
         enddo
         close(20)
      endif

      end








