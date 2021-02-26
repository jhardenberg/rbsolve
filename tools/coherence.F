#include "config.h"

      program main

#include "param.h"

      complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 v(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 w(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 t(nwax+1,2*nwaz+1,0:Nylmem+1)
c      real*8 to(Nx,Ny)
      real*8 ss(Nx),ssx(Nx),ssy(Nx)
      real*8 umax,umin,dk
      real*8 vmax,vmin
      real xmin,xmax
      integer istart,iend,istep,it,jmax,jmin,jj,iargc
      character*12 file0
      character*45 file

      character cc*1,str*30,cc2*1

      real*8 ptan,qtan
      real*8 DD,Lx,Lz,Re,invRe,dt,dx,dz
      real*8 kxmin,kzmin,Ra,Kscalar,Pr,sum

      integer*4 nsave,ttot,nn(Nx),nnx(Nx),nny(Nx),nmin,nmax
      real*8 alfa(3),beta(3),gamma(3),zeta(3)

      common/ALLPARA/ptan,qtan,DD,Re,Lx,Lz,invRe,dt,dx,dz,
     +   kxmin,kzmin,nsave,ttot,alfa,beta,gamma,zeta,Ra,Kscalar,Pr

      real*8 kxm(nwax+1),kzm(2*nwaz+1),phi(nwax+1,2*nwaz+1),phi_r(Nx)
      complex*16 mm(nwax+1),nnn(2*nwaz+1)
      complex*16 st1(nwax+1,2*nwaz+1), sv1(nwax+1,2*nwaz+1)
      complex*16 sm1(nwax+1,2*nwaz+1), st1r(Nx),sv1r(Nx),sm1r(Nx)

      real*8 K2(nwax+1,2*nwaz+1),dum,rt,rs,rr
      integer ik,ntim



#ifndef NOMPI 
      print *,"Please compile with the NOMPI option!"
      stop
#endif      

      iam=1

      if(iargc().lt.3) then
         print *, 'Usage: cospect cc1 cc2 y istart [iend istep]'
	 print *,'Computes the horizontal spectrum of variable yy'
	 print *, 'y<=0 means all levels'
         stop
      endif
      call slice()


      call getarg(1,cc)
      call getarg(2,cc2)

      call getarg(3,str)
      read(str,*) jj

      call getarg(4,str)
      read(str,*) istart

      if(iargc().gt.4) then
       call getarg(5,str)
       read(str,*) iend
       call getarg(6,str)
       read(str,*) istep
      else
        iend=istart
        istep=1
      endif

      call read_para(0)
      call Crea_wave(kxm,kzm,K2,mm,nnn,kxmin,kzmin)

c  use x scale as granularity	
	dk=kxm(2)
	print *,dk

	sum=0.

      do i=1,Nx
	ss(i)=0
        phi_r(i)=0.
        sm1r(i)=(0.d0,0.d0)
        sv1r(i)=(0.d0,0.d0)
        st1r(i)=(0.d0,0.d0)
      enddo
      do i=1,nwax+1
        do k=1,2*nwaz+1
           st1(i,k)=(0.d0,0.d0)
           sm1(i,k)=(0.d0,0.d0)
           sv1(i,k)=(0.d0,0.d0)
        enddo
      enddo  
      ntim=0
      do it=istart,iend,istep
        call readfield(t,cc,it)
        call readfield(v,cc2,it)
        ntim=ntim+1
        if(jj.le.0) then
          do j=-jj,Ny+1+jj
            do i=1,nwax+1
              do k=1,2*nwaz+1
                 st1(i,k)=t(i,k,j)*conjg(t(i,k,j))*2
                 sv1(i,k)=v(i,k,j)*conjg(v(i,k,j))*2
                 sm1(i,k)=t(i,k,j)*conjg(v(i,k,j))
     &                   +v(i,k,j)*conjg(t(i,k,j))
              enddo
            enddo
          enddo
        else
          do i=1,nwax+1
            do k=1,2*nwaz+1
               st1(i,k)=t(i,k,jj)*conjg(t(i,k,jj))*2
               sv1(i,k)=v(i,k,jj)*conjg(v(i,k,jj))*2
               sm1(i,k)=t(i,k,jj)*conjg(v(i,k,jj))
     &                 +v(i,k,jj)*conjg(t(i,k,jj))
            enddo
          enddo
        endif
      enddo
      do i=1,nwax+1
         do k=1,2*nwaz+1
           ik=sqrt(K2(i,k))/dk+1.5
           sm1r(ik)=sm1r(ik)+sm1(i,k)
           sv1r(ik)=sv1r(ik)+sv1(i,k)
           st1r(ik)=st1r(ik)+st1(i,k)
           nn(ik)=nn(ik)+1
          enddo
      enddo
      do ik=1,Nx/2+1
           phi_r(ik)=phi_r(ik)+abs(sm1r(ik))**2/(sv1r(ik)*st1r(ik))
      enddo

      call makenamebob(file0,'m',istart,'.cos')
      open(unit=20,file=file0,form='formatted')
      do i=1,Nx/2+1
           write(20,'(i5,4e16.8)') (i-1),phi_r(i),abs(sm1r(i))**2,
     &                            abs(sv1r(i)),abs(st1r(i))
      enddo
      close(20)

      end 


