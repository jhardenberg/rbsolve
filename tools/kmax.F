c Maximum K in spectrum as a function of time
#include "config.h"

      program main

#include "param.h"

      complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 v(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 w(nwax+1,2*nwaz+1,0:Nylmem+1)
      real*8 ss(Nx)
      integer*8 n(Nx),imode
      real xmin,xmax
	real*8 sum,sumk,suml,sump(Nx),sumpp,kmode
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
      real*8 K2(nwax+1,2*nwaz+1),lmax,ssmax,dk
	integer ik



#ifndef NOMPI 
      print *,"Please compile with the NOMPI option!"
      stop
#endif      

      iam=1

      if(iargc().lt.3) then
         print *, 'Usage: kmax cc y istart [iend istep]'
	 print *,'Computes the horizontal spectrum of variable yy'
	 print *, 'y<=0 means all levels'
         stop
      endif
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

      call read_para(0)
      call Crea_wave(kxm,kzm,K2,mm,nn,kxmin,kzmin)
c  use x scale as granularity   
        dk=kxm(2)



      open(unit=20,file='kmax.dat',form='formatted')
      do it=istart,iend,istep

      do i=1,Nx/2+1
         n(i)=0.	
         ss(i)=0	.
      enddo
        if(cc.eq.'e') then

             call readfield(u,'u',it)
             call readfield(v,'v',it)
             call readfield(w,'w',it)
	if(jj.le.0) then
c	Sum over shells

         do j=-jj,Ny+jj+1
           do i=1,nwax+1
              do k=1,2*nwaz+1
   	 	 ik=sqrt(K2(i,k))/dk+1.5      
	         ss(ik)=ss(ik)+dreal(dconjg(u(i,k,j))*u(i,k,j))+
     &	 dreal(dconjg(v(i,k,j))*v(i,k,j))+
     &	 dreal(dconjg(w(i,k,j))*w(i,k,j))
	         n(ik)=n(ik)+1
	      enddo
           enddo
	  enddo
	else
	    do i=1,nwax+1
              do k=1,2*nwaz+1
                 ik=sqrt(K2(i,k))/dk+1.5
	         ss(ik)=ss(ik)+dreal(dconjg(u(i,k,jj))*u(i,k,jj))+
     &	 dreal(dconjg(v(i,k,jj))*v(i,k,jj))+
     &	 dreal(dconjg(w(i,k,j))*w(i,k,j))
                 n(ik)=n(ik)+1
              enddo
           enddo
	endif
	else
         call readfield(u,cc,it)
        if(jj.le.0) then
c       Sum over shells
         do j=-jj,Ny+jj+1
           do i=1,nwax+1
              do k=1,2*nwaz+1
                 ik=sqrt(K2(i,k))/dk+1.5      
                 ss(ik)=ss(ik)+dreal(dconjg(u(i,k,j))*u(i,k,j))
                 n(ik)=n(ik)+1
              enddo
           enddo
          enddo
        else
            do i=1,nwax+1
              do k=1,2*nwaz+1
                 ik=sqrt(K2(i,k))/dk+1.5
                 ss(ik)=ss(ik)+dreal(dconjg(u(i,k,jj))*u(i,k,jj))
                 n(ik)=n(ik)+1
              enddo
           enddo
        endif

	endif

c  we want int(E(k)dk) = var
          do i=1,Nx/2+1
	     ss(i)=ss(i)/n(i)*(i-1)
	  enddo

	  sum=0
	  sumk=0
	  suml=0
	  lmax=0.
	  ssmax=0.

          imode=0
          do i=2,kx/2+1
            print *,i-1,ss(i)
            if(ss(i).gt.ssmax) then
                        ssmax=ss(i)
                        lmax=Lx/(i-1)
            endif
            suml=suml+ss(i)*Lx/(i-1)
            sumk=sumk+ss(i)*(i-1)
            sum=sum+ss(i)
          enddo
          sumpp=0.
          do i=2,kx/2+1
            sumpp=sumpp+ss(i)
	    sump(i)=sumpp
            if(sumpp.lt.sum*0.5) imode=i
          enddo
      kmode=(imode-1)+(sum*0.5-sump(imode))/(sump(imode+1)-sump(imode))

        write(20,'(i6,4e16.7)') it,lmax,suml/sum,Lx/sumk*sum,Lx/kmode

      enddo
      close(20)

      end



