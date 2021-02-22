C Horizontal Energy spectrum
#include "config.h"

      program main

#include "param.h"

      complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 v(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 w(nwax+1,2*nwaz+1,0:Nylmem+1)
      real*8 ss(Nx),ssx(Nx),ssy(Nx), ss3d(Nx/2+1,Ny)
      real*8 umax,umin,dk
      real*8 vmax,vmin
      real xmin,xmax
      integer istart,iend,istep,it,jmax,jmin,jj,iargc
      character*12 file0
      character*45 file

      character cc*1,str*30

      real*8 ptan,qtan
      real*8 DD,Lx,Lz,Re,invRe,dt,dx,dz
      real*8 kxmin,kzmin,Ra,Kscalar,Pr,sum

      integer*4 nsave,ttot,nn(Nx),nnx(Nx),nny(Nx),nmin,nmax
      real*8 alfa(3),beta(3),gamma(3),zeta(3)

      common/ALLPARA/ptan,qtan,DD,Re,Lx,Lz,invRe,dt,dx,dz,
     +   kxmin,kzmin,nsave,ttot,alfa,beta,gamma,zeta,Ra,Kscalar,Pr

      real*8 kxm(nwax+1),kzm(2*nwaz+1),pwr
      complex*16 mm(nwax+1),nnn(2*nwaz+1)
      real*8 K2(nwax+1,2*nwaz+1),dum,rt,rs,rr
      integer ik
      logical f3d

#ifndef NOMPI 
      print *,"Please compile with the NOMPI option!"
      stop
#endif      

      iam=1
      f3d=.false.

      if(iargc().lt.3) then
         print *, 'Usage: spectrumh cc y istart [iend istep]'
	 print *,'Computes the horizontal spectrum of variable yy'
	 print *, 'y<0 means: compute from level -y to Ny+1+y'
         print *, '(y=-1 means all levels)'
	 print *, 'y=0 means compute a 2D spectrum on each level'
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
      call Crea_wave(kxm,kzm,K2,mm,nnn,kxmin,kzmin)

c  use x scale as granularity	
	dk=kxm(2)
	print *,dk

      if((cc.eq.'b').or.(cc.eq.'f')) then
	print *,'Reading param0'
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

	sum=0.

      do i=1,Nx
	ss(i)=0
        nn(i)=0
      enddo
      do it=istart,iend,istep

       if(cc.eq.'e') then
             call readfield(u,'u',it)
             call readfield(w,'w',it)
             call readfield(v,'v',it)
         
	if(jj.lt.0) then
c	Sum over shells
	 do j=-jj,Ny+1+jj
           do i=1,nwax+1
              do k=1,2*nwaz+1
   	 	 ik=sqrt(K2(i,k))/dk+1.5
                 pwr=dreal(dconjg(u(i,k,j))*u(i,k,j))+
     &	 dreal(dconjg(w(i,k,j))*w(i,k,j))+
     &	 dreal(dconjg(v(i,k,j))*v(i,k,j))
                 ss(ik)=ss(ik)+pwr
                 nn(ik)=nn(ik)+1
	      enddo
           enddo
	sum=sum+1
	  enddo
        elseif(jj.eq.0) then
          f3d=.true.
          do j=1,Ny
           ss(:)=0.0
           nn(:)=0
           do i=1,nwax+1
              do k=1,2*nwaz+1
                 ik=sqrt(K2(i,k))/dk+1.5
                 pwr=dreal(dconjg(u(i,k,j))*u(i,k,j))+
     &  dreal(dconjg(w(i,k,j))*w(i,k,j))+
     &  dreal(dconjg(v(i,k,j))*v(i,k,j))
                 ss(ik)=ss(ik)+pwr
                 nn(ik)=nn(ik)+1
              enddo
           enddo
           do i=1,Nx/2+1
              ss3d(i,j)=ss(i)/nn(i)*(i-1)
           enddo
         enddo
	else
	    do i=1,nwax+1
              do k=1,2*nwaz+1
                 ik=sqrt(K2(i,k))/dk+1.5
	         pwr=dreal(dconjg(u(i,k,jj))*u(i,k,jj))+
     &	 dreal(dconjg(w(i,k,jj))*w(i,k,jj))+
     &	 dreal(dconjg(v(i,k,jj))*v(i,k,jj))
                 ss(ik)=ss(ik)+pwr
                 nn(ik)=nn(ik)+1
              enddo
           enddo
	   sum=sum+1
	endif
       else
c if cc ne e
        if(cc.eq.'b') then
        call readfield(u,'t',it)
        call readfield(v,'s',it)
             do j=1,nwax+1
              do k=1,2*nwaz+1
                 do i=1,Nx
                   u(i,k,j)=u(i,k,j)*rr-v(i,k,j)
                 enddo
              enddo
             enddo
        elseif(cc.eq.'f') then
        call readfield(u,'t',it)
        call readfield(v,'s',it)
        call readfield(w,'v',it)
             do j=1,nwax+1
              do k=1,2*nwaz+1
                 do i=1,Nx
                   u(i,k,j)=w(i,k,j)*(u(i,k,j)*rr-v(i,k,j))
                 enddo
              enddo
             enddo
	else
           call readfield(u,cc,it)
        endif

        if(jj.le.0) then
c       Sum over shells
         do j=-jj,Ny+1+jj
           do i=1,nwax+1
              do k=1,2*nwaz+1
                 ik=sqrt(K2(i,k))/dk+1.5
                 pwr=dreal(dconjg(u(i,k,j))*u(i,k,j))
                 ss(ik)=ss(ik)+pwr
                 nn(ik)=nn(ik)+1
              enddo
           enddo
        sum=sum+1
          enddo
        else
            do i=1,nwax+1
              do k=1,2*nwaz+1
                 ik=sqrt(K2(i,k))/dk+1.5
                 pwr=dreal(dconjg(u(i,k,jj))*u(i,k,jj))
                 ss(ik)=ss(ik)+pwr
                 nn(ik)=nn(ik)+1
              enddo
           enddo
        sum=sum+1
        endif
	endif

      enddo
      
      if (f3d) then
        call makenamebob(file0,cc,istart,'.s3d')
        open(unit=20,file=file0,form='formatted')
        do j=1,Ny
         do i=1,Nx/2+1
           write(20,'(i5,x,i5,x,2e16.8)') (i-1),j,ss3d(i,j)
         enddo
        enddo
        close(20)

      else
        call makenamebob(file0,cc,istart,'.spe')
        open(unit=20,file=file0,form='formatted')
        do i=2,Nx/2+1
  	  sum=sum+ss(i)
        enddo
        sum=sum*Nx*Nx
        do i=1,Nx/2+1
           write(20,'(i5,2e16.8)') (i-1),ss(i)/nn(i)*(i-1),sum
        enddo
        close(20)
      endif

      end



