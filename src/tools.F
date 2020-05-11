
#include "config.h"


c **********************************************************************

      subroutine salvat (ut,vt,wt,inum)
     
#include "param.h"
 
      complex*16 ut(nwax+1,2*nwaz+1,0:Nyl+1),vt(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 wt(nwax+1,2*nwaz+1,0:Nyl+1)
      real*8 um,wm,vm 
      integer*4 inum 
      
       do j=0,Nyl+1
       um=abs(ut(5,1,j)  ) 
       vm=abs(vt(5,1,j)  ) 
       wm=abs(wt(5,1,j)  ) 
           write(inum,*)  um,vm,wm       
       enddo

      return
      end

      subroutine salvatotv (ut,inum)
#include "param.h"
      real*8 ut(Nx,Nz,0:Nyl+1)
      real*4 a1(Nx,Nyl)
      integer inum
      k=Nz/2
      do j=1,Nyl
             do i=1,Nx
                a1(i,j)=ut(i,k,j)
             enddo
      enddo
      write(inum) a1
      return
      end

      subroutine salvatotr (ut,inum)
#include "param.h"
      real*8 ut(Nx,Nz,0:Nyl+1)
	real*4 a1(Nx,Nz),a2(Nx,Nz),a3(Nx,Nz)
      integer inum
      do k=1,Nz
          do i=1,Nx
            a1(i,k)= ut(i,k,10) 
            a2(i,k)= ut(i,k,65) 
            a3(i,k)= ut(i,k,120) 
          enddo
      enddo
      write(inum)  a1
      write(inum)  a2
      write(inum)  a3
      return
      end

      subroutine salvatot (ut,vt,wt,inum)
      
#include "param.h"
 
      complex*16 ut(nwax+1,2*nwaz+1,0:Nyl+1),vt(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 wt(nwax+1,2*nwaz+1,0:Nyl+1)
      integer*4 inum,jj
      
       do j=1,Nylmax
          jj=j+Nybase
          do k=1,nwaz*2+1
             do i=1,nwax+1
           write(inum,'(3(i3,x),3(e13.3,x))') 
     +     i,jj,k,dreal(ut(i,k,j)),dreal(vt(i,k,j)) ,dreal(wt(i,k,j) )     
       enddo
       enddo
       enddo
      return
      end

 
      subroutine slice()

*     Computes the vertical size Nyl of the domain associated
*     with this thread
*     Reads in cpu weights from file cpuweights

*     Computes also:
*     Nybase = the number of levels below this domain
*     Nylmin, Nylmax = inner levels of this domain (no ghost levels)
*     Nylmin0,Nylmax0 = the same but with ghost levels, except 
*                       for the lower and upper domain
*                       where special cases for the boundary 
*                       conditions are made
*     Nylmin1,Nylmax1 = inner levels except for 
*			lower where Nylmin1=1
*			upper one where Nylmax1=Nyl-1
*     Nyls(i)     = a list of all Nyl associated with node i

#include "param.h"
*	cauchy,lagr,navie,euler=1000,650,500,400

	integer il(MAXCPU),nsum2
        integer ierr
        real*4 ww(MAXCPU),wsum,rest

#ifndef NOMPI      
#ifdef AUTO_WEIGHTS
      open(unit=11,file='cpuweights')
      wsum=0
      do i=1,NPROC
            read(11,*) ww(i)
            wsum=wsum+ww(i)
      enddo
      close(11)
#else
      do i=1,NPROC
            ww(i)=1
            wsum=wsum+ww(i)
      enddo
#endif
      nsum2=0
      rest=0
      do i=1,NPROC
         il(i)=Ny*ww(i)/wsum
	 rest=rest+Ny*ww(i)/wsum-il(i)
        if(rest.gt.1) then
  	  il(i)=il(i)+1
	  rest=rest-1
        endif
        nsum2=nsum2+il(i)
      enddo
      il(1)=il(1)+Ny-nsum2
#else
      il(1)=Ny
#endif      
      do i=1,NPROC
         Nyls(i)=il(i)
      enddo
      Nyl=il(iam)

      Nybase=0
      do i=1,iam-1
	   Nybase=Nybase+Nyls(i)
      enddo

      if(Nyl.gt.Nylmem) then
	print *,'cpu',iam,' Nyl=',Nyl,' > Nylmem=',Nylmem
cmpi	call MPI_FINALIZE(ierr)
	stop
      endif

      Nylmax=Nyl
      Nylmin=1
      Nylmin0=0
      Nylmax0=Nyl+1
      Nylmin1=0
      Nylmax1=Nyl

#ifndef FREE_SLIP_BOTTOM
      if(iam.eq.1) then 
         Nylmin=2
         Nylmin0=1
	 Nylmin1=1
      endif
#endif
#ifndef FREE_SLIP_TOP
      if(iam.eq.NPROC) then
         Nylmax=Nyl-1
         Nylmax0=Nyl
	 Nylmax1=Nyl-1
      endif
#endif      

#ifdef TEMPERATURE
#ifndef TEMPERATURE_BOTTOM
      if(iam.eq.1) then 
         Nylmin=1
         Nylmin0=0
      endif
#endif
#ifndef TEMPERATURE_TOP
      if(iam.eq.NPROC) then
         Nylmax=Nyl
         Nylmax0=Nyl+1
      endif
#endif      
#endif


#ifdef SALINITY
#ifndef SCALAR_BOTTOM
      if(iam.eq.1) then 
         Nylmin=1
         Nylmin0=0
      endif
#endif
#ifndef SCALAR_TOP
      if(iam.eq.NPROC) then
         Nylmax=Nyl
         Nylmax0=Nyl+1
      endif
#endif      
#endif


#ifndef NOMPI
      print *,'Limits(',iam,'):  ',100*real(Nyl)/Ny,
     +          '%  size Nyl=',Nyl,' base=',Nybase,
     + ' Nylmin0=',Nylmin0,' Nylmin =',Nylmin,
     + ' Nylmax =',Nylmax, ' Nylmax0=',Nylmax0
#endif

      end

#ifndef NOMPI

*     igather is defined here and rgather in mpitools.F
*     to avoid a silly warning given by g77

      subroutine igather(value,array)
#include "param.h"
cmpi        include 'mpif.h'
      integer*4 value,array(MAXCPU)
      integer ierr
      array(iam)=value
      if(NPROC.ne.1) then
#ifndef NOMPI
        call MPI_ALLGATHER(value,1,MPI_INTEGER4,
     &              array,1,MPI_INTEGER4,MPI_COMM_WORLD,ierr)
#endif
      endif
      end

*     ibroadcast is defined here and rbroadcast in mpitools.F


      subroutine ibroadcast(value,iroot)
* Broadcasts a real value from root
#include "param.h"
cmpi        include 'mpif.h'
      integer*4 value
      integer ierr,iroot
      if(NPROC.ne.1) then
#ifndef NOMPI
      call MPI_BCAST (value,1,MPI_INTEGER4,iroot-1,MPI_COMM_WORLD,ierr) 
#endif
      endif
      end
#endif


      subroutine arraycopy (a,b)
*     copies a to b

#include "param.h"
     
      complex*16 a(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 b(nwax+1,2*nwaz+1,0:Nyl+1)

      do j=Nylmin0,Nylmax0
        do k=1,2*nwaz+1
 	 do i=1,nwax+1
            b(i,k,j)=a(i,k,j)
         enddo      
        enddo      
      enddo      

      end

#ifdef NAMELIST_PARAMS
      subroutine read_para(idisp)

c Read in the parameters
#include "param.h"
#include "allpara.h"

      real*8 kcx,kcz
      integer*4 idisp

      real*8 pi2,pi
      parameter (pi2=6.28318530717959d0)

      pi=pi2/2

      ptan=1.7
      qtan=1.5
      Re=300.d0
      DD=1.d0
      Lx=pi2
      Lz=pi2
      dt=1e-6
      ttot=50000
      nsave=500
      Ra=1e7
      Pr=0.71d0
      Rt=1.025e7
      Rs=1e7
      Le=3.d0
      Omegay=0.d0
      Omegaz=0.d0
      lapse=100.d0
      qrad=1.d0

      alfa(1)= 0.0d0
      alfa(2)=0.066666666666d0
      alfa(3)=0.166666666666d0
      beta(1)=0.533333333333d0
      beta(2)=0.066666666666d0
      beta(3)=0.166666666666d0
      gamma(1)=0.533333333333d0
      gamma(2)=0.416666666666d0
      gamma(3)=0.75d0
      zeta(1)=0.d0
      zeta(2)=-0.28333333333d0
      zeta(3)=-0.41666666666d0

      open(11,file='param_namelist')
      read (11, numeric)
      read (11, phys)
      close (11)

      invRe=1.d0/Re

      if(idisp.eq.1) then
        if(iam.eq.1) then
          write(*,*) 'Nx=',Nx
          write(*,*) 'Ny=',Ny
          write(*,*) 'Nz=',Nz
          write (*,*) 'ptan=',ptan
          write (*,*) 'qtan=',qtan
          write(*,*) 'DD=',DD
          write(*,*) 'Lx=',Lx
          write(*,*) 'Lz=',Lz
          write(*,*) 'dt=',dt
          write(*,*) 'ttot=',ttot
          write(*,*) 'nsave=',nsave
          write(*,*) 'Reynolds =', Re
          write(*,*) 'Rayleigh =',Ra
          write(*,*) 'Prandtl  =',Pr
          write(*,*) 'Rt',Rt
          write(*,*) 'Rs',Rs
          write(*,*) 'Le',Le
          write(*,*) 'Omegay',Omegay
          write(*,*) 'Omegaz',Omegaz
          write(*,*) 'lapse',lapse
        endif
#ifdef SALINITY
        write(*,*) 'SALT FINGERS'
#endif
      endif

      dx=Lx/dfloat(Nx)
      dz=Lz/dfloat(Nz)

      kcx=pi/dx
      kcz=pi/dz
      kxmin=2.d0*kcx/dfloat(Nx)
      kzmin=2.d0*kcz/dfloat(Nz)

      if(idisp.eq.1) then
c Write stuff to file
        open(12,file='tanh.par')
        write(12,*) 'ptan=',ptan,' ','qtan=',qtan
        write(12,*) 'Ny=',Ny,' ','Nx=',Nx,' ','Nz=',Nz
        write(12,*) 'Re=',Re,' ','Lx=',Lx,' ','Lz=',Lz
        write(12,*) 'dt=',dt,' ','ttot=',ttot,' ','nsave=',nsave
        write(12,*) 'dx=',dx,' ','dz=',dz
        close(12)
      endif

      end
#else
      subroutine read_para(idisp)

c Read in the parameters
#include "param.h"

#include "allpara.h"
      real*8 kcx,kcz
      integer*4 idisp

      real*8 pi2,pi
      parameter (pi2=6.28318530717959d0)

*      pi=asin(1.d0)*2.d0
*      pi2=pi*2.d0
      pi=pi2/2


      open(11,file='param0')
      read (11,*) ptan
      read (11,*) qtan
      read (11,*) Re
      read (11,*) DD
      read(11,*)  Lx
      read (11,*) Lz
      read (11,*) dt
      read (11,*) ttot
      read (11,*) nsave
      read (11,*)  Ra
      read (11,*)  Pr
      read (11,*) Rt
      read (11,*) Rs
      read (11,*) Le
      read (11,*) Omegay
      read (11,*) Omegaz
      close (11)

      invRe=1.d0/Re

	if(idisp.eq.1) then
      if(iam.eq.1) then
      write (*,*) 'ptan=',ptan
      write (*,*) 'qtan=',qtan
      write(*,*)'Reynolds =', Re
      write(*,*)'Rayleigh =',Ra
      write(*,*)'Prandtl  =',Pr
      write(*,*)'DD=',DD
      write(*,*)'Lx=',Lx
      write(*,*)'Lz=',Lz
      write(*,*)'dt=',dt
      write(*,*)'ttot=',ttot
      write(*,*)'nsave=',nsave
      write(*,*) 'Ny=',Ny
      write(*,*) 'Rt',Rt
      write(*,*) 'Rs',Rs
      write(*,*) 'Le',Le
      write(*,*) 'Omegay',Omegay
      write(*,*) 'Omegaz',Omegaz
      endif
#ifdef SALINITY
	write(*,*) 'SALT FINGERS'
#endif
	endif


c Read in coefficients for RK3

      open(12,file='param1')
      do k=1,3
       read(12,*) alfa(k)
       read(12,*) beta(k)
       read(12,*) gamma(k)
       read(12,*) zeta(k)
      end do
      close(12)

      dx=Lx/dfloat(Nx)
      dz=Lz/dfloat(Nz)

      kcx=pi/dx
      kcz=pi/dz
      kxmin=2.d0*kcx/dfloat(Nx)
      kzmin=2.d0*kcz/dfloat(Nz)

      if(idisp.eq.1) then
c Write stuff to file
        open(12,file='tanh.par')
        write(12,*) 'ptan=',ptan,' ','qtan=',qtan
        write(12,*) 'Ny=',Ny,' ','Nx=',Nx,' ','Nz=',Nz
        write(12,*) 'Re=',Re,' ','Lx=',Lx,' ','Lz=',Lz
        write(12,*) 'dt=',dt,' ','ttot=',ttot,' ','nsave=',nsave
        write(12,*) 'dx=',dx,' ','dz=',dz
        close(12)
      endif

      end
#endif


      subroutine crea_wave(kxm,kzm,K2,mm,nn,kxmin,kzmin)
#include "param.h"
      real*8 kxm(nwax+1),kzm(2*nwaz+1)
      complex*16 mm(nwax+1),nn(2*nwaz+1)
      real*8 K2(nwax+1,2*nwaz+1),kxmin,kzmin

c Computation of wave vector arrays

      kxm(1)=0.d0
      kzm(1)=0.d0
      do i=2,(nwax+1)
          kxm(i)=kxmin*dfloat(i-1)
      end do

      do k=2,(nwaz+1)
          kzm(k)=kzmin*dfloat(k-1)
      end do

      do k=nwaz+2,2*nwaz+1
        kzm(k)=-kzm(2*nwaz+3-k)
      end do

      do i=1,nwax+1
       do k=1,2*nwaz+1
        K2(i,k)=(kxm(i))**2+(kzm(k))**2
       end do
      end do

      do i=1,nwax+1
       mm(i)=(0.d0,1.d0)*kxm(i)
      end do

      do k=1,2*nwaz+1
       nn(k)=(0.d0,1.d0)*kzm(k)
      end do

      end
