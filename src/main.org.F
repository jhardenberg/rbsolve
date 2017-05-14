
c  As soon as we started programming, we found to our surprise that it
c  was not as easy to get programs right as we had thought.  Debugging had
c  to be discovered.  I can remember the exact instant when I realized
c  that a large part of my life from then on was going to be spent in
c  finding mistakes in my own programs.
c                 -- Maurice Wilkes discovers debugging, 1949


c   A debugged program is one for which you have
c   not yet found the conditions that make it fail.
c                -- Jerry Ogdin


#include "config.h"

C Integrates the Navier-Stokes equations using third order 
c Runge-Kutta (Le-Moin, 1981).
c Horizontal derivatives are spectral, vertical derivatives use finite differences.

C (c) 2003-2007 Rayleigh-Benard, Double diffusion, 
c MPI version by J. von Hardenberg (ISAC-CNR) and G. Passoni (Poli MI)

C (c) 2002 Navier-Stokes version by G. Passoni

c  if TEMPERATURE is defined, the code is convective
c  check config.h for the coefficients in the equations

c  if SALINITY is defined, an active scalar (salinity) is added
c  check config.h for the coefficients in the equations

c  Using RK the old fields are assumed not to be needed outside
c  the inner time loop. If you change the time-stepping scheme
c  you must consider this (example: changing
c  slice thicknesses will require also Lu_2 etc to be exchanged )



* This is a  MPI version
* define NOMPI in config.h if you want a scalar code 


c  No need to touch these
#define ADD1t
#define ADD1s
#ifdef TEMPERATURE
#undef ADD1t
#define ADD1t +1
#endif
#define ADD1s
#ifdef SALINITY
#undef ADD1s
#define ADD1s +1
#endif


      Program main  
#include "param.h"
cmpi      include 'mpif.h'
     
      integer ierr,npr

      iam=0
      npr=1

      write (*,*) NPROC
#ifndef NOMPI      
c     Initialize MPI stuff
  	print *,'This is the MPI version'

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, iam, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, npr, ierr )

*     my big MPI comm datatype
*     5(or 3) blocks of (nwax+1)*(2*nwaz+1) pointing at u,v,w,(t,s)
*     spaced 1*(nwax+1)*(2*nwaz+1)*(Nylmem+2)
*     myMPI_ALL_TYPE is ALL fields
*     myMPI_VEL_TYPE is velocity fields only

      call MPI_TYPE_VECTOR(3 ADD1t ADD1s ,(nwax+1)*(2*nwaz+1),
     + 1*(nwax+1)*(2*nwaz+1)*(Nylmem+2),
     + MPI_COMPLEX16,myMPI_ALL_TYPE,ierr)

      call MPI_TYPE_VECTOR(3 ,(nwax+1)*(2*nwaz+1),
     + 1*(nwax+1)*(2*nwaz+1)*(Nylmem+2),
     + MPI_COMPLEX16,myMPI_VEL_TYPE,ierr)

      CALL MPI_TYPE_COMMIT(myMPI_VEL_TYPE, ierr)
      CALL MPI_TYPE_COMMIT(myMPI_ALL_TYPE, ierr)

#endif

      if(npr.ne.NPROC) then
         print *,'Compiled for ',NPROC,' threads, using ',npr,' !'
!cmpi         call MPI_FINALIZE(ierr)
cmpi       call MPI_ABORT( MPI_COMM_WORLD ,6,ierr)
         stop         
      endif

      iam=iam+1   


c Initialize FFT
      if(iam.eq.1) print *,'Initializing FFT routines ...' 
      call initfour(1)

*     Compute the appropriate slice size Nyl

#ifdef AUTO_WEIGHTS
      call cpuweights()
#endif

      call slice()


#ifdef FREE_SLIP_TOP
      if(iam.eq.1) print *,'Top boundary   : free slip'
#else
      if(iam.eq.1) print *,'Top boundary   :   no slip'
#endif
#ifdef FREE_SLIP_BOTTOM
      if(iam.eq.1) print *,'Bottom boundary: free slip'
#else
      if(iam.eq.1) print *,'Bottom boundary:   no slip'
#endif

  
*     Call the main program
      call mainsub()

cmpi      call MPI_FINALIZE(ierr)

      end

      subroutine mainsub
#include "param.h"
      
      integer*4 l,conta
      integer*4 t,itime

#include "allpara.h"
	integer idbg

      real*8 kxm(nwax+1),kzm(2*nwaz+1) 
      complex*16 mm(nwax+1),nn(2*nwaz+1)
      real*8 K2(nwax+1,2*nwaz+1) 

      real*8 yp(0:Ny+1),yv(0:Ny+1) 
      real*8 qp(Ny),rp(Ny),sp(Ny),qv(Ny),rv(Ny),sv(Ny) 
      real*8 tauw

*     Physical  space fields
      real*8 uo(Nx,Nz,0:Nylmem+1)
      real*8 vo(Nx,Nz,0:Nylmem+1)
      real*8 wo(Nx,Nz,0:Nylmem+1)

      complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1,3 ADD1t ADD1s )

*     u(1,1,1,1) is u
*     u(1,1,1,2) is w
*     u(1,1,1,3) is v
*     u(1,1,1,4) is T (if it applies)
*     This structure is adopted to optimize MPI communications

      integer iu1,iut,iutemp,ilu1,ilu2,ierr,iv1,ivt,iw1,iwt
      parameter(iu1=1,iw1=2,iv1=3)
 
#define u_1 u(1,1,0,iu1)
#define v_1 u(1,1,0,iv1)
#define w_1 u(1,1,0,iw1)

#define Lu_1 Lu(1,1,0,ilu1)
#define Lv_1 Lv(1,1,0,ilu1)
#define Lw_1 Lw(1,1,0,ilu1)
#define Nu_1 Nu(1,1,0,ilu1)
#define Nv_1 Nv(1,1,0,ilu1)
#define Nw_1 Nw(1,1,0,ilu1)
#define Lu_2 Lu(1,1,0,ilu2)
#define Lv_2 Lv(1,1,0,ilu2)
#define Lw_2 Lw(1,1,0,ilu2)
#define Nu_2 Nu(1,1,0,ilu2)
#define Nv_2 Nv(1,1,0,ilu2)
#define Nw_2 Nw(1,1,0,ilu2)

      complex*16 Lu(nwax+1,2*nwaz+1,0:Nylmem+1,2)
      complex*16 Lv(nwax+1,2*nwaz+1,0:Nylmem+1,2)
      complex*16 Lw(nwax+1,2*nwaz+1,0:Nylmem+1,2)
      complex*16 Nu(nwax+1,2*nwaz+1,0:Nylmem+1,2)
      complex*16 Nv(nwax+1,2*nwaz+1,0:Nylmem+1,2)
      complex*16 Nw(nwax+1,2*nwaz+1,0:Nylmem+1,2)


#ifdef TEMPERATURE
      integer it1
      parameter (it1=4)

      complex*16 t_2(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 Lt(nwax+1,2*nwaz+1,0:Nylmem+1,2)
      complex*16 Nt(nwax+1,2*nwaz+1,0:Nylmem+1,2)
      real*8 to(Nx,Nz,0:Nylmem+1)

#define t_1 u(1,1,0,it1)
#define Lt_1 Lt(1,1,0,ilu1)
#define Lt_2 Lt(1,1,0,ilu2)
#define Nt_1 Nt(1,1,0,ilu1)
#define Nt_2 Nt(1,1,0,ilu2)

#endif

#ifdef SALINITY
      integer is1
#ifdef TEMPERATURE
      parameter (is1=5)
#else
      parameter (is1=4)
#endif

      complex*16 s_2(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 Ls(nwax+1,2*nwaz+1,0:Nylmem+1,2)
      complex*16 Ns(nwax+1,2*nwaz+1,0:Nylmem+1,2)
      real*8 so(Nx,Nz,0:Nylmem+1)

#define s_1 u(1,1,0,is1)
#define Ls_1 Ls(1,1,0,ilu1)
#define Ls_2 Ls(1,1,0,ilu2)
#define Ns_1 Ns(1,1,0,ilu1)
#define Ns_2 Ns(1,1,0,ilu2)

#endif

      complex*16 p(nwax+1,2*nwaz+1,0:Nylmem+1)

*     These arrays contain the principal diagonals of
*     the Helmholtz problem
      real*8 D(nwax+1,2*nwaz+1,Nylmem)
      real*8 E(nwax+1,2*nwaz+1,Nylmem)
      real*8 F(nwax+1,2*nwaz+1,Nylmem)

#define NX ((nwax+1)*(2*nwaz+1))

      real*8 xlh(NX,Nylmem)
      real*8 xuh(NX,Nylmem)

#define NNR (2*NPROC-2)
      real*8 reda(NX,0:NNR+1)
      real*8 redb(NX,0:NNR+1)
      real*8 redc(NX,0:NNR+1)
      complex*16 redr(NX,0:NNR+1)

      real*4 tdat0(2),tdat(2),gdat1(2),gdat0(2),udat0(2),udat(2)
      real*4 ttot1,ttot0,timu,timu0
      real*4 tglob1,tglob0
      real*8 usertot

      real etime

      integer idum
      real*8 ran2,rtime

	

      ilu1=1
      ilu2=2

c Read in initial condition

      itime=-1
      rtime=0.d0

      if(iam.eq.SAVENODE)  then
	print *,'reading nrec'
        open(unit=10,file='nrec.d')
        read(10,*,end=111) itime,rtime
        close(10)
 111    continue
#ifndef NOMPI
        call ibroadcast(itime,SAVENODE)
        call rbroadcast(rtime,SAVENODE)
      else
        call ibroadcast(itime,SAVENODE)    
        call rbroadcast(rtime,SAVENODE)
#endif
      endif

     
      if(itime.eq.-1) then
c-#define READINIFILE 
#ifdef READINIFILE   
      open(10,file='ini.dat',form='unformatted',status='unknown')
      do j=1,Nybase
      read(10) ((u(i,k,1,iu1),i=1,nwax+1),k=1,2*nwaz+1) 
      read(10) ((u(i,k,1,iv1),i=1,nwax+1),k=1,2*nwaz+1) 
      read(10) ((u(i,k,1,iw1),i=1,nwax+1),k=1,2*nwaz+1) 
      enddo
      do j=0,Nyl+1         
      read(10) ((u(i,k,j,iu1),i=1,nwax+1),k=1,2*nwaz+1) 
      read(10) ((u(i,k,j,iv1),i=1,nwax+1),k=1,2*nwaz+1) 
      read(10) ((u(i,k,j,iw1),i=1,nwax+1),k=1,2*nwaz+1) 
      enddo
      close(10) 
#else
*     Reads formatted temperature from fort.9
      do i=1,nwax+1
       do k=1,2*nwaz+1
        do j=0,Ny+1
          u(i,k,j,it1)=(0.0,0.0)
          u(i,k,j,iu1)=(0.0,0.0)
          u(i,k,j,iv1)=(0.0,0.0)
          u(i,k,j,iw1)=(0.0,0.0)
        enddo
       enddo
      enddo
*
c      open(12,file='u.50',form='unformatted',status='unknown')
c      open(13,file='v.50',form='unformatted',status='unknown')
c      open(14,file='w.50',form='unformatted',status='unknown')
c      open(15,file='t.50',form='unformatted',status='unknown')

c      read(12) (((u(i,k,j,iu1),j=0,Nyl+1),i=1,nwax+1),k=1,2*nwaz+1)
c      read(13) (((u(i,k,j,iv1),j=0,Nyl+1),i=1,nwax+1),k=1,2*nwaz+1)
c      read(14) (((u(i,k,j,iw1),j=0,Nyl+1),i=1,nwax+1),k=1,2*nwaz+1)
c      read(15) (((u(i,k,j,it1),j=0,Nyl+1),i=1,nwax+1),k=1,2*nwaz+1)

c      close(12)
c      close(13)
c      close(14)
c      close(15)


#endif
      itime=0

      else

c-#define READOLD
#ifdef READOLD
c     reads the old (antonio et al) format, works only on single cpu ....

       call readold(u_1,'u',itime)
       call readold(v_1,'v',itime)
       call readold(w_1,'w',itime)
#ifdef TEMPERATURE
       call readold(t_1,'T',itime)
#endif
#ifdef SALINITY
       call readold(s_1,'s',itime)
#endif

#else
       call readfield(u_1,'u',itime)
       call readfield(v_1,'v',itime)
       call readfield(w_1,'w',itime)
#ifdef TEMPERATURE
       call readfield(t_1,'t',itime)
#endif
#ifdef SALINITY
       call readfield(s_1,'s',itime)
#endif

#endif
      endif

c Controllo le condizioni iniziali di T

c	Diagnostics
#ifdef LINSTABDBG
      write(20,'(i5,4e16.7)') itime,
     + abs(  u(5,1,18,it1)),
     + abs(  u(5,1,18,is1)),
     + abs(  u(5,1,1,iu1)),
     + abs(  u(5,1,18,iv1))
#endif


      call dealias(u_1)
      call dealias(v_1)
      call dealias(w_1)
#ifdef TEMPERATURE
      call dealias(t_1)
#endif
#ifdef SALINITY
      call dealias(s_1)
#endif

c Read in parameters from files

      call read_para(1)
      ttot=ttot+itime
c Computation of wave vector arrays 

      call Crea_wave(kxm,kzm,K2,mm,nn,kxmin,kzmin)
      
c Computation of the y-grid parameters

      call punti (DD,ptan,qtan,yp,yv,qp,rp,sp,qv,rv,sv)

      
c Computation of the tridiagonal matrices for the pressure Helmoltz eq.
c These are not dependent on the wave number

      call Crea_DF (qp,sp,D,F)
      call Crea_E  (qp,rp,sp,K2,E)
  
* Compute static part of tridiagonal problem
* D and F are changed, uo and vo are used as temporary storage

      call tridlhs(D,E,F,reda,redb,redc,xlh,xuh,
     +             uo,vo,Nylmax1,NX,NPROC,iam)
 


c Mean pressure gradient
c      tauw=2.d0*invRe  (for laminar flow only with center velocity)
c  for turbulent flow with friction velocity:
      tauw= PRESSURE_GRADIENT   

***   Enforce vertical boundary conditions   ***

       call boundary(u_1,v_1,w_1)

#ifdef TEMPERATURE
       call boundary_temp(t_1,yv)
#endif
#ifdef SALINITY
       call boundary_scalar(s_1,yv)
#endif

****   Exchange velocities between domains   ***


       call mpighost(u_1)

***   Compute physical space velocities   ***

ctime       ttot0=etime(tdat0)    

#ifdef TEMPERATURE
       call Fourier1 (t_1,to)
#endif

#ifdef SALINITY
       call Fourier1 (s_1,so)
#endif
       
       call Fourier1 (u_1,uo) 
       call Fourier1 (v_1,vo) 
       call Fourier1 (w_1,wo)
	
       call check_CFL(uo,vo,wo,itime,dx,dz,yv,dt,rtime)

#ifdef JOST
       call savefield(u_1,'u',itime,rtime)
       call savefield(v_1,'v',itime,rtime)
       call savefield(w_1,'w',itime,rtime)

#ifdef TEMPERATURE
       call savefield(t_1,'t',itime,rtime)
#endif
#ifdef SALINITY
       call savefield(s_1,'s',itime,rtime)
#endif

#endif

***********   MAIN LOOP  **************************

ctime      tglob0=etime(gdat0)
      usertot=0.0


	idbg=0
      do while(itime.le.ttot) 

c      if(iam.eq.1) print *,' t=',itime

#ifdef AUTO_BALANCE 
ctime      timu0=etime(udat0)    
#endif

c Time-advancement over 3 substeps

      do 20 l=1,3

#ifdef SALINITY
       call diffus (s_1,SDIFF,qv,rv,sv,K2,Ls_1)
#endif

#ifdef TEMPERATURE
       call diffus (t_1,TDIFF,qv,rv,sv,K2,Lt_1)
#endif

***   Compute advective/convective terms   ***

ctime       ttot0=etime(tdat0)
*   Lu_1, Lv_1, Lw_1  are used as temporary workspace in conv

       call conv (uo,vo,wo,
     &   yv,mm,nn,Nu_1,Nv_1,Nw_1,
     &   Lu_1,Lv_1,Lw_1)


ctime       ttot1=etime(tdat)
ctime       print *,'Conv: t:',ttot1-ttot0,' u:',tdat(1)-tdat0(1),
ctime     &              ' s:',tdat(2)-tdat0(2)

***   Compute diffusive terms   ***

ctime       ttot0=etime(tdat0)

c VDIFF is defined as invRe and Pr ifdef TEMPERATURE
       call diffus (u_1,VDIFF,qv,rv,sv,K2,Lu_1)
       call diffus (v_1,VDIFF,qv,rv,sv,K2,Lv_1)
       call diffus (w_1,VDIFF,qv,rv,sv,K2,Lw_1)

       call addrotv (u_1,v_1,Omega,Lu_1,Lv_1)
c       Omega = Omega - dt*50000

***     Compute new temporary star velocities  
***     stored   u_1,v_1,w_1                    
***     starting from here the current velocities 
***     are only in uo,vo,wo 

       call build_star (u_1,dt,
     +    Lu_1,Lu_2,
     +    Nu_1,Nu_2,
     +    alfa,beta,gamma,zeta,l,u_1)
       call build_star (v_1,dt,
     +    Lv_1,Lv_2,
     +    Nv_1,Nv_2,
     +    alfa,beta,gamma,zeta,l,v_1)
       call build_star (w_1,dt,
     +    Lw_1,Lw_2,
     +    Nw_1,Nw_2,
     +    alfa,beta,gamma,zeta,l,w_1)
       call add_star_forcing(u_1,dt,
     +    alfa,beta,gamma,zeta,l,u_1,tauw)



#ifdef TEMPERATURE
*  Add temperature contribution to vertical velocity
       call add_star_temp(v_1,dt,
     +     t_1,t_2,
     +     alfa,beta,gamma,zeta,l,v_1,BUOYT)
c BUOYT is Ra*Pr for RB, 1 for Prandtl
c BUOYT is Pr*Le*Rt   for salt fingers
#endif
#ifdef SALINITY
*  Add salinity contribution to vertical velocity
       call add_star_temp(v_1,dt,
     +     s_1,s_2,
     +     alfa,beta,gamma,zeta,l,v_1,BUOYS)

c BUOYS is Pr*Le*Rs for the SF problem  
#endif

        

ctime       ttot1=etime(tdat)
ctime       print *,'Vstr: t:',ttot1-ttot0,' u:',tdat(1)-tdat0(1),
ctime     &              ' s:',tdat(2)-tdat0(2)

***   Enforce vertical boundary conditions   ***
       call boundary(u_1,v_1,w_1)

***   Propagate star velocities to domain below   ***

      call mpighostdown(u_1)

***   Solve the pressure Helmholtz eqs.  ***

ctime       ttot0=etime(tdat0)
       call Poiss_imp (u_1,v_1,w_1,
     +    D,E,F,yv,kxm,kzm,p,dt,alfa,beta,l,
     +     xlh,xuh,reda,redb,redc,redr )  

c      if(mod(itime,nsave).eq.0) then
c       call savefield(p,'p',itime,rtime)
c      endif

ctime       ttot1=etime(tdat)
ctime       print *,'Pres: t:',ttot1-ttot0,' u:',tdat(1)-tdat0(1),
ctime     &              ' s:',tdat(2)-tdat0(2)

***  Propagate pressure to domain above   ***

       call mpighostup(p)

***   Compute the solution  adding pressure term to the velocities ***

ctime       ttot0=etime(tdat0)


c       call savefield(v_1,'a',idbg)
c       call savefield(u_1,'c',idbg)

       call velocita (u_1,v_1,w_1,
     +    p,mm,nn,yp,yv,
     +    u_1,v_1,w_1,
     +    dt,alfa,beta,l)

c       call savefield(v_1,'b',idbg)
c       call savefield(u_1,'d',idbg)
c       call savefield(p,'p',idbg)

ctime       ttot1=etime(tdat)
ctime       print *,'Vel : t:',ttot1-ttot0,' u:',tdat(1)-tdat0(1),
ctime     &              ' s:',tdat(2)-tdat0(2)

***   Enforce vertical boundary conditions   ***
       call boundary(u_1,v_1,w_1)


#ifdef TEMPERATURE
*     For temperature we do not use the index switching trick
*     used for the diffusive and advective terms
*     because we chose to keep temperature in the same
*     big array with velocities and scalars 
*     (which is more efficient for communication)

       call arraycopy(t_1,t_2)
*      Nu_2 etc are now recycled as temporary storage ...
 
       call conv_scalar (uo,vo,wo,to,
     &   yv,mm,nn,Nt_1,
     &   Nu_2,Nv_2,Nw_2)

       call build_star (t_1,dt,
     +    Lt_1,Lt_2,
     +    Nt_1,Nt_2,
     +    alfa,beta,gamma,zeta,l,t_1)

       call boundary_temp(t_1,yv)


#endif

#ifdef SALINITY
      call arraycopy(s_1,s_2)
*    Nu_2 etc are now recycled as temporary storage ...

       call conv_scalar (uo,vo,wo,so,
     &   yv,mm,nn,Ns_1,
     &   Nu_2,Nv_2,Nw_2)

       call build_star (s_1,dt,
     +    Ls_1,Ls_2,
     +    Ns_1,Ns_2,
     +    alfa,beta,gamma,zeta,l,s_1)
c       call savefield(s_1,'a',idbg)
       call boundary_scalar(s_1,yv)
#endif

          idbg=idbg+1



***** This part used to be at the beginning of the loop
***** we do it here so that we end the loop with valid
***** physical velocities

****   Exchange velocities (and temp) between domains   ***

       call mpighost(u_1)

***   Compute physical space velocities   ***

ctime       ttot0=etime(tdat0)

#ifdef TEMPERATURE
       call Fourier1 (t_1,to)
#endif

#ifdef SALINITY
       call Fourier1 (s_1,so)
#endif
       
       call Fourier1 (u_1,uo) 
       call Fourier1 (v_1,vo) 
       call Fourier1 (w_1,wo)


ctime       ttot1=etime(tdat)
ctime       print *,'Four1: t:',ttot1-ttot0,' u:',tdat(1)-tdat0(1),
ctime     &              ' s:',tdat(2)-tdat0(2)

*****
*****
*****
c     Now Lu_1 becomes Lu_2 and so on ....
      iutemp=ilu2
      ilu2=ilu1
      ilu1=iutemp

20    continue 
c	stop

c End of ONE FULL TIME STEP (3 substeps)

      itime=itime+1
      rtime=rtime+dt

#ifdef LINSTABDBG
c	Lin stab. diagnostics
      write(20,'(i5,4e16.7)') itime,
     + abs(  u(5,1,18,it1)),
     + abs(  u(5,1,18,is1)),
     + abs(  u(5,1,1,iu1)),
     + abs(  u(5,1,18,iv1))
#endif

#ifdef AUTO_BALANCE

ctime       ttot0=etime(tdat0)
ctime      timu=etime(udat)
      usertot=usertot+udat(1)-udat0(1)
      if(mod(itime,NBALANCE).eq.0) then

	call autobalance(usertot,u_1)

*     Vertical slice sizes (Nyl) have changed so we need to
*     recompute these!
         call Crea_DF (qp,sp,D,F)
         call Crea_E  (qp,rp,sp,K2,E)
         call tridlhs(D,E,F,reda,redb,redc,xlh,xuh,
     +          uo,vo,Nylmax1,NX,NPROC,iam)
c     XXXX
         call Fourier1 (u_1,uo) 
         call Fourier1 (v_1,vo) 
         usertot=0.
      endif
ctime       ttot1=etime(tdat)
ctime       print *,'Bal. : t:',ttot1-ttot0,' u:',tdat(1)-tdat0(1),
ctime     &              ' s:',tdat(2)-tdat0(2)

#endif


      if(mod(itime,NCHECKCFL).eq.0) 
     &          call check_CFL(uo,vo,wo,itime,dx,dz,yv,dt,rtime)

c Save the spectral velocity fields


      if(mod(itime,nsave).eq.0) then


       call savefield(u_1,'u',itime,rtime)
       call savefield(v_1,'v',itime,rtime)
       call savefield(w_1,'w',itime,rtime)

#ifdef TEMPERATURE
       call savefield(t_1,'t',itime,rtime)
       
#endif
#ifdef SALINITY
       call savefield(s_1,'s',itime,rtime)
#endif

      end if       

      enddo
********* END OF MAIN LOOP **************************
	
ctime       tglob1=etime(gdat1)
       print *,'TOTAL (',iam,'): t:',tglob1-tglob0,
     &                         ' u:',gdat1(1)-gdat0(1),
     &                         ' s:',gdat1(2)-gdat0(2)

1200  end 


      subroutine check_CFL(uo,vo,wo,t,dx,dz,yv,dt,rtime)

#include "param.h"
#include "mpif.h"
      integer*4 t,jj
      real*8 uo(Nx,Nz,0:Nylmem+1)
      real*8 vo(Nx,Nz,0:Nylmem+1)
      real*8 wo(Nx,Nz,0:Nylmem+1)
      real*8 yv(0:Ny+1) 
      real*8 dt,dx,dz,cflx,cflz,cfly,cflymax,umax,wmax,dy,cdy
      real*8 enx,eny,enz,rtime
      real*8 cfl(3,MAXCPU),ic(3),ene(3,MAXCPU)
      integer*4 icoord(3,MAXCPU)
      integer ierr

*     Find local CFLs
*     this could also be done with an MPI reduce operation
*     but we want also a scalar NOMPI version

      wmax=0
      umax=0
      cflymax=0

      enx=0.
      eny=0.
      enz=0.
      do j=1,Nyl
         jj=j+Nybase
         dy=abs(yv(jj)-yv(jj-1))
         if(abs(yv(jj+1)-yv(jj)).lt.dy) dy=abs(yv(jj+1)-yv(jj))
         cdy=dt/dy
         do k=1,Nz
            do i=1,Nx   
               enx=enx+uo(i,k,j)*uo(i,k,j)
               eny=eny+vo(i,k,j)*vo(i,k,j)
               enz=enz+wo(i,k,j)*wo(i,k,j)
               if(abs(uo(i,k,j)).gt.umax) umax=abs(uo(i,k,j))
               if(abs(wo(i,k,j)).gt.wmax) wmax=abs(wo(i,k,j))
               cfly=abs(vo(i,k,j))*cdy
               if(cfly.gt.cflymax) cflymax=cfly
            enddo
         enddo
      enddo

      cfl(1,iam)=umax*dt/dx
      cfl(2,iam)=cflymax
      cfl(3,iam)=wmax*dt/dz

      ene(1,iam)=enx
      ene(2,iam)=eny
      ene(3,iam)=enz

*     Send them all to root
      call rgatherroot(cfl(1,iam),cfl)
      call rgatherroot(ene(1,iam),ene)
  
      if(iam.eq.1) then
         cflx=0
         cfly=0
         cflz=0
         enx=0
         eny=0
         enz=0
         do i=1,NPROC
            if(cfl(1,i).gt.cflx) cflx=cfl(1,i)
            if(cfl(2,i).gt.cfly) cfly=cfl(2,i)
            if(cfl(3,i).gt.cflz) cflz=cfl(3,i)
            enx=enx+ene(1,i)
            eny=eny+ene(2,i)
            enz=enz+ene(3,i)
         enddo
         enx=enx/(2*Nx*Ny*Nz)
         eny=eny/(2*Nx*Ny*Nz)
         enz=enz/(2*Nx*Ny*Nz)
       
         write(*,'(i7,x,f11.6,x,a4,x,3f9.6,x,a4,x,3e13.5)') 
     &         t,rtime,'CFL:',cflx,cfly,cflz,'ENE:',enx,eny,enz

	if(cfly.gt.0.5) then
	  write(*,*) 'CFLy exceeded 0.5, aborting ...'
cmpi      call MPI_ABORT( MPI_COMM_WORLD ,5,ierr)
!ccmpi      call MPI_FINALIZE(ierr)
	  stop
	endif
         
      endif
      end


#ifdef AUTO_BALANCE

        subroutine autobalance(usertime,uu)
*     Every NBALANCE steps call this subroutine
*     which checks how slow Iam has been in the past steps
*     and, if needed gives (acquires) one level from the fastest (slower)
*     neighbours

#include "param.h"
cmpi        include 'mpif.h'

        complex*16 uu(nwax+1,2*nwaz+1,0:Nylmem+1,6) 
*        complex*16 lu(nwax+1,2*nwaz+1,0:Nyl+1,6) 
*        complex*16 nu(nwax+1,2*nwaz+1,0:Nyl+1,6)
 
        real*8 usertime,utime(MAXCPU),utot,uratioup,uratiodown,thres
        integer ierr,l
        integer*4 isend(0:MAXCPU+1)
        real*8 ud(4)
        data ud/10,8,9,10/
*        usertime=ud(iam)

        if(iam.eq.1) print *,'AUTO BALANCING'

        do i=0,NPROC+1
           isend(i)=0
        enddo

*     Exchange all neighbouring levels
        call mpighost(uu)

*     Distribute speeds
        call rgather(usertime,utime)

        utot=0.0
        do i=1,NPROC
           utot=utot+utime(i)
           isend(i)=0
        enddo
        do i=1,NPROC
           utime(i)=utime(i)/utot
        enddo

        uratioup=0
        uratiodown=0

* This is the threshold if it is worth or not to exchange a level
        thres=1.+1./Nyls(iam)

        if(iam.ne.NPROC) uratioup=utime(iam)/utime(iam+1)
        if(iam.ne.1) uratiodown=utime(iam)/utime(iam-1)

        if(Nyl.gt.1) then
        if(uratioup.gt.uratiodown) then
           if(uratioup.gt.thres) then 
              if(Nyls(iam+1).lt.Nylmem-1) then
                isend(iam)=+1 
              else
	       print *,iam+1,' Nylmem too small!'
	      endif
           endif 
        else
           if(uratiodown.gt.thres) then 
              if(Nyls(iam-1).lt.Nylmem-1) then
  		isend(iam)=-1 
	      else
	       print *,iam-1,' Nylmem too small!'
	      endif
           endif 
        endif
        endif

*  Isend contains info on who wants to send a level where
*  Distribute it

      call igather(isend(iam),isend(1))

      if(isend(iam).eq.1) then
*     we want to send a level up
         Nyl=Nyl-1
      endif
      if(isend(iam).eq.-1) then
*     we want to send a level below
         do l=1,5,2
         do j=1,Nyl
         do k=1,2*nwaz+1
         do i=1,nwax+1
            uu(i,k,j,l)=uu(i,k,j+1,l)
         enddo
         enddo
         enddo
         enddo
         Nyl=Nyl-1
      endif
      if(isend(iam+1).eq.-1) then
*     we should acquire the level from above
         Nyl=Nyl+1
      endif
      if(isend(iam-1).eq.1) then
*     we should get the level from below
         do l=1,5,2
         do j=Nyl,0,-1
         do k=1,2*nwaz+1
         do i=1,nwax+1
            uu(i,k,j+1,l)=uu(i,k,j,l)
         enddo
         enddo
         enddo
         enddo
         Nyl=Nyl+1
      endif
    
      call igather(Nyl,Nyls)
    
      if(iam.eq.1) then

         print *,'% CPU user times      New slice thicknesses:'
         do i=1,NPROC
            print *,i,utime(i),'     ',Nyls(i)
         enddo
      endif
        
      Nybase=0
      do i=1,iam-1
	   Nybase=Nybase+Nyls(i)
      enddo

      Nylmax=Nyl
      Nylmin=1
      Nylmin0=0
      Nylmax0=Nyl+1

#ifndef FREE_SLIP_BOTTOM
      if(iam.eq.1) then 
         Nylmin=2
         Nylmin0=1
      endif
#endif

#ifndef FREE_SLIP_TOP
      if(iam.eq.NPROC) then
         Nylmax=Nyl-1
         Nylmax0=Nyl
      endif
#endif      
 
        end
#endif

#ifdef AUTO_WEIGHTS

      subroutine cpuweights()
*     Estimates cpu weights
*     writing them to the file cpuweights

#include "param.h"
cmpi      include 'mpif.h'


*	cauchy,lagr,navie,euler=1000,650,500,400

      real*8 uo(Nx,Nz,0:Nylmem+1),al(3)
      complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 Lu2(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 Lu1(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 uu(nwax+1,2*nwaz+1,0:Nylmem+1)

      real*4 tdat0(2),tdat(2)
      real*4 ttot1,ttot0
      
      real*8 utime(MAXCPU),utot
      integer ierr
      real etime

      Nyl=Ny/NPROC

      al(1)=1
      al(2)=1
      al(3)=1
      do i=1,nwax+1
       do k=1,2*nwaz+1
        do j=0,Nyl+1
         u(i,k,j)=(2.d0,1.d0)
         Lu1(i,k,j)=(1.d0,2.d0)
         Lu2(i,k,j)=(2.d0,1.d0)
        end do
       end do
      end do

      if(iam.eq.1) print *,'Testing CPU weights'

ctime      ttot0=etime(tdat0)
      do i=1,10
      call Fourier1 (u,uo)   
      call Fourier2 (uo,uo,uu,0,Nyl+1)
      call Fourier2 (uo,uo,uu,0,Nyl+1)
      call build_star (u,1.d0,Lu1,Lu2,Lu1,Lu2,al,al,al,al,2,u)
      call build_star (u,1.d0,Lu1,Lu2,Lu1,Lu2,al,al,al,al,2,u)
      enddo
ctime      ttot1=etime(tdat)
      utot=tdat(1)-tdat0(1)
      call rgather(utot,utime)

c      if(iam.eq.2) then
         open(unit=11,file='cpuweights')
         utot=1e6
         do i=1,NPROC
            if(utime(i).lt.utot) utot=utime(i)
         enddo
         do i=1,NPROC

      if (iam.eq.1) 
     + print *,'CPU(',i,'): estimated cpu weight=',utot/utime(i)*1000
            write(11,*) utot/utime(i)*1000
         enddo
         close(11)
c      endif


      end
#endif
