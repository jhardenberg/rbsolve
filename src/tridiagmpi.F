#include "config.h"


*     Parallel tridiagonal solver
*     J.v.H. & G.P. 2002

*     Based on Mattor et al., Parallel Computing 21 (11) 1769--1782 (1995)


#define NR (2*NPR-2)
#define REAL real*8
#define COMPLEX complex*16

*     Tridiagonal MPI solver

      subroutine tridlhs(a,b,c,reda,redb,redc,xlh,xuh,
     &                   wom,womd,M,N,NPR,ip)

*     Solves all that can be solved once for all at the beginning
*     if the LHS is constant
*     Input:  a,b,c
*     Output: b,c,reda,redb,redc,xlh,xuh
*     Temporary: wom,womd

      implicit none
      integer ip,M,N,NPR,i
      REAL reda(N,0:NR+1),redb(N,0:NR+1)
      REAL redc(N,0:NR+1)
      REAL xlh(N,*),xuh(N,*)
      REAL a(N,*),b(N,*),c(N,*),wom(N,*),womd(N,*)

      call trid1pass(a,b,c,xlh,xuh,wom,womd,M,N)

#ifndef NOMPI
      if(NPR.ne.1) then
        call buildr_lhs(reda,redb,redc,xlh,xuh,M,N,NPR,ip) 
      endif
#endif

      end

      subroutine tridmpi(a,b,c,r,reda,redb,redc,redr,
     &                   xlh,xuh,xr,M,N,NPR,ip)
*     Solves the tridiagonal system using the output of tridlhs
*     Input:  a,b,c,r,reda,redb,redc,xlh,xuh
*     Output: xr
*     Temporary workspace: redr

      implicit none
cmpi      include 'mpif.h'
      integer ip,M,N,j,i,NPR
      REAL reda(N,0:NR+1),redb(N,0:NR+1)
      REAL redc(N,0:NR+1)
      COMPLEX redr(N,0:NR+1)
      REAL xlh(N,*),xuh(N,*)
      COMPLEX xr(N,*)
      REAL a(N,*),b(N,*),c(N,*)
      COMPLEX r(N,*)
	integer ierr
      call trid2pass(a,b,c,r,xr,M,N)

#ifndef NOMPI
      if(NPR.ne.1) then

      call buildr_rhs(redr,xr,M,N,NPR,ip)

*     Here we pass the POINTERS to row 1 of the reduced
*     arrays using the trick of passing element (1,1)

c      call tridag(reda(1,1),redb(1,1),redc(1,1),redr(1,1),
c     &            redr(1,1),NR,N)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call tridag(reda,redb,redc,redr,redr,NR,N)
      call tridassemble(xr,xlh,xr,xuh,redr,
     &            M,N,NPR,ip)
      
      endif
#endif

      end

#ifndef NOMPI

      subroutine buildr_lhs(reda,redb,redc,xlh,xuh,M,N,NPR,ip)
*     Builds the LHS of the reduced system.
*     Communication intensive but can be run only once
*     if the LHS matrix is constant.

      implicit none
cmpi      include 'mpif.h'

      integer N,M,NPR,NRl
      REAL reda(N,0:NR+1),redb(N,0:NR+1)
      REAL redc(N,0:NR+1)
      REAL sreda(N,0:1),sredb(N,0:1),sredc(N,0:1)
      REAL xlh(N,*),xuh(N,*)
      integer ip,i,ierr

      if(ip.eq.1) then
         do i=1,N
            sredb(i,1)=xlh(i,M)
            sreda(i,1)=0
            sredc(i,1)=-1
         enddo
      else if (ip.eq.NPR) then
         do i=1,N
            sredb(i,0)=xuh(i,1)
            sreda(i,0)=-1
            sredc(i,0)=0
         enddo
      else
         do i=1,N
          sreda(i,0)=-1
          sredb(i,0)=xuh(i,1)
          sredc(i,0)=xlh(i,1)
          sreda(i,1)=xuh(i,M)
          sredb(i,1)=xlh(i,M)
          sredc(i,1)=-1
         enddo
      endif

      call MPI_ALLGATHER(sreda,2*N,MPI_REAL8,
     &              reda,2*N,MPI_REAL8,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHER(sredb,2*N,MPI_REAL8,
     &              redb,2*N,MPI_REAL8,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHER(sredc,2*N,MPI_REAL8,
     &              redc,2*N,MPI_REAL8,MPI_COMM_WORLD,ierr)

      end

      subroutine buildr_rhs(redr,xr,M,N,NPR,ip)
*     Builds the RHS of the reduced system.
*     This is where all communication occurs.

      implicit none
cmpi      include 'mpif.h'

      integer N,M,NPR
      COMPLEX redr(N,0:NR+1)
      COMPLEX sredr(N,0:1)
      COMPLEX xr(N,*)
      integer ip,i,ierr

      if(ip.eq.1) then
         do i=1,N
            sredr(i,1)=-xr(i,M)
         enddo
      else if (ip.eq.NPR) then
         do i=1,N
          sredr(i,0)=-xr(i,1)
         enddo
      else
         do i=1,N
            sredr(i,0)=-xr(i,1)
            sredr(i,1)=-xr(i,M)
         enddo
      endif

      call MPI_ALLGATHER(sredr,2*N,MPI_COMPLEX16,
     &              redr,2*N,MPI_COMPLEX16,MPI_COMM_WORLD,ierr)

      end

*       NOMPI
#endif

      subroutine tridassemble(xs,xlh,xr,xuh,reds,M,N,NPR,ip)
*     Assembles the final solution using the solution of the
*     reduced system

      implicit none

      integer i,j,M,ip,N,ip2m2,ip2m1,NPR
      COMPLEX xr(N,M)
      REAL xlh(N,M),xuh(N,M)
      COMPLEX reds(N,0:M),xs(N,M)
      COMPLEX uhcoeff(N),lhcoeff(N)

      ip2m1=ip*2-1
      ip2m2=ip*2-2

      if(ip.eq.1) then
         do i=1,N
            lhcoeff(i)=reds(i,ip2m1)
         enddo
         do j=1,M
            do i=1,N
               xs(i,j)=xr(i,j)+lhcoeff(i)*xlh(i,j)
            enddo
         enddo
      else if(ip.eq.NPR) then
         do i=1,N
            uhcoeff(i)=reds(i,ip2m2)
         enddo
         do j=1,M
            do i=1,N
               xs(i,j)=xr(i,j)+uhcoeff(i)*xuh(i,j)
            enddo
         enddo
      else
         do i=1,N
           uhcoeff(i)=reds(i,ip2m2)
           lhcoeff(i)=reds(i,ip2m1)
         enddo
         do j=1,M
            do i=1,N
               xs(i,j)=xr(i,j)+lhcoeff(i)*xlh(i,j)+uhcoeff(i)*xuh(i,j)
            enddo
         enddo
      endif
   
      end

      subroutine trid1pass(a,b,c,xlh,xuh,wom,womd,M,N)
*     Precomputes all arrays which can be 
*     computed once for all at the beginning
*     if the LHS is constant
*     omega is stored in b and the denominators in c
*     wom and womd are workspace arrays 
*     r and xr could be used for these ...

      implicit none

      integer N,M
      REAL a(N,M),b(N,M),c(N,M)
      REAL xlh(N,M),xuh(N,M),womd(N,M),wom(N,M)
      REAL denom
      integer i,j

*     Forward substitution
      do i=1,N
         womd(i,1)=1./b(i,1)
         wom(i,1)=c(i,1)*womd(i,1)
      enddo
      do j=2,M
         do i=1,N
           denom=b(i,j)-a(i,j)*wom(i,j-1)
*          if(denom.eq.0) call diverror()
           womd(i,j)=1./denom
           wom(i,j)=c(i,j)*womd(i,j)
         enddo
      enddo

*     Back substitution
      do i=1,N
         xlh(i,M)=-wom(i,M)
         xuh(i,M)=a(i,M)/b(i,M)
      enddo
      do j=M-1,1,-1
         do i=1,N
            xlh(i,j)=-wom(i,j)*xlh(i,j+1)
            denom=b(i,j)-c(i,j)*xuh(i,j+1)
*           if(denom.eq.0) diverror()
            xuh(i,j)=a(i,j)/denom
         enddo
      enddo
      
*     Forward again
      do i=1,N
         xuh(i,1)=-xuh(i,1)
         b(i,1)=wom(i,1)
         c(i,1)=womd(i,1)
      enddo
      do j=2,M
         do i=1,N
            xuh(i,j)=-xuh(i,j)*xuh(i,j-1)
            b(i,j)=wom(i,j)
            c(i,j)=womd(i,j)
         enddo
      enddo

      end

      subroutine trid2pass(a,om,omd,r,xr,M,N)
*     Computes the xr array
*     om and omd are usually stored in b and c respectively

      implicit none

      integer N,M

      REAL a(N,M)
      COMPLEX r(N,M),xr(N,M)
      REAL omd(N,M),om(N,M)
      integer i,j

      do i=1,N
         xr(i,1)=r(i,1)*omd(i,1)
      enddo
*     Forward sobstitution
      do j=2,M
         do i=1,N
            xr(i,j)=(r(i,j)-a(i,j)*xr(i,j-1))*omd(i,j)
         enddo
      enddo
      
*     Back substitution
      do j=M-1,1,-1
         do i=1,N
            xr(i,j)=xr(i,j)-om(i,j)*xr(i,j+1)
         enddo
      enddo
 
      end

#ifndef NOMPI
      subroutine TRIDAG(A,B,C,TN,U,Nn,N)
*     Solves N tridiagonal systems using the Thomas algorithm
*     Based on a modification of the NR routine.
*     Safety checks are commented out for greater speed
*     on pipelining CPUs.
*     REAL version


#include "param.h"
#define NX ((nwax+1)*(2*nwaz+1)) 
#define NNR (2*NPROC-2) 

      integer Nn,N
      REAL A(N,0:Nn+1),B(N,0:Nn+1),C(N,0:Nn+1)
      COMPLEX TN(N,0:Nn+1),U(N,0:Nn+1)
      REAL bet(NX)
      REAL	gam( NX  ,NNR )


      do i=1,N
c         if (B(i,1).eq.0.d0) then
c		 print *,i,'b!'
c		stop
c	 endif
         bet(i)=1./B(i,1)
         U(i,1)=TN(i,1)*bet(i)
      enddo

      do 11 j=2,Nn
         do i=1,N
          gam(i,j)=C(i,j-1)*bet(i)
c         if ((B(i,j)-A(i,j)).eq.0.d0) then
c		print *,'ab!'
c		stop
c	  endif
          bet(i)=1./(B(i,j)-A(i,j)*gam(i,j))
c          if (bet(i).eq.0.d0)  print *,'bet!'
          U(i,j)=(TN(i,j)-A(i,j)*U(i,j-1))*bet(i)
         enddo
11    continue

      do 12 j=Nn-1,1,-1
         do i=1,N
          U(i,j)=U(i,j)-gam(i,j+1)*U(i,j+1)
         enddo
12    continue

      return
      end
#endif
