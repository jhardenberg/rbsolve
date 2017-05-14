#include "config.h"

c **********************************************************************

      subroutine savefield (ut,ch,nsalv,rtime)
*    Saves one field at step nsalv (w/o ghost points)

#include "param.h"
 
      character*1 ch
      complex*16 ut(nwax+1,2*nwaz+1,0:Nylmem+1)
      integer*4 nsalv
      real*8    rtime

* Two options: distributed saving or gather everything to the root node first

#ifdef ROOT_SAVE
      call saveroot(ut,ch,nsalv)
#else
      call savelocal(ut,ch,nsalv)
#endif
      open(10,file='nrec.d',form='formatted',status='unknown')
      write(10,*) nsalv,rtime
      close(10)

      end


      subroutine readfield (ut,ch,nsalv)
*    Reads one field at step nsalv (w/o ghost points)

#include "param.h"
 
      character*1 ch
      complex*16 ut(nwax+1,2*nwaz+1,0:Nylmem+1)
      integer*4 nsalv

* Two options: distributed saving or gather everything to the root node first

#ifdef ROOT_SAVE
      call readroot(ut,ch,nsalv)
#else
      call readlocal(ut,ch,nsalv)
#endif
      end


#ifdef ROOT_SAVE
      subroutine saveroot (ut,ch,nsalv)
*     Gathers everything to the node SAVENODE
*     and saves

#include "param.h"
cmpi        include 'mpif.h'

      character*1 ch
      integer nsalv,count

      complex*16 ut(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 utemp(nwax+1,2*nwaz+1,Nylmem) 
      integer ierr,st(MPI_STATUS_SIZE)
      character*12 file0

      if(iam.eq.1) write(*,*)'Saving ',ch, ',  t=', nsalv

      if(iam.eq.SAVENODE) then
* I am SAVENODE, so I have to listen and save
         call makename(file0,ch,nsalv)
#ifdef FORMATTED_OUTPUT
         open(12,file=file0,form='formatted',status='unknown')         
#else
         open(12,file=file0,form='unformatted',status='unknown')         
#endif
         do i=1,NPROC
            if(i.ne.iam) then
*     Listen to the others
               count=(nwax+1)*(2*nwaz+1)*Nyls(i)
*	  print *,iam,' listening for ',i , 'to send'
               call MPI_RECV(utemp,count,MPI_COMPLEX16,   
     +  i-1, 70+i, MPI_COMM_WORLD, st, ierr)

               call saveslice4(utemp,Nyls(i))
            else
*     Use my own data
               call saveslice4(ut(1,1,1),Nyl)
            endif
         enddo
         close(12)
      else
* I am not SAVENODE: send my data to SAVENODE
        count=(nwax+1)*(2*nwaz+1)*Nyl

*	  print *,iam,' sending to ',SAVENODE, 'to send'
        call MPI_SEND(ut(1,1,1),count,MPI_COMPLEX16,   
     +  SAVENODE-1, 70+iam, MPI_COMM_WORLD,st, ierr)
       
      endif
   
      end
#endif

      subroutine saveghost (ut,ch,nsalv)
*    Saves one field at step nsalv (with ghost points)

#include "param.h"

      character*1 ch
      complex*16 ut(nwax+1,2*nwaz+1,0:Nylmem+1)
      integer*4 nsalv
      character*12 file0

      if(iam.eq.1) write(*,*)'Saving (',iam,') ',ch, ',  t=', nsalv

      call makename(file0,ch,nsalv)
      open(12,file=file0,form='formatted',status='unknown')
      do j=0,Nyl+1 
          write(12,*) j,real(ut(4,1,j)),real(ut(1,1,j))
      enddo
      close(12)
      return
      end

                  
      subroutine savelocal (ut,ch,nsalv)
*    Saves one field at step nsalv (w/out ghost points)

#include "param.h"
 
      character*1 ch
      complex*16 ut(nwax+1,2*nwaz+1,0:Nylmem+1)
      integer*4 nsalv
      character*12 file0
            
      if(iam.eq.1) write(*,*)'Saving (',iam,') ',ch, ',  t=', nsalv
      call makename(file0,ch,nsalv)
#ifdef FORMATTED_OUTPUT
      open(12,file=file0,form='formatted',status='unknown')         
#else
      open(12,file=file0,form='unformatted',status='unknown')         
#endif
      call saveslice4(ut(1,1,1),Nyl)
      close(12)
      return
      end


      subroutine readlocal (ut,ch,nsalv)
      
*    Reads one field at step nsalv (w/out ghost points)

#include "param.h"
 
      character*1 ch
      complex*16 ut(nwax+1,2*nwaz+1,0:Nylmem+1)
      integer*4 nsalv
      character*12 file0
            
      if(iam.eq.1) write(*,*)'Reading (',iam,') ',ch, ',  t=', nsalv
      call makename(file0,ch,nsalv)
#ifdef FORMATTED_OUTPUT
      open(12,file=file0,form='formatted',status='unknown') 
#else
      open(12,file=file0,form='unformatted',status='unknown')
#endif

      call readslice4(ut(1,1,1),Nyl)
      close(12)

      return
      end



      subroutine saveslice4 (ut,Nys)
#include "param.h"
      integer Nys
      complex*16 ut(nwax+1,2*nwaz+1,Nys)
      complex*8  ut4(nwax+1,2*nwaz+1)

#ifdef FORMATTED_OUTPUT
      do j=1,Nys
         do k=1,2*nwaz+1
            do i=1,nwax+1
               write(12,'(e16.7)') ut(i,k,j)
            enddo
         enddo
      enddo
#else
      do j=1,Nys
         do k=1,2*nwaz+1
            do i=1,nwax+1
               ut4(i,k)=ut(i,k,j)
            enddo
         enddo
         write(12) ((ut4(i,k),i=1,nwax+1),k=1,2*nwaz+1)
      enddo
#endif
      return
      end

      subroutine saveslice1 (ut)
#include "param.h"
      complex*16 ut(nwax+1,2*nwaz+1)
      complex*8  ut4(nwax+1,2*nwaz+1)

#ifdef FORMATTED_OUTPUT
         do k=1,2*nwaz+1
            do i=1,nwax+1
               write(12,'(e16.7)') ut(i,k,j)
            enddo
         enddo
#else
         do k=1,2*nwaz+1
            do i=1,nwax+1
               ut4(i,k)=ut(i,k)
            enddo
         enddo
         write(12) ((ut4(i,k),i=1,nwax+1),k=1,2*nwaz+1)
#endif
      return
      end


      subroutine readslice4 (ut,Nys)
#include "param.h"
      integer Nys
      complex*16 ut(nwax+1,2*nwaz+1,Nys)
      complex*8 ut4(nwax+1,2*nwaz+1)
#ifdef FORMATTED_OUTPUT
      do j=1,Nys
         read(12,*) ((ut(i,k,j),i=1,nwax+1),k=1,2*nwaz+1) 
      enddo
#else
      do j=1,Nys
         read(12) ((ut4(i,k),i=1,nwax+1),k=1,2*nwaz+1)
         do k=1,2*nwaz+1
            do i=1,nwax+1
               ut(i,k,j)=ut4(i,k)
            enddo
         enddo
      enddo
#endif
      return
      end


      subroutine readslice1 (ut,ifile)
#include "param.h"
      integer ifile
      complex*16 ut(nwax+1,2*nwaz+1)
      complex*8 ut4(nwax+1,2*nwaz+1)
#ifdef FORMATTED_OUTPUT
         read(ifile,*) ((ut(i,k),i=1,nwax+1),k=1,2*nwaz+1)
#else
         read(ifile) ((ut4(i,k),i=1,nwax+1),k=1,2*nwaz+1)
         do k=1,2*nwaz+1
            do i=1,nwax+1
               ut(i,k)=ut4(i,k)
            enddo
         enddo
#endif
      return
      end



#ifdef ROOT_SAVE
      subroutine readroot (ut,ch,nsalv)
*     Reads on node  SAVENODE
*     and distributes to the other nodes

#include "param.h"
cmpi        include 'mpif.h'

      character*1 ch
      integer nsalv,count

      complex*16 ut(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 utemp(nwax+1,2*nwaz+1,Nylmem) 
      integer ierr,st(MPI_STATUS_SIZE)
      character*12 file0

      if(iam.eq.1) write(*,*)'Reading ',ch, ',  t=', nsalv

      if(iam.eq.SAVENODE) then
* I am SAVENODE, so I have to read and send to the others
         call makename(file0,ch,nsalv)
#ifdef FORMATTED_OUTPUT
         open(12,file=file0,form='formatted',status='unknown')         
#else
         open(12,file=file0,form='unformatted',status='unknown')         
#endif
         do i=1,NPROC
            if(i.ne.iam) then
*     Read and send to the others
               call readslice4(utemp,Nyls(i))
               count=(nwax+1)*(2*nwaz+1)*Nyls(i)
               call MPI_SEND(utemp,count,MPI_COMPLEX16,   
     +  i-1, 50+i, MPI_COMM_WORLD, st,ierr)
            else
*     Read directly
               call readslice4(ut(1,1,1),Nyl)
            endif
         enddo
         close(12)
      else
* I am not SAVENODE: get my data from SAVENODE
        count=(nwax+1)*(2*nwaz+1)*Nyl
        call MPI_RECV(ut(1,1,1),count,MPI_COMPLEX16,   
     +  SAVENODE-1, 50+iam, MPI_COMM_WORLD, st,ierr)
      endif
      end
#endif

#ifdef FORMATTED_OUTPUT
#define EXTENSION  '.dat'
#define EXTENSION2  '.d'
#else
#define EXTENSION  '.unf'
#define EXTENSION2  '.u'
#endif

      subroutine makename (file0,ch,nsalv)      

 
      character ch,SYMBOL(12)
      integer*4 nsalv,i
      character*12 file0

#if defined( ROOT_SAVE ) || defined ( NOMPI )
      write(file0,'(A1,I7,A4)') ch,nsalv,EXTENSION
#else
      write(file0,'(A1,I7,A2,i2)') ch,nsalv,EXTENSION2,iam
#endif
      read(file0,'(12A1)') SYMBOL
      do i=1,12
        if(SYMBOL(i).eq.' ') SYMBOL(i)='0'  
        if(SYMBOL(i).eq.'*') SYMBOL(i)='X'
       enddo
      write(file0,'(12A1)') SYMBOL
      return
      end

      subroutine makenamebob (file0,ch,nsalv,ext)      

      character ch,SYMBOL(12)
      integer*4 nsalv,i
      character*12 file0
      character*4 ext


      write(file0,'(A1,I7,A4)') ch,nsalv,ext

      read(file0,'(12A1)') SYMBOL
      do i=1,12
        if(SYMBOL(i).eq.' ') SYMBOL(i)='0'  
        if(SYMBOL(i).eq.'*') SYMBOL(i)='X'
       enddo
      write(file0,'(12A1)') SYMBOL
      return
      end



      subroutine readold (ut,ch,nsalv)
      
*    Reads fields saved in the old format at step nsalv (w/out ghost points)

#include "param.h"
 
      character*1 ch
      complex*16 ut(nwax+1,2*nwaz+1,0:Nylmem+1)
      integer*4 nsalv
      character*5 file0
            
      write(*,*)'Reading (OLD) (',iam,') ',ch, ',  t=', nsalv

      if(nsalv.lt.10) then
         write(file0,'(a1,a1,i1)')ch,'.',nsalv
      elseif(nsalv.gt.9.and.nsalv.lt.100) then
         write(file0,'(a1,a1,i2)')ch,'.',nsalv
      elseif(nsalv.gt.99) then
         write(file0,'(a1,a1,i3)')ch,'.',nsalv
      endif

      open(12,file=file0,form='unformatted',status='unknown') 
      read(12) (((ut(i,k,j),j=0,Nyl+1),i=1,nwax+1),k=1,2*nwaz+1)
      close(12)

      return
      end   

/*
      subroutine savebob (uo,ch,nsalv,xmin,xmax)
      
#include "param.h"
 
      character*1 ch
      real*8 uo(Nx,Nz,0:Nylmem+1)
      real xmin,xmax,center,cff
      integer*4 nsalv,ix
      character*12 file0
      character*1 bob(Nx,Nz,Ny)
      integer ist


      call makenamebob(file0,ch,nsalv,'.bob')
      print *,'saving ', file0, ' ', xmin,xmax

      cff=255/(xmax-xmin)
      do j=1,Ny
        do k=1,Nz
         do i=1,Nx
          ix=int(cff*(uo(i,k,j)-xmin))
          if(ix.lt.0) ix=0
          if(ix.gt.255) ix=255
          bob(i,k,j)=char(ix)
         enddo
        enddo
       enddo
       
      open(unit=12,file=file0,access='DIRECT',recl=Nx*Ny*Nz/NWORD)
      write(12,rec=1,iostat=ist) bob
      close(12)
      if(ist.ne.0) write(*,*) 'I/O Error writing ',file0
      close(12)
c      open(unit=12,file=file0,access='append',form='formatted')
c      write(12,*) 'MIN',xmin,' MAX',xmax
c      close(12)

      return
      end
*/
