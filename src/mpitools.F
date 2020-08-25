
#include "config.h"

        subroutine mpighost(ff)
*     Exchanges ghost boundary conditions using MPI_SENDRCV
*     calls.  Uses the new datatype myMPI_ALL_TYPE
*     defined in the main, exchanging all 3 velocities plus
*     scalar and temperature (if they exist) at once

#include "param.h"
cmpi        include 'mpif.h'

        complex*16 ff(nwax+1,2*nwaz+1,0:Nyl+1) 

#ifndef NOMPI
        integer sbuf,ierr,st(MPI_STATUS_SIZE),datatype
        
        datatype=myMPI_ALL_TYPE
       
        if(NPROC.ne.1) then
        if(iam.eq.1) then
      call MPI_SENDRECV(ff(1,1,1),1,datatype,MPI_PROC_NULL,99,
     &                   ff(1,1,Nyl+1),1,datatype,iam,99,
     &        MPI_COMM_WORLD,st,ierr)
      call MPI_SENDRECV(ff(1,1,Nyl),1,datatype,iam,98,
     &                   ff(1,1,0),1,datatype,MPI_PROC_NULL,98,
     &        MPI_COMM_WORLD,st,ierr)
        else if(iam.eq.NPROC) then
      call MPI_SENDRECV(ff(1,1,1),1,datatype,iam-2,99,
     &        ff(1,1,Nyl+1),1,datatype,MPI_PROC_NULL,99,
     &        MPI_COMM_WORLD,st,ierr)
      call MPI_SENDRECV(ff(1,1,Nyl),1,datatype,MPI_PROC_NULL,98,
     &                   ff(1,1,0),1,datatype,iam-2,98,
     &        MPI_COMM_WORLD,st,ierr)
        else
      call MPI_SENDRECV(ff(1,1,1),1,datatype,iam-2,99,
     &                   ff(1,1,Nyl+1),1,datatype,iam,99,
     &        MPI_COMM_WORLD,st,ierr)
      call MPI_SENDRECV(ff(1,1,Nyl),1,datatype,iam,98,
     &                   ff(1,1,0),1,datatype,iam-2,98,
     &        MPI_COMM_WORLD,st,ierr)
        endif
        endif
#endif
        end

        subroutine mpighost3(ff)
*     Exchanges ghost boundary conditions using MPI_SENDRCV
*     calls.  Uses the new datatype myMPI_VEL_TYPE
*     defined in the main, exchanging all 3 velocities 

#include "param.h"
cmpi        include 'mpif.h'

        complex*16 ff(nwax+1,2*nwaz+1,0:Nyl+1)

#ifndef NOMPI
        integer sbuf,ierr,st(MPI_STATUS_SIZE),datatype
        

        datatype=myMPI_VEL_TYPE

        if(NPROC.ne.1) then
        if(iam.eq.1) then
      call MPI_SENDRECV(ff(1,1,1),1,datatype,MPI_PROC_NULL,95,
     &                   ff(1,1,Nyl+1),1,datatype,iam,95,
     &        MPI_COMM_WORLD,st,ierr)
      call MPI_SENDRECV(ff(1,1,Nyl),1,datatype,iam,94,
     &                   ff(1,1,0),1,datatype,MPI_PROC_NULL,94,
     &        MPI_COMM_WORLD,st,ierr)
        else if(iam.eq.NPROC) then
      call MPI_SENDRECV(ff(1,1,1),1,datatype,iam-2,95,
     &        ff(1,1,Nyl+1),1,datatype,MPI_PROC_NULL,95,
     &        MPI_COMM_WORLD,st,ierr)
      call MPI_SENDRECV(ff(1,1,Nyl),1,datatype,MPI_PROC_NULL,94,
     &                   ff(1,1,0),1,datatype,iam-2,94,
     &        MPI_COMM_WORLD,st,ierr)
        else
      call MPI_SENDRECV(ff(1,1,1),1,datatype,iam-2,95,
     &                   ff(1,1,Nyl+1),1,datatype,iam,95,
     &        MPI_COMM_WORLD,st,ierr)
      call MPI_SENDRECV(ff(1,1,Nyl),1,datatype,iam,94,
     &                   ff(1,1,0),1,datatype,iam-2,94,
     &        MPI_COMM_WORLD,st,ierr)
        endif
        endif
#endif
        end

        subroutine mpighost1(ff)
*     Exchanges ghost boundary conditions using MPI_SENDRCV
*     calls. Exchanges BCs for one single field

#include "param.h"
cmpi        include 'mpif.h'

        complex*16 ff(nwax+1,2*nwaz+1,0:Nyl+1)

#ifndef NOMPI
        integer sbuf,ierr,st(MPI_STATUS_SIZE),datatype

        datatype=MPI_COMPLEX16
        sbuf=(nwax+1)*(2*nwaz+1)

        if(NPROC.ne.1) then
        if(iam.eq.1) then
      call MPI_SENDRECV(ff(1,1,1),sbuf,datatype,MPI_PROC_NULL,93,
     &                   ff(1,1,Nyl+1),sbuf,datatype,iam,93,
     &        MPI_COMM_WORLD,st,ierr)
      call MPI_SENDRECV(ff(1,1,Nyl),sbuf,datatype,iam,92,
     &                   ff(1,1,0),sbuf,datatype,MPI_PROC_NULL,92,
     &        MPI_COMM_WORLD,st,ierr)
        else if(iam.eq.NPROC) then
      call MPI_SENDRECV(ff(1,1,1),sbuf,datatype,iam-2,93,
     &        ff(1,1,Nyl+1),sbuf,datatype,MPI_PROC_NULL,93,
     &        MPI_COMM_WORLD,st,ierr)
      call MPI_SENDRECV(ff(1,1,Nyl),sbuf,datatype,MPI_PROC_NULL,92,
     &                   ff(1,1,0),sbuf,datatype,iam-2,92,
     &        MPI_COMM_WORLD,st,ierr)
        else
      call MPI_SENDRECV(ff(1,1,1),sbuf,datatype,iam-2,93,
     &                   ff(1,1,Nyl+1),sbuf,datatype,iam,93,
     &        MPI_COMM_WORLD,st,ierr)
      call MPI_SENDRECV(ff(1,1,Nyl),sbuf,datatype,iam,92,
     &                   ff(1,1,0),sbuf,datatype,iam-2,92,
     &        MPI_COMM_WORLD,st,ierr)
        endif
        endif
#endif
        end




      subroutine mpighostup(ff)
*     Exchanges ghost boundary conditions using MPI_SENDRCV
*     calls.
*     This one sends up only (used for p)

#include "param.h"
cmpi        include 'mpif.h'

        complex*16 ff(nwax+1,2*nwaz+1,0:Nyl+1) 

#ifndef NOMPI
        integer sbuf,ierr,st(MPI_STATUS_SIZE)
        
        sbuf=(nwax+1)*(2*nwaz+1)
       
        if(NPROC.ne.1) then

        if(iam.eq.1) then

      call MPI_SENDRECV(ff(1,1,Nyl),sbuf,
     &        MPI_COMPLEX16,iam,98,
     &        ff(1,1,0),sbuf,MPI_COMPLEX16,MPI_PROC_NULL,98,
     &        MPI_COMM_WORLD,st,ierr)

        else if(iam.eq.NPROC) then
        
      call MPI_SENDRECV(ff(1,1,Nyl),sbuf,MPI_COMPLEX16,
     &        MPI_PROC_NULL,98,
     &        ff(1,1,0),sbuf,MPI_COMPLEX16,iam-2,98,
     &        MPI_COMM_WORLD,st,ierr)

        else
           
      call MPI_SENDRECV(ff(1,1,Nyl),sbuf,
     &        MPI_COMPLEX16,iam,98,
     &        ff(1,1,0),sbuf,MPI_COMPLEX16,iam-2,98,
     &        MPI_COMM_WORLD,st,ierr)

        endif

        endif

#endif
        
        end



        subroutine mpighostdown(ff)
*     Exchanges  ghost boundary conditions 
*     Uses the new datatype myMPI_VEL_TYPE
*     defined in the main, exchanging all three velocities at once
*     This one sends down only

#include "param.h"
cmpi        include 'mpif.h'

        complex*16 ff(nwax+1,2*nwaz+1,0:Nyl+1) 

#ifndef NOMPI
        integer sbuf,ierr,st(MPI_STATUS_SIZE),datatype

        datatype=myMPI_VEL_TYPE
        
        if(NPROC.ne.1) then

        if(iam.eq.1) then

      call MPI_SENDRECV(ff(1,1,1),1,datatype,MPI_PROC_NULL,97,
     &                   ff(1,1,Nyl+1),1,datatype,iam,97,
     &        MPI_COMM_WORLD,st,ierr)

        else if(iam.eq.NPROC) then
        
      call MPI_SENDRECV(ff(1,1,1),1,datatype,iam-2,97,
     &        ff(1,1,Nyl+1),1,datatype,MPI_PROC_NULL,97,
     &        MPI_COMM_WORLD,st,ierr)

        else
           
      call MPI_SENDRECV(ff(1,1,1),1,datatype,iam-2,97,
     &                   ff(1,1,Nyl+1),1,datatype,iam,97,
     &        MPI_COMM_WORLD,st,ierr)

        endif

        endif

#endif

        end



      subroutine rgatherroot(value,array)
#include "param.h"
cmpi        include 'mpif.h'
      real*8 value(3),array(3,MAXCPU)
      integer ierr

      array(1,iam)=value(1)
      array(2,iam)=value(2)
      array(3,iam)=value(3)
      if(NPROC.ne.1) then
#ifndef NOMPI
        call MPI_GATHER(value,3,MPI_REAL8,
     &              array,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
#endif
      endif
      end




      subroutine rgather(value,array)
#include "param.h"
cmpi        include 'mpif.h'
      real*8 value,array(MAXCPU)
      integer ierr
      array(iam)=value
      if(NPROC.ne.1) then
#ifndef NOMPI
        call MPI_ALLGATHER(value,1,MPI_REAL8,
     &              array,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
#endif
      endif
      end



      subroutine rbroadcast(value,iroot)
* Broadcasts a real value from root
#include "param.h"
cmpi        include 'mpif.h'
      real*8 value
      integer ierr,iroot
      if(NPROC.ne.1) then
#ifndef NOMPI
        call MPI_BCAST ( value,1,MPI_REAL8 ,iroot-1,MPI_COMM_WORLD,ierr) 
#endif
      endif
      end






        
