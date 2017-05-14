c Converts NS/RB files to bob format
#include "config.h"

      program main

#include "param.h"


      complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1)
      real*8 uo(Nx,Nz,0:Nylmem+1),umax,umin
      integer jmax,jmin,step,iargc
      character*12 file0
	character*45 file
      character cc*1,str*30,aa*30


#ifndef NOMPI 
      print *,"Please compile with the NOMPI option!"
      stop
#endif      
      
      iam=1

      if(iargc().lt.3) then
         print *, 'Usage: anomaly cc out time'
         stop
      endif
      call initfour(0)
      call slice()


      call getarg(1,cc)

      call getarg(2,aa)
      
      call getarg(3,str)
      read(str,*) step 

      if(iargc().gt.3) then
          print *, 'Usage: anomaly cc out time'
         stop
      endif


             call readfield(u,cc,step)

         do j=0,Ny+1
            u(1,1,j)=(0.d0,0.d0)
         enddo
 
             
        call savefield(u,aa,step)
     
      end








