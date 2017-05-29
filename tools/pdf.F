c  calcola pdf o

        
#include "config.h" 

        program main    

#include "param.h"

        integer nc,nn,npt,n
	real h(1000),aa,bb,xx,dx,xmin,xmax,xnorm,xmino,xmaxo
        character cc*1
        character*30 str
        integer istart,iend,istep,it,jmax,jmin,jj,iargc
        complex*16 u(nwax+1,2*nwaz+1,0:Nylmem+1)
        real*8 uo(Nx,Nz,0:Nylmem+1),umax,umin
        character*45 file  
        character*12 file0




#ifndef NOMPI 
      print *,"Please compile with the NOMPI option!"
      stop
#endif   
       iam=1

        if(iargc().lt.3) then
         print *, 'Usage: pdf cc jj nc xmino xmaxo istart [iend istep]'
         stop
      endif
   
       call initfour(0)
       call slice()
 
       call getarg(1,cc)
    
       call getarg(2,str)        
       read(str,*) jj
       
       call getarg(3,str)
       read(str,*) nc 

       call getarg(4,str)
       read(str,*) xmino 

       call getarg(5,str)
       read(str,*) xmaxo 

       call getarg(6,str)
       read(str,*) istart 

       

       if(iargc().gt.6) then
       call getarg(7,str)
       read(str,*) iend
       call getarg(8,str)
       read(str,*) istep
      else
        iend=istart
        istep=1
      endif

      do it=istart,iend,istep

         call makenamebob(file0,cc,it,'.ppp')
         call readfield(u,cc,it)
         call Fourier1 (u,uo) 

        do i=1,1000
         h(i)=0
        enddo

          xmin=1e10
          xmax=-1e10

  
        if(xmino.eq.0.and.xmaxo.eq.0) then    

               do k=1,Nz
                 do i=1,Nx
                   if(uo(i,k,jj).gt.xmax) then
                      xmax=uo(i,k,jj)
                   endif
                   if(uo(i,k,jj).lt.xmin) then
                      xmin=uo(i,k,jj)
                   endif
                 enddo
               enddo
        write(*,*) xmin
        write(*,*) xmax
          
         else
           xmin=xmino
           xmax=xmaxo

         endif
 
        dx=abs((xmax-xmin)/(real(nc)-1))
        
        npt=0
       



        do k=1,Nz
          do i=1,Nx
            xx=uo(i,k,jj)
                 n=(xx-xmin)/dx+1
                 if(n.le.nc) then
           
                   if(n.gt.0) then
             
                     h(n)=h(n)+1 
                     npt=npt+1
                 endif
               endif
     
         enddo
       enddo 

       xnorm=1.0/(npt*dx)
       
       write(*,*) 'N punti:',npt
       write(*,*) 'dx=',dx 
       file=file0


       open(unit=20,file=file,form='formatted')

        do i=1,nc
           write(20,*) xmin+dx*(i-1),h(i)*xnorm
        enddo
       close(20)
      enddo
      end
