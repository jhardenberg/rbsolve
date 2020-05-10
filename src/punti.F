#include "config.h"

      subroutine punti (DD,ptan,qtan,yp,yv,qp,rp,sp,qv,rv,sv)
       
#include "param.h"
      real*8 yp(0:Ny+1),yv(0:Ny+1)
      real*8 qp(Ny),rp(Ny),sp(Ny),qv(Ny),rv(Ny),sv(Ny)
      real*8 ptan,qtan,DD
      integer*4 nstrat

      do j=1,Ny+1
       yp(j)=0.d0
       yv(j)=0.d0
      end do
      do j=1,Ny-1
       qp(j)=0.d0
       rp(j)=0.d0
       sp(j)=0.d0
      end do
      do j=1,Ny
       qv(j)=0.d0
       rv(j)=0.d0
       sv(j)=0.d0
      end do


      call stretch(yv(1),Ny,-DD/2,DD/2,-1.d0,0.d0,qtan,ptan)
      call stretch(yv(1),Ny,-DD/2,DD/2,0.d0,1.d0,qtan,ptan)

c does not create the lower BL stretching
c stretch(yv(1),Ny,-DD/2,DD/2,-1.d0,0.d0,qtan,ptan)
c does not create the upper BL stretching
c stretch(yv(1),Ny,-DD/2,DD/2,0.d0,1.d0,qtan,ptan)
c does not create the lower BL stretching

#ifndef HTROPO
      call stretch(yv(1),Ny,-DD/2,DD/2,-1.d0,1.d0,qtan,ptan)
#else 
      nstrat=Ny-Ntropo+1
      call stretch(yv(1),Ntropo,0.d0,HTROPO,-1.d0,1.d0,qtan,ptan)
      call stretch(yv(Ntropo),nstrat,HTROPO,DD,-1.d0,1.d0,qtan,ptan)      
#endif      

      yv(0)=2.d0*yv(1)-yv(2)
      yv(Ny+1)=2.d0*yv(Ny)-yv(Ny-1)
      
      do j=0,Ny
         yp(j)=(yv(j)+yv(j+1))/2.d0
      end do

      do j=1,Ny
       qp(j)=1.d0/((yp(j+1)-yp(j-1))*(yp(j)-yp(j-1)))
       rp(j)=1.d0/((yp(j+1)-yp(j))*(yp(j)-yp(j-1)))
       sp(j)=1.d0/((yp(j+1)-yp(j-1))*(yp(j+1)-yp(j)))
      end do

      do j=1,Ny
       qv(j)=1.d0/((yv(j+1)-yv(j-1))*(yv(j)-yv(j-1)))
       rv(j)=1.d0/((yv(j+1)-yv(j))*(yv(j)-yv(j-1)))
       sv(j)=1.d0/((yv(j+1)-yv(j-1))*(yv(j+1)-yv(j)))
      end do
      
      open (13,file='coord')

      do j=0,Ny+1
       write(13,*) 'j=',j,' ','yv=',' ',yv(j)
      end do
      do j=0,Ny
       write(13,*) 'j=',j,' ','yp=',' ',yp(j)
      end do

      close(13)
      
      return 
      end


      subroutine stretch (yv,nn,ymin,ymax,yeqmin,yeqmax,qtan,ptan)
       
#include "param.h"
      integer nn
      real*8 yv(nn),yeq(nn),ymax,ymin,yeqmin,yeqmax
      real*8 ptan,qtan,dyeq,yv1,yvn

      dyeq=(yeqmax-yeqmin)/dfloat(nn-1)
     
      do j=1,nn
       yeq(j)=yeqmin+ (j-1)*dyeq 
      enddo

#ifdef OLDPUNTI
      do j=1,nn
        yv(j)=  ptan*yeq(j)+
     +   (1.d0-ptan)* ( dsign(1.d0,yeq(j))*
     +        (1.d0-(dtanh(qtan*(1.d0-
     +        abs(yeq(j))))/dtanh(1.d0*qtan)))  )
      end do
#else
      if(qtan.ne.0) then
       do j=1,nn
        yv(j)=  tanh(yeq(j)*qtan)/tanh(qtan)
       end do
      else
       do j=1,nn
        yv(j)= yeq(j)
       end do
      endif
#endif
  
      yv1=yv(1)
      do j=1,nn
         yv(j)=yv(j)-yv1
      end do

      yvn=yv(nn)
      do j=1,nn
         yv(j)=yv(j)/ yvn *( ymax-ymin ) + ymin 
      end do
      
      print *,yv(1),yv(2),yv(nn-1),yv(nn)

      return 
      end
