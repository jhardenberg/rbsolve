#include "config.h"

      subroutine build_star (vel_1,dt,diff_1,diff_2,
     +   conv_1,conv_2,alfa,beta,gamma,zeta,l,ve_star)

*     Assembles new star fields (timestepping)

#include "param.h"
     
      complex*16 diff_1(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 diff_2(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 conv_1(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 conv_2(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 vel_1(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 ve_star(nwax+1,2*nwaz+1,0:Nylmem+1)
      real*8 alfa(3),beta(3),gamma(3),zeta(3),dt
      real*8 cff1,cff2,cff3,cff4
 
      integer*4 l
     
         cff1=dt*alfa(l)
         cff2=dt*beta(l)
         cff3=dt*gamma(l)
         cff4=dt*zeta(l)

      if(l.eq.1) then

      do j=Nylmin,Nylmax
        do k=1,2*nwaz+1
 	 do i=1,nwax+1
          ve_star(i,k,j)=vel_1(i,k,j)+
     +    cff2*diff_1(i,k,j)-cff3*conv_1(i,k,j)
         enddo      
        enddo      
      enddo      

      else

      do j=Nylmin,Nylmax 
        do k=1,2*nwaz+1
 	 do i=1,nwax+1
          ve_star(i,k,j)=vel_1(i,k,j)+cff1*(diff_2(i,k,j)+
     +    diff_1(i,k,j))-cff3*conv_1(i,k,j)-cff4*conv_2(i,k,j)
         enddo      
        enddo      
      enddo      

      endif


      return
      end



      subroutine add_star_forcing (vel_1,dt,
     +   alfa,beta,gamma,zeta,l,ve_star,tauw)

*     Adds forcing

#include "param.h"
     
      complex*16 vel_1(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 ve_star(nwax+1,2*nwaz+1,0:Nylmem+1)
      real*8 alfa(3),beta(3),gamma(3),zeta(3),dt,tauw
      real*8 cff1,cff2,cff3,cff4
 
      integer*4 l
     
         cff1=dt*alfa(l)
         cff2=dt*beta(l)
         cff3=dt*gamma(l)
         cff4=dt*zeta(l)


        cff1=tauw*(alfa(l)+beta(l))*dt
        do j=Nylmin,Nylmax
         ve_star(1,1,j)=vel_1(1,1,j)+cff1
        enddo

 
      return
      end

#ifdef TEMPERATURE

      subroutine add_star_temp (vel_1,dt,t_1,t_2,
     +  alfa,beta,gamma,zeta,l,ve_star,rapr)
     
#include "param.h"
     
      complex*16 t_1(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 t_2(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 vel_1(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 ve_star(nwax+1,2*nwaz+1,0:Nylmem+1)
      real*8 alfa(3),beta(3),gamma(3),zeta(3),dt,tauw
      real*8 cff1,cff2,cff3,cff4,rapr,cff1t,cff2t
 
      integer*4 l

      cff1=dt*alfa(l)*rapr
      cff2=dt*beta(l)*rapr

* XXXXXXXXXXXX

      if(l.eq.1) then

      do j=Nylmin,Nylmax
        do k=1,2*nwaz+1
 	 do i=1,nwax+1
          ve_star(i,k,j)=vel_1(i,k,j)+cff2*t_1(i,k,j)
         enddo      
        enddo      
      enddo      

      else

       do j=Nylmin,Nylmax 

        do k=1,2*nwaz+1
 	 do i=1,nwax+1
          ve_star(i,k,j)=vel_1(i,k,j)+cff1*(t_2(i,k,j)+
     +    t_1(i,k,j))
         enddo      
        enddo      
       enddo      

      endif

      return
      end



#endif
