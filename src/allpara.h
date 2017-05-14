
      real*8 ptan,qtan 
      real*8 DD,Lx,Lz,Re,invRe,dt,dx,dz 
      real*8 kxmin,kzmin,Ra,Kscalar,Pr,Rt,Rs,Le,Omega

      real*8 alfa(3),beta(3),gamma(3),zeta(3) 
      integer*4 nsave,ttot

      common/ALLPARA/ptan,qtan,DD,Re,Lx,Lz,invRe,dt,dx,dz,
     +   kxmin,kzmin,nsave,ttot,alfa,beta,gamma,zeta,Ra,Kscalar,Pr,
     +   Rt,Rs,Le,Omega
