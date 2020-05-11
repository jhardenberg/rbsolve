
      real*8 ptan,qtan 
      real*8 DD,Lx,Lz,Re,invRe,dt,dx,dz 
      real*8 kxmin,kzmin,Ra,Pr,Rt,Rs,Le,Omegay,Omegaz,lapse,qrad

      real*8 alfa(3),beta(3),gamma(3),zeta(3) 
      integer*4 nsave,ndt

      common/ALLPARA/ptan,qtan,DD,Re,Lx,Lz,invRe,dt,dx,dz,
     +   kxmin,kzmin,nsave,ndt,alfa,beta,gamma,zeta,Ra,Pr,
     +   Rt,Rs,Le,Omegay,Omegaz,lapse,qrad

      namelist /numeric/ptan,qtan,dt,ndt,nsave,alfa,beta,gamma,zeta
      namelist /phys/Lx,Lz,DD,Ra,Re,Pr,Rt,Rs,Le,Omegay,Omegaz,lapse,qrad
