#include "config.h"
#ifdef PRINF
c ______________________________________________________________________
c Pressure Helmholtz eq. for velocity
c ______________________________________________________________________

      subroutine prinf_helmu (p,t,D,E,F,yp,yv,
     +     kxm,kzm,u,v,w,xlh,xuh,reda,redb,redc,redr,Ra )
 
#include "param.h"
      integer*4 l,jj
      real*8 D(nwax+1,2*nwaz+1,Nyl), F(nwax+1,2*nwaz+1,Nyl)
      real*8 E(nwax+1,2*nwaz+1,Nyl), yp(0:Ny+1),yv(0:Ny+1)
      complex*16 u(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 v(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 w(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 p(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 t(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 TNU(nwax+1,2*nwaz+1,Nylmem)
      complex*16 TNV(nwax+1,2*nwaz+1,Nylmem)
      complex*16 TNW(nwax+1,2*nwaz+1,Nylmem)
      complex*16 a,mm,nn

      real*8 kxm(nwax+1),kzm(2*nwaz+1),Ra,cff

#define NX ((nwax+1)*(2*nwaz+1))

      real*8 xlh(NX,Nyl)
      real*8 xuh(NX,Nyl)

#define NNR (2*NPROC-2)

      real*8 reda(NX,0:NNR+1)
      real*8 redb(NX,0:NNR+1)
      real*8 redc(NX,0:NNR+1)
      complex*16 redr(NX,0:NNR+1)

      u(:,:,:)=1e23
      v(:,:,:)=1e23
      w(:,:,:)=1e23
      call TN_prinf_helmu (p,t,yp,yv,kxm,kzm,TNU,TNV,TNW,Ra)
      call tridmpi(D,E,F,TNU,reda,redb,redc,redr,
     &     xlh,xuh,u(1,1,1),
     &     Nyl,NX,NPROC,iam)
      call tridmpi(D,E,F,TNW,reda,redb,redc,redr,
     &     xlh,xuh,w(1,1,1),
     &     Nyl,NX,NPROC,iam)
      call tridmpi(D,E,F,TNV,reda,redb,redc,redr,
     &     xlh,xuh,v(1,1,1),
     &     Nyl,NX,NPROC,iam)
!      print *,'v=',v(1,1,1),v(1,1,2),v(1,1,Ny),v(1,1,Ny-1)
!      print *,'v=',v(5,1,1),v(5,1,2),v(5,1,Ny),v(5,1,Ny-1)
!      print *,v(1,1,1),v(1,1,2),v(1,1,3)
!      print *,v(2,1,1),v(2,1,2),v(2,1,3)
!      print *,v(1,1,Nyl),v(1,1,Nyl-1),v(1,1,Nyl-2)
!      print *,v(2,1,Nyl),v(2,1,Nyl-1),v(2,1,Nyl-2)

!      print *,'test p'
!      print *,p(1,1,0),p(1,1,1),p(1,1,2),p(5,1,0),p(5,1,1),p(5,1,2)
!      print *,p(1,1,Nyl),p(1,1,Nyl-1),p(1,1,Nyl-2),p(5,1,Nyl),
!     &       p(5,1,Nyl-1),p(5,1,Nyl-2)
!      print *,'test u'
!      print *,u(1,1,0),u(1,1,1),u(1,1,2),u(5,1,0),u(5,1,1),u(5,1,2)
!      print *,u(1,1,Nyl+1),u(1,1,Nyl),u(1,1,Nyl-1),u(5,1,Nyl+1),
!     &       u(5,1,Nyl),u(5,1,Nyl-1)
!      print *,'test v'
!!      print *,v(1,1,0),v(1,1,1),v(1,1,2),v(5,1,0),v(5,1,1),v(5,1,2)
!      print *,v(1,1,Nyl+1),v(1,1,Nyl),v(1,1,Nyl-1),v(5,1,Nyl+1),
!     &       v(5,1,Nyl),v(5,1,Nyl-1)
!      stop
      v(1,1,:)=0.0
      if(iam.eq.1) then 
         u(:,:,1)=0.0
         v(:,:,1)=0.0
         w(:,:,1)=0.0
      endif
      if(iam.eq.NPROC) then
         u(:,:,Nyl)=0.0
         v(:,:,Nyl)=0.0
         w(:,:,Nyl)=0.0
      endif

      return
      end

      subroutine prinf_helmv (p,t,D,E,F,yp,yv,
     +     kxm,kzm,u,v,w,xlh,xuh,reda,redb,redc,redr,Ra )
 
#include "param.h"
      integer*4 jj
      real*8 D(nwax+1,2*nwaz+1,Nyl),F(nwax+1,2*nwaz+1,Nyl)
      real*8 E(nwax+1,2*nwaz+1,Nyl),yp(0:Ny+1),yv(0:Ny+1)
      complex*16 u(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 v(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 w(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 p(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 t(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 TNU(nwax+1,2*nwaz+1,Nylmem)
      complex*16 TNV(nwax+1,2*nwaz+1,Nylmem)
      complex*16 TNW(nwax+1,2*nwaz+1,Nylmem)
      complex*16 a,mm,nn

      real*8 kxm(nwax+1),kzm(2*nwaz+1),Ra,cff

#define NX ((nwax+1)*(2*nwaz+1))

      real*8 xlh(NX,Nyl)
      real*8 xuh(NX,Nyl)

#define NNR (2*NPROC-2)

      real*8 reda(NX,0:NNR+1)
      real*8 redb(NX,0:NNR+1)
      real*8 redc(NX,0:NNR+1)
      complex*16 redr(NX,0:NNR+1)

      v(:,:,:)=1e23
      call TN_prinf_helmv(p,t,yp,yv,kxm,kzm,TNU,TNV,TNW,Ra)
      call tridmpi(D,E,F,TNU,reda,redb,redc,redr,
     &     xlh,xuh,u(1,1,1),
     &     Nyl,NX,NPROC,iam)
      call tridmpi(D,E,F,TNW,reda,redb,redc,redr,
     &     xlh,xuh,w(1,1,1),
     &     Nyl,NX,NPROC,iam)
      call tridmpi(D,E,F,TNV,reda,redb,redc,redr,
     &     xlh,xuh,v(1,1,1),
     &     Nyl,NX,NPROC,iam)
      print *,'v=',v(1,1,1),v(1,1,2),v(1,1,Ny),v(1,1,Ny-1)
      print *,'v=',v(5,1,1),v(5,1,2),v(5,1,Ny),v(5,1,Ny-1)
      print *,'u=',u(1,1,1),u(1,1,2),u(1,1,Ny),u(1,1,Ny-1)
      print *,'u=',u(5,1,1),u(5,1,2),u(5,1,Ny),u(5,1,Ny-1)
      v(1,1,:)=0.0
      if(iam.eq.1) then
         u(:,:,1)=0.0
         v(:,:,1)=0.0
         w(:,:,1)=0.0
      endif
      if(iam.eq.NPROC) then
         u(:,:,Nyl)=0.0
         v(:,:,Nyl)=0.0
         w(:,:,Nyl)=0.0
      endif

      return
      end
#endif
