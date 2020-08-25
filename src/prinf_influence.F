#include "config.h"

      subroutine prinf_influence (u,v,w,p,u1,u2,v1,v2,w1,w2,p1,p2,
     &                            yp,yv,mm,nn)
 
#include "param.h"
cmpi        include 'mpif.h'

      complex*16 u(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 v(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 w(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 p(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 u1(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 u2(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 v1(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 v2(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 w1(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 w2(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 p1(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 p2(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 vt(nwax+1,2*nwaz+1)
      complex*16 v1t(nwax+1,2*nwaz+1)
      complex*16 v2t(nwax+1,2*nwaz+1)
      complex*16 vb(nwax+1,2*nwaz+1)
      complex*16 v1b(nwax+1,2*nwaz+1)
      complex*16 v2b(nwax+1,2*nwaz+1)
      complex*16 l1(nwax+1,2*nwaz+1)
      complex*16 l2(nwax+1,2*nwaz+1)
      complex*16 mm(nwax+1),nn(2*nwaz+1),nn1
      real*8 cff,cff1,cff2,yp(0:Ny+1),yv(0:Ny+1)
      integer jj

#ifndef NOMPI
      complex*16 vbuf(nwax+1,2*nwaz+1,3)
      integer sbuf,ierr,st(MPI_STATUS_SIZE)

      sbuf=(nwax+1)*(2*nwaz+1)*3
      if(iam.eq.NPROC) then
         jj=Ny-1
         j=Nyl-1
         cff1=(yv(jj)-yv(jj-1))/( (yv(jj+1)-yv(jj))*(yv(jj+1)-yv(jj-1)))
         cff2=(yv(jj+1)-yv(jj))/( (yv(jj)-yv(jj-1))*(yv(jj+1)-yv(jj-1)))
         vt(:,:)=(v(:,:,j+1)-v(:,:,j))*cff1+
     &           (v(:,:,j)-v(:,:,j-1))*cff2
         v1t(:,:)=(v1(:,:,j+1)-v1(:,:,j))*cff1+
     &            (v1(:,:,j)-v1(:,:,j-1))*cff2
         v2t(:,:)=(v2(:,:,j+1)-v2(:,:,j))*cff1+
     &            (v2(:,:,j)-v2(:,:,j-1))*cff2
         do k=1,2*nwaz+1
           nn1=nn(k)
           do i=1,nwax+1
             vt(i,k)=vt(i,k)+mm(i)*u(i,k,j)+nn1*w(i,k,j)
             v1t(i,k)=v1t(i,k)+mm(i)*u1(i,k,j)+nn1*w1(i,k,j)
             v2t(i,k)=v2t(i,k)+mm(i)*u2(i,k,j)+nn1*w2(i,k,j)
           end do
         end do
         vbuf(:,:,1)=vt(:,:) 
         vbuf(:,:,2)=v1t(:,:) 
         vbuf(:,:,3)=v2t(:,:) 
         call MPI_SEND(vbuf(1,1,1),sbuf,MPI_COMPLEX16,0,96,
     &                 MPI_COMM_WORLD,ierr)
      elseif(iam.eq.1) then
         jj=2
         j=2
         cff1=(yv(jj)-yv(jj-1))/( (yv(jj+1)-yv(jj))*(yv(jj+1)-yv(jj-1)))
         cff2=(yv(jj+1)-yv(jj))/( (yv(jj)-yv(jj-1))*(yv(jj+1)-yv(jj-1)))
         vb(:,:)=(v(:,:,j+1)-v(:,:,j))*cff1+(v(:,:,j)-v(:,:,j-1))*cff2
         v1b(:,:)=(v1(:,:,j+1)-v1(:,:,j))*cff1+
     &            (v1(:,:,j)-v1(:,:,j-1))*cff2
         v2b(:,:)=(v2(:,:,j+1)-v2(:,:,j))*cff1+
     &            (v2(:,:,j)-v2(:,:,j-1))*cff2
         do k=1,2*nwaz+1
           nn1=nn(k)
           do i=1,nwax+1
             vb(i,k)=vb(i,k)+mm(i)*u(i,k,j)+nn1*w(i,k,j)
             v1b(i,k)=v1b(i,k)+mm(i)*u1(i,k,j)+nn1*w1(i,k,j)
             v2b(i,k)=v2b(i,k)+mm(i)*u2(i,k,j)+nn1*w2(i,k,j)
           end do
         end do
         call MPI_RECV(vbuf(1,1,1),sbuf,MPI_COMPLEX16,NPROC-1,96,
     &                 MPI_COMM_WORLD,st,ierr)
         vt(:,:)=vbuf(:,:,1)
         v1t(:,:)=vbuf(:,:,2)
         v2t(:,:)=vbuf(:,:,3)

         l1(:,:)=(vb(:,:)*v2t(:,:)-vt(:,:)*v2b(:,:))/
     &        (v1t(:,:)*v2b(:,:)-v1b(:,:)*v2t(:,:))
         l2(:,:)=(vt(:,:)*v1b(:,:)-vb(:,:)*v1t(:,:))/
     &        (v1t(:,:)*v2b(:,:)-v1b(:,:)*v2t(:,:))
         l1(1,1)=0.d0
         l2(1,1)=0.d0
      endif
      sbuf=(nwax+1)*(2*nwaz+1)
      call MPI_BCAST(l1(1,1),sbuf,MPI_COMPLEX16,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(l2(1,1),sbuf,MPI_COMPLEX16,0,MPI_COMM_WORLD,ierr)
#else
      jj=Ny-1
      j=Ny-1
      cff1=(yv(jj)-yv(jj-1))/( (yv(jj+1)-yv(jj))*(yv(jj+1)-yv(jj-1)))
      cff2=(yv(jj+1)-yv(jj))/( (yv(jj)-yv(jj-1))*(yv(jj+1)-yv(jj-1)))
      vt(:,:)=(v(:,:,j+1)-v(:,:,j))*cff1+(v(:,:,j)-v(:,:,j-1))*cff2
      v1t(:,:)=(v1(:,:,j+1)-v1(:,:,j))*cff1+(v1(:,:,j)-v1(:,:,j-1))*cff2
      v2t(:,:)=(v2(:,:,j+1)-v2(:,:,j))*cff1+(v2(:,:,j)-v2(:,:,j-1))*cff2
      do k=1,2*nwaz+1
       nn1=nn(k)
       do i=1,nwax+1
        vt(i,k)=vt(i,k)+mm(i)*u(i,k,j)+nn1*w(i,k,j)
        v1t(i,k)=v1t(i,k)+mm(i)*u1(i,k,j)+nn1*w1(i,k,j)
        v2t(i,k)=v2t(i,k)+mm(i)*u2(i,k,j)+nn1*w2(i,k,j)
       end do
      end do

      jj=2
      j=2
      cff1=(yv(jj)-yv(jj-1))/( (yv(jj+1)-yv(jj))*(yv(jj+1)-yv(jj-1)))
      cff2=(yv(jj+1)-yv(jj))/( (yv(jj)-yv(jj-1))*(yv(jj+1)-yv(jj-1)))
      vb(:,:)=(v(:,:,j+1)-v(:,:,j))*cff1+(v(:,:,j)-v(:,:,j-1))*cff2
      v1b(:,:)=(v1(:,:,j+1)-v1(:,:,j))*cff1+(v1(:,:,j)-v1(:,:,j-1))*cff2
      v2b(:,:)=(v2(:,:,j+1)-v2(:,:,j))*cff1+(v2(:,:,j)-v2(:,:,j-1))*cff2
      do k=1,2*nwaz+1
       nn1=nn(k)
       do i=1,nwax+1
        vb(i,k)=vb(i,k)+mm(i)*u(i,k,j)+nn1*w(i,k,j)
        v1b(i,k)=v1b(i,k)+mm(i)*u1(i,k,j)+nn1*w1(i,k,j)
        v2b(i,k)=v2b(i,k)+mm(i)*u2(i,k,j)+nn1*w2(i,k,j)
       end do
      end do

      l1(:,:)=(vb(:,:)*v2t(:,:)-vt(:,:)*v2b(:,:))/
     &        (v1t(:,:)*v2b(:,:)-v1b(:,:)*v2t(:,:))
      l2(:,:)=(vt(:,:)*v1b(:,:)-vb(:,:)*v1t(:,:))/
     &        (v1t(:,:)*v2b(:,:)-v1b(:,:)*v2t(:,:))
      l1(1,1)=0.d0
      l2(1,1)=0.d0
#endif
      do j=Nylmin,Nylmax
        u(:,:,j)=u(:,:,j)+l1(:,:)*u1(:,:,j)+l2(:,:)*u2(:,:,j)
        v(:,:,j)=v(:,:,j)+l1(:,:)*v1(:,:,j)+l2(:,:)*v2(:,:,j)
        w(:,:,j)=w(:,:,j)+l1(:,:)*w1(:,:,j)+l2(:,:)*w2(:,:,j)
      enddo
      do j=1,Nylmax1
        p(:,:,j)=p(:,:,j)+l1(:,:)*p1(:,:,j)+l2(:,:)*p2(:,:,j)
      enddo

      return
      end

      subroutine prinf_checkdiv (u,v,w,uv,yp,yv,mm,nn)
 
#include "param.h"
      complex*16 u(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 v(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 w(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 uv(nwax+1,2*nwaz+1,0:Nyl+1)
      complex*16 mm(nwax+1),nn(2*nwaz+1),nn1
      real*8 cff,cff1,cff2,yp(0:Ny+1),yv(0:Ny+1)
      integer jj

      do j=2,Nyl-2
       jj=j+Nybase
       cff1=(yv(jj)-yv(jj-1))/( (yv(jj+1)-yv(jj))*(yv(jj+1)-yv(jj-1)))
       cff2=(yv(jj+1)-yv(jj))/( (yv(jj)-yv(jj-1))*(yv(jj+1)-yv(jj-1)))
       do k=1,2*nwaz+1
        nn1=nn(k)
        do i=1,nwax+1
         uv(i,k,j)=mm(i)*u(i,k,j)+nn1*w(i,k,j)+
     +    (v(i,k,j+1)-v(i,k,j))*cff1+(v(i,k,j)-v(i,k,j-1))*cff2
        end do
       end do
      end do
      return
      end

