      implicit none
      
      integer*4 kx,kz
      integer*4 Nx,Ny,Nz
      integer*4 i,j,k

#ifdef HTROPO
	integer*4 ntropo

	parameter (Ntropo=129)
#endif
	
!       parameter (Nx=384,Ny=385,Nz=384)
!       parameter (kx=336	,kz=336)

!       parameter (Nx=256,Ny=192,Nz=256)
!       parameter (kx=170,kz=170)
!       parameter (kx=204,kz=204)

!       parameter (Nx=384,Ny=257,Nz=384)
!       parameter (kx=256,kz=256)

!       parameter (Nx=64,Ny=129,Nz=64)
!       parameter (kx=56	,kz=56)
!       parameter (Nx=64,Ny=97,Nz=64)
!       parameter (kx=56,kz=56)
!       parameter (Nx=128,Ny=97,Nz=128)
!       parameter (kx=102,kz=102)

!       parameter (Nx=96,Ny=129,Nz=64)
!       parameter (kx=84	,kz=56)
!       parameter (Nx=96,Ny=65,Nz=96)
!       parameter (kx=84	,kz=84)
!       parameter (Nx=64,Ny=97,Nz=64)
!       parameter (kx=56,kz=56)
!       parameter (Nx=16,Ny=97,Nz=16)
!       parameter (kx=14,kz=14)

!       parameter (Nx=256,Ny=3072,Nz=256)
!       parameter (kx=204,kz=204)
!       parameter (Nx=256,Ny=192,Nz=256)
!       parameter (kx=204,kz=204)
!       parameter (Nx=768,Ny=5120,Nz=768)

        parameter (Nx=192,Ny=129,Nz=1)
        parameter (kx=154,kz=0)

!       parameter (Nx=192,Ny=257,Nz=192)
!       parameter (kx=154,kz=154)

!       parameter (Nx=256,Ny=1025,Nz=256)
!       parameter (kx=204	,kz=204)

!       parameter (Nx=192,Ny=129,Nz=2)
!       parameter (kx=154,kz=2)

!       parameter (Nx=64,Ny=129,Nz=2)
!       parameter (kx=51,kz=2)

!       parameter (Nx=128,Ny=97,Nz=128)
!       parameter (kx=102,kz=102)

!       parameter (Nx=192,Ny=257,Nz=192)
!       parameter (kx=154,kz=154)

!       parameter (Nx=512,Ny=129,Nz=512)
!       parameter (kx=410,kz=410)

!       parameter (Nx=768,Ny=65,Nz=768)
!       parameter (kx=614,kz=614)
!       parameter (Nx=192,Ny=192,Nz=192)     
!       parameter (kx=168 ,kz=168)   
!
!       parameter (Nx=256,Ny=68,Nz=256)
!       parameter (kx=224	,kz=224)
!       parameter (Nx=512,Ny=50,Nz=512)
!       parameter (kx=448	,kz=448)

!       parameter (Nx=256,Ny=65,Nz=256)
!       parameter (kx=224,kz=224)
	
!       parameter (Nx=768,Ny=64,Nz=192)
!       parameter (kx=614	,kz=154)
!       parameter (Nx=768,Ny=5120,Nz=192)
!       parameter (Nx=768,Ny=5120,Nz=192)
!       parameter (kx=614,kz=154)

!c kx,kz must be smaller than Nx,Nz and even
      integer*4 nwax,nwaz
!      parameter (nwax=kx/2, nwaz=kz/2)
      parameter (nwax=kx/2, nwaz=kz/2)
	 
      integer Nyl,NPROC,iam,Nylmin,Nylmax,Nybase,Nylmem,Nylmin0,Nylmax0
      integer Nylmin1,Nylmax1
!c	Nylmin1 and Nylmax1 are 1-Nyl for all nodes except the first and
!c	last where they are 2 and Ny-1

      integer Nyls(MAXCPU)
      integer myMPI_VEL_TYPE,myMPI_ALL_TYPE
      integer*8 f_plan,  i_plan


      common /MPIGLOBAL/ f_plan,i_plan, Nyl,iam,Nylmin,
     &          Nylmax,Nybase,Nylmin0,Nylmax0,
     &          Nylmax1,Nylmin1,
     &          Nyls,myMPI_ALL_TYPE,myMPI_VEL_TYPE

!* This is the maximum  slice thickness 
#ifdef NOMPI
      parameter (Nylmem=Ny)
      parameter  (NPROC=1)
#else
      parameter (Nylmem=33)
      parameter  (NPROC=4)
#endif

