#include "config.h"

#ifdef NETCDF_OUTPUT

c **********************************************************************
c  NetCDF4 save/read routines for rbsolve  --  PHYSICAL SPACE version
c
c  Fields are transformed to physical space (inverse FFT via Fourier1)
c  before writing. NetCDF dimensions: x (Nx), z (Nz), y (Ny), time.
c  Variables are single real (NF90_FLOAT), no _re/_im split.
c
c  File naming: rb_NNNNNNN.nc  (first itime in that file)
c  Compression: NetCDF4/HDF5, deflate level 4, chunked per (x,z) slab
c **********************************************************************

c ------------------------------------------------------------------
c  nc_check: error handler (no param.h / no implicit typing needed)
c ------------------------------------------------------------------
      subroutine nc_check(ret, context)
      use netcdf
      implicit none
      integer,          intent(in) :: ret
      character(len=*), intent(in) :: context
      if(ret.ne.NF90_NOERR) then
         write(*,*) 'NetCDF error in ',trim(context),': ',
     +              trim(nf90_strerror(ret))
         stop 1
      endif
      return
      end

c ------------------------------------------------------------------
c  nc_defvar: define a NF90_FLOAT variable with chunk+deflate
c ------------------------------------------------------------------
      subroutine nc_defvar(ncid, vname, lname, dims4, chunk4,
     +                     dlev, vid)
      use netcdf
      implicit none
      integer,          intent(in)  :: ncid
      integer,          intent(in)  :: dims4(4), chunk4(4), dlev
      character(len=*), intent(in)  :: vname, lname
      integer,          intent(out) :: vid
      integer                       :: ret

      ret = nf90_def_var(ncid, vname, NF90_FLOAT, dims4, vid)
      call nc_check(ret, 'def_var '//trim(vname))
      ret = nf90_def_var_deflate(ncid, vid,
     +      shuffle=1, deflate=1, deflate_level=dlev)
      ret = nf90_def_var_chunking(ncid, vid, NF90_CHUNKED, chunk4)
      ret = nf90_put_att(ncid, vid, 'long_name', trim(lname))
      return
      end

c ------------------------------------------------------------------
c  nc_open_or_create: open (append) or create a NetCDF4 file.
c  Only SAVENODE does real work; others just get ncstep=1.
c ------------------------------------------------------------------
      subroutine nc_open_or_create(itime0, ncid,
     +   vid_time,
     +   vid_u, vid_v, vid_w,
#ifdef TEMPERATURE
     +   vid_t,
#endif
#ifdef SALINITY
     +   vid_s,
#endif
#ifdef SAVE_PRESSURE
     +   vid_p,
#endif
     +   ncstep)

      use netcdf
#include "param.h"
#include "allpara.h"

      integer itime0, ncid
      integer vid_time
      integer vid_u, vid_v, vid_w
#ifdef TEMPERATURE
      integer vid_t
#endif
#ifdef SALINITY
      integer vid_s
#endif
#ifdef SAVE_PRESSURE
      integer vid_p
#endif
      integer ncstep

      character*20 fname
      integer ret
      integer dim_time, dim_x, dim_z, dim_y
      integer dims4(4), chunk4(4)
      integer deflate_level, existing_len, did_time
      parameter(deflate_level=4)

      write(fname,'(a3,i7.7,a3)') 'rb_',itime0,'.nc'

      if(iam.ne.SAVENODE) then
         ncstep = 1
         return
      endif

c Try to open existing file for append
      ret = nf90_open(trim(fname), NF90_WRITE, ncid)
      if(ret.eq.NF90_NOERR) then
         ret = nf90_inq_varid(ncid, 'time', vid_time)
         ret = nf90_inq_varid(ncid, 'u',    vid_u)
         ret = nf90_inq_varid(ncid, 'v',    vid_v)
         ret = nf90_inq_varid(ncid, 'w',    vid_w)
#ifdef TEMPERATURE
         ret = nf90_inq_varid(ncid, 't',    vid_t)
#endif
#ifdef SALINITY
         ret = nf90_inq_varid(ncid, 's',    vid_s)
#endif
#ifdef SAVE_PRESSURE
         ret = nf90_inq_varid(ncid, 'p',    vid_p)
#endif
         ret = nf90_inq_dimid(ncid, 'time', did_time)
         call nc_check(ret,'inq_dimid time')
         ret = nf90_inquire_dimension(ncid, did_time,
     +                                len=existing_len)
         call nc_check(ret,'inquire_dim time')
         ncstep = existing_len + 1
      else
c Create new file (NetCDF4/HDF5)
         ret = nf90_create(trim(fname),
     +         ior(NF90_NETCDF4,NF90_CLOBBER), ncid)
         call nc_check(ret,'create '//trim(fname))

c Global attributes
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'Nx', Nx)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'Ny', Ny)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'Nz', Nz)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'Ra',   Ra)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'Pr',   Pr)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'Lx',   Lx)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'Lz',   Lz)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'itime0', itime0)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'nsave',  nsave)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'nsnc',   nsnc)
         ret = nf90_put_att(ncid, NF90_GLOBAL,
     +      'description', 'rbsolve fields in physical space')

c Define dimensions  (Fortran col-major: x runs fastest in memory)
         ret = nf90_def_dim(ncid,'time', NF90_UNLIMITED, dim_time)
         call nc_check(ret,'def_dim time')
         ret = nf90_def_dim(ncid,'y',  Ny,  dim_y)
         call nc_check(ret,'def_dim y')
         net_z: if(Nz.gt.1) then
            ret = nf90_def_dim(ncid,'z',  Nz,  dim_z)
            call nc_check(ret,'def_dim z')
         else
            ret = nf90_def_dim(ncid,'z',  1,   dim_z)
            call nc_check(ret,'def_dim z (scalar)')
         end if net_z
         ret = nf90_def_dim(ncid,'x',  Nx,  dim_x)
         call nc_check(ret,'def_dim x')

c Dim/chunk order: x,z,y,time  (Fortran storage order)
         dims4(1) = dim_x
         dims4(2) = dim_z
         dims4(3) = dim_y
         dims4(4) = dim_time
         chunk4(1) = Nx
         chunk4(2) = max(1,Nz)
         chunk4(3) = 1
         chunk4(4) = 1

c Time coordinate
         ret = nf90_def_var(ncid, 'time', NF90_DOUBLE,
     +                      (/dim_time/), vid_time)
         call nc_check(ret,'def_var time')
         ret = nf90_put_att(ncid, vid_time, 'long_name','time')
         ret = nf90_put_att(ncid, vid_time, 'units','seconds')

c Field variables
         call nc_defvar(ncid,'u','x-velocity',
     +      dims4,chunk4,deflate_level,vid_u)
         call nc_defvar(ncid,'v','y-velocity',
     +      dims4,chunk4,deflate_level,vid_v)
         call nc_defvar(ncid,'w','z-velocity',
     +      dims4,chunk4,deflate_level,vid_w)
#ifdef TEMPERATURE
         call nc_defvar(ncid,'t','temperature',
     +      dims4,chunk4,deflate_level,vid_t)
#endif
#ifdef SALINITY
         call nc_defvar(ncid,'s','salinity',
     +      dims4,chunk4,deflate_level,vid_s)
#endif
#ifdef SAVE_PRESSURE
         call nc_defvar(ncid,'p','pressure',
     +      dims4,chunk4,deflate_level,vid_p)
#endif
         ret = nf90_enddef(ncid)
         call nc_check(ret,'enddef')
         ncstep = 1
      endif

      return
      end


c ------------------------------------------------------------------
c  nc_save_all: transform fields to physical space and save.
c
c  Called from main.F when mod(itime,nsave)==0.
c  Applies Fourier1 (inverse FFT) to each spectral field first.
c  Then gathers each y-slab from all MPI ranks to SAVENODE and writes.
c ------------------------------------------------------------------
      subroutine nc_save_all(
     +   uu, vv, ww,
#ifdef TEMPERATURE
     +   tt,
#endif
#ifdef SALINITY
     +   ss,
#endif
#ifdef SAVE_PRESSURE
     +   pp,
#endif
     +   itime, rtime)

      use netcdf
#include "param.h"
#include "allpara.h"
cmpi      include 'mpif.h'

      complex*16 uu(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 vv(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 ww(nwax+1,2*nwaz+1,0:Nylmem+1)
#ifdef TEMPERATURE
      complex*16 tt(nwax+1,2*nwaz+1,0:Nylmem+1)
#endif
#ifdef SALINITY
      complex*16 ss(nwax+1,2*nwaz+1,0:Nylmem+1)
#endif
#ifdef SAVE_PRESSURE
      complex*16 pp(nwax+1,2*nwaz+1,0:Nylmem+1)
#endif
      integer*4 itime
      real*8    rtime

c Physical-space work buffer (inverse FFT output)
      real*8    uo(Nx,Nz,0:Nylmem+1)

c Persistent file state
      integer, save :: ncid_sv   = -1
      integer, save :: ncstep_sv =  0
      integer, save :: ifile_t0  = -1
      integer, save :: vid_time_sv
      integer, save :: vid_u_sv, vid_v_sv, vid_w_sv
#ifdef TEMPERATURE
      integer, save :: vid_t_sv
#endif
#ifdef SALINITY
      integer, save :: vid_s_sv
#endif
#ifdef SAVE_PRESSURE
      integer, save :: vid_p_sv
#endif

      integer ifile_new, ret
      real*8  rtime_arr(1)

      ifile_new = (itime / nsnc) * nsnc

c Rotate to a new file if needed
      if(ncid_sv.lt.0 .or. ifile_new.ne.ifile_t0) then
         if(ncid_sv.ge.0 .and. iam.eq.SAVENODE) then
            ret = nf90_close(ncid_sv)
            call nc_check(ret,'close on rotate')
         endif
         ifile_t0 = ifile_new
         call nc_open_or_create(ifile_t0, ncid_sv,
     +      vid_time_sv,
     +      vid_u_sv, vid_v_sv, vid_w_sv,
#ifdef TEMPERATURE
     +      vid_t_sv,
#endif
#ifdef SALINITY
     +      vid_s_sv,
#endif
#ifdef SAVE_PRESSURE
     +      vid_p_sv,
#endif
     +      ncstep_sv)
      endif

      if(iam.eq.1) write(*,'(a,i7,a,i7,a)')
     +   ' NC saving step ',itime,' (record ',ncstep_sv,' in file)'

c Transform and write each field
      call Fourier1(uu, uo)
      call nc_gather_write(uo, ncid_sv, vid_u_sv, ncstep_sv)

      call Fourier1(vv, uo)
      call nc_gather_write(uo, ncid_sv, vid_v_sv, ncstep_sv)

      call Fourier1(ww, uo)
      call nc_gather_write(uo, ncid_sv, vid_w_sv, ncstep_sv)

#ifdef TEMPERATURE
      call Fourier1(tt, uo)
      call nc_gather_write(uo, ncid_sv, vid_t_sv, ncstep_sv)
#endif
#ifdef SALINITY
      call Fourier1(ss, uo)
      call nc_gather_write(uo, ncid_sv, vid_s_sv, ncstep_sv)
#endif
#ifdef SAVE_PRESSURE
      call Fourier1(pp, uo)
      call nc_gather_write(uo, ncid_sv, vid_p_sv, ncstep_sv)
#endif

c Write time coordinate and flush
      if(iam.eq.SAVENODE) then
         rtime_arr(1) = rtime
         ret = nf90_put_var(ncid_sv, vid_time_sv, rtime_arr,
     +                      start=(/ncstep_sv/), count=(/1/))
         call nc_check(ret,'put time')
         ret = nf90_sync(ncid_sv)
         call nc_check(ret,'sync')
      endif

      ncstep_sv = ncstep_sv + 1

c Update nrec.d
      if(iam.eq.SAVENODE) then
         open(10,file='nrec.d',form='formatted',status='unknown')
         write(10,*) itime, rtime
         close(10)
      endif

      return
      end


c ------------------------------------------------------------------
c  nc_gather_write: gather physical field to SAVENODE, write to NC.
c  field  : real*8(Nx,Nz,0:Nylmem+1), local y-slice for this rank
c  Nyls is in COMMON /MPIGLOBAL/ via param.h -- not re-declared here.
c ------------------------------------------------------------------
      subroutine nc_gather_write(field, ncid, vid, irec)

      use netcdf
#include "param.h"
cmpi      include 'mpif.h'

      real*8     field(Nx,Nz,0:Nylmem+1)
      integer    ncid, vid, irec

      real*8     ubuf(Nx,Nz,Nylmem)
      real*4     slab(Nx,max(1,Nz))

      integer nc_start(4), nc_count(4)
      integer cnt, ret, ierr, iproc, ybase_i
cmpi      integer st(MPI_STATUS_SIZE)

      nc_count(1) = Nx
      nc_count(2) = max(1,Nz)
      nc_count(3) = 1
      nc_count(4) = 1

      if(iam.eq.SAVENODE) then
         ybase_i = 0
         do iproc = 1, NPROC
            if(iproc.ne.iam) then
               cnt = Nx*max(1,Nz)*Nyls(iproc)
cmpi               call MPI_RECV(ubuf, cnt, MPI_DOUBLE_PRECISION,
cmpi     +            iproc-1, 80+iproc, MPI_COMM_WORLD, st, ierr)
            else
               do j=1,Nyls(iproc)
                  ubuf(1:Nx,1:max(1,Nz),j) =
     +               field(1:Nx,1:max(1,Nz),j)
               enddo
            endif
            do j=1,Nyls(iproc)
               do i=1,Nx
                  slab(i,1:max(1,Nz)) = real(ubuf(i,1:max(1,Nz),j))
               enddo
               nc_start(1) = 1
               nc_start(2) = 1
               nc_start(3) = ybase_i + j
               nc_start(4) = irec
               ret = nf90_put_var(ncid, vid, slab,
     +                            start=nc_start, count=nc_count)
               call nc_check(ret,'put_var slab')
            enddo
            ybase_i = ybase_i + Nyls(iproc)
         enddo
      else
         cnt = Nx*max(1,Nz)*Nyl
cmpi         call MPI_SEND(field(1,1,1), cnt, MPI_DOUBLE_PRECISION,
cmpi     +      SAVENODE-1, 80+iam, MPI_COMM_WORLD, st, ierr)
      endif

      return
      end


c ------------------------------------------------------------------
c  nc_read_all: read physical fields back from NetCDF for restart,
c  then apply forward FFT (Fourier2) to get spectral fields.
c
c  IMPORTANT: the caller must call initfour() before nc_read_all.
c  Fourier2 (r2c) is the forward transform that converts
c  physical real*8(Nx,Nz,y) -> spectral complex*16(nwax+1,2*nwaz+1,y).
c ------------------------------------------------------------------
      subroutine nc_read_all(
     +   uu, vv, ww,
#ifdef TEMPERATURE
     +   tt,
#endif
#ifdef SALINITY
     +   ss,
#endif
     +   itime)

      use netcdf
#include "param.h"
#include "allpara.h"
cmpi      include 'mpif.h'

      complex*16 uu(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 vv(nwax+1,2*nwaz+1,0:Nylmem+1)
      complex*16 ww(nwax+1,2*nwaz+1,0:Nylmem+1)
#ifdef TEMPERATURE
      complex*16 tt(nwax+1,2*nwaz+1,0:Nylmem+1)
#endif
#ifdef SALINITY
      complex*16 ss(nwax+1,2*nwaz+1,0:Nylmem+1)
#endif
      integer*4 itime

c Physical work buffer
      real*8    uo(Nx,Nz,0:Nylmem+1)
c Ones buffer so Fourier2(uo,ones,field,1,Ny) = FFT(uo*1) = FFT(uo)
      real*8    ones(Nx,Nz,0:Nylmem+1)

      integer ncid_r, ncstep_r
      integer vid_time_r
      integer vid_u_r, vid_v_r, vid_w_r
#ifdef TEMPERATURE
      integer vid_t_r
#endif
#ifdef SALINITY
      integer vid_s_r
#endif
      integer dummy1

      integer ifile_t0_r, ret, ntimes, irec, ierr, did_time
cmpi      integer st(MPI_STATUS_SIZE)

      ifile_t0_r = (itime / nsnc) * nsnc

      if(iam.eq.1) write(*,*) 'NC reading step ',itime,
     +   ' from file starting at ',ifile_t0_r

      irec = -1

      if(iam.eq.SAVENODE) then
         call nc_open_or_create(ifile_t0_r, ncid_r,
     +      vid_time_r,
     +      vid_u_r, vid_v_r, vid_w_r,
#ifdef TEMPERATURE
     +      vid_t_r,
#endif
#ifdef SALINITY
     +      vid_s_r,
#endif
#ifdef SAVE_PRESSURE
     +      dummy1,
#endif
     +      ncstep_r)

         ret = nf90_inq_dimid(ncid_r, 'time', did_time)
         call nc_check(ret,'inq_dimid time on read')
         ret = nf90_inquire_dimension(ncid_r, did_time, len=ntimes)
         call nc_check(ret,'inquire_dim time on read')

c Record j corresponds to step: ifile_t0_r + (j-1)*nsave
         do j = 1, ntimes
            if(ifile_t0_r + (j-1)*nsave .eq. itime) irec = j
         enddo

         if(irec.lt.0) then
            write(*,*)'nc_read_all: cannot find record for step',itime
            stop 1
         endif
      endif

c Initialise ones buffer (used to forward-FFT a single field)
      do j=0,Nylmem+1
        ones(1:Nx,1:Nz,j) = 1.d0
      enddo

cmpi      call MPI_BCAST(irec, 1, MPI_INTEGER4,
cmpi     +               SAVENODE-1, MPI_COMM_WORLD, ierr)

c Read physical slabs then forward-FFT to spectral
c  Fourier2(uo, ones, field, 1, Ny) computes FFT(uo*1) = FFT(uo)
      call nc_read_scatter(uo, ncid_r, vid_u_r, irec)
      call Fourier2(uo, ones, uu, 1, Ny)

      call nc_read_scatter(uo, ncid_r, vid_v_r, irec)
      call Fourier2(uo, ones, vv, 1, Ny)

      call nc_read_scatter(uo, ncid_r, vid_w_r, irec)
      call Fourier2(uo, ones, ww, 1, Ny)

#ifdef TEMPERATURE
      call nc_read_scatter(uo, ncid_r, vid_t_r, irec)
      call Fourier2(uo, ones, tt, 1, Ny)
#endif
#ifdef SALINITY
      call nc_read_scatter(uo, ncid_r, vid_s_r, irec)
      call Fourier2(uo, ones, ss, 1, Ny)
#endif

      if(iam.eq.SAVENODE) then
         ret = nf90_close(ncid_r)
         call nc_check(ret,'close on read')
      endif

      return
      end


c ------------------------------------------------------------------
c  nc_read_scatter: SAVENODE reads physical slabs -> distributes
c  to all MPI ranks as real*8(Nx,Nz,0:Nylmem+1).
c ------------------------------------------------------------------
      subroutine nc_read_scatter(field, ncid, vid, irec)

      use netcdf
#include "param.h"
cmpi      include 'mpif.h'

      real*8  field(Nx,Nz,0:Nylmem+1)
      integer ncid, vid, irec

      real*8  ubuf(Nx,Nz,Nylmem)
      real*4  slab(Nx,max(1,Nz))

      integer nc_start(4), nc_count(4)
      integer cnt, ret, ierr, iproc, ybase_i
cmpi      integer st(MPI_STATUS_SIZE)

      nc_count(1) = Nx
      nc_count(2) = max(1,Nz)
      nc_count(3) = 1
      nc_count(4) = 1

      if(iam.eq.SAVENODE) then
         ybase_i = 0
         do iproc = 1, NPROC
            do j = 1, Nyls(iproc)
               nc_start(1) = 1
               nc_start(2) = 1
               nc_start(3) = ybase_i + j
               nc_start(4) = irec
               ret = nf90_get_var(ncid, vid, slab,
     +                            start=nc_start, count=nc_count)
               call nc_check(ret,'get_var slab')
               do i = 1, Nx
                  ubuf(i, 1:max(1,Nz), j) =
     +               dble(slab(i, 1:max(1,Nz)))
               enddo
            enddo
            if(iproc.ne.iam) then
               cnt = Nx*max(1,Nz)*Nyls(iproc)
cmpi               call MPI_SEND(ubuf(1,1,1), cnt, MPI_DOUBLE_PRECISION,
cmpi     +            iproc-1, 90+iproc, MPI_COMM_WORLD, st, ierr)
            else
               do j=1,Nyls(iproc)
                  field(1:Nx,1:max(1,Nz),j) =
     +               ubuf(1:Nx,1:max(1,Nz),j)
               enddo
            endif
            ybase_i = ybase_i + Nyls(iproc)
         enddo
      else
         cnt = Nx*max(1,Nz)*Nyl
cmpi         call MPI_RECV(field(1,1,1), cnt, MPI_DOUBLE_PRECISION,
cmpi     +      SAVENODE-1, 90+iam, MPI_COMM_WORLD, st, ierr)
      endif

      return
      end

#endif /* NETCDF_OUTPUT */
