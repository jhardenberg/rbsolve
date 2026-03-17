#include "config.h"

#ifdef NETCDF_OUTPUT

c **********************************************************************
c  NetCDF4 save/read routines for rbsolve
c
c  - File naming: rb_NNNNNNN.nc  (first itime stored in that file)
c  - Dims: kx (nwax+1), kz (2*nwaz+1), yi (Ny), time (unlimited)
c  - All complex*16 fields stored as real*4 _re / _im pairs
c  - NetCDF4/HDF5 with deflate level 4 compression
c  - ROOT_SAVE only: SAVENODE gathers and writes
c
c  NOTE: param.h declares  i,j,k  as integer*4 and sets implicit none.
c        Do NOT re-declare i,j,k in any subroutine that includes param.h.
c **********************************************************************

c ------------------------------------------------------------------
c  nc_check: standalone error handler (no param.h needed)
c ------------------------------------------------------------------
      subroutine nc_check(ret, context)
      use netcdf
      implicit none
      integer, intent(in)          :: ret
      character(len=*), intent(in) :: context
      if(ret.ne.NF90_NOERR) then
         write(*,*) 'NetCDF error in ',trim(context),': ',
     +              trim(nf90_strerror(ret))
         stop 1
      endif
      return
      end

c ------------------------------------------------------------------
c  nc_defvar: define a NF90_FLOAT var with chunking and compression
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
c  nc_open_or_create
c ------------------------------------------------------------------
      subroutine nc_open_or_create(itime0, ncid,
     +   vid_time,
     +   vid_ure, vid_uim,
     +   vid_vre, vid_vim,
     +   vid_wre, vid_wim,
#ifdef TEMPERATURE
     +   vid_tre, vid_tim,
#endif
#ifdef SALINITY
     +   vid_sre, vid_sim,
#endif
#ifdef SAVE_PRESSURE
     +   vid_pre, vid_pim,
#endif
     +   ncstep)

      use netcdf
#include "param.h"
#include "allpara.h"

      integer itime0, ncid
      integer vid_time
      integer vid_ure, vid_uim
      integer vid_vre, vid_vim
      integer vid_wre, vid_wim
#ifdef TEMPERATURE
      integer vid_tre, vid_tim
#endif
#ifdef SALINITY
      integer vid_sre, vid_sim
#endif
#ifdef SAVE_PRESSURE
      integer vid_pre, vid_pim
#endif
      integer ncstep

      character*20 fname
      integer ret
      integer dim_time, dim_kx, dim_kz, dim_yi
      integer dims4(4)
      integer chunk4(4)
      integer deflate_level, existing_len, did_time
      parameter(deflate_level=4)

      write(fname,'(a3,i7.7,a3)') 'rb_',itime0,'.nc'

      if(iam.ne.SAVENODE) then
         ncstep=1
         return
      endif

c Try to open existing file (for append on restart)
      ret = nf90_open(trim(fname), NF90_WRITE, ncid)
      if(ret.eq.NF90_NOERR) then
         ret = nf90_inq_varid(ncid, 'time',  vid_time)
         ret = nf90_inq_varid(ncid, 'u_re',  vid_ure)
         ret = nf90_inq_varid(ncid, 'u_im',  vid_uim)
         ret = nf90_inq_varid(ncid, 'v_re',  vid_vre)
         ret = nf90_inq_varid(ncid, 'v_im',  vid_vim)
         ret = nf90_inq_varid(ncid, 'w_re',  vid_wre)
         ret = nf90_inq_varid(ncid, 'w_im',  vid_wim)
#ifdef TEMPERATURE
         ret = nf90_inq_varid(ncid, 't_re',  vid_tre)
         ret = nf90_inq_varid(ncid, 't_im',  vid_tim)
#endif
#ifdef SALINITY
         ret = nf90_inq_varid(ncid, 's_re',  vid_sre)
         ret = nf90_inq_varid(ncid, 's_im',  vid_sim)
#endif
#ifdef SAVE_PRESSURE
         ret = nf90_inq_varid(ncid, 'p_re',  vid_pre)
         ret = nf90_inq_varid(ncid, 'p_im',  vid_pim)
#endif
         ret = nf90_inq_dimid(ncid, 'time', did_time)
         call nc_check(ret,'inq_dimid time')
         ret = nf90_inquire_dimension(ncid, did_time,
     +                                len=existing_len)
         call nc_check(ret,'inquire_dim time')
         ncstep = existing_len + 1
      else
c Create new file
         ret = nf90_create(trim(fname),
     +         ior(NF90_NETCDF4,NF90_CLOBBER), ncid)
         call nc_check(ret,'create '//trim(fname))

         ret = nf90_put_att(ncid, NF90_GLOBAL, 'Nx', Nx)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'Ny', Ny)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'Nz', Nz)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'nwax', nwax)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'nwaz', nwaz)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'Ra',   Ra)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'Pr',   Pr)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'itime0', itime0)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'nsave',  nsave)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'nsnc',   nsnc)
         ret = nf90_put_att(ncid, NF90_GLOBAL, 'description',
     +         'rbsolve spectral fields (Fourier space)')

         ret = nf90_def_dim(ncid,'time', NF90_UNLIMITED, dim_time)
         call nc_check(ret,'def_dim time')
         ret = nf90_def_dim(ncid,'yi',   Ny,             dim_yi)
         call nc_check(ret,'def_dim yi')
         ret = nf90_def_dim(ncid,'kz',  2*nwaz+1,       dim_kz)
         call nc_check(ret,'def_dim kz')
         ret = nf90_def_dim(ncid,'kx',   nwax+1,         dim_kx)
         call nc_check(ret,'def_dim kx')

         dims4(1) = dim_kx
         dims4(2) = dim_kz
         dims4(3) = dim_yi
         dims4(4) = dim_time
         chunk4(1) = nwax+1
         chunk4(2) = 2*nwaz+1
         chunk4(3) = 1
         chunk4(4) = 1

         ret = nf90_def_var(ncid,'time', NF90_DOUBLE,
     +                      (/dim_time/), vid_time)
         call nc_check(ret,'def_var time')
         ret = nf90_put_att(ncid, vid_time, 'long_name', 'time')
         ret = nf90_put_att(ncid, vid_time, 'units', 'seconds')

         call nc_defvar(ncid,'u_re','x-velocity real part',
     +    dims4,chunk4,deflate_level,vid_ure)
         call nc_defvar(ncid,'u_im','x-velocity imag part',
     +    dims4,chunk4,deflate_level,vid_uim)
         call nc_defvar(ncid,'v_re','y-velocity real part',
     +    dims4,chunk4,deflate_level,vid_vre)
         call nc_defvar(ncid,'v_im','y-velocity imag part',
     +    dims4,chunk4,deflate_level,vid_vim)
         call nc_defvar(ncid,'w_re','z-velocity real part',
     +    dims4,chunk4,deflate_level,vid_wre)
         call nc_defvar(ncid,'w_im','z-velocity imag part',
     +    dims4,chunk4,deflate_level,vid_wim)
#ifdef TEMPERATURE
         call nc_defvar(ncid,'t_re','temperature real part',
     +    dims4,chunk4,deflate_level,vid_tre)
         call nc_defvar(ncid,'t_im','temperature imag part',
     +    dims4,chunk4,deflate_level,vid_tim)
#endif
#ifdef SALINITY
         call nc_defvar(ncid,'s_re','salinity real part',
     +    dims4,chunk4,deflate_level,vid_sre)
         call nc_defvar(ncid,'s_im','salinity imag part',
     +    dims4,chunk4,deflate_level,vid_sim)
#endif
#ifdef SAVE_PRESSURE
         call nc_defvar(ncid,'p_re','pressure real part',
     +    dims4,chunk4,deflate_level,vid_pre)
         call nc_defvar(ncid,'p_im','pressure imag part',
     +    dims4,chunk4,deflate_level,vid_pim)
#endif

         ret = nf90_enddef(ncid)
         call nc_check(ret,'enddef')
         ncstep = 1
      endif

      return
      end


c ------------------------------------------------------------------
c  nc_save_all: save all fields for one timestep
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

c Persistent state (one file open at a time)
      integer, save :: ncid_sv   = -1
      integer, save :: ncstep_sv =  0
      integer, save :: ifile_t0  = -1
      integer, save :: vid_time_sv
      integer, save :: vid_ure_sv, vid_uim_sv
      integer, save :: vid_vre_sv, vid_vim_sv
      integer, save :: vid_wre_sv, vid_wim_sv
#ifdef TEMPERATURE
      integer, save :: vid_tre_sv, vid_tim_sv
#endif
#ifdef SALINITY
      integer, save :: vid_sre_sv, vid_sim_sv
#endif
#ifdef SAVE_PRESSURE
      integer, save :: vid_pre_sv, vid_pim_sv
#endif

      integer ifile_new, ret
      real*8 rtime_arr(1)

      ifile_new = (itime / nsnc) * nsnc

      if(ncid_sv.lt.0 .or. ifile_new.ne.ifile_t0) then
         if(ncid_sv.ge.0 .and. iam.eq.SAVENODE) then
            ret = nf90_close(ncid_sv)
            call nc_check(ret,'close on rotate')
         endif
         ifile_t0 = ifile_new
         call nc_open_or_create(ifile_t0, ncid_sv,
     +      vid_time_sv,
     +      vid_ure_sv, vid_uim_sv,
     +      vid_vre_sv, vid_vim_sv,
     +      vid_wre_sv, vid_wim_sv,
#ifdef TEMPERATURE
     +      vid_tre_sv, vid_tim_sv,
#endif
#ifdef SALINITY
     +      vid_sre_sv, vid_sim_sv,
#endif
#ifdef SAVE_PRESSURE
     +      vid_pre_sv, vid_pim_sv,
#endif
     +      ncstep_sv)
      endif

      if(iam.eq.1) write(*,'(a,i7,a,i7,a)')
     +   ' NC saving step ',itime,' (record ',ncstep_sv,' in file)'

      call nc_gather_write(uu, ncid_sv, vid_ure_sv, vid_uim_sv,
     +                     ncstep_sv)
      call nc_gather_write(vv, ncid_sv, vid_vre_sv, vid_vim_sv,
     +                     ncstep_sv)
      call nc_gather_write(ww, ncid_sv, vid_wre_sv, vid_wim_sv,
     +                     ncstep_sv)
#ifdef TEMPERATURE
      call nc_gather_write(tt, ncid_sv, vid_tre_sv, vid_tim_sv,
     +                     ncstep_sv)
#endif
#ifdef SALINITY
      call nc_gather_write(ss, ncid_sv, vid_sre_sv, vid_sim_sv,
     +                     ncstep_sv)
#endif
#ifdef SAVE_PRESSURE
      call nc_gather_write(pp, ncid_sv, vid_pre_sv, vid_pim_sv,
     +                     ncstep_sv)
#endif

      if(iam.eq.SAVENODE) then
         rtime_arr(1) = rtime
         ret = nf90_put_var(ncid_sv, vid_time_sv, rtime_arr,
     +                      start=(/ncstep_sv/), count=(/1/))
         call nc_check(ret,'put time')
         ret = nf90_sync(ncid_sv)
         call nc_check(ret,'sync')
      endif

      ncstep_sv = ncstep_sv + 1

      if(iam.eq.SAVENODE) then
         open(10,file='nrec.d',form='formatted',status='unknown')
         write(10,*) itime, rtime
         close(10)
      endif

      return
      end


c ------------------------------------------------------------------
c  nc_gather_write: MPI-gather one field and write to NetCDF
c  (Nyls is in COMMON /MPIGLOBAL/ from param.h - not a dummy arg)
c ------------------------------------------------------------------
      subroutine nc_gather_write(field, ncid, vid_re, vid_im, irec)

      use netcdf
#include "param.h"
cmpi      include 'mpif.h'

      complex*16 field(nwax+1,2*nwaz+1,0:Nylmem+1)
      integer    ncid, vid_re, vid_im, irec

      complex*16 ubuf(nwax+1,2*nwaz+1,Nylmem)
      real*4     slab_re(nwax+1,2*nwaz+1)
      real*4     slab_im(nwax+1,2*nwaz+1)

c  Note: i,j,k already declared by param.h  -- do not re-declare them
      integer nc_start(4), nc_count(4)
      integer cnt, ret, ierr, iproc, ybase_i
cmpi      integer st(MPI_STATUS_SIZE)

      nc_count(1) = nwax+1
      nc_count(2) = 2*nwaz+1
      nc_count(3) = 1
      nc_count(4) = 1

      if(iam.eq.SAVENODE) then
         ybase_i = 0
         do iproc = 1, NPROC
            if(iproc.ne.iam) then
               cnt = (nwax+1)*(2*nwaz+1)*Nyls(iproc)
cmpi               call MPI_RECV(ubuf, cnt, MPI_COMPLEX16,
cmpi     +            iproc-1, 80+iproc, MPI_COMM_WORLD, st, ierr)
            else
               do j=1,Nyls(iproc)
                  ubuf(1:nwax+1, 1:2*nwaz+1, j) =
     +               field(1:nwax+1, 1:2*nwaz+1, j)
               enddo
            endif
            do j=1,Nyls(iproc)
               do i=1,nwax+1
                  slab_re(i,1:2*nwaz+1) = real(ubuf(i,1:2*nwaz+1,j))
                  slab_im(i,1:2*nwaz+1) = aimag(ubuf(i,1:2*nwaz+1,j))
               enddo
               nc_start(1) = 1
               nc_start(2) = 1
               nc_start(3) = ybase_i + j
               nc_start(4) = irec
               ret = nf90_put_var(ncid, vid_re, slab_re,
     +                            start=nc_start, count=nc_count)
               call nc_check(ret,'put_var re')
               ret = nf90_put_var(ncid, vid_im, slab_im,
     +                            start=nc_start, count=nc_count)
               call nc_check(ret,'put_var im')
            enddo
            ybase_i = ybase_i + Nyls(iproc)
         enddo
      else
         cnt = (nwax+1)*(2*nwaz+1)*Nyl
cmpi         call MPI_SEND(field(1,1,1), cnt, MPI_COMPLEX16,
cmpi     +      SAVENODE-1, 80+iam, MPI_COMM_WORLD, st, ierr)
      endif

      return
      end


c ------------------------------------------------------------------
c  nc_read_all: read all fields for restart
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

      integer ncid_r, ncstep_r
      integer vid_time_r
      integer vid_ure_r, vid_uim_r
      integer vid_vre_r, vid_vim_r
      integer vid_wre_r, vid_wim_r
#ifdef TEMPERATURE
      integer vid_tre_r, vid_tim_r
#endif
#ifdef SALINITY
      integer vid_sre_r, vid_sim_r
#endif
      integer dummy1, dummy2

      integer ifile_t0_r, ret, ntimes, irec, ierr, did_time
cmpi      integer st(MPI_STATUS_SIZE)

      ifile_t0_r = (itime / nsnc) * nsnc

      if(iam.eq.1) write(*,*) 'NC reading step ',itime,
     +   ' from file starting at ',ifile_t0_r

      irec = -1

      if(iam.eq.SAVENODE) then
         call nc_open_or_create(ifile_t0_r, ncid_r,
     +      vid_time_r,
     +      vid_ure_r, vid_uim_r,
     +      vid_vre_r, vid_vim_r,
     +      vid_wre_r, vid_wim_r,
#ifdef TEMPERATURE
     +      vid_tre_r, vid_tim_r,
#endif
#ifdef SALINITY
     +      vid_sre_r, vid_sim_r,
#endif
#ifdef SAVE_PRESSURE
     +      dummy1, dummy2,
#endif
     +      ncstep_r)

         ret = nf90_inq_dimid(ncid_r, 'time', did_time)
         call nc_check(ret,'inq_dimid time on read')
         ret = nf90_inquire_dimension(ncid_r, did_time, len=ntimes)
         call nc_check(ret,'inquire_dim time on read')

c Find record for itime: record j holds step ifile_t0_r + (j-1)*nsave
         do j = 1, ntimes
            if(ifile_t0_r + (j-1)*nsave .eq. itime) irec = j
         enddo

         if(irec.lt.0) then
            write(*,*)'nc_read_all: cannot find record for step',itime
            stop 1
         endif
      endif

cmpi      call MPI_BCAST(irec, 1, MPI_INTEGER4,
cmpi     +               SAVENODE-1, MPI_COMM_WORLD, ierr)

      call nc_read_scatter(uu, ncid_r, vid_ure_r, vid_uim_r, irec)
      call nc_read_scatter(vv, ncid_r, vid_vre_r, vid_vim_r, irec)
      call nc_read_scatter(ww, ncid_r, vid_wre_r, vid_wim_r, irec)
#ifdef TEMPERATURE
      call nc_read_scatter(tt, ncid_r, vid_tre_r, vid_tim_r, irec)
#endif
#ifdef SALINITY
      call nc_read_scatter(ss, ncid_r, vid_sre_r, vid_sim_r, irec)
#endif

      if(iam.eq.SAVENODE) then
         ret = nf90_close(ncid_r)
         call nc_check(ret,'close on read')
      endif

      return
      end


c ------------------------------------------------------------------
c  nc_read_scatter: SAVENODE reads one field, distributes to all
c ------------------------------------------------------------------
      subroutine nc_read_scatter(field, ncid, vid_re, vid_im, irec)

      use netcdf
#include "param.h"
cmpi      include 'mpif.h'

      complex*16 field(nwax+1,2*nwaz+1,0:Nylmem+1)
      integer    ncid, vid_re, vid_im, irec

      complex*16 ubuf(nwax+1,2*nwaz+1,Nylmem)
      real*4     slab_re(nwax+1,2*nwaz+1)
      real*4     slab_im(nwax+1,2*nwaz+1)

      integer nc_start(4), nc_count(4)
      integer cnt, ret, ierr, iproc, ybase_i
cmpi      integer st(MPI_STATUS_SIZE)

      nc_count(1) = nwax+1
      nc_count(2) = 2*nwaz+1
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
               ret = nf90_get_var(ncid, vid_re, slab_re,
     +                            start=nc_start, count=nc_count)
               call nc_check(ret,'get_var re')
               ret = nf90_get_var(ncid, vid_im, slab_im,
     +                            start=nc_start, count=nc_count)
               call nc_check(ret,'get_var im')
               do i = 1, nwax+1
                  ubuf(i, 1:2*nwaz+1, j) =
     +               cmplx(slab_re(i, 1:2*nwaz+1),
     +                     slab_im(i, 1:2*nwaz+1), kind=8)
               enddo
            enddo
            if(iproc.ne.iam) then
               cnt = (nwax+1)*(2*nwaz+1)*Nyls(iproc)
cmpi               call MPI_SEND(ubuf(1,1,1), cnt, MPI_COMPLEX16,
cmpi     +            iproc-1, 90+iproc, MPI_COMM_WORLD, st, ierr)
            else
               do j=1,Nyls(iproc)
                  field(1:nwax+1, 1:2*nwaz+1, j) =
     +               ubuf(1:nwax+1, 1:2*nwaz+1, j)
               enddo
            endif
            ybase_i = ybase_i + Nyls(iproc)
         enddo
      else
         cnt = (nwax+1)*(2*nwaz+1)*Nyl
cmpi         call MPI_RECV(field(1,1,1), cnt, MPI_COMPLEX16,
cmpi     +      SAVENODE-1, 90+iam, MPI_COMM_WORLD, st, ierr)
      endif

      return
      end

#endif /* NETCDF_OUTPUT */
