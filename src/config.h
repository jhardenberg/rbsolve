#define ONLY2D
!* This version of the code has an active scalar S */
!* if used for salt fingers define Rt,Rs,Le,Pr (Ra is not used) */
!**********************************
!* User defined options            *
!* Comment out what does not apply *
!**********************************/
!**************************
!*    NUMERICS            *
!**************************

!-#define OLDPUNTI

!**************************
!*   PHYSICS              *
!**************************

!* Convective code *
#define TEMPERATURE

!-#define HTROPO (1.d0)
!* problem with lapse rate

!* Use also an additional field S *
!-#define SALINITY
!-#define FINGER

!* Lagrangian tracers *
!-#define LAGRANGIAN

!* The mean pressure gradient *
#define PRESSURE_GRADIENT 0.d0

!***** BOUNDARIES ******

!* Velocities: NO SLIP is the default, 
!*  uncomment these for free slip instead   
#define FREE_SLIP_TOP
#define FREE_SLIP_BOTTOM

!* Dirichlet boundary conditions. Comment to leave free *

#define TEMPERATURE_BOTTOM 0.5d0
#define TEMPERATURE_TOP -0.5d0

!* Temperature Flux boundary conditions. *

!-#define FLUXT_BOTTOM (Re*Pr) 
!-#define FLUXT_BOTTOM (1.d0)
!-#define NOFLUXT_BOTTOM

!-#define FLUXT_TOP (1.d0)
!-#define FLUXT_TOP (0.d0)
!-#define NOFLUXT_TOP

#define SCALAR_BOTTOM -0.5d0
#define SCALAR_TOP 0.5d0

!-#define FLUXS_BOTTOM (1.d0)
!-#define NOFLUXS_BOTTOM
!-#define FLUXS_TOP (1.d0)
!-#define NOFLUXS_TOP

!* Non-homogeneous Temperature boundary conditions (Dirichlet or flux)
!-#define PATTERNT_BOTTOM
!-#define PATTERNT_TOP

!**************************
!* Technical code options *
!*************************

!* Do you want a scalar version w/out MPI? uncomment this!  *
! Actually not needed anymore ... ignore #define NOMPI

!* Use the FFTW libraries, else fall back to the Temperton FFTs.
!*  If you have these libraries use them! If not, consider installing 
!*  them! (http://www.fftw.org)  *
#define HAVE_FFTW

!* Save only on one node (SAVENODE). If not, every node saves separately *
#define ROOT_SAVE
#define SAVENODE (1)

!* FORMATTED or UNFORMATTED output (.dat and .unf files respectively) *
!-#define FORMATTED_OUTPUT

!* Let the MPI code estimate by itself CPU weights at startup *
!-#define AUTO_WEIGHTS

!* Check load balance every NBALANCE steps and adjust it.
!*  This option activates ROOT_SAVE by default */
!-#define AUTO_BALANCE
#define NBALANCE (5)

!* Uncomment this if you want to do some profiling on the code *
!-#define LOGTIME

!* How often to check CFL, energy etc *
#define NCHECKCFL (10)

!* Define the coefficients in the equations

#ifndef TEMPERATURE
!* Navier-Stokes code
#define VDIFF invRe

#else
!* convective code

#ifdef FINGER
!* Salt fingers

#define VDIFF (Pr*Le)
#define BUOYT (+Pr*Le*Rt)
#define BUOYS (-Pr*Le*Rs)
#define TDIFF Le
#define SDIFF (1.d0)

#else
!* Rayleigh-Benard
#define VDIFF (Pr)
#define BUOYT (Ra*Pr)
#define TDIFF (1.d0)

#endif
#endif


!******* Do not bother to look below this line **************

!* Maximum number of threads *
#define MAXCPU (1000)

#ifndef NOMPI
#define cmpi 
#endif
#ifdef LOGTIME
#define ctime
#endif


#ifdef AUTO_BALANCE 
#ifndef ROOT_SAVE
#define ROOT_SAVE
#endif
#endif

#ifdef HAVE_FFTW
#define INITFOUR_FFTW initfour
#define FOURIER1_FFTW fourier1
#define FOURIER2_FFTW fourier2
#define FOURIER2_FFTW_SINGLE2D fourier2_single2d
#else
#define INITFOUR_TMPT initfour
#define FOURIER1_TMPT fourier1
#define FOURIER2_TMPT fourier2
#define FOURIER2_TMPT_SINGLE2D fourier2_single2d
#endif


!* If no MPI remove useless options *
#ifdef NOMPI
#ifdef AUTO_BALANCE
#undef AUTO_BALANCE
#endif
#ifdef ROOT_SAVE
#undef ROOT_SAVE
#endif
#undef SAVENODE
#define SAVENODE (1)
#ifdef AUTO_WEIGHTS
#undef AUTO_WEIGHTS
#endif
#endif

#ifndef NWORD
#define NWORD (4)
#endif

#ifdef FFTWUNDER
#define MPI_INIT mpi_init_
#define MPI_COMM_RANK mpi_comm_rank_
#define MPI_COMM_SIZE mpi_comm_size_
#define MPI_TYPE_VECTOR mpi_type_vector_
#define MPI_FINALIZE mpi_finalize_
#define MPI_SENDRECV MPI_SENDRECV_
#define MPI_GATHER mpi_gather_
#define MPI_ALLGATHER mpi_allgather_
#define MPI_BCAST mpi_bcast_
#define MPI_RECV mpi_recv_
#define MPI_SEND mpi_send_
#define MPI_BARRIER mpi_barrier_
#define MPI_COMMIT mpi_commit_
#endif
