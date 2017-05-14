# Makefile for RB/Navier Stokes code

#Comment this line for MPI
NOMPI = 1

MPIF77COMP = mpif77
OMPI_F77 =      gfortran

#Gfortran

F77COMP	=	gfortran
LINKER =	gfortran
CFLAGS =       -O3 -ffixed-line-length-none
#LDFLAGS        =   -static -lrfftw -lfftw   -lm
LDFLAGS        =  -L/opt/local/lib  -lfftw3 -lm


# Intel ?86/g77 options:
#F77COMP	=	f77 
#LINKER =	f77   
#CFLAGS =       -O6 -mcpu=pentiumpro -malign-double -fno-schedule-insns \
               -fschedule-insns2 -fomit-frame-pointer \
               -ffortran-bounds-check
#LDFLAGS        =        -lfftw3  -lm
#LDFLAGS        =   -static -L/home/jost/local/lib  -lrfftw -lfftw   -lm

FFTWDIR  = ${FFTW_LIBDIR}

# Intel Fortran Compiler

#F77COMP =       ifort
#LINKER =        ifort
#CFLAGS =        -O3   -no-prec-div -mcmodel=medium -shared-intel
#LDFLAGS =  -lfftw3  -lm farg.o  -L ${FFTWDIR}
#LDFLAGS =  -lfftw3   -L ${FFTWDIR}


# Intel Fortran Compiler OSX
#F77COMP =       ifort
#LINKER =        ifort
#CFLAGS =	-nowarn -O3 -ftz -ip -ipo  -us -I/usr/local/include
#LDFLAGS        =  -L/opt/local/lib  -lfftw3 -lm

#Digital g77
#F77COMP = g77
#LINKER =	g77
#CFLAGS =       -O6
#LDFLAGS        =   -L.    -lrfftw -lfftw -lm



ifdef NOMPI
COMP	= ${F77COMP}
else
MPICH_F77 = ${F77COMP}
MPICH_F77LINKER = ${F77COMP}
COMP   = ${MPIF77COMP}    
LINKER = ${MPIF77COMP} 
export MPICH_F77
export MPICH_FLINKER
endif 


NAVIER_OBJ	=	boundary.o build_star.o conv.o crea_df.o \
			crea_e.o diffus.o fourier1.o  \
			fourier2.o  \
			main.o mpitools.o poiss_im.o punti.o salva.o tn_poiss.o \
			tools.o tridiagmpi.o velocita.o farg.o

NAVIER_REDOM_OBJ	=	boundary.o build_star.o conv.o crea_df.o \
			crea_e.o diffus.o fourier1.o  \
			fourier2.o  \
			main_redom.o mpitools.o poiss_im.o punti.o salva.o tn_poiss.o \
			tools.o tridiagmpi.o velocita.o farg.o

DEPEND= makedepend 


all:	mpi

clean:
	-rm *.o core 

depend: 
	$(DEPEND) *.F

mpi: $(NAVIER_OBJ)
	$(LINKER)  $(CFLAGS) -o rb $(NAVIER_OBJ)  $(LDFLAGS) 

redom: $(NAVIER_OBJ)
	$(LINKER)  $(CFLAGS) -o rb_redom $(NAVIER_REDOM_OBJ)  $(LDFLAGS) 


NUSSELTSALT_OBJ    =       nusseltsalt.o fourier1.o fourier1_tmpt.o tools.o salva.o
nusseltsalt: $(NUSSELTSALT_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o nusseltsalt $(NUSSELTSALT_OBJ) $(LDFLAGS)

NUSSELT_OBJ    =       nusselt.o fourier1.o fourier1_tmpt.o tools.o salva.o punti.o
nusselt: $(NUSSELT_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o nusselt $(NUSSELT_OBJ) $(LDFLAGS)

SECTIONVSER_OBJ    =       sectionvser.o fourier1ser.o tools.o salva.o
sectionvser: $(SECTIONVSER_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o sectionvser $(SECTIONVSER_OBJ) $(LDFLAGS)

CONVERTBOVSER_OBJ    =       convertbovser.o fourier1ser.o tools.o salva.o
convertbovser: $(CONVERTBOVSER_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o convertbovser $(CONVERTBOVSER_OBJ) $(LDFLAGS)

SECTIONV_OBJ    =       sectionv.o fourier1.o fourier1_tmpt.o tools.o salva.o
sectionv: $(SECTIONV_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o sectionv $(SECTIONV_OBJ) $(LDFLAGS)

SECTIONVANOM_OBJ    =       sectionvanom.o fourier1.o fourier1_tmpt.o tools.o salva.o
sectionvanom: $(SECTIONVANOM_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o sectionvanom $(SECTIONVANOM_OBJ) $(LDFLAGS)


MEANFLUXH_OBJ    =       meanfluxh.o fourier1.o fourier1_tmpt.o tools.o salva.o
meanfluxh: $(MEANFLUXH_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o meanfluxh $(MEANFLUXH_OBJ) $(LDFLAGS)

MEANH_OBJ    =       meanh.o fourier1.o fourier1_tmpt.o tools.o salva.o
meanh: $(MEANH_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o meanh $(MEANH_OBJ) $(LDFLAGS)

SECTIONH_OBJ    =       sectionh.o fourier1.o fourier1_tmpt.o tools.o salva.o
sectionh: $(SECTIONH_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o sectionh $(SECTIONH_OBJ) $(LDFLAGS)

FLUX_OBJ    =       flux.o fourier1.o fourier1_tmpt.o tools.o salva.o
flux: $(FLUX_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o flux $(FLUX_OBJ) $(LDFLAGS)

EXTR_OBJ    =       extr.o fourier1.o fourier1_tmpt.o tools.o salva.o
extr: $(FLUX_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o extr.x $(EXTR_OBJ) $(LDFLAGS)


FLUXV_OBJ    =       fluxv.o fourier1.o fourier1_tmpt.o tools.o salva.o
fluxv: $(FLUXV_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o fluxv $(FLUXV_OBJ) $(LDFLAGS)

SPECTRUMH_OBJ    =       spectrumh.o tools.o salva.o
spectrumh: $(SPECTRUMH_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o spectrumh $(SPECTRUMH_OBJ) $(LDFLAGS)

SPECTRUMHBAND_OBJ    =       spectrumhband.o  fourier1ser.o tools.o salva.o
spectrumhband: $(SPECTRUMHBAND_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o spectrumhband $(SPECTRUMHBAND_OBJ) $(LDFLAGS)

SPECTRUMX_OBJ    =       spectrumx.o fourier1.o  tools.o salva.o
spectrumx: $(SPECTRUMX_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o spectrumx $(SPECTRUMX_OBJ) $(LDFLAGS)

VORT_OBJ    =       vort.o  fourier1.o fourier1_tmpt.o tools.o salva.o 
vort: $(KMAX_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o vort $(VORT_OBJ) $(LDFLAGS)

KMAXZ_OBJ    =       kmaxz.o tools.o salva.o
kmaxz: $(KMAXZ_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o kmaxz $(KMAXZ_OBJ) $(LDFLAGS)

KMAXT_OBJ    =       kmaxt.o tools.o salva.o
kmaxt: $(KMAXT_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o kmaxt $(KMAXT_OBJ) $(LDFLAGS)

KMAX1D_OBJ    =       kmax1d.o tools.o salva.o
kmax1d: $(KMAX1D_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o kmax1d $(KMAX1D_OBJ) $(LDFLAGS)

KMAX_OBJ    =       kmax.o tools.o salva.o
kmax: $(KMAX_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o kmax $(KMAX_OBJ) $(LDFLAGS)

CONVERTBOB_OBJ	=	convertbob.o fourier1.o fourier1_tmpt.o tools.o salva.o
convertbob: $(CONVERTBOB_OBJ)	
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o convertbob $(CONVERTBOB_OBJ) $(LDFLAGS)

CONVERTBOV_OBJ	=	convertbov.o fourier1.o fourier1_tmpt.o tools.o salva.o
convertbov: $(CONVERTBOV_OBJ)	
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o convertbov $(CONVERTBOV_OBJ) $(LDFLAGS)

CENSUS_OBJ	=	census.o fourier1.o fourier1_tmpt.o tools.o salva.o
census: $(CENSUS_OBJ)	
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o census $(CENSUS_OBJ) $(LDFLAGS)

PROFILES_OBJ	=	prof.o fourier1.o fourier1_tmpt.o tools.o salva.o farg.o
prof: $(PROFILES_OBJ)	
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o prof $(PROFILES_OBJ) $(LDFLAGS)

PROFVAR_OBJ	=	profvar.o fourier1.o fourier1_tmpt.o tools.o salva.o farg.o
profvar: $(PROFVAR_OBJ)	
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o profvar $(PROFVAR_OBJ) $(LDFLAGS)

PROFILESSER_OBJ    =       profser.o fourier1ser.o tools.o salva.o
profser: $(PROFILESSER_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o profser $(PROFILESSER_OBJ) $(LDFLAGS)

PDFJOINTBF_OBJ    =       pdfjointbf.o fourier1ser.o tools.o salva.o
pdfjointbf: $(PDFJOINTBF_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o pdfjointbf $(PDFJOINTBF_OBJ) $(LDFLAGS)



INICONDS_OBJ    =       iniconds.o tools.o salva.o ran2.o punti.o boundary.F farg.o
iniconds: $(INICONDS_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o iniconds $(INICONDS_OBJ) $(LDFLAGS)

INICONDSALTSER_OBJ    =       inicondsaltser.o tools.o salva1.o ran2.o punti.o farg.o
inicondsaltser: $(INICONDSALTSER_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o inicondsaltser $(INICONDSALTSER_OBJ) $(LDFLAGS)

INICONDSALTSER2_OBJ    =       inicondsaltser2.o tools.o salva1.o ran2.o punti.o farg.o
inicondsaltser2: $(INICONDSALTSER2_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o inicondsaltser2 $(INICONDSALTSER2_OBJ) $(LDFLAGS)


INICONDSALT_OBJ    =       inicondsalt.o tools.o salva.o ran2.o punti.o boundary.F farg.o
inicondsalt: $(INICONDSALT_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o inicondsalt $(INICONDSALT_OBJ) $(LDFLAGS)

INICOND_OBJ    =       inicond.o tools.o salva.o ran2.o punti.o boundary.F
inicond: $(INICOND_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o inicond $(INICOND_OBJ) $(LDFLAGS)

ANOMALY_OBJ    =       anomaly.o tools.o salva.o punti.o fourier1.o \
                       fourier1_tmpt.o
anomaly: $(ANOMALY_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o anomaly $(ANOMALY_OBJ) $(LDFLAGS)

FLUXCONV_OBJ    =       fluxconv.o fourier1.o fourier1_tmpt.o tools.o \
                        salva.o fourier2.o
fluxconv: $(FLUXCONV_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o fluxconv $(FLUXCONV_OBJ) $(LDFLAGS)

PDF_OBJ    =            pdf.o fourier1.o fourier1_tmpt.o tools.o \
                        salva.o 
pdf: $(PDF_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o pdf $(PDF_OBJ) $(LDFLAGS)

CFL_OBJ    =            cfl.o fourier1.o fourier1_tmpt.o tools.o \
                        salva.o punti.o
cfl: $(CFL_OBJ)
	$(COMP)  $(CFLAGS) -D"NOMPI"  -o cfl $(CFL_OBJ) $(LDFLAGS)



ifdef NOMPI
%.o : %.F
	$(COMP) $(CFLAGS) -D"NOMPI" -c $< -o $@
census.o : census.F
	$(COMP) $(CFLAGS) -D"NOMPI" -c   $< -o $@
else
%.o : %.F
	$(COMP) $(CFLAGS) -c $< -o $@
endif

# DO NOT DELETE

anomaly.o: config.h param.h
boundary.o: config.h param.h
build_star.o: config.h param.h
cfl.o: config.h param.h
convertbob.o: config.h param.h
convertbov.o: config.h param.h
conv.o: config.h param.h
crea_df.o: config.h param.h
crea_e.o: config.h param.h
diffus.o: config.h param.h
div.o: config.h param.h
extr.o: config.h param.h
flux.o: config.h param.h
fluxconv.o: config.h param.h 
fluxv.o: config.h param.h
fourier1.o: config.h
fourier1_tmpt.o: config.h param.h
fourier2.o: config.h
fourier2_tmpt.o: config.h param.h
fourier2_tst.o: config.h param.h
inicond.o: config.h param.h
kmax.o: config.h param.h
kmaxt.o: config.h param.h
kmaxv.o: config.h param.h
main.o: config.h param.h
main_redom.o: config.h param.h
maino.o: config.h param.h
mpitools.o: config.h param.h
pdf.o: config.h param.h
poiss_im.o: config.h param.h
prof.o: config.h param.h
profser.o: config.h param.h
punti.o: config.h param.h
salva.o: config.h param.h
sectionh.o: config.h param.h
meanh.o: config.h param.h
meanfluxh.o: config.h param.h
sectionv.o: config.h param.h
sectionz.o: config.h param.h
slice.o: config.h param.h
spectrumh.o: config.h param.h
spectrumx.o: config.h param.h
spectrumt.o: config.h param.h
tn_poiss.o: config.h param.h
tools.o: config.h param.h
tridiagmpi.o: config.h
velocita.o: config.h param.h
vort.o: config.h param.h



