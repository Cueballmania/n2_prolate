BINDIR = ..

FC = gfortran
FFLAGS = -Wall -std=f95    -fPIC

LIB = -Wl,-framework -Wl,veclib ~/ARPACK/libarpack_MACOSX.a

.SUFFIXES : .o .f90 .f95 .f08

OBJS = mainprolate.o SVD_orthogonalization.o gaussquad.f08 diagwrap.o hamiltonian.o\
         stos.o V_potential.o index_2d.o xi_dvr.o eta_dvr.o SVD_insert.o svd_lagpot.o\
         gaussmat.o xform.o readmesa.o SVD_eigen.o arnoldi.o

.f08.o:
	$(FC) -c $(FFLAGS) $<

.f90.o:
	$(FC) -c $(FFLAGS) $<

.f95.o:
	$(FC) -c $(FFLAGS) $<

xprolate: $(OBJS) gaussquad.mod
	$(FC) $(OBJS) -o $(BINDIR)/$@ $(LIB)
