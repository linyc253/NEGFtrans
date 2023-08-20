FFLAGS=-O3
FC=gfortran
MPIFC=mpifort
SOURCEDIR=../src
TARGET=NEGFtrans

# For the file 'fftw3.f03'
INCLUDE= -I/usr/local/include

LIBS=-lfftw3 -llapack -lblas

OBJS=grid.o math_kernel.o negf.o global.o

# Check tab by: cat -t makefile
all: $(TARGET).x

$(TARGET).x: $(SOURCEDIR)/main.f90 $(OBJS)
	$(MPIFC) $(FFLAGS) $(INCLUDE) $(SOURCEDIR)/main.f90 -o $(TARGET).x $(OBJS) $(LIBS)


%.o:$(SOURCEDIR)/%.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c $< -o $@

clean:
	rm *.o *.mod $(TARGET).x

#=========== Deal with dependency here============#
negf.o:grid.o math_kernel.o global.o
grid.o:global.o
math_kernel.o:global.o
