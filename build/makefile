FFLAGS=-O3
FC=gfortran
SOURCEDIR=../src
TARGET=NEGFtrans
INCLUDE= -I/usr/local/include
LIBS=-lfftw3 -llapack -lblas

OBJS=fxc.o anderson_mixing.o grid.o negf.o

# Check tab by: cat -t makefile
all: $(TARGET).x

$(TARGET).x: $(SOURCEDIR)/main.f90 $(OBJS)
	$(FC) $(FFLAGS) $(INCLUDE) $(SOURCEDIR)/main.f90 -o $(TARGET).x $(OBJS) $(LIBS)


%.o:$(SOURCEDIR)/%.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c $< -o $@

clean:
	rm *.o *.mod $(TARGET).x

#=========== Deal with dependency here============#
#test.o:test2.o
