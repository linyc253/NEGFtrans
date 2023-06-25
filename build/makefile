FFLAGS=-O3
FC=gfortran
SOURCEDIR=../src
TARGET=NEGFtrans

OBJS=test.o test2.o

# Check tab by: cat -t makefile
all: $(TARGET).x

$(TARGET).x: $(SOURCEDIR)/main.f90 $(OBJS)
	$(FC) $(FFLAGS) $(SOURCEDIR)/main.f90 -o $(TARGET).x


%.o:$(SOURCEDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm *.o *.mod $(TARGET).x

#=========== Deal with dependency============#
test.o:test2.o
