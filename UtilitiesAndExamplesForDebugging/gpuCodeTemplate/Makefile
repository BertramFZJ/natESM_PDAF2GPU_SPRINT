# ARCH = cc80
ARCH = sm_86

CUDAVER = cuda12.3
LTOFLAG = -gpu=lto

CUBLAS_LIB = -L/opt/nvidia/hpc_sdk/Linux_x86_64/23.11/math_libs/lib64

# CXX = /sw/spack-levante/openmpi-4.1.4-3qb4sy/bin/mpic++
CXX = mpic++
CXX_COMPILER_OPTIONS = -mp -acc=gpu -gpu=ptxinfo -gpu=$(ARCH) -gpu=$(CUDAVER) $(LTOFLAG) -cuda -Minfo=accel -O2 -Mnoautoinline
SOURCES_CXX =
OBJECTS_CXX = $(patsubst %.cpp, %.o, $(SOURCES_CXX))

# FC = /sw/spack-levante/openmpi-4.1.4-3qb4sy/bin/mpif90
FC = mpif90
FORTRAN_COMPILER_OPTIONS = -mp -acc=gpu -gpu=ptxinfo -gpu=$(ARCH) -gpu=$(CUDAVER) $(LTOFLAG) -cuda -cpp -Minfo=accel -O2 -Mnoautoinline
FORTRAN_LINKER_OPTIONS = $(FORTRAN_COMPILER_OPTIONS)
FORTRAN_LINKER_LIBS = $(CUBLAS_LIB) -cudalib=cublas -cudalib=cusolver -lblas -llapack
SOURCES_FORTRAN = linearAlgebraModule.F90 main.F90
OBJECTS_FORTRAN = $(patsubst %.F90, %.o, $(SOURCES_FORTRAN))

# CC = /sw/spack-levante/openmpi-4.1.4-3qb4sy/bin/mpicc
CC = mpicc
C_COMPILER_OPTIONS = -D_GNU_SOURCE -mp -acc=gpu -gpu=ptxinfo -gpu=$(ARCH) $(LTOFLAG) -gpu=$(CUDAVER) -cuda -Minfo=accel -O2 -Mnoautoinline
SOURCES_C =
OBJECTS_C = $(patsubst %.c, %.o, $(SOURCES_C))

TARGET = gpuCodeTemplate

%.o : %.c
	@echo "compiling $^ ..."
	@$(CC) $(C_COMPILER_OPTIONS) -c $^ -o $@

%.o : %.cpp
	@echo "compiling $^ ..."
	@$(CXX) $(CXX_COMPILER_OPTIONS) -c $^ -o $@

%.o : %.F90
	@echo "compiling $^ ..."
	@$(FC) $(FORTRAN_COMPILER_OPTIONS) -c $^ -o $@

all : $(OBJECTS_C) $(OBJECTS_CXX) $(OBJECTS_FORTRAN)
	@echo "building $(TARGET).exe ..."
	@$(FC) $(OBJECTS_C) $(OBJECTS_CXX) $(OBJECTS_FORTRAN) $(FORTRAN_LINKER_OPTIONS) $(FORTRAN_LINKER_LIBS) -o $(TARGET).exe
#	@echo "copying $(TARGET).exe to bin directory ..."
#	@mv $(TARGET).exe ./bin/

clean :
	@rm -f *.mod ./bin/$(TARGET).exe $(TARGET).exe $(OBJECTS_C) $(OBJECTS_CXX) $(OBJECTS_FORTRAN)

