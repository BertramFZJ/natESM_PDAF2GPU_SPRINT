######################################################
# Include file with machine-specific definitions     #
# for building PDAF.                                 #
#                                                    #
# Variant for Levante@DKRZ                           #
# (Intel Fortran Compiler and Intel MPI)             #
######################################################

# Compiler, Linker, and Archiver
FC = /sw/spack-levante/openmpi-4.1.4-3qb4sy/bin/mpif90 
LD = $(FC)
AR = ar
RANLIB = ranlib 

# C preprocessor
# (only required, if preprocessing is not performed via the compiler)
# CPP = /usr/bin/cpp NATESM

# C compiler
CC = /sw/spack-levante/openmpi-4.1.4-3qb4sy/bin/mpicc
OPTCC= -D_GNU_SOURCE -mp -O2 -acc=gpu -gpu=ptxinfo -gpu=cc80 -gpu=cuda11.7 -Minfo=accel

# Definitions for CPP
# Define USE_PDAF to include PDAF
# (if the compiler does not support get_command_argument()
# from Fortran 2003 you should define F77 here.)
CPP_DEFS = -DUSE_PDAF

# Optimization specs for compiler
#   (You should explicitly define double precision for floating point
#   variables in the compilation)  
OPT= -mp -g -O2 -r8 -cpp -acc=gpu -gpu=ptxinfo -gpu=cc80 -gpu=cuda11.7 -Minfo=accel

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
# from email by Mathis and https://docs.dkrz.de/doc/levante/code-development/compiling-and-linking.html for netcdf-fortran
LINK_LIBS = -Wl,-rpath,/sw/spack-levante/netcdf-fortran-4.5.4-syv4qr/lib \
            -cudalib=cublas -cudalib=cusolver -lblas -llapack

# Specifications for the archiver
AR_SPEC = 

# Specifications for ranlib
RAN_SPEC =

# Specification for directory holding modules (-module for Intel, -J for GNU)
MODULEOPT = -module

# Include path for MPI header file
MPI_INC =            # dummy MPI for PDAF versions < 2

# Object for nullMPI - if compiled without MPI library
OBJ_MPI =            # needed when running without MPI for PDAF versions < 2

# NetCDF 
# June 2022 on levante
NC_LIB   = -L/sw/spack-levante/netcdf-fortran-4.5.4-syv4qr/lib -lnetcdff
NC_INC   = -I/sw/spack-levante/netcdf-fortran-4.5.4-syv4qr/include
