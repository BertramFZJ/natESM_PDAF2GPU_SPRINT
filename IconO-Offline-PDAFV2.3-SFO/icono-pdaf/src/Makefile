# $Id: Makefile 746 2009-08-04 12:16:28Z lnerger $
#######################################################
# Generic Makefile for to build PDAF with dummy model #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

# User specifications
# 1. Set BASEDIR, the directory where the PDAF package resides
# 2. Set PDAF_ARCH to include compile definitions
#    (See directory BASEDIR/make.arch for files. PDAF_ARCH is filename without .h)

# Root directory of PDAF package
BASEDIR = ../../PDAF

# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include $(BASEDIR)/make.arch/$(PDAF_ARCH).h

# Name of executable
EXE = icono_offline

# End of user specifications
######################################################

.SUFFIXES: .F90 .o 

# Modules used for the model part
MODULES =  	mod_parallel.o \
		timer.o \
		mod_memcount.o \
		mod_model.o \
		parser_mpi.o \
		mod_statevector_pdaf.o \
		mod_io_pdaf.o

# Module required for assimilation
MOD_ASSIM = 	mod_assimilation.o 

# Model routines
OBJ_MODEL =  	main_offline.o \
		initialize.o

# Routines of observation handling (PDAF-OMI)
OBJ_USER_PDAFOMI = obs_EN4_tho_pdafomi.o \
		obs_EN4_sao_pdafomi.o \
		callback_obs_pdafomi.o

# Interface to PDAF - model sided
OBJ_PDAF_INT =  init_parallel_pdaf.o \
		init_pdaf_offline.o \
		init_pdaf_parse.o \
		init_pdaf_info.o \
		assimilation_pdaf_offline.o

# Generic user-supplied routines
OBJ_USER_GEN  = init_ens_offline.o \
		collect_state_pdaf_offline.o \
		prepoststep_ens_offline.o #\
#		read_dims.o \
#		read_ens.o \
#		write_ens.o

# User-supplied routines for localized analysis
OBJ_USER_LOCAL = init_n_domains_pdaf.o \
		init_dim_l_pdaf.o \
		g2l_state_pdaf.o \
		l2g_state_pdaf.o

# Full list of user-supplied routines for online modes
OBJ_PDAF_USER = $(OBJ_USER_PDAFOMI) $(OBJ_USER_GEN) $(OBJ_USER_LOCAL)

######################################################

$(EXE) : $(MOD_ASSIM) $(MODULES) $(OBJ_MODEL) $(OBJ_PDAF_USER) $(OBJ_PDAF_INT)
	$(LD)  $(OPT_LNK)  -o $@  \
	$(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_PDAF_USER) $(OBJ_PDAF_INT) $(OBJ_MPI) \
	-L$(BASEDIR)/lib -lpdaf-d \
	$(NC_LIB) $(LINK_LIBS) 
	@echo "++++++ Done ++++++"

icono_pdaf-ol.a : $(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER)
	@echo "++++++ Generate archive for model offline_1D ++++++"
	$(AR) -r $(AR_SPEC) $@ $(MODULES) $(MOD_ASSIM) $(OBJ_MODEL) $(OBJ_PDAF_INT) $(OBJ_PDAF_USER)
	@echo "++++++ Done ++++++"

.F90.o :
	$(FC) $(OPT) $(CPP_DEFS) -I$(BASEDIR)/include $(MPI_INC) $(NC_INC) -c $*.F90

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

clean :
	rm -f *.o *.mod *.a $(EXE)
