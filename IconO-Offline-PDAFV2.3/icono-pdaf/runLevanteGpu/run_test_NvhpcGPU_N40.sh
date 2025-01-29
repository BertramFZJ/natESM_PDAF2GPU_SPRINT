#!/bin/bash

# Run script to perform a single offline analysis with icono-pdaf

### Levante
#SBATCH --job-name=pdaf_OpenMP_GPU # Specify job name

#SBATCH --partition=gpu    # Specify partition name
#SBATCH --nodes=1             # Specify number of nodes
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16    # Specify number of CPUs per task
#SBATCH --time=00:30:00        # Set a limit on the total run time

#SBATCH --output=my_job.o%j    # File name for standard output
#SBATCH --error=my_job.e%j     # File name for standard error output

#SBATCH --mem=0                # Request full memory of a node (needed for scalability tests)
#SBATCH --account=ka1298      # Charge resources on this project account
#SBATCH --mail-type=NONE

##SBATCH --gpus=4
#SBATCH --gpus-per-node=4
#SBATCH --gpus-per-task=1
#SBATCH --constraint=a100_80

export LD_LIBRARY_PATH=/sw/spack-levante/nvhpc-22.5-v4oky3/Linux_x86_64/22.5/math_libs/lib64

### Bind your OpenMP threads
export OMP_NUM_THREADS=16
# export OMP_SCHEDULE="auto"
export OMP_SCHEDULE="dynamic,100"
export OMP_STACKSIZE=128M

## export KMP_AFFINITY="verbose,granularity=fine,scatter"
## export KMP_LIBRARY="turnaround"

set -e
ulimit -s 204800
ulimit -c 0

export OMPI_MCA_osc="ucx"
export OMPI_MCA_pml="ucx"
export OMPI_MCA_btl="self"
export UCX_HANDLE_ERRORS="bt"
export OMPI_MCA_pml_ucx_opal_mem_hooks=1

export OMPI_MCA_io="romio321"          # basic optimisation of I/O
export UCX_TLS="shm,rc_mlx5,rc_x,self" # for jobs using LESS than 150 nodes

# export UCX_TLS="shm,dc_mlx5,dc_x,self" # for jobs using MORE than 150 nodes
# export UCX_UNIFIED_MODE="y"            # JUST for homogeneous jobs on CPUs, do not use for GPU nodes

####################################################################

set -x

# first test setup for ICONO R2B6, 21.01.2020
y1=1960
m1=01
y2=1959
m2=12
echo $y1 $m1 $y2 $m2

# Paths and file names
basedir=/work/ab0995/a270007/ICONO-PDAF
#pdaf_exe=$basedir'/icono_offline_mpi_PDAF221/src/icono_offline'
pdaf_exe='./icono_offline'
indir=$basedir'/Pdaf_restart/Unfltd'
prefix=hpo0011_restart_oce_r
suffix='_'$y1$m1'01T000000Z'
obs_dir=$basedir'/EN4'
obs_file='en4_'$y2$m2'_1800x900_R2B6.nc'
grid_dir=$basedir'/grid'
grid_file='icon_grid_0031_R02B06_O.nc'
outdir='.'
#pdaf_fltd=$basedir'/assim_test_mpi'
pdaf_prefix='ana'


# Configuration for PDAF
dim_ens=40          # Ensemble size
multivar='T'        # T: multivariate DA; F: univariate 2-step DA
assim_tho='T'       # T: Assimilate EN4 temperature data
assim_sao='T'       # T: Assimilate EN4 salinity data
rms_tho=1.0         # observation error (0= overwrite with obs, 100= keep model state)
rms_sao=1.0         # observation error (0= overwrite with obs, 100= keep model state)
rms_min=0.001       # Minimum limit for observation error
forget=1.0          # Forgetting factor
filtertype=7        # Filter type: (7) LESTKF
cradius=5.0         # horizontal radius in coords units 
                    #    (for disttype=1: degree, 0=only at obs grid cell, 360= all over globe ="global EnKF")
                    #    (for disttype>2: km, 0=only at obs grid cell, 360= all over globe ="global EnKF")
cradius_z=1.0       # vertical radius in levels (0= only in obs level, 1=in obs level and 1 above and 1 below, 64= all levels)
disttype=11         # Type of distances: (1) Cartesian, (2) simplified geographic, (3) haversine function
                    #    (11) Cartesian 2+1D factorized
locweight=2         # 5th order polynomial Gaspari and Cohn
writeens='T'        # Write full analysis ensemble
writemean='T'       # Write analysis ensemble mean state

# Replace slashes by colons for command-line parsing
indirTemp=`echo ${indir} | sed 's/\//:/g'`
outdirTemp=`echo ${outdir} | sed 's/\//:/g'`
obs_dirTemp=`echo ${obs_dir} | sed 's/\//:/g'`
grid_dirTemp=`echo ${grid_dir} | sed 's/\//:/g'`

#
echo ${pdaf_exe} -indir ${indirTemp} -inprefix ${prefix} -insuffix ${suffix} \
	    -obsdir ${obs_dirTemp} -obsfile ${obs_file} \
	    -griddir ${grid_dirTemp} -gridfile ${grid_file} \
	    -outdir ${outdirTemp} -outprefix ${pdaf_prefix}_ \
	    -dim_ens ${dim_ens} -multivar ${multivar} \
	    -forget ${forget} -filtertype ${filtertype} \
            -cradius ${cradius} -cradius_z ${cradius_z} -locweight ${locweight} \
	    -disttype ${disttype} -rms_min ${rms_min} \
	    -assim_EN4_tho ${assim_tho} -assim_EN4_sao ${assim_sao} \
	    -rms_obs_EN4_tho ${rms_tho} -rms_obs_EN4_sao ${rms_sao} \
	    -writeens ${writeens} -writemean ${writemean}
	    
# srun -l --cpu_bind=verbose --cpus-per-task=4 --hint=nomultithread \
#  --distribution=block:cyclic:block ./myprog

## export NVCOMPILER_ACC_NOTIFY=3

srun -l  --cpu_bind=verbose --cpus-per-task=${OMP_NUM_THREADS} --hint=nomultithread --distribution=block:block:cyclic \
	    ${pdaf_exe} -indir ${indirTemp} -inprefix ${prefix} -insuffix ${suffix} \
	    -obsdir ${obs_dirTemp} -obsfile ${obs_file} \
	    -griddir ${grid_dirTemp} -gridfile ${grid_file} \
	    -outdir ${outdirTemp} -outprefix ${pdaf_prefix}_ \
	    -dim_ens ${dim_ens} -multivar ${multivar} \
	    -forget ${forget} -filtertype ${filtertype} \
            -cradius ${cradius} -cradius_z ${cradius_z} -locweight ${locweight}\
	    -disttype ${disttype} -rms_min ${rms_min} \
	    -assim_EN4_tho ${assim_tho} -assim_EN4_sao ${assim_sao} \
	    -rms_obs_EN4_tho ${rms_tho} -rms_obs_EN4_sao ${rms_sao} \
            -writeens ${writeens} -writemean ${writemean} \
            1>out.icono-pdaf 2>err.icono-pdaf

