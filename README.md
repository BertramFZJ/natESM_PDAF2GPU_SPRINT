<h1 style="text-align: center; color: blue; font-size: 28pt;">Sprint #16: Performance evaluation, optimization, and strategy development for porting PDAF to GPUs</h1>

This GIT repository contains modified and added during the sprint source code files of the PDAF (The Parallel Data Assimilation Framework) and
the standalone offline assimilation program PDAF-ICONO (PDAF implementation used as assimilation method for …), as well as configuration files,
makefiles and runscripts for compiling and running PDAF-ICONO executables on the compute and gpu partitions of the Levante HPC System (DKRZ).
The three main directories with source code files relate to three different issues:
<table>
<tr> <td>1.</td> <td>IconO-Offline-PDAFV2.3</td>
<td>Porting the <code>PDAF_lestkf_analysis.F90::PDAF_lestkf_analysis</code> subroutine to GPU</td>
</tr>
<tr> <td>2.</td> <td>IconO-Offline-PDAFV2.3-SFO</td>
<td> Porting to GPU the subroutine<br><code>PDAFomi_dim_obs_l.F90::PDAFomi_check_dist2_noniso_loop</code>,
which is part of the search for local observations module</td>
</tr>
<tr> <td>3.</td> <td>IconO-Offline-PDAFV2.3-CPU</td>
<td> Linking of the module <code>mo_mem_workspaces.F90::mo_mem_workspaces</code> for allocation of temporary arrays 
from the heap withing <code>PDAF_lestkf_analysis.F90::PDAF_lestkf_analysis</code> subroutine</td>
</tr>
</table>
<i>Note:</i> The initial commit was made by uploading already modified files.

<h2 style="text-align: center; font-size: 24pt">Workspace setup</h2>

The workspace setup procedure includes creating a working directory and subsequently cloning three Git repositories into it:<br>
<b>A1.</b> Create the working directory <code>{work_dir}</code> and navigate into it<br>
<code>mkdir {work_dir}
cd {work_dir}</code>                    
<b>A2.</b> Clone the PDAF repository into the working directory and switch to branch 2.3<br>
<code>git clone git@github.com:PDAF/PDAF.git <br>
cd PDAF <br>
git checkout PDAF_V2.3 <br>
cd ..</code><br>
<b>A3.</b> Clone the PDAF-ICONO repository and switch to branch for-natESM<br>
<code>git clone git@gitlab.dkrz.de:a270007/icono-pdaf.git <br>
cd icono-pdaf <br>
git checkout for-natESM <br>
cd ..</code><br>
<b>A4.</b> Clone this repository into the working directory<br>
<code>git clone git@github.com:BertramFZJ/natESM_PDAF2GPU_SPRINT.git</code><br>

<h2 style="text-align: center; font-size: 24pt">Compilation and launch of the standalone offline data assimilation program PDAF-ICONO with external library linking on Levante (DKRZ)</h2>

<b>B1.</b> Copy the libraries added during the sprint to the folder containing the source files of the PDAF framework<br>
<code>cp -r ./IconO-Offline-PDAFV2.3/PDAF/external/natEsmP2G ./PDAF/external/</code><br>
<b>B2.</b> Copy the updated header files required for compilation<br>
<code>cp ./IconO-Offline-PDAFV2.3/PDAF/make.arch/* ./PDAF/make.arch/<br></code>
<br><b>B3.</b> Copy the updated Makefiles for the framework and the standalone program<br>
<code>cp ./IconO-Offline-PDAFV2.3/PDAF/src/Makefile ./PDAF/src/<br>
cp ./IconO-Offline-PDAFV2.3/icono-pdaf/src/Makefile ./icono-pdaf/src/</code><br>
<b>B4.</b> Create the required additional subdirectories in the root directory of the framework (this step is performed only once)<br>
<code>mkdir ./PDAF/build ./PDAF/include ./PDAF/lib</code><br>
<b>B5.</b> Compile the framework using one of the three available configurations<br>
<code>cd ./PDAF/src/<br>
make PDAF_ARCH=Levante_mpiifort_impi<br>
make PDAF_ARCH=Levante_nvhpc_openmpi<br>
make PDAF_ARCH=Levante_nvhpc_gpu_openmpi<br>
cd ../../</code><br>
Here, the name of the header file containing the compiler name and settings is specified as an argument to the make command:<br>
<table>
<tr> <td>Levante_mpiifort_impi</td> <td>Intel compiler & Intel MPI</td> </tr>
<tr> <td>Levante_nvhpc_openmpi</td> <td>Nvidia compiler (cpu executable) & OpenMPI</td> </tr>
<tr> <td>Levante_nvhpc_gpu_openmpi</td> <td>Nvidia compiler (gpu executable) & OpenMPI</td> </tr>
</table>
To compile using the Intel compiler, two modules need to be preloaded<br>
<code>module load intel-oneapi-compilers<br>
module load intel-oneapi-mpi/2021.5.0-gcc-11.2.0</code><br>
To compile using the Nvidia compiler (version 22.5), no additional modules need to be loaded.<br>
<b>B6.</b> Compile the standalone program code using the same settings as specified in section <b>B5</b><br>
<code>cd ./icono-pdaf/src/<br>
make PDAF_ARCH=Levante_mpiifort_impi<br>
make PDAF_ARCH=Levante_nvhpc_openmpi<br>
make PDAF_ARCH=Levante_nvhpc_gpu_openmpi<br>
cd ../../</code><br>
If the compilation is successful, an executable file named <code>icono_offline</code> will be created in the <code>./icono-pdaf/src/</code> directory.<br>
<b>B7.</b> Create a directory for running the experiment and move the executable file into it<br>
<code>mkdir ./icono-pdaf/runExecutable<br>
mv ./icono-pdaf/src/icono_offline ./icono-pdaf/runExecutable/</code><br>
<b>B8.</b> Copy the updated runscripts into the directory created in section <b>B7</b><br>
<code>cp ./natESM_PDAF2GPU_SPRINT/IconO-Offline-PDAFV2.3/icono-pdaf/runLevanteCompute/* ./icono-pdaf/runExecutable/<br>
cp ./natESM_PDAF2GPU_SPRINT/IconO-Offline-PDAFV2.3/icono-pdaf/runLevanteGpu/* ./icono-pdaf/runExecutable/</code><br>
<b>B9.</b> Launch the runscript that corresponds to the compilation settings<br>
<code>sbatch run_test_IntelCPU_N40.sh<br>
sbatch run_test_NvhpcCPU_N40.sh<br>
sbatch run_test_NvhpcGPU_N40.sh</code><br>
<i>Note:</i> Don't forget to update the "account" field in the SLURM configuration settings.<br>
In the latter case, a CPU executable will actually be built and after that run on the GPU partition of Levante. 
However, the added external libraries will have OpenACC directives enabled, and the CUDA libraries will be linked during compilation.<br>

<h2 style="text-align: center; font-size: 24pt">Enabling some useful program features in the CPU code</h2>

To enable one or more features, the corresponding source code files must be copied, followed 
by recompilation and program launch (<b>B5</b> - <b>B9</b>).<br>
<b>C1.</b> Displaying the active binding of MPI processes and OpenMP threads to CPU cores.
An option for manual binding of processes to cores can also be enabled here. Its implementation is included in the external module.
Additionally, heap memory allocation for temporary arrays can be activated via compiler directives<br>
<code>cp ./natESM_PDAF2GPU_SPRINT/IconO-Offline-PDAFV2.3-SFO/icono-pdaf/src/main_offline.F90 ./icono-pdaf/src/</code><br>
<b>C2.</b> Adding execution timers for main subroutines on each MPI-process<br>
<code>cp ./natESM_PDAF2GPU_SPRINT/IconO-Offline-PDAFV2.3-SFO/PDAF/src/PDAF_lestkf_update.F90 ./PDAF/src/</code><br>
<b>C3.</b> Substituting the CPU version of the <code>PDAFomi_check_dist2_noniso_loop</code> subroutine with an optimized variant<br>
<code>cp ./natESM_PDAF2GPU_SPRINT/IconO-Offline-PDAFV2.3-SFO/PDAF/src/cpuCodeOptimization/PDAFomi_dim_obs_l.F90 ./PDAF/src/</code><br>
<i>Note:</i> Performance improvement is observed only when using the Nvidia compiler.<br>

