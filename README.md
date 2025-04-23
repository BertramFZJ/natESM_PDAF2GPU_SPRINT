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

