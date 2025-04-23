<h1 style="text-align: center; color: blue;">Sprint #16: Performance evaluation, optimization, and strategy development for porting PDAF to GPUs</h1>

This GIT repository contains modified and added during the sprint source code files of the PDAF (The Parallel Data Assimilation Framework) and
the stand-alone offline assimilation program PDAF-ICONO (PDAF implementation used as assimilation method for …), as well as configuration files,
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

