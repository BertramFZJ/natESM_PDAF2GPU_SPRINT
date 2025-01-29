#ifndef _TALCC_HEADER_
#define _TALCC_HEADER_

#define _TALCC_USE_MPI_PLUGINS_

int talcc_getmasksize(void);

int talcc_getaffinitycore(int coreid);

void talcc_setaffinitycore(int coreid);

void talcc_setaffinitysinglecore(int coreid);

void talcc_unsetaffinitycore(int coreid);

void talcc_printaffinitymaskstringstdout(int mpiid, int ompid);

#ifdef _TALCC_USE_MPI_PLUGINS_
int talcc_getnummpinodesplugin(void);
void talcc_setprogramaffinityplugin(int numcorespernode, int numopenmpthreads);
#endif

#endif
