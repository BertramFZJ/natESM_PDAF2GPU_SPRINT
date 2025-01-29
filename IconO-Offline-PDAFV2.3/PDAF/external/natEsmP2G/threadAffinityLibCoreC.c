#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sched.h>
#include <sys/types.h>

#include "threadAffinityLibCoreC.h"

#define ERR_EXIT_LIB_MACRO exit(0)

int talcc_getmasksize(void)
{				
    int maskSize = (int)CPU_SETSIZE;   

    return maskSize;
}

int talcc_getaffinitycore(int coreid)
{
    cpu_set_t mask;			
    int ret = -1;
    int get;

    get = sched_getaffinity(0, sizeof(cpu_set_t), &mask);
    if(get != 0) {
        fprintf(stderr, "talcc_getaffinitycore: sched_getaffinity returned an error code %d\n", get); fflush(stderr);
        ERR_EXIT_LIB_MACRO; }

    if( (coreid < 0) || (coreid >= CPU_SETSIZE) ) {
        fprintf(stderr, "talcc_getaffinitycore: invalid coreid = %d value\n", coreid); fflush(stderr);
        ERR_EXIT_LIB_MACRO; }

    if(CPU_ISSET(coreid, &mask) != 0) ret = 1;
                                 else ret = 0;

    return ret;
}

void talcc_setaffinitycore(int coreid)
{
    cpu_set_t mask;
    int get;

    get = sched_getaffinity(0, sizeof(cpu_set_t), &mask);
    if(get != 0) {
        fprintf(stderr, "talcc_setaffinitycore: sched_getaffinity returned an error code %d\n", get); fflush(stderr);
        ERR_EXIT_LIB_MACRO; }

    if( (coreid < 0) || (coreid >= CPU_SETSIZE) ) {
        fprintf(stderr, "talcc_setaffinitycore: invalid coreid = %d value\n", coreid); fflush(stderr);
        ERR_EXIT_LIB_MACRO; }

    if(CPU_ISSET(coreid, &mask) == 0)
    {
        CPU_SET(coreid, &mask);
		get = sched_setaffinity(0, sizeof(cpu_set_t), &mask);
        if(get != 0) {
            fprintf(stderr, "talcc_setaffinitycore: sched_setaffinity returned an error code %d\n", get); fflush(stderr);
            ERR_EXIT_LIB_MACRO; }

        CPU_ZERO(&mask);
        get = sched_getaffinity(0, sizeof(cpu_set_t), &mask);
        if(get != 0) {
            fprintf(stderr, "talcc_setaffinitycore: sched_getaffinity returned an error code %d\n", get); fflush(stderr);
            ERR_EXIT_LIB_MACRO; }

        if(CPU_ISSET(coreid, &mask) == 0) {
            fprintf(stderr, "talcc_setaffinitycore: failed to enable core %d\n", coreid); fflush(stderr);
            ERR_EXIT_LIB_MACRO; }
    }

}

void talcc_setaffinitysinglecore(int coreid)
{
    cpu_set_t mask;
    int i, get;

    if( (coreid < 0) || (coreid >= CPU_SETSIZE) ) {
        fprintf(stderr, "talcc_setaffinitysinglecore: invalid coreid = %d value\n", coreid); fflush(stderr);
        ERR_EXIT_LIB_MACRO; }

    CPU_ZERO(&mask);
    CPU_SET(coreid, &mask);

    get = sched_setaffinity(0, sizeof(cpu_set_t), &mask);
    if(get != 0) {
        fprintf(stderr, "talcc_setaffinitysinglecore: sched_setaffinity returned an error code %d\n", get); fflush(stderr);
        ERR_EXIT_LIB_MACRO; }

    CPU_ZERO(&mask);
    get = sched_getaffinity(0, sizeof(cpu_set_t), &mask);
    if(get != 0) {
        fprintf(stderr, "talcc_setaffinitysinglecore: sched_getaffinity returned an error code %d\n", get); fflush(stderr);
        ERR_EXIT_LIB_MACRO; }

    for(i=0; i<CPU_SETSIZE; i++) {
        if( ( (i != coreid) && (CPU_ISSET(i, &mask) != 0) ) || ( (i == coreid) && (CPU_ISSET(i, &mask) == 0) ) ) {
            fprintf(stderr, "talcc_setaffinitysinglecore: failed to set single core %d mask\n", coreid); fflush(stderr);
            ERR_EXIT_LIB_MACRO; }        
    }

}

void talcc_unsetaffinitycore(int coreid)
{
    cpu_set_t mask;
    int get;

    get = sched_getaffinity(0, sizeof(cpu_set_t), &mask);
    if(get != 0) {
        fprintf(stderr, "talcc_unsetaffinitycore: sched_getaffinity returned an error code %d\n", get); fflush(stderr);
        ERR_EXIT_LIB_MACRO; }

    if( (coreid < 0) || (coreid >= CPU_SETSIZE) ) {
        fprintf(stderr, "talcc_unsetaffinitycore: invalid coreid = %d value\n", coreid); fflush(stderr);
        ERR_EXIT_LIB_MACRO; }

    if(CPU_ISSET(coreid, &mask) != 0)
    {
        CPU_CLR(coreid, &mask);
		get = sched_setaffinity(0, sizeof(cpu_set_t), &mask);
        if(get != 0) {
            fprintf(stderr, "talcc_unsetaffinitycore: sched_setaffinity returned an error code %d\n", get); fflush(stderr);
            ERR_EXIT_LIB_MACRO; }

        CPU_ZERO(&mask);
        get = sched_getaffinity(0, sizeof(cpu_set_t), &mask);
        if(get != 0) {
            fprintf(stderr, "talcc_unsetaffinitycore: sched_getaffinity returned an error code %d\n", get); fflush(stderr);
            ERR_EXIT_LIB_MACRO; }

        if(CPU_ISSET(coreid, &mask) != 0) {
            fprintf(stderr, "talcc_unsetaffinitycore: failed to enable core %d\n", coreid); fflush(stderr);
            ERR_EXIT_LIB_MACRO; }
    }

}

#define TALCC_MAX_STRING_LENGTH 2048
void talcc_printaffinitymaskstringstdout(int mpiid, int ompid)
{
    cpu_set_t mask;
    char stdString[TALCC_MAX_STRING_LENGTH];
    int pos, activeCoreNum, begin, end;
    int get, i;

    for(i=0; i<TALCC_MAX_STRING_LENGTH; i++) stdString[i] = '\0';
    
    get = sched_getaffinity(0, sizeof(cpu_set_t), &mask);
    if(get != 0) {
        fprintf(stderr, "talcc_printaffinitymaskstringstdout: sched_getaffinity returned an error code %d\n", get); fflush(stderr);
        ERR_EXIT_LIB_MACRO; }

    if( (mpiid <  0) && (ompid <  0) ) sprintf(stdString, "MASK: ");
    if( (mpiid >= 0) && (ompid <  0) ) sprintf(stdString, "MPIID %6d MASK: ", mpiid);
    if( (mpiid >= 0) && (ompid >= 0) ) sprintf(stdString, "MPIID %6d OMPTHREAD %3d MASK: ", mpiid, ompid);

    pos = activeCoreNum = 0;
    while(pos < CPU_SETSIZE) {

        begin = end = -1;

        for(; pos<CPU_SETSIZE; pos++) {
            if(CPU_ISSET(pos, &mask) != 0) {
                begin = end = pos;
                break;
            } // if
        } // for

        if(begin == -1) break;

        for(; end+1<CPU_SETSIZE; end++)
            if(CPU_ISSET(end+1, &mask) == 0)
                break;

        pos = end + 1;

        if(activeCoreNum != 0) sprintf(stdString + strlen(stdString), ", ");
        if(begin != end) sprintf(stdString + strlen(stdString), "%d-%d", begin, end);
                    else sprintf(stdString + strlen(stdString), "%d", begin);

        activeCoreNum += end - begin + 1;
    } // while
    
    if(activeCoreNum == 0) sprintf(stdString + strlen(stdString), "list is empty");
                      else sprintf(stdString + strlen(stdString), " => total %d core(s)", activeCoreNum);
    sprintf(stdString + strlen(stdString), "\n");
    printf("%s", stdString); fflush(stdout);

}
#undef TALCC_MAX_STRING_LENGTH

#ifdef _TALCC_USE_MPI_PLUGINS_
#include <mpi.h>
#include <omp.h>

int talcc_getnummpinodesplugin(void)
{
    char *listOfNodeNames = NULL;

    int mpiRank, mpiSize;
    int numNodes = 0;

    int ret, i, j;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    if(mpiRank == 0) {
        listOfNodeNames = (char *)malloc(mpiSize * MPI_MAX_PROCESSOR_NAME * sizeof(char));
        if(listOfNodeNames == NULL) {
            fprintf(stderr, "talcc_getnummpinodesplugin: can't allocate memory\n"); fflush(stderr);
            ERR_EXIT_LIB_MACRO; }
        for(i=0; i<mpiSize * MPI_MAX_PROCESSOR_NAME; i++) listOfNodeNames[i] = '\n';

        MPI_Get_processor_name(listOfNodeNames, &ret);

        for(i=1; i<mpiSize; i++) {
            MPI_Recv(listOfNodeNames + i*MPI_MAX_PROCESSOR_NAME, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE); }

        numNodes = 1;
        for(i=1; i<mpiSize; i++) {

            for(j=0; j<numNodes; j++) {
                ret = strcmp(listOfNodeNames + j*MPI_MAX_PROCESSOR_NAME, listOfNodeNames + i*MPI_MAX_PROCESSOR_NAME);
                if(ret == 0) break;
            } // for j

            if(j == numNodes) {
                if(i != j) memcpy(listOfNodeNames + j*MPI_MAX_PROCESSOR_NAME, 
                                  listOfNodeNames + i*MPI_MAX_PROCESSOR_NAME, 
                                  MPI_MAX_PROCESSOR_NAME * sizeof(char));
                numNodes += 1;
            } // if
        } // for i

#if 1
        printf("******** List of nodes ********\n");
        for(i=0; i<numNodes; i++)
            printf("NODE %5d => NAME: %s\n", i, listOfNodeNames + i*MPI_MAX_PROCESSOR_NAME);
        printf("*******************************\n"); fflush(stdout);
#endif

        free(listOfNodeNames);
        MPI_Bcast(&numNodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    else {
        listOfNodeNames = (char *)malloc(          MPI_MAX_PROCESSOR_NAME * sizeof(char));
        if(listOfNodeNames == NULL) {
            fprintf(stderr, "talcc_getnummpinodesplugin: can't allocate memory\n"); fflush(stderr);
            ERR_EXIT_LIB_MACRO; }
        for(i=0; i<          MPI_MAX_PROCESSOR_NAME; i++) listOfNodeNames[i] = '\n';

        MPI_Get_processor_name(listOfNodeNames, &ret);

        MPI_Send(listOfNodeNames, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, mpiRank, MPI_COMM_WORLD);

        free(listOfNodeNames);
        MPI_Bcast(&numNodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    return numNodes;
}

void talcc_setprogramaffinityplugin(int numcorespernode, int numopenmpthreads)
{
    int mpiRank, mpiSize;
    int numNodes, nOmpThreads;
    int numMpiProcPerNode, numCoresPerMpiProc, numCoresPerOpenMPThread;
    int i, get, coreId = -1;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    MPI_Barrier(MPI_COMM_WORLD);

    if(numopenmpthreads > 0) nOmpThreads = numopenmpthreads;
    else {
        #pragma omp parallel shared(nOmpThreads)
        {
            #pragma omp single
            nOmpThreads = omp_get_num_threads();
            #pragma omp barrier
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    numNodes = talcc_getnummpinodesplugin();
    numMpiProcPerNode = mpiSize / numNodes;
    numCoresPerMpiProc = numcorespernode / numMpiProcPerNode;
    numCoresPerOpenMPThread = numCoresPerMpiProc / nOmpThreads;
    MPI_Barrier(MPI_COMM_WORLD);

#if 1
    if(mpiRank == 0) {
        printf("**************** CONFIGURATION ****************\n");
        printf("                 Number of mpi processes: %d\n", mpiSize);
        printf("                     Number of HPC nodes: %d\n", numNodes);
        printf("                Number of cores per node: %d\n", numcorespernode);
        printf("          Number of mpi process per node: %d\n", numMpiProcPerNode);
        printf("         Number of cores per mpi process: %d\n", numCoresPerMpiProc);
        printf("Number of OpenMP threads per mpi process: %d\n", nOmpThreads);
        printf("       Number of cores per OpenMP thread: %d\n", numCoresPerOpenMPThread);
        printf("***********************************************\n"); fflush(stdout);
    } // if
    MPI_Barrier(MPI_COMM_WORLD);

    get = 0;
    if(numMpiProcPerNode > numcorespernode) get = 1;
    if(numCoresPerMpiProc < nOmpThreads) get = 1;
    MPI_Barrier(MPI_COMM_WORLD);

    if(get == 1) { fprintf(stderr, "[%d] talcc_setprogramaffinityplugin: wrong configuration parameters\n", mpiRank); fflush(stderr); }
    MPI_Barrier(MPI_COMM_WORLD);
    if(get == 1) ERR_EXIT_LIB_MACRO;
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    #pragma omp parallel num_threads(nOmpThreads) private(i, coreId)
    {
#if 1
        coreId  = (mpiRank % numMpiProcPerNode) * numCoresPerMpiProc;
        coreId += omp_get_thread_num() * numCoresPerOpenMPThread;
        for(i=0; i<numCoresPerOpenMPThread; i++) {
            if(i == 0) talcc_setaffinitysinglecore(coreId);
                  else talcc_setaffinitycore(coreId + i);
        } // for i
#else
        // coreId = 16 + 32 * mpiRank + omp_get_thread_num();
        coreId = 32 * mpiRank + 2 * omp_get_thread_num();
        talcc_setaffinitysinglecore(coreId);        
#endif
    }
    MPI_Barrier(MPI_COMM_WORLD);

    #pragma omp parallel num_threads(nOmpThreads)
    {
        talcc_printaffinitymaskstringstdout(mpiRank, omp_get_thread_num()); 
    }
    MPI_Barrier(MPI_COMM_WORLD);

}
#endif
