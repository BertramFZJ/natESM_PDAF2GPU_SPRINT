!>  Main driver for PDAF offline tutorial
!!
!! This is the main program for an example implementation of
!! PDAF with domain-decomposition and offline configuration.
!!
!! In the offline mode, we assume that the ensemble
!! integrations are performed in a separate program (model)
!! and the forecasted ensemble can be read from files. After
!! initializing the ensemble information by reading model
!! outputs, a single analysis step is performed. Subsequently,
!! the analysis ensemble can be written to files that can be 
!! used to initialize another ensemble forecast.
!!
!! Using PDAF for domain-decomposition, the offline
!! mode can be used to perform assimilation with domain-
!! decomposed models. If the models write results for each 
!! sub-domain, these can be read here using the same 
!! parallelization. Then, the filter analysis can be 
!! performed utilizing this parallelization. If the files
!! contain the full model state, PDAF in offline mode
!! can be used either on a single processor, or the 
!! fields can be distributed in this program to utilize
!! the parallelization of the filters.
!!
!! Parameters can be set in the code, or - preferably -
!! by command line arguments that are parsed by the 
!! module PARSER. The format for this is
!! EXECUTABLE -HANDLE1 VALUE1 -HANDLE2 VALUE2 ...
!! The handles are defined in the code before the calls
!! to the routine PARSE.
!!
!! __Revision history:__
!! * 2008-07 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
!! adapted for mpiom pdaf in offline mode;
!! nothing done
!! S.Brune, 12.11.2012
!!
PROGRAM MAIN

! USE mpi
  USE omp_lib
#ifdef _OPENACC
  use openacc
  use cudafor
#endif

  USE mod_parallel, &     ! Parallelization variables
       ONLY: MPIerr, npes_world, mype_world, MPI_COMM_WORLD, &
       init_parallel, finalize_parallel
  USE timer, &            ! Timing
       ONLY: timeit, time_tot
  USE mod_memcount, &     ! Counting allocated memory
       ONLY: memcount, memcount_ini, memcount_get

  USE talcc_c_bindings, only: printAffinityBindingC, setAffinityBindingC
  USE mo_mem_worspaces, only: memInitLibrary, memSetOmpNumThreads, memAllocateAnalysisHeap

  IMPLICIT NONE

! Local variables
  INTEGER :: i                 ! Counter
  
  INTEGER :: numCores = 128 


! **********************
! *** Initialize MPI ***
! **********************

  CALL init_parallel() ! initializes MPI

! ***************************
! *** Initialize AFFINITY ***
! ***************************

#if 1
  !$OMP PARALLEL SHARED(mype_world)
  call printAffinityBindingC(mype_world, omp_get_thread_num())
  !$OMP END PARALLEL

  IF( .FALSE. ) THEN
    call setAffinityBindingC(numCores, -1)

    ! !$OMP PARALLEL SHARED(mype_world)
    ! call printAffinityBindingC(mype_world, omp_get_thread_num())
    ! !$OMP END PARALLEL
  ENDIF
#endif  

#ifdef _OPENACC
  IF (mype_world == 0) THEN
    WRITE(*,*) 'NATESM: OPENACC ENABLED ==> ', _OPENACC, ' DEVICE COUNT = ', acc_get_num_devices(acc_device_nvidia)
  ENDIF
#endif

#if 1
  CALL memInitLibrary()
#ifdef _OPENMP
  !$OMP PARALLEL
  !$OMP MASTER
  CALL memSetOmpNumThreads( OMP_GET_NUM_THREADS() )
  !$OMP END MASTER
  !$OMP END PARALLEL
#else
  CALL memSetOmpNumThreads( 1 )
#endif
  CALL memAllocateAnalysisHeap( 65536, numaInitFlag = .TRUE., accEnterDataFlag = .TRUE.)
#endif

! ********************************
! ***      INITIALIZATION      ***
! ********************************

! *** Initial Screen output ***
  initscreen: IF (mype_world == 0) THEN

     WRITE (*, '(/8x, a/)') '+++++ PDAF for ICON-O  - offline mode +++++'
     WRITE (*, '(16x, a)') 'Data assimilation with PDAF'

     IF (npes_world > 1) THEN
        WRITE (*, '(/21x, a, i3, a/)') 'Running on ', npes_world, ' PEs'
     ELSE
        WRITE (*, '(/21x, a/)') 'Running on 1 PE'
     END IF
     
  END IF initscreen

! *** set number of timers ***
  CALL timeit(6, 'ini')

! *** set first timer ***
  CALL timeit(1, 'new')

! *** set number of memory counters ***
  CALL memcount_ini(4)

  
! *** Initialize MPI communicators for PDAF (model and filter) ***
! *** NOTE: It is always n_modeltasks=1 for offline mode       ***
  CALL init_parallel_pdaf(0, 1)

! *** Initialize model information ***
! *** This should only be information on the model dimension
! *** Generally, this could be joined with init_pdaf.

  CALL timeit(2, 'new')
  CALL initialize()
  CALL timeit(2, 'old')


! *******************************
! ***      ASSIMILATION       ***
! *******************************

  ! *** Initialize PDAF ***

  CALL timeit(4, 'new')
  CALL init_pdaf()
  CALL timeit(4, 'old')

  ! *** Perform analysis ***

  CALL timeit(3, 'new')
  IF (mype_world == 0) &
       WRITE (*, '(/2x, a)') 'PDAF for ICON-O - offline mode: START ASSIMILATION'
  CALL assimilation_pdaf_offline()

  ! Syncronize at barrier for exit
  CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr) 
  !WRITE (*,*) 'model PE exited: mype ', mype_world

  CALL timeit(3, 'old')


! ********************
! *** Finishing up ***
! ********************

  CALL timeit(1, 'old')

! *** Final screen output ***
  screen3: IF (mype_world == 0) THEN
     WRITE (*, '(/1x, a)') 'PDAF for ICON-O - offline mode: EXITED ASSIMILATION'

     ! *** Show allocated memory for the model ***
     WRITE (*, '(/18x, a)') 'Model - Memory overview'
     WRITE (*, '(10x, 45a)') ('-', i=1, 45)
     WRITE (*, '(21x, a, f10.3, a)') 'Allocated memory  (MB)'
     WRITE (*, '(14x, a, f10.3, a)') &
          'Model field:', memcount_get(1, 'M'), ' MB (persistent)'
     WRITE (*, '(12x, a, f10.3, a)') &
          'ensemble init:', memcount_get(2, 'M'), ' MB (temporary)'
     WRITE (*, '(13x, a, f10.3, a)') &
          'Pre-Poststep:', memcount_get(3, 'M'), ' MB (temporary)'
     WRITE (*, '(13x, a, f10.3, a)') &
          'Observations:', memcount_get(4, 'M'), ' MB (temporary)'

     ! Show allocated memory for PDAF
     CALL PDAF_print_info(10)

     ! *** Print timings onto screen ***

     ! Show timings for PDAF
!     CALL PDAF_print_info(1)  ! normal info
     CALL PDAF_print_info(3)  ! info on call-back routines
     CALL PDAF_print_info(5)  ! Detailed info

     WRITE (*, '(/17x, a)') 'Offline - Timing information'
     WRITE (*, '(10x, 45a)') ('-', i=1, 45)
     ! Timing summary for assimilation
     WRITE (*, '(19x, a, F11.3, 1x, a)') 'initialize model:', time_tot(2), 's'
     WRITE (*, '(18x, a, F11.3, 1x, a)') 'initialize filter:', time_tot(4), 's'
     WRITE (*, '(23x, a, F11.3, 1x, a)') 'assimilation:', time_tot(3), 's'
     WRITE (*, '(19x, a, F11.3, 1x, a)') 'total run time:', time_tot(1), 's'

     WRITE (*, '(/1x, a)') 'PDAF for ICON-O - offline mode: END'
  END IF screen3

! *** deallocate timers ***
  CALL timeit(6, 'fin')

!! *** deallocate PDAF arrays ***
  CALL PDAF_deallocate()

! *** Terminate MPI
  CALL finalize_parallel()

END PROGRAM MAIN
