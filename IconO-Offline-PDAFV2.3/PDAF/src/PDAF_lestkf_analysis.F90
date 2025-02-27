! #define _RSE_USE_CUBLAS_DGEMM
#define _RSE_ALLOCATE_SUBROUTINE_PARAMETERS_

#define _RSE_CHECK_DnDsyevd_BUFFER_SIZE
#ifndef _RSE_CHECK_DnDsyevd_BUFFER_SIZE
#define _RSE_DnDsyevd_BUFFER_SIZE 39616
#endif

! Copyright (c) 2004-2024 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id: PDAF_lestkf_analysis.F90 1147 2023-03-12 16:14:34Z lnerger $
!BOP
!
! !ROUTINE: PDAF_lestkf_analysis --- LESTKF analysis/transformation
!
! !INTERFACE:
SUBROUTINE PDAF_lestkf_analysis(domain_p, step, dim_l, dim_obs_f, dim_obs_l, &
     dim_ens, rank, state_l, Ainv_l, ens_l, HX_f, &
     HXbar_f, state_inc_l, OmegaT_in, forget, U_g2l_obs, &
     U_init_obs_l, U_prodRinvA_l, U_init_obsvar_l, U_init_n_domains_p, &
     screen, incremental, type_forget, type_sqrt, TA, flag)

! !DESCRIPTION:
! Analysis step of the LESTKF filter with direct
! transformation of the forecast into the 
! analysis ensemble. This variant does not
! compute the analysis state, but only the
! analysis ensemble, whose mean is the analysis
! state.
! Supported is also the adaptive forgetting factor.
!
! Variant for domain decomposed states.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2011-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype
  USE PDAF_mod_filter, &
       ONLY: obs_member, debug
  USE PDAFomi, &
       ONLY: omi_omit_obs => omit_obs
#if defined (_OPENMP)
  USE omp_lib, &
       ONLY: omp_get_num_threads, omp_get_thread_num
#endif

#ifdef _OPENACC
  use openacc
  use cudafor
  use cublas_v2
  use cusolverDn
#endif

  IMPLICIT NONE

! !ARGUMENTS:
! ! Variable naming scheme:
! !   suffix _p: Denotes a full variable on the PE-local domain
! !   suffix _l: Denotes a local variable on the current analysis domain
  INTEGER, INTENT(in) :: domain_p    ! Current local analysis domain
  INTEGER, INTENT(in) :: step        ! Current time step
  INTEGER, INTENT(in) :: dim_l       ! State dimension on local analysis domain
  INTEGER, INTENT(in) :: dim_obs_f   ! PE-local dimension of full observation vector
  INTEGER, INTENT(in) :: dim_obs_l   ! Size of obs. vector on local ana. domain
  INTEGER, INTENT(in) :: dim_ens     ! Size of ensemble 
  INTEGER, INTENT(in) :: rank        ! Rank of initial covariance matrix
  REAL, INTENT(inout) :: state_l(dim_l)        ! on exit: state on local analysis domain (RSE: INTENT(IN) ?)
  REAL, INTENT(inout) :: Ainv_l(rank, rank)    ! Inverse of matrix U - temporary use only
  REAL, INTENT(inout) :: ens_l(dim_l, dim_ens) ! Local state ensemble
  REAL, INTENT(in) :: HX_f(dim_obs_f, dim_ens) ! PE-local full observed state ens.
  REAL, INTENT(in) :: HXbar_f(dim_obs_f)       ! PE-local full observed ens. mean
  REAL, INTENT(in) :: state_inc_l(dim_l)       ! Local state increment
  REAL, INTENT(in) :: OmegaT_in(rank, dim_ens) ! Matrix Omega
  REAL, INTENT(inout) :: forget      ! Forgetting factor
  INTEGER, INTENT(in) :: screen      ! Verbosity flag
  INTEGER, INTENT(in) :: incremental ! Control incremental updating
  INTEGER, INTENT(in) :: type_forget ! Type of forgetting factor
  INTEGER, INTENT(in) :: type_sqrt   ! Type of square-root of A
                                     ! (0): symmetric sqrt; (1): Cholesky decomposition
  REAL, INTENT(inout) :: TA(dim_ens, dim_ens)  ! Ensemble transformation matrix
  INTEGER, INTENT(inout) :: flag     ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_g2l_obs, &   ! Restrict full obs. vector to local analysis domain
       U_init_obs_l, &       ! Init. observation vector on local analysis domain
       U_init_obsvar_l, &    ! Initialize local mean observation error variance
       U_init_n_domains_p, & ! Provide PE-local number of local analysis domains
       U_prodRinvA_l         ! Provide product R^-1 A for local analysis domain

! !CALLING SEQUENCE:
! Called by: PDAF_lestkf_update
! Calls: U_g2l_obs
! Calls: U_init_obs_l
! Calls: U_prodRinvA_l
! Calls: PDAF_set_forget_local
! Calls: PDAF_estkf_AOmega
! Calls: PDAF_estkf_OmegaA
! Calls: gemmTYPE (BLAS; dgemm or sgemm dependent on precision)
! Calls: gemvTYPE (BLAS; dgemv or sgemv dependent on precision)
! Calls: gesvTYPE (LAPACK; dgesv or sgesv dependent on precision)
! Calls: syevTYPE (LAPACK; dsyev or ssyev dependent on precision)
! Calls: potrfTYPE (LAPACK; dpotrf or spotrf dependent on precision)
! Calls: trtrsTYPE (LAPACK; dtrtrs or strtrs dependent on precision)
!EOP
       
! *** local variables ***
  INTEGER :: i, j, member, col, row    ! Counters
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER :: lib_info                  ! Status flag for LAPACK calls
  INTEGER :: ldwork                    ! Size of work array for SYEVTYPE
  INTEGER :: maxblksize, blkupper, blklower  ! Variables for blocked ensemble update
  REAL    :: fac                       ! Temporary variable sqrt(dim_ens) or sqrt(rank)
  INTEGER, SAVE :: lastdomain = -1     ! store domain index
  LOGICAL, SAVE :: screenout = .true.  ! Whether to print information to stdout
  REAL, ALLOCATABLE :: HL_l(:,:)       ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: RiHL_l(:,:)     ! Temporary matrices for analysis
  REAL, ALLOCATABLE :: resid_l(:)      ! observation residual
  REAL, ALLOCATABLE :: obs_l(:)        ! local observation vector
  REAL, ALLOCATABLE :: HXbar_l(:)      ! state projected onto obs. space
  REAL, ALLOCATABLE :: RiHLd_l(:)      ! local RiHLd
  REAL, ALLOCATABLE :: VRiHLd_l(:)     ! Temporary vector for analysis
  REAL, ALLOCATABLE :: tmp_Ainv_l(:,:) ! Temporary storage of Ainv
  REAL, ALLOCATABLE :: OmegaT(:,:)   ! Transpose of Omega
  REAL, ALLOCATABLE :: ens_blk(:,:)  ! Temporary blocked state ensemble
  REAL, ALLOCATABLE :: svals(:)      ! Singular values of Ainv
  REAL, ALLOCATABLE :: work(:)       ! Work array for syevTYPE
  INTEGER, ALLOCATABLE :: ipiv(:)    ! vector of pivot indices for GESVTYPE
  INTEGER, SAVE :: mythread, nthreads  ! Thread variables for OpenMP
  REAL :: state_inc_l_dummy(1)       ! Dummy variable to avoid compiler warning

#ifdef _OPENACC
  TYPE(cublasHandle), SAVE :: handleCuBlas
  TYPE(cublasHandle), SAVE :: handleCuBlasDevice
  type(cusolverDnHandle), SAVE :: handleCuSolverDn
  INTEGER(kind=cuda_stream_kind), SAVE :: cudaStream
  INTEGER, SAVE :: accStreamId = -1
  INTEGER, SAVE :: handleCuBlasFlag = 0
  INTEGER :: ierror, devInfo, cuBlasHandleMode

  INTEGER :: iRow, jCol, kCell, iBlk
  INTEGER, SAVE :: accMemRefCounter = 0
#else
  INTEGER :: iBlk
#endif

!$OMP THREADPRIVATE(mythread, nthreads, lastdomain, allocflag, screenout)
#ifdef _OPENACC
!$OMP THREADPRIVATE(accStreamId, handleCuBlas, handleCuBlasFlag, handleCuSolverDn, cudaStream)
!$OMP THREADPRIVATE(accMemRefCounter)
#endif

! ***********************
! *** GPU Preparation ***
! ***********************
#ifdef _OPENACC

  IF(handleCuBlasFlag .NE. 1) THEN

  ! >>>>> !$OMP CRITICAL INITGPUIS

    ierror = cudaSetDevice( 0 )
    IF(ierror .NE. 0) THEN
      WRITE(0,*) "cudaSetDevice returned err code ", ierror
      STOP
    ENDIF

    ierror = cublasCreate( handleCuBlas )
    IF(ierror .NE. 0) THEN
      WRITE(0,*) "cublasCreate returned err code ", ierror
      STOP
    ENDIF

    ierror = cublasCreate( handleCuBlasDevice )
    IF(ierror .NE. 0) THEN
      WRITE(0,*) "cublasCreate returned err code ", ierror
      STOP
    ENDIF

    ierror = cusolverDnCreate( handleCuSolverDn )
    IF(ierror .NE. 0) THEN
      WRITE(0,*) "cusolverDnCreate returned err code ", ierror
      STOP
    ENDIF

#if(1)
    accStreamId = omp_get_thread_num() + 1
    ierror = cudaStreamCreate( cudaStream )
    IF(ierror .NE. 0) THEN
      WRITE(0,*) "cudaStreamCreate returned err code ", ierror
      STOP
    ENDIF
    CALL acc_set_cuda_stream( accStreamId, cudaStream )
#else
    ! accStreamId = MOD(omp_get_thread_num(), 4) + 1
    accStreamId = 1
    ierror = cudaStreamCreate( cudaStream )
    IF(ierror .NE. 0) THEN
      WRITE(0,*) "cudaStreamCreate returned err code ", ierror
      STOP
    ENDIF
#endif

    ierror = cublasSetStream( handleCuBlas, cudaStream )
    IF(ierror .NE. 0) THEN
      WRITE(0,*) "cublasSetStream returned err code ", ierror
      STOP
    ENDIF

    ierror = cublasSetStream( handleCuBlasDevice, cudaStream )
    IF(ierror .NE. 0) THEN
      WRITE(0,*) "cublasSetStream returned err code ", ierror
      STOP
    ENDIF

    ierror = cusolverDnSetStream( handleCuSolverDn, cudaStream )
    IF(ierror .NE. 0) THEN
      WRITE(0,*) "cusolverDnSetStream returned err code ", ierror
      STOP
    ENDIF

    ierror = cublasSetPointerMode( handleCuBlasDevice, CUBLAS_POINTER_MODE_DEVICE)
    IF(ierror .NE. 0) THEN
      WRITE(0,*) "cublasSetPointerMode returned err code ", ierror
      STOP
    ENDIF

    ierror = cublasGetPointerMode( handleCuBlas, cuBlasHandleMode)
    IF(ierror .NE. 0) THEN
      WRITE(0,*) "cublasGetPointerMode #1 returned err code ", ierror
      STOP
    ENDIF
    IF(cuBlasHandleMode .NE. CUBLAS_POINTER_MODE_HOST) THEN
      WRITE(0,*) "WRONG handleCuBlas mode"
      STOP
    ENDIF

    ierror = cublasGetPointerMode( handleCuBlasDevice, cuBlasHandleMode)
    IF(ierror .NE. 0) THEN
      WRITE(0,*) "cublasGetPointerMode #2 returned err code ", ierror
      STOP
    ENDIF
    IF(cuBlasHandleMode .NE. CUBLAS_POINTER_MODE_DEVICE) THEN
      WRITE(0,*) "WRONG handleCuBlasDevice mode"
      STOP
    ENDIF
    
    handleCuBlasFlag = 1

  ! >>>>> !$OMP END CRITICAL INITGPUIS

  ENDIF

#endif

  ! Allocating GPU memory for arrays passed to the subroutine as parameters
#ifdef _OPENACC
  IF(accMemRefCounter .NE. 0) THEN
    STOP "RSE: ERROR:: ACC MEMORY REFERENCE COUNTER = 1"
  ENDIF
  accMemRefCounter = 1

#ifdef _RSE_ALLOCATE_SUBROUTINE_PARAMETERS_
  !$ACC ENTER DATA CREATE(state_l(:), Ainv_l(:,:), ens_l(:,:), TA(:,:)) &
  !$ACC            CREATE(OmegaT_in(:,:)) ASYNC(accStreamId)
  ! Currently not used on GPU :: HX_f(:,:) + HXbar_f(:) + state_inc_l(:) ==> IN
  !$ACC UPDATE DEVICE(OmegaT_in(:,:)) ASYNC(accStreamId)
  !$ACC WAIT(accStreamId)
#endif
#endif
  ! Allocating GPU memory for arrays passed to the subroutine as parameters

! *******************
! *** Preparation ***
! *******************

  CALL PDAF_timeit(51, 'new')

#if defined (_OPENMP)
  nthreads = omp_get_num_threads()
  mythread = omp_get_thread_num()
#else
  nthreads = 1
  mythread = 0
#endif

  ! Initialize variable to prevent compiler warning
  state_inc_l_dummy(1) = state_inc_l(1)

  ! Control screen output
  IF (lastdomain<domain_p .AND. lastdomain>-1) THEN
     screenout = .false.
  ELSE
     screenout = .true.

     ! In case of OpenMP, let only thread 0 write output to the screen
     IF (mythread>0) screenout = .false.

     ! Output, only in case of OpenMP parallelization
#if defined (_OPENMP)
     IF (screenout) THEN
        WRITE (*,'(a, 5x, a, i5, a)') &
             'PDAF','--- Use OpenMP parallelization with ', nthreads, ' threads'
     END IF
#endif
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lestkf_analysis -- START'

  CALL PDAF_timeit(51, 'old')


! ************************
! *** Compute residual ***
! ***   d = y - H x    ***
! ************************

  CALL PDAF_timeit(12, 'new')

  haveobsB: IF (dim_obs_l > 0) THEN
     ! *** The residual only exists for domains with observations ***

     ! WRITE(0,*) 'RSE: haveobsB #1'

     ALLOCATE(resid_l(dim_obs_l), obs_l(dim_obs_l), HXbar_l(dim_obs_l))
     !$ACC ENTER DATA CREATE(resid_l(:), obs_l(:), HXbar_l(:)) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)

     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 3 * dim_obs_l)

     ! Restrict mean obs. state onto local observation space
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lestkf_analysis -- call g2l_obs for mean'

     CALL PDAF_timeit(46, 'new')
     obs_member = 0
     ! WRITE(0,*) '**************** RSE: CALL U_g2l_obs #1'
     CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HXbar_f, HXbar_l)
     ! U_g2l_obs == PDAFomi_g2l_obs_cb [PDAFomi_callback.F90]
     !           L1 PDAFomi_g2l_obs [PDAFomi_obs_l.F90]
     !           L2 PDAFomi_g2l_obs_internal [PDAFomi_obs_l.F90]
     ! WRITE(0,*) '**************** RSE: EXIT U_g2l_obs #1'
     CALL PDAF_timeit(46, 'old')

     ! get local observation vector
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lestkf_analysis -- call init_obs_l'

     CALL PDAF_timeit(47, 'new')
     ! WRITE(0,*) '**************** RSE: CALL U_init_obs_l #1'
     CALL U_init_obs_l(domain_p, step, dim_obs_l, obs_l)
     ! U_init_obs_l == PDAFomi_init_obs_l_cb [PDAFomi_callback.F90]
     !              L1 PDAFomi_init_obs_l [PDAFomi_obs_l.F90]
     !              L2 PDAFomi_g2l_obs_internal [PDAFomi_obs_l.F90]
     ! ?? Internal memory allocation for global variables
     ! WRITE(0,*) '**************** RSE: EXIT U_init_obs_l #1'
     CALL PDAF_timeit(47, 'old')

     ! Get residual as difference of observation and observed state
     CALL PDAF_timeit(51, 'new')
#ifndef _OPENACC

     resid_l = obs_l - HXbar_l

#else

     !$ACC UPDATE DEVICE(obs_l(:), HXbar_l(:)) ASYNC(accStreamId)
     !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(accStreamId)
     DO kCell = 1, dim_obs_l
        resid_l(kCell) = obs_l(kCell) - HXbar_l(kCell)
     ENDDO
     !$ACC END PARALLEL LOOP     
     !$ACC WAIT(accStreamId)

#endif
     CALL PDAF_timeit(51, 'old')

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  innovation d_l', resid_l

  END IF haveobsB

  ! Omit observations with too high innovation
  IF (omi_omit_obs .and. incremental /= 2) THEN
    WRITE(0,*) '**************** RSE: CALL PDAFomi_omit_by_inno_l_cb #1'
    CALL PDAFomi_omit_by_inno_l_cb(domain_p, dim_obs_l, resid_l, obs_l)
    WRITE(0,*) '**************** RSE: EXIT PDAFomi_omit_by_inno_l_cb #1'
#ifdef _OPENACC
    STOP "RSE: NOT PORTED CODE BLOCK #1"
#endif
  END IF

  CALL PDAF_timeit(12, 'old')


! *************************************************
! ***   Compute analyzed matrix Ainv            ***
! ***                                           ***
! ***  -1                       T  -1           ***
! *** A  = forget*(N-1) I + (HL)  R  (HL)       ***
! ***                                           ***
! *************************************************

  CALL PDAF_timeit(10, 'new')

  haveobsA: IF (dim_obs_l > 0) THEN
     ! *** The contribution of observation matrix ist only ***
     ! *** computed for domains with observations          ***

     ! WRITE(0,*) '************** RSE: ENTER haveobsA #1'

     CALL PDAF_timeit(30, 'new')

     ! *** Compute HL = [Hx_1 ... Hx_N] Omega
     ALLOCATE(HL_l(dim_obs_l, dim_ens))
     !$ACC ENTER DATA CREATE(HL_l(:,:)) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l * dim_ens)

     CALL PDAF_timeit(46, 'new')

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lestkf_analysis -- call g2l_obs', dim_ens, 'times'

     ! WRITE(0,*) 'RSE: dim_ens = ', dim_ens
     ENS: DO member = 1, dim_ens ! dim_ens = 40
        ! Store member index
        obs_member = member

        ! [Hx_1 ... Hx_(r+1)] for local analysis domain
        CALL U_g2l_obs(domain_p, step, dim_obs_f, dim_obs_l, HX_f(:, member), &
             HL_l(:, member))
     END DO ENS

     CALL PDAF_timeit(46, 'old')

     ! *** Set the value of the forgetting factor  ***
     ! *** Inserted here, because HL_l is required ***
     CALL PDAF_timeit(51, 'new')
     IF (type_forget == 2) THEN
        WRITE(0,*) '**************** RSE: CALL PDAF_set_forget_local #1'
#ifdef _OPENACC
        STOP "RSE: NOT PORTED CODE BLOCK #2"    
#endif
        CALL PDAF_set_forget_local(domain_p, step, dim_obs_l, dim_ens, HL_l, &
             HXbar_l, resid_l, obs_l, U_init_n_domains_p, U_init_obsvar_l, &
             forget)
        WRITE(0,*) '**************** RSE: EXIT PDAF_set_forget_local #1'
     ENDIF
     
     !$ACC EXIT DATA DELETE(HXbar_l) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)
     DEALLOCATE(HXbar_l)

     ! Complete HL = [Hx_1 ... Hx_N] Omega
     ! WRITE(0,*) '**************** RSE: CALL PDAF_estkf_AOmega #1'
#ifndef _OPENACC
     CALL PDAF_estkf_AOmega(dim_obs_l, dim_ens, HL_l)
#else
    !$ACC UPDATE DEVICE(HL_l(:,:)) ASYNC(accStreamId)
    !$ACC WAIT(accStreamId)
    CALL PDAF_estkf_AOmega(dim_obs_l, dim_ens, HL_l, accStreamId)
    !$ACC UPDATE HOST(HL_l(:,:)) ASYNC(accStreamId)    
    !$ACC WAIT(accStreamId)
#endif
     ! WRITE(0,*) '**************** RSE: EXIT PDAF_estkf_AOmega #1'

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  HXT_l', HL_l(:, 1:dim_ens-1)

     CALL PDAF_timeit(51, 'old')
     CALL PDAF_timeit(30, 'old')
     CALL PDAF_timeit(31, 'new')


     ! ***                RiHL = Rinv HL                 ***
     ! *** this is implemented as a subroutine thus that ***
     ! *** Rinv does not need to be allocated explicitly ***
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lestkf_analysis -- call prodRinvA_l'

     ALLOCATE(RiHL_l(dim_obs_l, rank))
     !$ACC ENTER DATA CREATE(RiHL_l(:,:)) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_l * rank)

     CALL PDAF_timeit(48, 'new')
     CALL U_prodRinvA_l(domain_p, step, dim_obs_l, rank, obs_l, HL_l, RiHL_l)
     CALL PDAF_timeit(48, 'old')

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  R^-1(HXT_l)', RiHL_l

     !$ACC EXIT DATA DELETE(obs_l) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)
     DEALLOCATE(obs_l)
 
     CALL PDAF_timeit(51, 'new')

     ! *** Initialize Ainv = (N-1) I ***
#ifndef _OPENACC

     Ainv_l = 0.0
     DO i = 1, rank
        Ainv_l(i,i) = REAL(dim_ens - 1)
     END DO

#else

     !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(accStreamId)
     DO jCol = 1, rank
        DO iRow = 1, rank
            IF(jCol .EQ. iRow) THEN
                Ainv_l(iRow, jCol) = REAL(dim_ens - 1)
            ELSE
                Ainv_l(iRow, jCol) = 0.0
            END IF
        END DO
     END DO
     !$ACC END PARALLEL LOOP
     !$ACC WAIT(accStreamId)

#endif

     ! ***             T        ***
     ! ***  Compute  HL  RiHL   ***
     ALLOCATE(tmp_Ainv_l(rank, rank))
     !$ACC ENTER DATA CREATE(tmp_Ainv_l(:,:)) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank**2)
     
#ifndef _OPENACC

     CALL gemmTYPE('t', 'n', rank, rank, dim_obs_l, &
          1.0, HL_l, dim_obs_l, RiHL_l, dim_obs_l, &
          0.0, tmp_Ainv_l, rank)

#else

#ifdef _RSE_USE_CUBLAS_DGEMM
     !$ACC UPDATE DEVICE(HL_l(:,:), RiHL_l(:,:)) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)

     !$ACC HOST_DATA USE_DEVICE(HL_l, RiHL_l, tmp_Ainv_l)
     ierror = cublasDgemm(handleCuBlas, CUBLAS_OP_T, CUBLAS_OP_N, rank, rank, dim_obs_l, &
                          1.0d0, HL_l, dim_obs_l, &
                          RiHL_l, dim_obs_l, &
                          0.0d0, tmp_Ainv_l, rank)
     !$ACC END HOST_DATA     
     ierror = cudaStreamSynchronize( cudaStream )
     
#else
     !$ACC UPDATE DEVICE(HL_l(:,:), RiHL_l(:,:)) ASYNC(accStreamId)
    
     !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) PRIVATE(kCell) ASYNC(accStreamId)
     DO jCol = 1, rank
        DO iRow = 1, rank

            tmp_Ainv_l(iRow, jCol) = 0.0d0
            !$ACC LOOP SEQ
            DO kCell = 1, dim_obs_l
                tmp_Ainv_l(iRow, jCol) = tmp_Ainv_l(iRow, jCol) + HL_l(kCell, iRow) * RiHL_l(kCell, jCol)
            ENDDO
            !$ACC END LOOP

        ENDDO
     ENDDO
     !$ACC END PARALLEL LOOP     
     !$ACC WAIT(accStreamId)
#endif

#endif

     !$ACC EXIT DATA DELETE(HL_l) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)
     DEALLOCATE(HL_l)
     CALL PDAF_timeit(51, 'old')

     ! WRITE(0,*) '************** RSE: EXIT haveobsA #1'

  ELSE haveobsA
     ! *** For domains with dim_obs_l=0 there is no ***
     ! *** direct observation-contribution to Ainv  ***

     WRITE(0,*) '************** RSE: ENTER ELSE haveobsA #1'
#ifdef _OPENACC
     STOP "RSE: NOT PORTED CODE BLOCK #3"
#endif

     CALL PDAF_timeit(31, 'new')
     CALL PDAF_timeit(51, 'new')

     ! *** Initialize Ainv = (N-1) I ***
     Ainv_l = 0.0
     DO i = 1, rank
        Ainv_l(i,i) = REAL(dim_ens - 1)
     END DO

     ! No observation-contribution to Ainv from this domain
     ALLOCATE(tmp_Ainv_l(rank, rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank**2)

     tmp_Ainv_l = 0.0

     CALL PDAF_timeit(51, 'old')

     WRITE(0,*) '************** RSE: EXIT ELSE haveobsA #1'

  END IF haveobsA

  ! *** Complete computation of Ainv  ***
  ! ***   -1                T         ***
  ! ***  A  = forget I  + HL RiHL     ***
  CALL PDAF_timeit(51, 'new')
#ifndef _OPENACC
  Ainv_l = forget * Ainv_l + tmp_Ainv_l
#else  
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(accStreamId)
  DO jCol = 1, rank
    DO iRow = 1, rank
        Ainv_l(iRow, jCol) = forget * Ainv_l(iRow, jCol) + tmp_Ainv_l(iRow, jCol)
    END DO
  END DO
  !$ACC END PARALLEL LOOP  
  !$ACC WAIT(accStreamId)
#endif
  CALL PDAF_timeit(51, 'old')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  A^-1_l', Ainv_l

  CALL PDAF_timeit(31, 'old')
  CALL PDAF_timeit(10, 'old')


! ***********************************************
! *** Compute weight for model state update   ***
! ***                                         ***
! ***              T                    f     ***
! ***    w = A RiHL d  with d = (y - H x )    ***
! ***********************************************

  CALL PDAF_timeit(13, 'new')
  CALL PDAF_timeit(51, 'new')

  ! *** Compute RiHLd = RiHL^T d ***
  ALLOCATE(RiHLd_l(rank))
  !$ACC ENTER DATA CREATE(RiHLd_l(:)) ASYNC(accStreamId)     
  !$ACC WAIT(accStreamId)
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank)

  haveobsC: IF (dim_obs_l > 0) THEN
     ! *** RiHLd_p/=0 only with observations ***

#ifndef _OPENACC

     CALL gemvTYPE('t', dim_obs_l, rank, 1.0, RiHL_l, &
          dim_obs_l, resid_l, 1, 0.0, RiHLd_l, 1)

#else     

     !$ACC HOST_DATA USE_DEVICE(RiHL_l, resid_l, RiHLd_l)
     ierror = cublasDgemv(handleCuBlas, CUBLAS_OP_T, dim_obs_l, rank, 1.0, &
                          RiHL_l, dim_obs_l, resid_l, 1, 0.0d0, RiHLd_l, 1)
     !$ACC END HOST_DATA
     ierror = cudaStreamSynchronize( cudaStream )
     
#endif

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  (HXT_l R^-1)^T d_l', RiHLd_l

     !$ACC EXIT DATA DELETE(RiHL_l, resid_l) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)
     DEALLOCATE(RiHL_l, resid_l)

  ELSE haveobsC

#ifdef _OPENACC
    STOP "RSE: NOT PORTED CODE BLOCK #4"
#endif

     RiHLd_l = 0.0

  END IF haveobsC


  ! *** Compute weight vector for state analysis:     ***
  ! ***          w = A RiHLd                          ***
  ! *** For this, two variants are implemented:       ***
  ! *** 1. solve for w in:                            ***
  ! ***           -1                                  ***
  ! ***          A  w = RiHLd                         ***
  ! ***   We use the LAPACK routine gesvTYPE          ***
  ! *** 2. Compute singular value decomposition       ***
  ! ***   of Ainv: Ainv = USV^T                       ***
  ! ***   Then: A = U S^(-1) V^T                      ***
  ! ***   This is combined with a symmetric           ***
  ! ***   square-root for the ensemble transformation ***

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lestkf_analysis -- type_sqrt', type_sqrt
  typeainv1: IF (type_sqrt==1) THEN
     ! *** Variant 1: Solve Ainv w= RiHLd for w

#ifdef _OPENACC
    STOP "RSE: NOT PORTED CODE BLOCK #5"
#endif

     ALLOCATE(ipiv(rank))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'i', rank)

     ! save matrix Ainv
     tmp_Ainv_l = Ainv_l

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, &
          '  Invert A^-1_l using solver GESV'

     ! call solver (gesvTYPE - LU solver)
     WRITE(0,*) 'RSE: GESV #1'
     CALL gesvTYPE(rank, 1, tmp_Ainv_l, rank, ipiv, &
          RiHLd_l, rank, lib_info)
     DEALLOCATE(ipiv)

  ELSE typeainv1
     ! *** Variant 2: Invert Ainv using SVD

     ALLOCATE(svals(rank))
     !$ACC ENTER DATA CREATE(svals(:)) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)
     
#ifndef _OPENACC

     ALLOCATE(work(3 * rank))
     ldwork = 3 * rank
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 4 * rank)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, &
          '  Compute eigenvalue decomposition of A^-1_l'

     ! Compute SVD of Ainv     
     CALL syevTYPE('v', 'l', rank, Ainv_l, rank, svals, work, ldwork, lib_info)

#else

     ! !$ACC UPDATE DEVICE(svals(:)) ASYNC(accStreamId)
     ! !$ACC WAIT(accStreamId)

#ifdef _RSE_CHECK_DnDsyevd_BUFFER_SIZE
     !$ACC HOST_DATA USE_DEVICE(Ainv_l, svals)
     ! jobz = CUSOLVER_EIG_MODE_NOVECTOR || CUSOLVER_EIG_MODE_VECTOR
     ! uplo = CUBLAS_FILL_MODE_UPPER || CUBLAS_FILL_MODE_LOWER
     ierror = cusolverDnDsyevd_buffersize(handleCuSolverDn, &
                                          CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, &
                                          rank, Ainv_l, rank, svals, ldwork)
     !$ACC END HOST_DATA
     ierror = cudaStreamSynchronize( cudaStream )
     ! WRITE(0,*) 'NATSEM LDWORK: ', ldwork
#else
     ldwork = _RSE_DnDsyevd_BUFFER_SIZE
#endif
     
     ALLOCATE( work(ldwork) )
     !$ACC ENTER DATA CREATE( work(:) , lib_info) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 4 * rank)     

     ! Compute SVD of Ainv
     !$ACC HOST_DATA USE_DEVICE(Ainv_l, svals, work, lib_info)
     ! jobz = CUSOLVER_EIG_MODE_NOVECTOR || CUSOLVER_EIG_MODE_VECTOR
     ! uplo = CUBLAS_FILL_MODE_UPPER || CUBLAS_FILL_MODE_LOWER
     ierror = cusolverDnDsyevd(handleCuSolverDn, &
                               CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, &
                               rank, Ainv_l, rank, svals, work, ldwork, lib_info)
     !$ACC END HOST_DATA
     ierror = cudaStreamSynchronize( cudaStream )

     !$ACC UPDATE HOST(lib_info) ASYNC(accStreamId)
     !$ACC EXIT DATA DELETE(work(:), lib_info) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)

#endif

     DEALLOCATE(work)

     ! Compute product A RiHLd
     IF (lib_info==0) THEN
        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  eigenvalues', svals

        ALLOCATE(VRiHLd_l(rank))
        !$ACC ENTER DATA CREATE(VRiHLd_l(:)) ASYNC(accStreamId)     
        !$ACC WAIT(accStreamId)
        IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank)

#ifndef _OPENACC

        CALL gemvTYPE('t', rank, rank, 1.0, Ainv_l, &
             rank, RiHLd_l, 1, 0.0, VRiHLd_l, 1)
     
        DO row = 1,rank
           VRiHLd_l(row) = VRiHLd_l(row) / svals(row)
        END DO
  
        CALL gemvTYPE('n', rank, rank, 1.0, Ainv_l, &
             rank, VRiHLd_l, 1, 0.0, RiHLd_l, 1)

#else        

        !$ACC HOST_DATA USE_DEVICE(Ainv_l, RiHLd_l, VRiHLd_l)
        ierror = cublasDgemv(handleCuBlas, CUBLAS_OP_T, rank, rank, 1.0, &
                             Ainv_l, rank, RiHLd_l, 1, 0.0d0, VRiHLd_l, 1)
        !$ACC END HOST_DATA
        ierror = cudaStreamSynchronize( cudaStream )

        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(accStreamId)
        DO row = 1,rank
            VRiHLd_l(row) = VRiHLd_l(row) / svals(row)
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC WAIT(accStreamId)

        !$ACC HOST_DATA USE_DEVICE(Ainv_l, RiHLd_l, VRiHLd_l)
        ierror = cublasDgemv(handleCuBlas, CUBLAS_OP_N, rank, rank, 1.0, &
                             Ainv_l, rank, VRiHLd_l, 1, 0.0d0, RiHLd_l, 1)
        !$ACC END HOST_DATA
        ierror = cudaStreamSynchronize( cudaStream )        

#endif

        !$ACC EXIT DATA DELETE(VRiHLd_l(:)) ASYNC(accStreamId)     
        !$ACC WAIT(accStreamId)
        DEALLOCATE(VRiHLd_l)
     END IF
  END IF typeainv1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  wbar_l', RiHLd_l

  CALL PDAF_timeit(51, 'old')
  CALL PDAF_timeit(13, 'old')

  ! *** check if SVD was successful
  IF (lib_info == 0) THEN
     flag = 0
  ELSE
     WRITE (*, '(/5x, a, i10, a/)') &
          'PDAF-ERROR(1): Domain ', domain_p, ' Problem in computation of analysis weights !!!'
     flag = 1
  END IF


! ****************************************************************
! ***     Transform state ensemble                             ***
! ***              a   _f   f                                  ***
! ***             X  = X + X  W                                ***
! *** with                                   -T      T         ***
! ***          W = Omega (RiHLd + sqrt(N-1) C   Omega )        ***
! ****************************************************************

! *** Prepare weight matrix for ensemble transformation ***

  CALL PDAF_timeit(51, 'new')
  check1: IF (flag == 0) THEN

     ! allocate fields
     ALLOCATE(OmegaT(rank, dim_ens))
     !$ACC ENTER DATA CREATE(OmegaT(:,:)) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', rank * dim_ens)

     IF (mype == 0 .AND. screen > 0 .AND. screenout) THEN
        WRITE (*, '(a, 5x, a)') 'PDAF', 'Transform state ensemble'
        IF (type_sqrt == 1) THEN
           WRITE (*, '(a, 5x, a)') 'PDAF', '--- use Cholesky square-root of A'
        ELSE
           WRITE (*, '(a, 5x, a)') 'PDAF', '--- use symmetric square-root of A'
        END IF
     END IF

     CALL PDAF_timeit(20, 'new')
     CALL PDAF_timeit(32, 'new')

     ! Part 1: square-root of A
     typeainv2: IF (type_sqrt == 1) THEN
        ! Variant, if Ainv has been inverted above by solving

#ifdef _OPENACC
        STOP "RSE: NOT PORTED CODE BLOCK #6"
#endif

        ! Store Ainv for temporary use
        tmp_Ainv_l(:, :) = Ainv_l(:, :)

        ! Cholesky decomposition of tmp_Ainv_l = C C^T
        WRITE(0,*) 'RSE: DPOTRF #1'
        CALL potrfTYPE('l', rank, tmp_Ainv_l, rank, lib_info)

        ! check if Cholesky decomposition was successful
        CholeskyOK: IF (lib_info == 0) THEN
           ! Decomposition OK, continue
           flag = 0
        ELSE
           ! Decomposition failed
           WRITE (*, '(/5x, a, i8, a/)') &
                '!!! PDAF-ERROR(3): Problem with Cholesky decomposition of Ainv - domain ', &
                domain_p, ' !!!'
           flag = 3
        ENDIF CholeskyOK

     ELSE typeainv2
        ! Variant, if SVD inversion of Ainv has been performed

        ! Use OmegaT as temporary array        
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) COLLAPSE(2) ASYNC(accStreamId)
        DO col = 1, rank
           DO row = 1, rank
              OmegaT(row, col) = Ainv_l(row, col) / SQRT(svals(col))
           END DO
        END DO
        !$ACC END PARALLEL LOOP        
        !$ACC WAIT(accStreamId)

#ifndef _OPENACC
        CALL gemmTYPE('n', 't', rank, rank, rank, &
             1.0, OmegaT, rank, Ainv_l, rank, &
             0.0, tmp_Ainv_l, rank)
#else
        !$ACC HOST_DATA USE_DEVICE(OmegaT, Ainv_l, tmp_Ainv_l)
        ierror = cublasDgemm(handleCuBlas, CUBLAS_OP_N, CUBLAS_OP_T, rank, rank, rank, &
                             1.0d0, OmegaT, rank, &
                             Ainv_l, rank, &
                             0.0d0, tmp_Ainv_l, rank)
        !$ACC END HOST_DATA

        ierror = cudaStreamSynchronize( cudaStream )        
#endif        
        
        !$ACC EXIT DATA DELETE(svals(:)) ASYNC(accStreamId)
        !$ACC WAIT(accStreamId)
        DEALLOCATE(svals)

        ! Set flag
        flag = 0

     END IF typeainv2

     CALL PDAF_timeit(32, 'old')

  END IF check1

  check2: IF (flag == 0) THEN
    
     ! *** Part 2: Product A^(1/2) Omega ***
    
     CALL PDAF_timeit(34, 'new')
     IF (type_sqrt == 1) THEN

#ifdef _OPENACC
        STOP "RSE: NOT PORTED CODE BLOCK #7"
#endif
        
        ! Initialize the matrix Omega from argument OmegaT_in
        OmegaT = OmegaT_in

        ! A = (Omega C^(-1)) by solving Ct A = OmegaT for A
        WRITE(0,*) 'RSE: DTRTRS #1'
        CALL trtrsTYPE('L', 'T', 'N', rank, dim_ens, &
             tmp_Ainv_l, rank, OmegaT, rank, lib_info)
     ELSE
        ! TEMP_AINV already contains matrix C (no more inversion)

#ifndef _OPENACC
        CALL gemmTYPE('n', 'n', rank, dim_ens, rank, &
             1.0, tmp_Ainv_l, rank, OmegaT_in, rank, &
             0.0, OmegaT, rank)
#else        
        !$ACC HOST_DATA USE_DEVICE(tmp_Ainv_l, OmegaT_in, OmegaT)
        ierror = cublasDgemm(handleCuBlas, CUBLAS_OP_N, CUBLAS_OP_N, rank, dim_ens, rank, &
                             1.0d0, tmp_Ainv_l, rank, &
                             OmegaT_in, rank, &
                             0.0d0, OmegaT, rank)
        !$ACC END HOST_DATA

        ierror = cudaStreamSynchronize( cudaStream )        
#endif

        lib_info = 0

     END IF
     CALL PDAF_timeit(34, 'old')

     ! check if solve was successful
     solveOK: IF (lib_info == 0) THEN
        ! Solve for A OK, continue
        flag = 0
     ELSE
        ! Solve for A failed
        WRITE (*, '(/5x, a, i8, a/)') &
             '!!! PDAF-ERROR(2): Problem in computation of transformation matrix - domain ', &
             domain_p, ' !!!'
        flag = 2

        CALL PDAF_timeit(20, 'old')
     END IF solveOK

  END IF check2

  check3: IF (flag == 0) THEN

     CALL PDAF_timeit(35, 'new')

     ! *** Part 4: Add RiHLd and multiply by scaling factor

     ! *** Add RiHLd to A^T stored in OmegaT
     fac = SQRT(REAL(dim_ens - 1))

     IF (incremental /= 2) THEN
        
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(accStreamId)
        DO j = 1, dim_ens
           DO i = 1, rank
              OmegaT(i, j) = fac * OmegaT(i, j) + RiHLd_l(i)
           END DO
        END DO
        !$ACC END PARALLEL LOOP        
        !$ACC WAIT(accStreamId)

     ELSE

#ifdef _OPENACC
        STOP "RSE: NOT PORTED CODE BLOCK #8"
#endif

        ! For ensemble 3D-Var update only ensemble perturbations
        DO j = 1, dim_ens
           DO i = 1, rank
              OmegaT(i, j) = fac * OmegaT(i, j)
           END DO
        END DO

     END IF

     !$ACC EXIT DATA DELETE(RiHLd_l) ASYNC(accStreamId)     
     !$ACC WAIT(accStreamId)
     DEALLOCATE(RiHLd_l)
      
     ! *** Omega A^T (A^T stored in OmegaT_l) ***
     ! WRITE(0,*) '**************** RSE: CALL PDAF_estkf_OmegaA #1'
#ifndef _OPENACC
     CALL PDAF_estkf_OmegaA(rank, dim_ens, OmegaT, TA)
#else
     ! RUN ON GPU
     CALL PDAF_estkf_OmegaA(rank, dim_ens, OmegaT, TA, accStreamId)
     ! RUN ON GPU     
#endif
     ! WRITE(0,*) '**************** RSE: EXIT PDAF_estkf_OmegaA #1'

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_lestkf_analysis:', debug, '  transform', TA

     CALL PDAF_timeit(35, 'old')
     CALL PDAF_timeit(20, 'old')


! *** Perform ensemble transformation ***

     ! Use block formulation for transformation
     maxblksize = 200
     IF (mype == 0 .AND. screen > 0 .AND. screenout) &
          WRITE (*, '(a, 5x, a, i5)') &
          'PDAF', '--- use blocking with size ', maxblksize

     ALLOCATE(ens_blk(maxblksize, dim_ens))
     !$ACC ENTER DATA CREATE(ens_blk(:,:)) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', maxblksize * dim_ens)

     ! WRITE(0,*) "RSE: blocking: blklower = 1,", dim_l, ",", maxblksize
     ! blocking: DO blklower = 1, dim_l(2), maxblksize(200)
     !$ACC UPDATE DEVICE(ens_l(:,:), state_l(:)) ASYNC(accStreamId)
     blocking: DO blklower = 1, dim_l, maxblksize

        blkupper = MIN(blklower + maxblksize - 1, dim_l)

        ! Store old state ensemble
        CALL PDAF_timeit(21, 'new')        

        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(accStreamId)
        DO col = 1, dim_ens
            DO iBlk = 1, blkupper - blklower + 1
                ens_blk(iBlk, col) = ens_l(blklower + iBlk - 1, col)
            END DO
        END DO
        !$ACC END PARALLEL LOOP

        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(accStreamId)
        DO col = 1, dim_ens
            DO iBlk = blklower, blkupper
                ens_l(iBlk, col) = state_l(iBlk)
            END DO
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC WAIT(accStreamId)

        CALL PDAF_timeit(21, 'old')

        !                        a  _f   f    T
        ! Transform ensemble:   X = X + X  T(A )
        CALL PDAF_timeit(22, 'new')

#ifndef _OPENACC

        CALL gemmTYPE('n', 'n', blkupper - blklower + 1, dim_ens, dim_ens, &
             1.0, ens_blk, maxblksize, TA, dim_ens, &
             1.0, ens_l(blklower, 1), dim_l)

#else

        !$ACC HOST_DATA USE_DEVICE(ens_blk, TA, ens_l)
        ierror = cublasDgemm(handleCuBlas, CUBLAS_OP_N, CUBLAS_OP_N, blkupper - blklower + 1, dim_ens, dim_ens, &
                             1.0d0, ens_blk, maxblksize, &
                             TA, dim_ens, &
                             1.0d0, ens_l(blklower, 1), dim_l)
        !$ACC END HOST_DATA

        ierror = cudaStreamSynchronize( cudaStream )
#endif

        CALL PDAF_timeit(22, 'old')        

     END DO blocking

     !$ACC UPDATE HOST(ens_l(:,:)) ASYNC(accStreamId)        
     
     !$ACC EXIT DATA DELETE(ens_blk(:,:)) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)
     DEALLOCATE(ens_blk)

     !$ACC EXIT DATA DELETE(OmegaT(:,:)) ASYNC(accStreamId)
     !$ACC WAIT(accStreamId)
     DEALLOCATE(OmegaT)

  END IF check3

  CALL PDAF_timeit(51, 'old')

! ************************
! *** GPU Finishing up ***
! ************************
#ifdef _OPENACC
  
  IF( .FALSE. ) THEN
    ierror = cublasDestroy(handleCuBlas)
    IF(ierror .NE. 0) THEN
      WRITE(*,*) "cublasDestroy returned err code ", ierror
      STOP
    ENDIF
  ENDIF

#endif

! ********************
! *** Finishing up ***
! ********************

  !$ACC EXIT DATA DELETE(tmp_Ainv_l) ASYNC(accStreamId)
  !$ACC WAIT(accStreamId)
  DEALLOCATE(tmp_Ainv_l)

  IF (allocflag == 0) allocflag = 1

  ! Store domain index
  lastdomain = domain_p

  ! Deallocating GPU memory for arrays passed to the subroutine as parameters
#ifdef _OPENACC
  accMemRefCounter = 0

#ifdef _RSE_ALLOCATE_SUBROUTINE_PARAMETERS_
  !$ACC UPDATE HOST(Ainv_l(:,:), TA(:,:)) ASYNC(accStreamId)
  !$ACC EXIT DATA DELETE(state_l(:), Ainv_l(:,:), ens_l(:,:), TA(:,:)) &
  !$ACC           DELETE(OmegaT_in(:,:)) ASYNC(accStreamId)
  ! Currently not used on GPU :: HX_f(:,:) + HXbar_f(:) + state_inc_l(:) ==> IN
  !$ACC WAIT(accStreamId)
#endif
#endif
  ! Deallocating GPU memory for arrays passed to the subroutine as parameters

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_lestkf_analysis -- END'

END SUBROUTINE PDAF_lestkf_analysis
