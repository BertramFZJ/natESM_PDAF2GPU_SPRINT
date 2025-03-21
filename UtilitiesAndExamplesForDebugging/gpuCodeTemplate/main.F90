PROGRAM MAIN

    USE OMP_LIB
    USE OPENACC
    USE CUDAFOR
    USE CUBLAS_V2    
    
    USE linearAlgebraGpuTemplate, ONLY: type_laTask, accMatrixMatrixMultVectorLevel, &
                                        cudaProcessSingleTask, cudaResetMatrixR

    IMPLICIT NONE

    INTEGER, PARAMETER :: numberOfOpenmpThreadsCPU = 4
        
    INTEGER, PARAMETER             :: numberOfCudaStreams = 8
    INTEGER(kind=cuda_stream_kind) :: cudaStream(1:numberOfCudaStreams)
    TYPE(cublasHandle)             :: handleCuBlas(1:numberOfCudaStreams)
    INTEGER                        :: ierror
    INTEGER                        :: idCudaStream

    INTEGER :: numberOfTasks ! Total number of individual tasks (virtual grid points)
    TYPE(type_laTask), ALLOCATABLE :: tasks(:) ! An array of structures with data for processing individual tasks

    INTEGER :: i, im, jm, km

    DOUBLE PRECISION :: alpha = 1.0d0
    DOUBLE PRECISION :: beta  = 0.0d0

    DOUBLE PRECISION :: deviation, ldeviation
    DOUBLE PRECISION :: execTime    
    
    DOUBLE PRECISION, POINTER :: ptrA(:,:), ptrB(:,:), ptrR(:,:)    
    
    WRITE(6,'(1x, A50)') "START  PROGRAM"
    FLUSH(6)

! #########################################################################################################################################################
! #########################################################################################################################################################
! #########################################################################################################################################################

    ! Memory allocation and initialization of an array of structures with parameters of individual tasks (grid points) on the CPU and GPU.
    ! Individual task - multiplication of square matrices matrixR := alpha * matrixA * matrixB + beta * matrixR
    ! The matrix MatrixEX stores the exact solution - the result of calculations on the CPU. It is used to check the correctness
    ! of the result of calculations on the GPU.
    ! The size of the matrices Matrix<A/B/C>(1:taskSize,1:taskSize) in different tasks may differ.

    numberOfTasks = 2048 * 16    
    ALLOCATE(tasks(1:numberOfTasks))

    DO i = 1, numberOfTasks
        
        tasks(i)%taskIndex = i
        tasks(i)%taskSize  = 64

        ALLOCATE( tasks(i)%matrixA (1:tasks(i)%taskSize, 1:tasks(i)%taskSize) )
        ALLOCATE( tasks(i)%matrixB (1:tasks(i)%taskSize, 1:tasks(i)%taskSize) )
        ALLOCATE( tasks(i)%matrixR (1:tasks(i)%taskSize, 1:tasks(i)%taskSize) )
        ALLOCATE( tasks(i)%matrixEX(1:tasks(i)%taskSize, 1:tasks(i)%taskSize) )

        tasks(i)%matrixA  = 1.0 + 0.001 * REAL(i, kind = 8)
        tasks(i)%matrixB  = 1.0 - 0.001 * REAL(i, kind = 8)
        tasks(i)%matrixR  = 0.0
        tasks(i)%matrixEX = 0.0

    END DO

    !$ACC ENTER DATA CREATE( tasks(1:numberOfTasks) ) ASYNC(1)
    !$ACC UPDATE DEVICE( tasks(1:numberOfTasks) ) ASYNC(1)
    DO i = 1, numberOfTasks
        !$ACC ENTER DATA CREATE( tasks(i)%matrixA(:,:), tasks(i)%matrixB(:,:) ) &
        !$ACC            CREATE( tasks(i)%matrixR(:,:) ) ASYNC(1)
        !$ACC UPDATE DEVICE( tasks(i)%matrixA(:,:), tasks(i)%matrixB(:,:) ) &
        !$ACC        DEVICE( tasks(i)%matrixR(:,:) ) ASYNC(1)
    END DO
    !$ACC WAIT(1)
    WRITE(6,'(1x, A50)') "CRATE DATA ON CPU & GPU"
    FLUSH(6)

! #########################################################################################################################################################
! #########################################################################################################################################################
! #########################################################################################################################################################

    ! Matrix multiplication on CPU using dgemm subroutine from BLAS library.
    
    execTime = omp_get_wtime()
    !$OMP PARALLEL DO NUM_THREADS(numberOfOpenmpThreadsCPU)
    DO i = 1, numberOfTasks
        CALL dgemm('n', 'n', tasks(i)%taskSize, tasks(i)%taskSize, tasks(i)%taskSize, &
                   alpha, tasks(i)%matrixA , tasks(i)%taskSize, &
                          tasks(i)%matrixB , tasks(i)%taskSize, &
                   beta,  tasks(i)%matrixEX, tasks(i)%taskSize)
    END DO
    !$OMP END PARALLEL DO
    execTime = omp_get_wtime() - execTime
    WRITE(6,'(1x, A50, 1x, F20.4, 1x, A, 1x, I3, 1x, A)') "[CPU] BLAS dgemm EXECUTION TIME:", execTime, "SEC FOR ", numberOfOpenmpThreadsCPU, "OpenMP THREADS"
    FLUSH(6)

! #########################################################################################################################################################
! #########################################################################################################################################################
! #########################################################################################################################################################

    ! Collective matrix multiplication on GPU using GANG VECTOR parallelization scheme

    execTime = omp_get_wtime()
    !$ACC PARALLEL DEFAULT(PRESENT) FIRSTPRIVATE(alpha, beta) VECTOR_LENGTH(256) ASYNC(1)
    !$ACC LOOP GANG(STATIC:1) PRIVATE(im, jm, km)
    DO i = 1, numberOfTasks

        !$ACC LOOP VECTOR COLLAPSE(2)
        DO jm = 1, tasks(i)%taskSize
            DO im = 1, tasks(i)%taskSize

                tasks(i)%matrixR(im, jm) = tasks(i)%matrixR(im, jm) * beta
                !$ACC LOOP SEQ
                DO km = 1, tasks(i)%taskSize
                    tasks(i)%matrixR(im, jm) = tasks(i)%matrixR(im, jm) &
                    + alpha * tasks(i)%matrixA(im, km) * tasks(i)%matrixB(km, jm)
                END DO
                !$END ACC LOOP

            END DO
        END DO
        !$ACC END LOOP

    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL
    !$ACC WAIT(1)
    execTime = omp_get_wtime() - execTime
    WRITE(6,'(1x, A50, 1x, F20.4, 1x, A)') "[GPU] GANG VECTOR LOOP EXECUTION TIME:", execTime, "SEC"
    FLUSH(6)

    CALL calculateDeviatoions(tasks, numberOfTasks, .TRUE.)

! #########################################################################################################################################################
! #########################################################################################################################################################
! #########################################################################################################################################################

    ! Call "ACC ROUTINE VECTOR" subroutine for individual tasks within a parallel GANG loop over grid nodes

    CALL resetMatrixR(tasks, numberOfTasks)

    execTime = omp_get_wtime()
    ! In this particular case, the compiler ignores the vector size settings (vector_length clause) and sets the vector length to 32.
    ! The number of gangs is equal to the number of grid nodes (individual tasks).
    !$ACC PARALLEL LOOP GANG(STATIC:1) DEFAULT(PRESENT) FIRSTPRIVATE(alpha, beta) ASYNC(1)
    DO i = 1, numberOfTasks

        ! If the subroutine is declared in an external module/file, compilation must be performed with the "-gpu=lto" flag        
        CALL accMatrixMatrixMultVectorLevel(tasks(i)%matrixA(:,:), tasks(i)%matrixB, tasks(i)%matrixR, &
                                            alpha, beta, tasks(i)%taskSize)        

    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)    
    execTime = omp_get_wtime() - execTime
    WRITE(6,'(1x, A50, 1x, F20.4, 1x, A)') "[GPU] GANG LOOP + VECTOR SUBROUTINE ET:", execTime, "SEC"
    FLUSH(6)

    CALL calculateDeviatoions(tasks, numberOfTasks, .TRUE.)

! #########################################################################################################################################################
! #########################################################################################################################################################
! #########################################################################################################################################################

    ! Concurrent execution of the cublasDgemm subroutine for individual tasks in different CUDA streams

    CALL resetMatrixR(tasks, numberOfTasks)

    ierror = cudaSetDevice( 0 )
    IF(ierror .NE. 0) THEN
      WRITE(0,*) "cudaSetDevice returned err code ", ierror; STOP
    ENDIF

    DO i = 1, numberOfCudaStreams
        ierror = cublasCreate( handleCuBlas(i) )
        IF(ierror .NE. 0) THEN
            WRITE(0,*) "cublasCreate returned err code ", ierror; STOP
        ENDIF
    END DO

    DO i = 1, numberOfCudaStreams
        ierror = cudaStreamCreate( cudaStream(i) )
        IF(ierror .NE. 0) THEN
            WRITE(0,*) "cudaStreamCreate returned err code ", ierror; STOP
        ENDIF
    END DO

    DO i = 1, numberOfCudaStreams
        ierror = cublasSetStream( handleCuBlas(i), cudaStream(i) )
        IF(ierror .NE. 0) THEN
            WRITE(0,*) "cublasSetStream returned err code ", ierror; STOP
        ENDIF
    END DO    

    execTime = omp_get_wtime()    
    !$OMP PARALLEL NUM_THREADS(numberOfCudaStreams) PRIVATE(idCudaStream, ierror, ptrA, ptrB, ptrR)

    idCudaStream = omp_get_thread_num() + 1

    !$OMP DO SCHEDULE(STATIC)
    DO i = 1, numberOfTasks

        ! Passing direct pointers to structure members (tasks(i)%matrixA) as cublasDgemm parameters results in an error:
        ! NVFORTRAN-S-0155-Could not resolve generic procedure cublasdgemm (main.F90: XXX)
        ptrA => tasks(i)%matrixA
        ptrB => tasks(i)%matrixB
        ptrR => tasks(i)%matrixR

        !$ACC HOST_DATA USE_DEVICE(ptrA, ptrB, ptrR)
        ierror = cublasDgemm(handleCuBlas(idCudaStream), CUBLAS_OP_N, CUBLAS_OP_N, &
                             tasks(i)%taskSize, tasks(i)%taskSize, tasks(i)%taskSize, &
                             alpha, ptrA, tasks(i)%taskSize, &
                             ptrB, tasks(i)%taskSize, &
                             beta, ptrR, tasks(i)%taskSize)
        !$ACC END HOST_DATA

        ! cublasDgemm is an asynchronous subroutine. Therefore, theoretically, checking for completion of tasks
        ! in the stream can be moved outside the loop
        ierror = cudaStreamSynchronize( cudaStream(idCudaStream) )
    END DO
    !$OMP END DO
    ! Host/device synchronization after all kernels are queued in the CUDA stream
    ! ierror = cudaStreamSynchronize( cudaStream(idCudaStream) )

    !$OMP END PARALLEL
    execTime = omp_get_wtime() - execTime
    WRITE(6,'(1x, A50, 1x, F20.4, 1x, A, 1x, I3, 1x, A)') "[CPU] CuBLAS cublasDgemm EXECUTION TIME:", execTime, "SEC FOR ", &
          numberOfCudaStreams, "CUDA STREAMS / OpenMP THREADS"
    FLUSH(6)

    CALL calculateDeviatoions(tasks, numberOfTasks, .TRUE.)

! #########################################################################################################################################################
! #########################################################################################################################################################
! #########################################################################################################################################################

    CALL resetMatrixR(tasks, numberOfTasks)

    execTime = omp_get_wtime()
#if 0
    DO i = 1, numberOfTasks
        !$ACC HOST_DATA USE_DEVICE(tasks(i))
        CALL cudaProcessSingleTask<<<4,128>>>(tasks(i), alpha, beta)
        !$ACC END HOST_DATA
    END DO
    ierror = cudaDeviceSynchronize()
#else
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) DEVICEPTR(tasks)
    DO i = 1, numberOfTasks
        ! !$ACC HOST_DATA USE_DEVICE(tasks(i))
        CALL cudaProcessSingleTask<<<4,128>>>(tasks(i), alpha, beta)
        ! !$ACC END HOST_DATA
    END DO
    !$ACC END PARALLEL LOOP
#endif
    execTime = omp_get_wtime() - execTime
    WRITE(6,'(1x, A50, 1x, F20.4, 1x, A)') "[GPU] cudaProcessSingleTask LOOP EXECUTION TIME:", execTime, "SEC"
    

    CALL calculateDeviatoions(tasks, numberOfTasks, .TRUE.)
    ! WRITE(6,*) tasks(1)%matrixR(:,:)

    
    
    WRITE(6,*) "FINISH PROGRAM"

    CONTAINS

#if 0
    subroutine callCuBlasDgemm(handle, matrixDim, MA, MB, MR, alpha, beta)

        TYPE(cublasHandle), INTENT(in) :: handle
        INTEGER, INTENT(in) :: matrixDim
        DOUBLE PRECISION, INTENT(in) :: MA(:,:), MB(:,:)
        DOUBLE PRECISION, INTENT(inout) :: MR(:,:)
        DOUBLE PRECISION, INTENT(in) :: alpha, beta

        INTEGER :: ierror

    #if 1
        DOUBLE PRECISION, ALLOCATABLE :: CMA(:,:), CMB(:,:), CMR(:,:)
        ALLOCATE(CMA(1:matrixDim,1:matrixDim), CMB(1:matrixDim,1:matrixDim), CMR(1:matrixDim,1:matrixDim))

        CMA(1:matrixDim,1:matrixDim) = MA(1:matrixDim,1:matrixDim)
        CMB(1:matrixDim,1:matrixDim) = MB(1:matrixDim,1:matrixDim)
        CMR(1:matrixDim,1:matrixDim) = MR(1:matrixDim,1:matrixDim)

        !$ACC DATA COPYIN(CMA(:,:), CMB(:,:)) COPYOUT(CMR(:,:))

        !$ACC HOST_DATA USE_DEVICE(CMA, CMB, CMR)
        ierror = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, &
                             matrixDim, matrixDim, matrixDim, &
                             alpha, CMA, matrixDim, &
                             CMB, matrixDim, &
                             beta, CMR, matrixDim)
        !$ACC END HOST_DATA
        ierror = cudaDeviceSynchronize()        
        
        !$ACC END DATA

        MR(1:matrixDim,1:matrixDim) = CMR(1:matrixDim,1:matrixDim)

        DEALLOCATE(CMA, CMB, CMR)
    #else
        !$ACC DATA COPYIN(MA(:,:), MB(:,:)) COPYOUT(MR(:,:))

        !$ACC HOST_DATA USE_DEVICE(MA, MB, MR)
        ierror = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, &
                             matrixDim, matrixDim, matrixDim, &
                             alpha, MA, matrixDim, &
                             MR, matrixDim, &
                             beta, MR, matrixDim)
        !$ACC END HOST_DATA
        ierror = cudaDeviceSynchronize()        
        
        !$ACC END DATA        
    #endif

    end subroutine callCuBlasDgemm
#endif

    subroutine compareMatrixes(matrixA, matrixB, nRow, nCol, deviation)

        double precision, intent(in) :: matrixA(nRow, nCol)
        double precision, intent(in) :: matrixB(nRow, nCol)
        integer, intent(in) :: nRow, nCol
        double precision, optional, intent(out) :: deviation

        integer :: iterRow, iterCol
        double precision :: subDeviation, itDeviation

        subDeviation = 0.0
        DO iterCol = 1, nCol
            DO iterRow = 1, nRow
                itDeviation = ABS(matrixB(iterRow, iterCol) - matrixA(iterRow, iterCol))
                IF(itDeviation .GT. subDeviation) subDeviation = itDeviation
            ENDDO
        ENDDO

        IF( PRESENT(deviation) ) THEN
            deviation = subDeviation
        ELSE
            WRITE(*,*) 'DEVIATION: ', subDeviation
        ENDIF        

    end subroutine compareMatrixes

    subroutine calculateDeviatoions(tasks, N, updateHost)

        type(type_laTask), intent(in) :: tasks(1:N)
        integer, intent(in)           :: N
        logical, intent(in)           :: updateHost

        integer                       :: i
        double precision              :: subDeviation, taskDeviation

        IF(updateHost) THEN            
            DO i = 1, N
                !$ACC UPDATE HOST( tasks(i)%matrixR(:,:) )
            END DO        
        END IF

        subDeviation = 0.0d0
        DO i = 1, N            
            CALL compareMatrixes(tasks(i)%matrixEX, tasks(i)%matrixR, &
                                 tasks(i)%taskSize, tasks(i)%taskSize, taskDeviation)
            IF (subDeviation < taskDeviation) subDeviation = taskDeviation
        END DO
        
        WRITE(6,'(1x, A50, 1x, E20.8)') "DEVIATION:", subDeviation
        FLUSH(6)

    end subroutine calculateDeviatoions

    subroutine resetMatrixR(tasks, N)

        type(type_laTask), intent(in) :: tasks(1:N)
        integer, intent(in)           :: N        

        integer                       :: i

        DO i = 1, N
            tasks(i)%matrixR  = 0.0
            !$ACC UPDATE DEVICE( tasks(i)%matrixR ) ASYNC(1)
        END DO
        !$ACC WAIT(1)

    end subroutine resetMatrixR


END PROGRAM MAIN
