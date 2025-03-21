MODULE linearAlgebraGpuTemplate

    USE CUDAFOR
    USE CUBLAS_V2

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: type_laTask
    PUBLIC :: accMatrixMatrixMultVectorLevel
    PUBLIC :: callCuBlasDgemmWrapper
    PUBLIC :: cudaProcessSingleTask
    PUBLIC :: cudaProcessTasksBatchedDP
    PUBLIC :: cudaProcessTasksBatched

    TYPE type_laTask

    INTEGER :: taskIndex ! Comment
    INTEGER :: taskSize ! Comment

    ! Про указатель вместо ALLOCATABLE
    DOUBLE PRECISION, POINTER :: matrixA(:,:)
    DOUBLE PRECISION, POINTER :: matrixB(:,:)
    DOUBLE PRECISION, POINTER :: matrixR(:,:)    
    DOUBLE PRECISION, POINTER :: matrixEX(:,:)

    END TYPE type_laTask

    CONTAINS

    SUBROUTINE accMatrixMatrixMultVectorLevel(matrixA, matrixB, matrixR, alpha, beta, mSize)
    !$ACC ROUTINE VECTOR

        DOUBLE PRECISION, INTENT(in   ) :: matrixA(:,:)
        DOUBLE PRECISION, INTENT(in   ) :: matrixB(:,:)
        DOUBLE PRECISION, INTENT(inout) :: matrixR(:,:)
        DOUBLE PRECISION, INTENT(in) :: alpha, beta
        INTEGER, INTENT(in) :: mSize        

        INTEGER :: im, jm, km 

        !$ACC LOOP VECTOR COLLAPSE(2)
        DO jm = 1, mSize
            DO im = 1, mSize

                matrixR(im, jm) = matrixR(im, jm) * beta
                !$ACC LOOP SEQ
                DO km = 1, mSize
                    matrixR(im, jm) = matrixR(im, jm) + alpha * matrixA(im, km) * matrixB(km, jm)
                END DO
                !$END ACC LOOP

            END DO
        END DO
        !$ACC END LOOP

    END SUBROUTINE accMatrixMatrixMultVectorLevel

    SUBROUTINE callCuBlasDgemmWrapper(handle, matrixDim, MA, MB, MR, alpha, beta)

        TYPE(cublasHandle), INTENT(in) :: handle
        INTEGER, INTENT(in) :: matrixDim
        DOUBLE PRECISION, INTENT(in) :: MA(:,:), MB(:,:)
        DOUBLE PRECISION, INTENT(inout) :: MR(:,:)
        DOUBLE PRECISION, INTENT(in) :: alpha, beta

        INTEGER :: ierror       

        !$ACC HOST_DATA USE_DEVICE(MA, MB, MR)
        ierror = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, &
                             matrixDim, matrixDim, matrixDim, &
                             alpha, MA, matrixDim, &
                             MR, matrixDim, &
                             beta, MR, matrixDim)
        !$ACC END HOST_DATA        

    END SUBROUTINE callCuBlasDgemmWrapper

    ATTRIBUTES(global) SUBROUTINE cudaResetMatrixR(tasks, N)
        TYPE(type_laTask), DEVICE, INTENT(inout) :: tasks(1:N)
        INTEGER, VALUE, INTENT(in)       :: N

        INTEGER :: idTask

        idTask = threadidx%x

        DO WHILE(idTask .LE. N)
            
            tasks(idTask)%MatrixR = 7.0d0
            idTask = idTask + blockdim%x

        END DO
        
    END SUBROUTINE cudaResetMatrixR

    ATTRIBUTES(global) SUBROUTINE cudaProcessSingleTask(task, alpha, beta)
        TYPE(type_laTask), DEVICE, INTENT(inout) :: task
        DOUBLE PRECISION, VALUE, INTENT(in)      :: alpha, beta
        
        INTEGER :: numCells, idCell, idCol, idRow
        INTEGER :: km

        numCells = task%taskSize * task%taskSize
        idCell = threadidx%x + (blockidx%x - 1) * blockdim%x - 1

        DO WHILE(idCell .LT. numCells)
            
            idCol =     idCell / task%taskSize  + 1
            idRow = MOD(idCell,  task%taskSize) + 1

            task%matrixR(idRow, idCol) = task%matrixR(idRow, idCol) * beta
            DO km = 1, task%taskSize
                task%matrixR(idRow, idCol) = task%matrixR(idRow, idCol) &
                + alpha * task%matrixA(idRow, km) * task%matrixB(km, idCol)
            END DO

            idCell = idCell + griddim%x * blockdim%x

        END DO
        
    END SUBROUTINE cudaProcessSingleTask

    ATTRIBUTES(global) SUBROUTINE cudaProcessTasksBatchedDP(tasks, N, alpha, beta)
        TYPE(type_laTask), DEVICE, INTENT(inout) :: tasks(1:N)
        INTEGER, VALUE, INTENT(IN)               :: N
        DOUBLE PRECISION, VALUE, INTENT(in)      :: alpha, beta

        INTEGER :: idTask
        INTEGER :: gridSize, blockSize

        idTask = threadidx%x + (blockidx%x - 1) * blockdim%x ! - 1

        DO WHILE(idTask .LE. N)

            blockSize = 128
#if 1
            gridSize = tasks(idTask)%taskSize * tasks(idTask)%taskSize / blockSize
            IF(MOD(tasks(idTask)%taskSize * tasks(idTask)%taskSize, blockSize) /= 0) THEN
                gridSize = gridSize + 1
            END IF
#else
            gridSize = 1
#endif

            CALL cudaProcessSingleTask<<<gridSize, blockSize>>>(tasks(idTask), alpha, beta)

            idTask = idTask + griddim%x * blockdim%x

        END DO        
        
    END SUBROUTINE cudaProcessTasksBatchedDP

    ATTRIBUTES(global) SUBROUTINE cudaProcessTasksBatched(tasks, N, alpha, beta)
        TYPE(type_laTask), DEVICE, INTENT(inout) :: tasks(1:N)
        INTEGER, VALUE, INTENT(IN)               :: N
        DOUBLE PRECISION, VALUE, INTENT(in)      :: alpha, beta

        INTEGER :: idTask, numCells, idCell, idCol, idRow
        INTEGER :: km

        idTask = blockidx%x

        DO WHILE(idTask .LE. N)

            numCells = tasks(idTask)%taskSize * tasks(idTask)%taskSize
            idCell = threadidx%x - 1

            DO WHILE(idCell .LT. numCells)
            
                idCol =     idCell / tasks(idTask)%taskSize  + 1
                idRow = MOD(idCell,  tasks(idTask)%taskSize) + 1
    
                tasks(idTask)%matrixR(idRow, idCol) = tasks(idTask)%matrixR(idRow, idCol) * beta
                DO km = 1, tasks(idTask)%taskSize
                    tasks(idTask)%matrixR(idRow, idCol) = tasks(idTask)%matrixR(idRow, idCol) &
                    + alpha * tasks(idTask)%matrixA(idRow, km) * tasks(idTask)%matrixB(km, idCol)
                END DO
    
                idCell = idCell + blockdim%x
    
            END DO

            idTask = idTask + griddim%x

        END DO
        
    END SUBROUTINE cudaProcessTasksBatched

END MODULE linearAlgebraGpuTemplate
