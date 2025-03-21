MODULE linearAlgebraGpuTemplate

    USE CUDAFOR

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: type_laTask
    PUBLIC :: accMatrixMatrixMultVectorLevel
    PUBLIC :: cudaProcessSingleTask

    PUBLIC :: cudaResetMatrixR

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

#if 0
            task%matrixR(idRow, idCol) = REAL(idCell) + 0.1d0 * REAL(idRow) + 0.01d0 * REAL(idCol)
#else
            task%matrixR(idRow, idCol) = task%matrixR(idRow, idCol) * beta
            DO km = 1, task%taskSize
                task%matrixR(idRow, idCol) = task%matrixR(idRow, idCol) &
                + alpha * task%matrixA(idRow, km) * task%matrixB(km, idCol)
            END DO
#endif

            idCell = idCell + griddim%x * blockdim%x

        END DO
        
    END SUBROUTINE cudaProcessSingleTask

END MODULE linearAlgebraGpuTemplate
