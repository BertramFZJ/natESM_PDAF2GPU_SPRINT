MODULE linearAlgebraGpuTemplate

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: type_laTask
    PUBLIC :: accMatrixMatrixMultVectorLevel

    TYPE type_laTask

    INTEGER :: taskIndex
    INTEGER :: taskSize

    ! DOUBLE PRECISION, ALLOCATABLE :: matrixA(:,:)
    ! DOUBLE PRECISION, ALLOCATABLE :: matrixB(:,:)
    ! DOUBLE PRECISION, ALLOCATABLE :: matrixR(:,:)

    DOUBLE PRECISION, POINTER :: matrixA(:,:)
    DOUBLE PRECISION, POINTER :: matrixB(:,:)
    DOUBLE PRECISION, POINTER :: matrixR(:,:)
    
    DOUBLE PRECISION, ALLOCATABLE :: matrixEX(:,:)

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


END MODULE linearAlgebraGpuTemplate
