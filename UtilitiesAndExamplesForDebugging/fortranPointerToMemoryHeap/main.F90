program main

    USE OMP_LIB

    USE mo_mem_worspaces, ONLY: memInitLibrary, memSetOmpNumThreads, memAllocateAnalysisHeap, &
                                memDeallocateAnalysisHeap, memGetPointerAnalysisHeap, &
                                memGetThreadPointerAnalysisHeap

    IMPLICIT NONE

    INTEGER, PARAMETER :: threadGroupSize = 8 ! Number of OpenMP threads in a group
    INTEGER, PARAMETER :: heapSize = 17 ! Memory block size for a single thread

    INTEGER, POINTER :: heap(:,:) ! A heap of memory

    ! Temporary arrays
    INTEGER, POINTER :: array2D(:,:) => null() ! array2D(1:nRow, 1:nCol)
    INTEGER, POINTER :: array1D(:) => null() ! array1D(1:nCell)
    INTEGER, POINTER :: threadHeap(:) => null()
    INTEGER :: nCol, nRow, nCell

    INTEGER :: col, row, cell, iThread, threadId, groupSize
    INTEGER :: heapOffset
    9997 FORMAT( 24(:,1X,I4) )

    ! Allocating memory for the heap
    ! The heap has the form of a two-dimensional matrix. Each OpenMP thread uses one of the matrix columns.
    CALL memInitLibrary()
    CALL memSetOmpNumThreads( threadGroupSize )
    CALL memAllocateAnalysisHeap( heapSize = heapSize, numaInitFlag = .TRUE., accEnterDataFlag = .TRUE.)
    CALL memGetPointerAnalysisHeap( heap )
    
    !$ACC KERNELS DEFAULT(PRESENT)
    heap = - 1
    !$ACC END KERNELS
    !$ACC UPDATE HOST(heap(:,:))

    !$OMP PARALLEL NUM_THREADS(threadGroupSize) PRIVATE(col, nCol, row, nRow, cell, nCell, iThread, threadId, groupSize) &
    !$OMP PRIVATE(array2D, array1D, threadHeap, heapOffset) SHARED(heap)

    threadId = OMP_GET_THREAD_NUM()
    groupSize = OMP_GET_NUM_THREADS()
    
    !$OMP SINGLE
    WRITE(*,*) "OMP THREAD NUMBER: ", groupSize
    WRITE(*,*)
    !$OMP END SINGLE
    !$OMP BARRIER

    ! Initializing the sizes of temporary arrays
    nCol = 3;
    nRow = 2 + MOD(threadId, 3)
    nCell = 3

    DO iThread = 0, groupSize
        IF(iThread .EQ. threadId) THEN
            WRITE(*,*) iThread, ":: MATRIX SIZE: ", nRow, nCol
        END IF
        !$OMP BARRIER
    END DO
    !$OMP SINGLE
    WRITE(*,*)
    !$OMP END SINGLE

    ! Allocating memory for temporary arrays
    ! The shift parameter stores the size of the memory block already occupied by temporary arrays.
#if 1
    CALL memGetThreadPointerAnalysisHeap(threadHeap, threadId)
    heapOffset = 0
    array2D(1:nRow, 1:nCol) => threadHeap(heapOffset + 1:heapOffset + nCol * nRow)
    heapOffset = heapOffset + nRow * nCol
    array1D(1:nCell) => threadHeap(heapOffset + 1:heapOffset + nCell)
    heapOffset = heapOffset + nCell
#else
    heapOffset = 0
    array2D(1:nRow, 1:nCol) => heap(heapOffset + 1:heapOffset + nCol * nRow, threadId + 1)
    heapOffset = heapOffset + nRow * nCol
    array1D(1:nCell) => heap(heapOffset + 1:heapOffset + nCell, threadId + 1)
    heapOffset = heapOffset + nCell
#endif

    !$OMP BARRIER
    DO iThread = 0, groupSize
        IF(iThread .EQ. threadId) THEN
            WRITE(*,*) iThread, ":: HEAP OFFSET: ", heapOffset
        END IF
        !$OMP BARRIER
    END DO
    !$OMP SINGLE
    WRITE(*,*)
    !$OMP END SINGLE

    DO col = 1, nCol
        DO row = 1, nRow
            array2D(row, col) = threadId*100 + col*10 + row
        END DO
    END DO

    DO cell = 1, nCell
        array1D(cell) = threadId*200 + cell
    END DO
    
    !$OMP BARRIER
    DO iThread = 0, groupSize - 1
        IF(iThread .EQ. threadId) THEN
            WRITE(*,*)
            WRITE(*,*) "********************", iThread, "********************"
            CALL PRINT_INTEGER_MATRIX_2D("matrix", nRow, nCell, array2D, nRow )
            WRITE(*,*) "vector"
            WRITE(*, 9997) (array1D(cell), cell=1,nCell)
        END IF
        !$OMP BARRIER
    END DO

    !$OMP BARRIER
    !$OMP SINGLE
    WRITE(*,*)
    CALL PRINT_INTEGER_MATRIX_2D("heap", heapSize, threadGroupSize, heap, heapSize )    
    !$OMP END SINGLE
    !$OMP BARRIER

#ifdef _OPENACC

    !$ACC UPDATE DEVICE(array2D(:,:), array1D(:))

    !$OMP BARRIER
    !$OMP SINGLE
    heap = -1
    !$ACC UPDATE HOST(heap(:,:))
    WRITE(*,*)
    CALL PRINT_INTEGER_MATRIX_2D("heap CPU ==> GPU", heapSize, threadGroupSize, heap, heapSize )
    
    heap = -1
    !$ACC UPDATE DEVICE(heap(:,:))
    !$OMP END SINGLE    
    !$OMP BARRIER

    !$ACC PARALLEL LOOP GANG VECTOR NUM_GANGS(1) VECTOR_LENGTH(64) DEFAULT(PRESENT)
    DO col = 1, nCol
        DO row = 1, nRow
            array2D(row, col) = threadId*100 + col*10 + row
        END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP GANG VECTOR NUM_GANGS(1) VECTOR_LENGTH(32) DEFAULT(PRESENT)
    DO cell = 1, nCell
        array1D(cell) = threadId*200 + cell
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC UPDATE HOST(array2D(:,:), array1D(:))
    !$OMP BARRIER

    !$OMP SINGLE
    WRITE(*,*)
    CALL PRINT_INTEGER_MATRIX_2D("heap GPU ==> CPU", heapSize, threadGroupSize, heap, heapSize )    
    !$OMP END SINGLE    
    !$OMP BARRIER

#endif
    
    !$OMP END PARALLEL

    CALL memDeallocateAnalysisHeap()

    CONTAINS

    SUBROUTINE PRINT_INTEGER_MATRIX_2D( DESC, M, N, A, LDA )
        CHARACTER(*)    DESC
        INTEGER          M, N, LDA
        INTEGER          A( LDA, * )

        INTEGER          I, J

        ! WRITE(*,*)
        WRITE(*,*) DESC
        DO I = 1, M
           WRITE(*,9998) ( A( I, J ), J = 1, N )
        END DO

    9998 FORMAT( 24(:,1X,I4) )
        RETURN
    END SUBROUTINE PRINT_INTEGER_MATRIX_2D

end program main
