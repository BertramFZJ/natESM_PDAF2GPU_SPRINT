#define DATA_TYPE INTEGER
#define INIT_CONSTANT_VALUE -1

MODULE mo_mem_worspaces

    USE OMP_LIB

#ifdef _OPENACC
    USE OPENACC
#endif

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: memInitLibrary
    PUBLIC :: memSetOmpNumThreads
    PUBLIC :: memAllocateAnalysisHeap
    PUBLIC :: memGetPointerAnalysisHeap
    PUBLIC :: memDeallocateAnalysisHeap
    PUBLIC :: memGetThreadPointerAnalysisHeap    

    INTEGER :: numOmpThreads

    INTEGER :: heapAnalysisSize
    DATA_TYPE, ALLOCATABLE, TARGET :: heapAnalysis(:,:)

    CONTAINS

    SUBROUTINE memInitLibrary()

        numOmpThreads = -1

        heapAnalysisSize = 0

    END SUBROUTINE memInitLibrary

    SUBROUTINE memSetOmpNumThreads( OmpNumThreads )

        INTEGER, INTENT(IN) :: OmpNumThreads

        numOmpThreads = OmpNumThreads        

    END SUBROUTINE memSetOmpNumThreads

    SUBROUTINE memAllocateAnalysisHeap( heapSize, heapPointer, numaInitFlag, accEnterDataFlag )

        INTEGER, INTENT(IN) :: heapSize
        DATA_TYPE, OPTIONAL, POINTER, INTENT(OUT) :: heapPointer(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: numaInitFlag
        LOGICAL, OPTIONAL, INTENT(IN) :: accEnterDataFlag

        heapAnalysisSize = heapSize
        ALLOCATE(heapAnalysis(1:heapAnalysisSize, 1:numOmpThreads))

        IF( PRESENT(numaInitFlag) .AND. (numaInitFlag .EQV. .TRUE.) ) THEN

            !$OMP PARALLEL NUM_THREADS(numOmpThreads)
            heapAnalysis(:, OMP_GET_THREAD_NUM() + 1) = INIT_CONSTANT_VALUE
            !$OMP END PARALLEL

        END IF

#ifdef _OPENACC
        IF( PRESENT(accEnterDataFlag) .AND. (accEnterDataFlag .EQV. .TRUE.) ) THEN

            !$ACC ENTER DATA CREATE(heapAnalysis(:,:))
    
        END IF
#endif

        IF( PRESENT(heapPointer) ) THEN

            heapPointer => heapAnalysis
        
        END IF

    END SUBROUTINE memAllocateAnalysisHeap

    SUBROUTINE memGetPointerAnalysisHeap( heapPointer )

        DATA_TYPE, POINTER, INTENT(OUT) :: heapPointer(:,:)
        
        heapPointer => heapAnalysis        

    END SUBROUTINE memGetPointerAnalysisHeap

    SUBROUTINE memGetThreadPointerAnalysisHeap( heapThreadPointer, threadId )

        DATA_TYPE, POINTER, INTENT(OUT) :: heapThreadPointer(:)
        INTEGER, OPTIONAL, INTENT(IN) :: threadId

        INTEGER :: colHeap

        IF( PRESENT(threadId) ) THEN
            colHeap = threadId + 1
        ELSE
            colHeap = OMP_GET_THREAD_NUM() + 1
        END IF
        
        heapThreadPointer => heapAnalysis(1:heapAnalysisSize, colHeap)   

    END SUBROUTINE memGetThreadPointerAnalysisHeap

    SUBROUTINE memDeallocateAnalysisHeap()

#ifdef _OPENACC
        IF( acc_is_present( heapAnalysis(:,:) ) ) THEN
            !$ACC EXIT DATA DELETE(heapAnalysis(:,:))    
        END IF
#endif

        DEALLOCATE(heapAnalysis)

    END SUBROUTINE memDeallocateAnalysisHeap

END MODULE mo_mem_worspaces

#undef DATA_TYPE
#undef INIT_CONSTANT_VALUE
