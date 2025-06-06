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
!$Id$
!BOP
!
! !ROUTINE: PDAF_estkf_AOmega --- Operate matrix Omega on A as AOmega
!
! !INTERFACE:
SUBROUTINE PDAF_estkf_AOmega(dim, dim_ens, A, accStreamIndex)

! !DESCRIPTION:
! Operate matrix Omega on another matrix as
!         $A = A Omega$.
!
! Omega is a dim_ens x (dim_ens-1) matrix with matximum
! rank and zero column sums. It is computed by the 
! Householder reflection associate with the vector
! (N-1)^(-1) (1,...,1)^T.
!
! The values of Omega are
!    1 - 1 / (N (1/sqrt(N) + 1)) for i=j
!    - 1 / (N (1/sqrt(N) + 1)) for i/=j, i<N
!    - 1 / sqrt(N) for i=N
!
! In this routine the product A Omega is implemented as
! operations:
! 1. Compute the row sums of A
! 2. Normalize row sums by 1/(sqrt(N) + N)
! 3. Subtract value of last column multiplied by 1/(1+sqrt(N))
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2011-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim               ! dimension of states
  INTEGER, INTENT(in) :: dim_ens           ! Size of ensemble
  REAL, INTENT(inout) :: A(dim, dim_ens)   ! Input/output matrix
  INTEGER, INTENT(in), OPTIONAL :: accStreamIndex

! !CALLING SEQUENCE:
! Called by: PDAF_estkf_analysis
! Calls PDAF_memcount
!EOP
  
! *** local variables ***
  INTEGER :: row, col  ! Counters
  REAL :: normsum      ! Normalization for row sum
  REAL :: normlast     ! Normalization for last column
  REAL :: val          ! Temporary variable
  INTEGER, SAVE :: allocflag = 0  ! Flag for dynamic allocation
  REAL, ALLOCATABLE :: rowsums(:) ! Row sums of A

!$OMP threadprivate(allocflag)


! **********************
! *** INITIALIZATION ***
! **********************

  ALLOCATE(rowsums(dim))
  !$ACC ENTER DATA CREATE(rowsums(:)) ASYNC(accStreamIndex)

  ! !$ACC ENTER DATA CREATE(A(:,:))
  ! !$ACC UPDATE DEVICE(A(:,:))

  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL PDAF_memcount(3, 'r', dim)
     allocflag = 1
  END IF

  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(accStreamIndex)
  DO row = 1, dim
    rowsums(row)   = 0.0
  END DO
  !$ACC END PARALLEL LOOP

  !$ACC WAIT(accStreamIndex)
  ! Initialize normalization values
  normsum = 1.0 / REAL(dim_ens) / (1.0/SQRT(REAL(dim_ens))+1.0)
  normlast = 1.0 / (1.0 + SQRT(REAL(dim_ens)))


  ! *** Compute row sums of A ***
#ifndef _OPENACC
  DO col = 1, dim_ens     
     DO row = 1, dim
        rowsums(row) = rowsums(row) + A(row, col)
     END DO
  END DO
#else
  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(col) ASYNC(accStreamIndex)
  DO row = 1, dim
     !$ACC LOOP SEQ
     DO col = 1, dim_ens     
        rowsums(row) = rowsums(row) + A(row, col)
     END DO
  END DO
  !$ACC END PARALLEL LOOP
#endif

  ! Scale by NORMSUM
  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(accStreamIndex)
  DO row = 1, dim
    rowsums(row) = normsum * rowsums(row)
  END DO
  !$ACC END PARALLEL LOOP

  ! Substract scale value for last column
  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(val) ASYNC(accStreamIndex)
  DO row = 1, dim
     val = A(row, dim_ens) * normlast
     rowsums(row) = rowsums(row) + val
  END DO
  !$ACC END PARALLEL LOOP


! **********************************************
! ***  Operate Omega on A                    ***
! **********************************************

  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(accStreamIndex)
  DO col = 1, dim_ens - 1
     DO row = 1, dim
        A(row, col) = A(row, col) - rowsums(row)
     END DO
  END DO
  !$ACC END PARALLEL LOOP
  
  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(accStreamIndex)
  DO row = 1, dim
     A(row, dim_ens) = 0.0
  END DO
  !$ACC END PARALLEL LOOP


! ********************
! *** FINISHING UP ***
! ********************

  ! !$ACC UPDATE HOST(A(:,:))
  ! !$ACC EXIT DATA DELETE(A(:,:))  

  !$ACC EXIT DATA DELETE(rowsums(:)) ASYNC(accStreamIndex)
  !$ACC WAIT(accStreamIndex)
  DEALLOCATE(rowsums)

END SUBROUTINE PDAF_estkf_AOmega
