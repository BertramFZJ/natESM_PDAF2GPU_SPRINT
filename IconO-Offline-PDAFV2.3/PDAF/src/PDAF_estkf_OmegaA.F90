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
! !ROUTINE: PDAF_estkf_OmegaA --- Operate matrix Omega on some matrix
!
! !INTERFACE:
SUBROUTINE PDAF_estkf_OmegaA(rank, dim_col, A, B, accStreamIndex)

! !DESCRIPTION:
! Operate matrix Omega on another matrix as
!                 B = Omega A\\
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
! 1. Compute the column sums of A
! 2. Normalize column sums by 1/(sqrt(N) + N)
! 3. Subtract value of last row multiplied by 1/(1+sqrt(N))
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
  INTEGER, INTENT(in) :: rank               ! Rank of initial covariance matrix
  INTEGER, INTENT(in) :: dim_col            ! Number of columns in A and B
  REAL, INTENT(in)    :: A(rank, dim_col)   ! Input matrix
  REAL, INTENT(out)   :: B(rank+1, dim_col) ! Output matrix (TA)
  INTEGER, INTENT(in), OPTIONAL :: accStreamIndex

! !CALLING SEQUENCE:
! Called by: PDAF_estkf_analysis
! Calls: PDAF_memcount
!EOP
  
! *** local variables ***
  INTEGER :: row, col  ! counters
  REAL :: normsum      ! Normalization for row sum
  REAL :: normlast     ! Normalization for last column
  INTEGER, SAVE :: allocflag = 0  ! Flag for dynamic allocation
  REAL, ALLOCATABLE :: colsums(:) ! Mean values of columns of A

!$OMP threadprivate(allocflag)

! WRITE(0,*) '**************** RSE: CALL PDAF_estkf_OmegaA #1 +++'

! **********************
! *** INITIALIZATION ***
! **********************

  ALLOCATE(colsums(dim_col))
  IF (allocflag == 0) THEN
     ! count allocated memory
     CALL PDAF_memcount(3, 'r', dim_col)
     allocflag = 1
  END IF
  
  ! Initialize normalization values
  normsum = 1.0 / REAL(rank+1) / (1.0/SQRT(REAL(rank+1))+1.0)
  normlast = - 1.0 / SQRT(REAL(rank+1))

  !$ACC ENTER DATA CREATE(colsums(:)) ASYNC(accStreamIndex)
  ! !$ACC ENTER DATA CREATE(A(:,:), B(:,:)) ASYNC(accStreamIndex)
  ! !$ACC UPDATE DEVICE(A(:,:)) ASYNC(accStreamIndex)

  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(accStreamIndex)
  DO col = 1, dim_col
    colsums(col)   = 0.0
  END DO
  !$ACC END PARALLEL LOOP  

! *** Compute column sums of A ***
  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(row) FIRSTPRIVATE(rank) ASYNC(accStreamIndex)
  DO col = 1, dim_col
     !$ACC LOOP SEQ
     DO row = 1, rank
        colsums(col) = colsums(col) + A(row, col)
     END DO
  END DO
  !$ACC END PARALLEL LOOP

! ****************************************************
! ***  Operate Omega on A                          ***
! ****************************************************

  ! Initialize last row of B  
  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) FIRSTPRIVATE(rank) ASYNC(accStreamIndex)
  DO col = 1, dim_col
    B(rank+1, col) = colsums(col) * normlast    
  END DO
  !$ACC END PARALLEL LOOP  

  ! Scale by NORMSUM
  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(accStreamIndex)
  DO col = 1, dim_col
    colsums(col) = normsum * colsums(col)
  END DO
  !$ACC END PARALLEL LOOP

  ! first rank rows
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(accStreamIndex)
  DO col = 1, dim_col
     DO row = 1, rank
        B(row, col) = A(row, col) - colsums(col)
     END DO
  END DO
  !$ACC END PARALLEL LOOP

! ********************
! *** FINISHING UP ***
! ********************

  ! !$ACC UPDATE HOST(B(:,:)) ASYNC(accStreamIndex)
  ! !$ACC EXIT DATA DELETE(A(:,:), B(:,:)) ASYNC(accStreamIndex)
  
  !$ACC EXIT DATA DELETE(colsums(:)) ASYNC(accStreamIndex)
  !$ACC WAIT(accStreamIndex)
  DEALLOCATE(colsums)

END SUBROUTINE PDAF_estkf_OmegaA
