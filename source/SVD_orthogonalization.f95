SUBROUTINE SVD_Ortho(numsize, matrixin, matinverse, orthomat, northo, SVD_tol)
!
! This subroutine orthogonalizes a real symmetric matrix and gives back an inverse
! along with the orthogonalizing vectors.
!
! This subroutine uses SVD
! S = U * Sigma * V^T
!   = U * Sigma * U^T
! since S is real symmetric.
!
IMPLICIT NONE

! Define Double Precision
INTEGER, PARAMETER :: DBL = Selected_Real_Kind(p=13,r=200)

! Input variables
INTEGER, INTENT(IN) :: numsize                                ! Size of inputed square matrix
REAL(KIND=DBL), INTENT(IN) :: SVD_tol                         ! Tolerance value for singular values
REAL(KIND=DBL), INTENT(IN) :: matrixin(1:numsize,1:numsize)   ! Input matrix

! Output variables
INTEGER, INTENT(OUT) :: northo                                 ! Number of SVD values above tolerance
REAL(KIND=DBL), INTENT(OUT) :: matinverse(1:numsize,1:numsize) ! Matrix inverse with Singular values set to 0
REAL(KIND=DBL), INTENT(OUT) :: orthomat(1:numsize,1:numsize)   ! Orthonormalizating matrix U

! SVD variables
CHARACTER(len=1) :: jobu = "A", jobvt = "A"               ! Flags to output all of the SVD
INTEGER :: lwork                                          ! Size of work array (10*numsize)
INTEGER :: info                                           ! Error variable
REAL(KIND=DBL) :: work(1:10*numsize)                      ! Work array
REAL(KIND=DBL) :: sing_values(1:numsize)                  ! Singular values in decreasing order
REAL(KIND=DBL) :: umat(1:numsize,1:numsize)               ! Singular value decomposition matrices
REAL(KIND=DBL) :: vtmat(1:numsize,1:numsize)
REAL(KIND=DBL) :: matcopy(1:numsize,1:numsize)            ! Copy of the matrixin since SVD destroys it

! Loop variables and error
INTEGER :: i, j, k, ierror
REAL(KIND=DBL) :: sum

! DEBUG VARIABLES
REAL(KIND=DBL) :: tempmat(1:numsize,1:numsize) 

! Define the SVD variables and make a copy of the input matrix
lwork = 10*numsize
matcopy = matrixin

! CALL SVD LAPACK routine
CALL DGESVD(jobu, jobvt, numsize, numsize, matcopy, numsize, sing_values, umat, numsize, vtmat, numsize, work, lwork, info)

! IF SVD failes, return with error
IF(info /= 0) THEN
   WRITE(6, '(" SVD failed, info = ", I4)') info
   RETURN

! IF SVD succeeds, write information to file
ELSE
   OPEN(UNIT=20, FILE='SVD.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)

   ! If the opening of SVD.out resulted in an error
   IF(ierror /= 0) THEN
      WRITE(6, '("Error opening SVD.out  ierror= ",I4 )') ierror
      RETURN
   ENDIF
   
   ! Write out the Singular values
   WRITE(20,'(" SVD computed")')
   WRITE(20,'(" Singular values  ",/,(g21.10))') (sing_values(i),i=1,numsize)

   ! Create orthonormal orbitals
   WRITE(20,*) "Dropping vectors with singular value less than or equal to ",SVD_tol
   northo = 0
   DO j=1,numsize
      IF(sing_values(j) > SVD_tol) THEN
         northo = northo + 1
         DO i=1,numsize
            orthomat(i,j) = umat(i,j)/SQRT(sing_values(j))
         ENDDO
      ENDIF
   ENDDO
   WRITE(20,*) " constructed", northo," orthogonal orbitals"
   
   ! Create the inverse matrix  
   DO i=1, numsize
      DO j=1, numsize
         sum=0.0d0
         DO k=1, numsize
            IF(sing_values(k) > SVD_tol) THEN
               sum = sum+umat(i,k)*vtmat(k,j)/sing_values(k)
            ENDIF
         ENDDO
         matinverse(i,j) = sum
      ENDDO
   ENDDO
  
   CLOSE(20)

!*****************************DEBUGS******************************
! CHECK that the SVD gives back the input matrix
OPEN(UNIT=21, FILE='SVD_redo.out', STATUS='UNKNOWN', ACTION='WRITE')
tempmat = 0.0d0
DO i=1, numsize
   DO j=1, numsize
      sum=0.0d0
      DO k=1, numsize
         sum = sum + umat(i,k)*sing_values(k)*vtmat(k,j)
      ENDDO
   tempmat(i,j) = sum
   ENDDO
ENDDO
DO i=1, numsize
   WRITE(21, '(1x,100ES15.7)') (tempmat(i,j),j=1,numsize)
ENDDO
CLOSE(21)

!!$! CHECK that Umat is unitary gives back the input matrix
!!$OPEN(UNIT=22, FILE='SVD_umat.out', STATUS='UNKNOWN', ACTION='WRITE')
!!$DO i=1, numsize
!!$   WRITE(22, '(1x,100ES15.7)') (umat(i,j),j=1,numsize)
!!$ENDDO
!!$CLOSE(22)
!!$
!!$OPEN(UNIT=23, FILE='SVD_umat1.out', STATUS='UNKNOWN', ACTION='WRITE')
!!$tempmat = MATMUL(Transpose(umat),umat)
!!$DO i=1, numsize
!!$   WRITE(23, '(1x,100ES15.7)') (tempmat(i,j),j=1,numsize)
!!$ENDDO
!!$CLOSE(23)
!!$
!!$OPEN(UNIT=24, FILE='SVD_orthomat1.out', STATUS='UNKNOWN', ACTION='WRITE')
!!$tempmat = MATMUL(TRANSPOSE(orthomat),orthomat)
!!$DO i=1, numsize
!!$   WRITE(24, '(1x,100ES15.7)') (tempmat(i,j),j=1,numsize)
!!$ENDDO
!!$CLOSE(24)
!!$
!!$! Check that the u-matrix is orthonormalized
!!$OPEN(UNIT=25, FILE='SVD_orthomatunitary.out', STATUS='UNKNOWN', ACTION='WRITE')
!!$tempmat = MATMUL(TRANSPOSE(orthomat),MATMUL(matrixin,orthomat))
!!$DO i=1, numsize
!!$   WRITE(25, '(1x,100ES15.7)') (tempmat(i,j),j=1,numsize)
!!$ENDDO
!!$CLOSE(25)
!!$
!!$! Check the psuedoinverse
!!$OPEN(UNIT=26, FILE='SVD_inverse.out', STATUS='UNKNOWN', ACTION='WRITE')
!!$tempmat = MATMUL(matinverse,matrixin)
!!$DO i=1, numsize
!!$   WRITE(26, '(1x,100ES15.7)') (tempmat(i,j),j=1,numsize)
!!$ENDDO
!!$CLOSE(26)
!!$
!!$!*****************************DEBUGS******************************

ENDIF
END SUBROUTINE


