SUBROUTINE SVD_eigen(matrixin, numsize, bigeigen)
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
REAL(KIND=DBL), INTENT(IN) :: matrixin(1:numsize,1:numsize)   ! Input matrix

! Output variables
REAL(KIND=DBL), INTENT(OUT) :: bigeigen   ! Largest eigenvalue

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
INTEGER :: i, ierror

! DEBUG VARIABLES
REAL(KIND=DBL) :: tempmat(1:numsize,1:numsize) 

! Define the SVD variables and make a copy of the input matrix
lwork = 10*numsize
matcopy = matrixin

! CALL SVD LAPACK routine
CALL DGESVD(jobu, jobvt, numsize, numsize, matcopy, numsize, sing_values, umat, numsize, vtmat, numsize, work, lwork, info)

! IF SVD failes, return with error
IF(info /= 0) THEN
   WRITE(6, '(" SVD failed in SVD_eigen, info = ", I4)') info
   RETURN

! IF SVD succeeds, write information to file
ELSE
   OPEN(UNIT=20, FILE='SVD_eigen.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)

   ! If the opening of SVD.out resulted in an error
   IF(ierror /= 0) THEN
      WRITE(6, '("Error opening SVD_eigen.out  ierror= ",I4 )') ierror
      RETURN
   ENDIF
   
   ! Write out the Singular values
   WRITE(20,'(" SVD computed")')
   WRITE(20,'(" Singular values  ",/,(g21.10))') (sing_values(i),i=1,numsize)

   bigeigen = sing_values(1)
 
   CLOSE(20)

ENDIF
END SUBROUTINE


