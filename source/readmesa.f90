SUBROUTINE readmesa(filename, ndim, matrix)
IMPLICIT NONE
!
CHARACTER(LEN=12) :: filename
INTEGER :: ndim, hold
INTEGER :: i, j 
REAL(KIND=8) :: matrix(1:ndim,1:ndim)
!
OPEN(unit=10, file=filename,status='old')
!
DO i=1, ndim
   READ(10,*) (matrix(i,j),j=1,ndim)
ENDDO
!
CLOSE(10)
!
ENDSUBROUTINE
