SUBROUTINE xform(nbf, nprim, xformmat)
IMPLICIT NONE
!
CHARACTER :: temp
INTEGER :: i, j, nprim, nbf, sum, icount
INTEGER, ALLOCATABLE :: primfunctions(:)
REAL(KIND=8), ALLOCATABLE :: primcoeff(:)
REAL(KIND=8) :: xformmat(1:nbf,1:nprim)
!
OPEN(unit=9,file='primcoeff.dat',status='old')
!
READ(9,*) temp, nprim
READ(9,*) temp, nbf
!
ALLOCATE(primfunctions(1:nbf), primcoeff(1:nprim))
!
DO i=1, nprim
   READ(9,*) primcoeff(i)
ENDDO
CLOSE(9)
!
OPEN(unit=10,file='contracted.dat',status='old')
!
DO i=1, nbf
   READ(10,*) primfunctions(i)
ENDDO
CLOSE(10)
!
xformmat(:,:) = 0.d0
sum=0
icount=0
DO i=1, nbf
   j=primfunctions(i)
   IF (j==3) THEN
      sum=sum+1
      xformmat(i,sum) = primcoeff(sum-icount)
      xformmat(i,sum+3) = primcoeff(sum-icount+1)
      xformmat(i,sum+6) = primcoeff(sum-icount+2)
      icount=icount+1
   ELSE
      DO WHILE(j>0)
         sum=sum+1
         xformmat(i,sum) = primcoeff(sum)
         j=j-1
      ENDDO
   ENDIF
   IF (icount == 3) THEN
      icount=0
      sum=sum+6
   ENDIF
ENDDO
!
OPEN(unit=11,file='xform.dat',status='unknown')
DO i=1, nbf
   WRITE(11,*) (xformmat(i,j),j=1,nprim)
ENDDO
!
DEALLOCATE(primfunctions, primcoeff)
END SUBROUTINE
