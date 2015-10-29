!This subroutine constructs the matrix of the Laguerre potential energy
!from the old code's lower diagonal matrix
SUBROUTINE svd_lagpot(numgauss,overlaps,matrixin,tswitch,vswitch)
IMPLICIT NONE

INTEGER, PARAMETER :: DBL =  SELECTED_REAL_KIND(13,200)         ! Define Double Precision (precision,range)
INTEGER, INTENT(IN) :: tswitch

INTEGER ::numgauss, vswitch
REAL(KIND=DBL) :: matrixin(1:numgauss, 1:numgauss)
REAL(KIND=DBL) :: direct(1:numgauss, 1:numgauss)
REAL(KIND=DBL) :: exchange(1:numgauss, 1:numgauss)
REAL(KIND=DBL) :: vnuc(1:numgauss, 1:numgauss)
REAL(KIND=DBL) :: kinetic(1:numgauss, 1:numgauss)
REAL(KIND=DBL) :: overlaps(1:numgauss, 1:numgauss)

!!$! Diagonal debug
!!$INTEGER :: lwork, info, ierror
!!$REAL(KIND=DBL) :: w(1:numgauss), work(1:10*numgauss)
!!$lwork=10*numgauss

kinetic = 0.0d0
IF(tswitch == 2) THEN
   CALL readmesa('kineticx.dat',numgauss,kinetic)
ENDIF

CALL readmesa('overlaps.dat',numgauss,overlaps)

CALL readmesa('directxx.dat',numgauss,direct)

CALL readmesa('exchange.dat',numgauss,exchange)

IF(vswitch == 2) THEN
   CALL readmesa('vnucxxxx.dat',numgauss,vnuc)
ELSE IF(vswitch == 3) THEN
   CALL readmesa('vnucbeck.dat',numgauss,vnuc)
ENDIF

matrixin = kinetic + vnuc + direct - 0.5d0*exchange

!OPEN(UNIT=9, FILE='overlap.out', STATUS='UNKNOWN', ACTION='WRITE')
!   DO i=1, numgauss
!      WRITE(9,'(1x,800ES17.9)') (overlaps(i,j),j=1,numgauss)
!   ENDDO
!CLOSE(9)

!!$OPEN(UNIT=9, FILE='diagdebug.dat', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
!!$CALL DSYGV(1, 'N', 'L', numgauss, matrixin, numgauss, overlaps, numgauss, w, work, lwork, info)
!!$DO i=1, numgauss
!!$   WRITE(9,'(ES17.9)') w(i)
!!$ENDDO
!!$
!!$STOP

ENDSUBROUTINE
