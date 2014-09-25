SUBROUTINE diagwrap(nmax, matin, wrbig)
!
! This subroutine is a wrapper for ZGEEV
!

! Variables for Hamiltonian diagonalization
INTEGER, PARAMETER :: DBL =  SELECTED_REAL_KIND(13,200)         ! Define Double Precision (precision,range)
INTEGER, INTENT(IN) :: nmax
COMPLEX(KIND=DBL), INTENT(IN) :: matin(1:nmax,1:nmax)
COMPLEX(KIND=DBL), INTENT(OUT) :: wrbig(1:nmax)

INTEGER :: ldvr=1, ldvl=1, info
INTEGER :: lwork
REAL(KIND=DBL):: rwork(1:2*nmax)
COMPLEX(KIND=DBL) :: workspcbig(1:10*nmax), VL(1:nmax,1:nmax), VR(1:nmax,1:nmax)
lwork = 10*nmax


CALL ZGEEV('N','N',nmax,matin,nmax,wrbig,VL,ldvl,VR,ldvr,workspcbig,lwork,rwork,info)

IF(info == 0) THEN
   OPEN(UNIT=10, FILE='wrapeigen.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)  
   IF(ierror == 0) THEN
      DO i=1, nmax
         WRITE(10,'(1x,2ES19.9E4)') wrbig(i)
      ENDDO
   ELSE
      WRITE(6,*) " Error opening wrapeigen, ierror = ", ierror
   ENDIF
   CLOSE(10)
ELSE
   WRITE(*,*) "Wrapper for ZGEEV failed", info
ENDIF

END SUBROUTINE diagwrap
