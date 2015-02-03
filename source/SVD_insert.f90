!This subroutine makes the N^2 seperable potential evaluation 
!from an insertion of the Gaussian Basis potential
SUBROUTINE svd_insert(xi_pts,xi_wts,nbas_xi,eta_pts,eta_wts,nbas_eta,insertpot)
USE inputvars
IMPLICIT NONE

! Input variables
COMPLEX(KIND=DBL), INTENT(IN) :: xi_pts(1:n_xi)            ! xi quadrature points
COMPLEX(KIND=DBL), INTENT(IN) :: xi_wts(1:n_xi)            ! xi quadrature weights
COMPLEX(KIND=DBL), INTENT(IN) :: eta_pts(1:n_eta)          ! eta quadrature points
COMPLEX(KIND=DBL), INTENT(IN) :: eta_wts(1:n_eta)          ! eta quadrature weights
INTEGER, INTENT(IN) :: nbas_xi, nbas_eta                   ! Size of basis in both coordinates

! Output variable
COMPLEX(KIND=DBL), INTENT(OUT) :: insertpot(1:ntot,1:ntot) ! insertion potential

! Local variables
INTEGER :: i, j, k                             ! Loop variables
INTEGER :: ierror                              ! Error variable
INTEGER :: mfac=1

! Read-in matrices from Gaussian code
REAL(KIND=DBL) :: gausspotin(1:numgauss,1:numgauss) ! N^2 Seperable potential in Gaussians
REAL(KIND=DBL) :: overlaps(1:numgauss,1:numgauss)   ! Gaussian overlap matrix

! DVR - Primative Gaussian matrix
COMPLEX(KIND=DBL) :: gaussdvr(1:ntot,1:numprimg)
COMPLEX(KIND=DBL) :: conjgaussdvr(1:ntot,1:numprimg)

! DVR - contracted *****DEBUG ONLY*****
COMPLEX(KIND=DBL) :: basisdvr(1:ntot,1:numgauss)

! Transformation matrix to make contractions
REAL(KIND=DBL) :: xformmat(1:numgauss,1:numprimg)

! SVD variables
INTEGER :: norbits                             ! number of orthogonal orbitals
REAL(KIND=DBL) :: orthorbitals(1:numgauss,1:numgauss)
COMPLEX(KIND=DBL) :: unity(1:numgauss,1:numgauss)
REAL(KIND=DBL) :: inverse_overlaps(1:numgauss,1:numgauss)
REAL(KIND=DBL) :: diffoverlaps(1:numgauss,1:numgauss)
REAL(KIND=DBL) :: lgeigen

! Temp matrix
REAL(KIND=DBL) :: fnorm = 0.0d0
COMPLEX(KIND=DBL) :: smalltemp(1:numgauss,1:numgauss)
COMPLEX(KIND=DBL) :: largetemp(1:numprimg,1:numprimg)
COMPLEX(KIND=DBL) :: wrbig(1:ntot)
!!$!ZGEEV
!!$INTEGER :: ldvr=1, ldvl=1, lwork, info
!!$REAL(KIND=DBL) :: rwork(1:2*ntot)
!!$COMPLEX(KIND=DBL) :: wrbig(1:ntot), workspcbig(1:10*ntot), VL(1:ntot,1:ntot), VR(1:ntot,1:ntot)

! Read in gaussian exponents, overlaps and the potential energy in Gaussians
CALL svd_lagpot(numgauss,overlaps,gausspotin,tswitch)

! Create the Gaussian matrix evaluated in dvrs with the appropriate weights and factors
CALL gaussmat(xi_pts,xi_wts,nbas_xi,eta_pts,eta_wts,nbas_eta,gaussdvr,mfac)
WRITE(*,*) "Gaussmat done"

! Create the conjugate Gaussian matrix by complex conjugating the phi variable
mfac=-1
CALL gaussmat(xi_pts,xi_wts,nbas_xi,eta_pts,eta_wts,nbas_eta,conjgaussdvr,mfac)
WRITE(*,*) "Gaussmat conjugate done"

! Check to see if there is a contraction matrix, if so, read it in
!   else, set to the identity matrix
IF(numgauss == numprimg) THEN
   xformmat = 0.0d0
   DO i=1,numgauss
      xformmat(i,i) = 1.0d0
   ENDDO

ELSE
   OPEN(UNIT=95, FILE='xform.dat', STATUS='OLD', ACTION='READ', IOSTAT=ierror)
   IF(ierror /= 0) THEN
      WRITE(*,*) "Error reading xform.dat -- did you forget it? ierror = ", ierror
      STOP
   ELSE
      DO i=1, numgauss
         READ(95, *) (xformmat(i,j),j=1,numprimg)
      ENDDO
      CLOSE(95)
   ENDIF
ENDIF

WRITE(*,*) "xform done"

! CALL SVD routine to calculate the orthonormal orbitals or the pseudoinverse
CALL SVD_Ortho(numgauss, overlaps, inverse_overlaps, orthorbitals, norbits, SVD_tol)
WRITE(6,*) "SVD completed"


! Use SVD orbitals or inverse for the insertion?
!    SVDswitch == 1: orbitals
!    SVDswitch == 0: inverse
insertinverse: IF(svdswitch .EQ. 1) THEN
   WRITE(*,'(1x,"Using SVD orthonormal orbitals for the potential with a tolerance of ", ES15.6)')  SVD_tol
   WRITE(*,'(1x,"There are ", I5, " orthonormal orbitals")') norbits

   ! Evaluate the insertion multiplication
   smalltemp = MATMUL(TRANSPOSE(orthorbitals),MATMUL(gausspotin,orthorbitals))
   smalltemp = MATMUL(orthorbitals,MATMUL(smalltemp,TRANSPOSE(orthorbitals)))
   largetemp = MATMUL(TRANSPOSE(xformmat),MATMUL(smalltemp,xformmat))
   insertpot = MATMUL(conjgaussdvr,MATMUL(largetemp,TRANSPOSE(gaussdvr)))

!   unity = MATMUL(xformmat,MATMUL(TRANSPOSE(conjgaussdvr),MATMUL(gaussdvr,TRANSPOSE(xformmat))))

!   DO i=1, numgauss
!      WRITE(6362,'(5000ES16.6)')(REAL(unity(i,j)),j=1,numgauss)
!   ENDDO

!   DO i=1, numgauss
!      DO j=1, numgauss
!         diffoverlaps(i,j) = ABS(REAL(unity(i,j)) - overlaps(i,j))
!      ENDDO
!   ENDDO

!   CALL SVD_eigen(diffoverlaps,numgauss,lgeigen)

!   WRITE(*,'(" The diffence in overlaps is: ", ES16.6)') lgeigen

!   STOP
!!$   largetemp = MATMUL(TRANSPOSE(gaussdvr), conjgaussdvr)
!!$   smalltemp = MATMUL(xformmat,MATMUL(largetemp,TRANSPOSE(xformmat)))
!!$
!!$   DO i=1, numgauss
!!$      DO j=1, numgauss
!!$         fnorm = fnorm + (DREAL(smalltemp(i,J)) - overlaps(i,j))**2
!!$      ENDDO
!!$   ENDDO
!!$
!!$   WRITE(*,*) "Frobenius Norm = ", SQRT(fnorm)

!   smalltemp = MATMUL(inverse_overlaps,smalltemp)
!   smalltemp = MATMUL(orthorbitals,MATMUL(smalltemp,TRANSPOSE(orthorbitals)))

   ! Write out the orthonormal orbitals
   OPEN(UNIT=1015, FILE='orthorbitals.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)  
   IF(ierror == 0) THEN
      DO i=1, numgauss
         WRITE(1015,'(1x,800ES17.9)') (orthorbitals(i,k),k=1,numgauss)
      ENDDO
   ELSE
      WRITE(6,*) " Error opening orthorbital file to write, ierror = ", ierror
   ENDIF
   CLOSE(1015)

ELSE IF(svdswitch .EQ. 0) THEN insertinverse
   WRITE(*,*) "Using the SVD pseudoinverse for the potential"
   WRITE(*,*) "SVD tolerance is: ", SVD_tol

   !Create the inserted potential evaluated on the DVR
   smalltemp = MATMUL(inverse_overlaps,MATMUL(gausspotin,inverse_overlaps))
   largetemp = MATMUL(TRANSPOSE(xformmat),MATMUL(smalltemp,xformmat))
   insertpot = MATMUL(conjgaussdvr,MATMUL(largetemp,TRANSPOSE(gaussdvr)))


ELSE insertinverse
   WRITE(*,*) "Insertion potential is undefined"
   STOP
ENDIF insertinverse

! Insertion sucessful!
WRITE(*,*) "Potential Energy Inserted"

!!$! Write out the inserted potential
!!$OPEN(UNIT=1005, FILE='potentinsert.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
!!$IF(ierror == 0) THEN
!!$   DO i=1, ntot
!!$      WRITE(1005,'(1x,800ES17.9)') (insertpot(i,k),k=1,ntot)
!!$   ENDDO
!!$ELSE
!!$     WRITE(6,*) " Error opening orthorbital file to write, ierror = ", ierror
!!$ENDIF
!!$CLOSE(1005)


!!$lwork = 8*ntot
!!$CALL ZGEEV('N','N',ntot,insertpot,ntot,wrbig,VL,ldvl,VR,ldvr,workspcbig,lwork,rwork,info)
!!$
!!$ If there is an error
!!$IF(info .NE. 0) WRITE(6,921) info
!!$921 FORMAT(' problem in zgeev: info =',i5)
!!$
!!$OPEN(UNIT=1995, FILE='fullinsert.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
!!$IF(ierror == 0) THEN
!!$   DO i=1, ntot
!!$      WRITE(1995,'(1x,2ES17.9)') wrbig(i)
!!$   ENDDO
!!$ELSE
!!$     WRITE(6,*) " Error opening full insertion file to write, ierror = ", ierror
!!$ENDIF
!!$CLOSE(1995)
!!$STOP


!!$lwork = 8*ntot
!!$CALL ZGEEV('N','N',numgauss,smalltemp,numgauss,wrbig,VL,ldvl,VR,ldvr,workspcbig,lwork,rwork,info)
!!$
!!$! If there is an error
!!$IF(info .NE. 0) WRITE(6,921) info
!!$921 FORMAT(' problem in zgeev: info =',i5)
!!$
!!$OPEN(UNIT=1995, FILE='smallinsert.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
!!$IF(ierror == 0) THEN
!!$   DO i=1, numgauss
!!$      WRITE(1995,'(1x,2ES17.9)') wrbig(i)
!!$   ENDDO
!!$ELSE
!!$     WRITE(6,*) " Error opening full insertion file to write, ierror = ", ierror
!!$ENDIF
!!$CLOSE(1995)
!!$STOP

!!$! *********DEBUG**********
!!$
!!$OPEN(UNIT=1050, FILE='smalltemp.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
!!$IF(ierror == 0) THEN
!!$   DO i=1, numgauss
!!$      WRITE(1050,'(1x,800ES17.9)') (DREAL(smalltemp(i,k)),k=1,numgauss)
!!$   ENDDO
!!$ELSE
!!$     WRITE(6,*) " Error opening smalltemp file to write, ierror = ", ierror
!!$ENDIF
!!$CLOSE(1050)
!!$STOP
!!$
!!$OPEN(UNIT=1060, FILE='largetemp.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
!!$IF(ierror == 0) THEN
!!$   DO i=1, numprimg
!!$      WRITE(1060,'(1x,800ES17.9)') (largetemp(i,k),k=1,numprimg)
!!$   ENDDO
!!$ELSE
!!$     WRITE(6,*) " Error opening largetemp file to write, ierror = ", ierror
!!$ENDIF
!!$CLOSE(1060)
!!$
!!$OPEN(UNIT=1030, FILE='overlapsinverse.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
!!$IF(ierror == 0) THEN
!!$   DO i=1, numgauss
!!$      WRITE(1030,'(1x,800ES17.9)') (overlaps(i,k),k=1,numgauss)
!!$   ENDDO
!!$ELSE
!!$     WRITE(6,*) " Error opening overlapsinverse file to write, ierror = ", ierror
!!$ENDIF
!!$CLOSE(1030)
!!$
!!$OPEN(UNIT=1040, FILE='gpotin.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
!!$IF(ierror == 0) THEN
!!$   DO i=1, numgauss
!!$      WRITE(1040,'(1x,800ES17.9)') (gausspotin(i,k),k=1,numgauss)
!!$   ENDDO
!!$ELSE
!!$     WRITE(6,*) " Error opening gausspotin file to write, ierror = ", ierror
!!$ENDIF
!!$CLOSE(1040)
!!$
!!$OPEN(UNIT=1015, FILE='xformtest.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
!!$IF(ierror == 0) THEN
!!$   DO i=1, numgauss
!!$      WRITE(1015,'(1x,800ES17.9)') (xformmat(i,k),k=1,numprimg)
!!$   ENDDO
!!$ELSE
!!$     WRITE(6,*) " Error opening xformtest file to write, ierror = ", ierror
!!$ENDIF
!!$CLOSE(1015)
!!$
!!$OPEN(UNIT=1025, FILE='gaussdvr.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
!!$IF(ierror == 0) THEN
!!$   DO i=1, ntot
!!$      WRITE(1025,'(1x,800ES17.9)') (gaussdvr(i,k),k=1,numprimg)
!!$   ENDDO
!!$ELSE
!!$     WRITE(6,*) " Error opening gaussdvr file to write, ierror = ", ierror
!!$ENDIF
!!$CLOSE(1025)
!!$
!!$OPEN(UNIT=1026, FILE='congaussdvr.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
!!$IF(ierror == 0) THEN
!!$   DO i=1, ntot
!!$      WRITE(1026,'(1x,800ES17.9)') (conjgaussdvr(i,k),k=1,numprimg)
!!$   ENDDO
!!$ELSE
!!$     WRITE(6,*) " Error opening conjgaussdvr file to write, ierror = ", ierror
!!$ENDIF
!!$CLOSE(1026)
!!$
!!$basisdvr = MATMUL(gaussdvr,TRANSPOSE(xformmat))
!!$OPEN(UNIT=1010, FILE='basisdvr.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
!!$IF(ierror == 0) THEN
!!$   DO i=1, ntot
!!$      WRITE(1010,'(1x,800ES17.9)') (basisdvr(i,k),k=1,numgauss)
!!$   ENDDO
!!$ELSE
!!$     WRITE(6,*) " Error opening basisdvr file to write, ierror = ", ierror
!!$ENDIF
!!$CLOSE(1010)
!!$! *********DEBUG**********
ENDSUBROUTINE
