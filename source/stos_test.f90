SUBROUTINE stos(xi, eta, a, numgauss, valarray, m)
! This subroutine evaluates the overlap of a DVR function in prolate spherodial
! coordinates with a primative Gaussian function read from a file: expos.dat.
!
! This 3-D integral is done via the quadrature rule in both xi and eta while
! the integral in phi is done explicitly.
! The phi integral depends on m, and the complex conjugation needs to be done in this
! step because with a complex scaled xi results in a Gaussian function that DOES NOT
! become conjugated (that is, the Hamiltonian is not Hermitian usinc complex scaling) 
! except for the phi integral which is done exactly.
!
IMPLICIT NONE

! Parameters
INTEGER, PARAMETER :: DBL =  SELECTED_REAL_KIND(13,200)         ! Define Double Precision (precision,range)
REAL(KIND=DBL), PARAMETER :: pi = 3.14159265358979
REAL(KIND=DBL),PARAMETER :: sqtpi = 1.77245385090552
COMPLEX(KIND=DBL), PARAMETER :: ci = (0.0d0,1.0d0)

! Input variables
INTEGER, INTENT(IN) :: numgauss                          ! Number of primative Gaussians
COMPLEX(KIND=DBL), INTENT(IN) :: xi                      ! Quadrature point in xi
COMPLEX(KIND=DBL), INTENT(IN) :: eta                     ! Quadrature point in eta
REAL(KIND=DBL), INTENT(IN) :: a                          ! Internuclear distance/2
INTEGER, INTENT(IN) :: m                       ! Angular quantum number

! Output variable
COMPLEX(KIND=DBL), INTENT(OUT) :: valarray(1:numgauss)   ! Array with overlap of a DVR function
                                                         ! with each primative Gaussian

! Coordinate Variables
COMPLEX(KIND=DBL) :: xy                       ! Variable for x and y without phi
COMPLEX(KIND=DBL) :: rsq_2d                   ! Variable for x^2, y^2 and x^2 + y^2 without phi
COMPLEX(KIND=DBL) :: z                        ! Variable for z


! Read variables
INTEGER :: ierror, iread                      ! Read error variables
REAL(KIND=DBL) :: zoff                        ! The z-offset of the input Gaussian
REAL(KIND=DBL) :: expo                        ! Exponent of the Gaussian
CHARACTER(LEN=2) :: sym                       ! Read in the 2 character symmetry

! Local variables
INTEGER :: i                                  ! Loop variable
COMPLEX(KIND=DBL) :: evalue                   ! Temp value of the integral under quadrature

! Calculate the Cartesian variables
xy = a*SQRT((xi*xi-1.0d0)*(1.0d0-eta*eta))
rsq_2d=a*a*(xi*xi-1.0d0)*(1.0d0-eta*eta)
z = a*xi*eta

! Open the file!
OPEN(UNIT=10, FILE="expos.dat", STATUS="OLD", ACTION="READ", IOSTAT=ierror)
fileopen: IF (ierror == 0) THEN

   readeval:DO i=1, numgauss
      
      ! Read the ith entry in expos.dat
      READ(10,*,IOSTAT=iread) zoff, sym, expo

      ! Print message if there is a read error
      IF (iread /= 0) THEN
         WRITE(*,100) i, REAL(xi), AIMAG(xi), REAL(eta)
100      FORMAT (1x,"Error reading entry: ", I3, " from expos.dat for xi= ", 2ES15.6, " and eta= ", ES15.6)
         WRITE(*,101) zoff, sym, expo, iread
101      FORMAT (1x,"zoff= ", ES15.6, " sym= ", A, " expo= ", ES15.6, " iread= ", I10)
      ENDIF
      
      ! Evaluate the Gaussian expoential
      evalue = (2.0d0/pi)**(0.75)*EXP(-expo*(rsq_2d+(z-zoff)*(z-zoff)))

      ! Cases for the angular parts
      IF(sym .EQ. 'ss' .AND. m .EQ. 0) THEN
         evalue = evalue*SQRT(2.0d0*pi)

      ELSEIF(sym .EQ. 'px' .AND. m .EQ. 1) THEN
         evalue = xy*evalue*SQRT(0.5d0*pi)*2.0d0*expo**(1.25)
      ELSEIF(sym .EQ. 'px' .AND. m .EQ. -1) THEN
         evalue = xy*evalue*SQRT(0.5d0*pi)*2.0d0*expo**(1.25)      
      ELSEIF(sym .EQ. 'py' .AND. m .EQ. 1) THEN
         evalue = xy*evalue*ci*SQRT(0.5d0*pi)*2.0d0*expo**(1.25)  
      ELSEIF(sym .EQ. 'py' .AND. m .EQ. -1) THEN
         evalue = -xy*evalue*ci*SQRT(0.5d0*pi)*2.0d0*expo**(1.25)  
      ELSEIF(sym .EQ. 'pz' .AND. m .EQ. 0) THEN
         evalue = (z-zoff)*evalue*SQRT(2.0d0*pi)*2.0d0*expo**(1.25)

      ELSEIF(sym .EQ. 'xx' .AND. m .EQ. 0) THEN
         evalue = rsq_2d*evalue*SQRT(0.5d0*pi)*4.0d0*expo**(1.75)
      ELSEIF(sym .EQ. 'xx' .AND. m .EQ. 2) THEN
         evalue = rsq_2d*evalue*SQRT(0.5d0*pi)*2.0d0*expo**(1.75)
      ELSEIF(sym .EQ. 'xx' .AND. m .EQ. -2) THEN
         evalue = rsq_2d*evalue*SQRT(0.5d0*pi)*2.0d0*expo**(1.75)

      ELSEIF(sym .EQ. 'yy' .AND. m .EQ. 0) THEN
         evalue = rsq_2d*evalue*SQRT(0.5d0*pi)*4.0d0*expo**(1.75)
      ELSEIF(sym .EQ. 'yy' .AND. m .EQ. 2) THEN
         evalue = -rsq_2d*evalue*SQRT(0.5d0*pi)*2.0d0*expo**(1.75)
      ELSEIF(sym .EQ. 'yy' .AND. m .EQ. -2) THEN
         evalue = -rsq_2d*evalue*SQRT(0.5d0*pi)*2.0d0*expo**(1.75)

      ELSEIF(sym .EQ. 'zz' .AND. m .EQ. 0) THEN
         evalue = (z-zoff)*(z-zoff)*evalue*SQRT(2.0d0*pi)*4.0d0*expo**(1.75)

      ELSEIF(sym .EQ. 'xy' .AND. m .EQ. 2) THEN
         evalue = rsq_2d*evalue*ci*SQRT(0.5d0*pi)*2.0d0*expo**(1.75)
      ELSEIF(sym .EQ. 'xy' .AND. m .EQ. -2) THEN
         evalue = -rsq_2d*evalue*ci*SQRT(0.5d0*pi)*2.0d0*expo**(1.75)
      ELSEIF(sym .EQ. 'xz' .AND. m .EQ. 1) THEN
         evalue = xy*(z-zoff)*evalue*SQRT(0.5d0*pi)*4.0d0*expo**(1.75)
      ELSEIF(sym .EQ. 'xz' .AND. m .EQ. -1) THEN
         evalue = xy*(z-zoff)*evalue*SQRT(0.5d0*pi)*4.0d0*expo**(1.75)
      ELSEIF(sym .EQ. 'yz' .AND. m .EQ. 1) THEN
         evalue = xy*(z-zoff)*evalue*ci*SQRT(0.5d0*pi)*4.0d0*expo**(1.75)
      ELSEIF(sym .EQ. 'yz' .AND. m .EQ. -1) THEN
         evalue = -xy*(z-zoff)*evalue*ci*SQRT(0.5d0*pi)*4.0d0*expo**(1.75)
      ELSE
         evalue = (0.0d0,0.0d0)
      ENDIF

      valarray(i) = evalue
!      WRITE(1909,'(1x,I3,6ES16.7E2,1x,A3,5(ES15.6E3,1x))') m,xi,eta,zoff,a,sym,expo,evalue
   ENDDO readeval

! If there was a read error, let me know!
ELSE fileopen
   WRITE(*,*) " Error opening expos.dat for xi= ", REAL(xi), AIMAG(xi), " and eta= ", REAL(eta)
ENDIF fileopen
CLOSE(10)
ENDSUBROUTINE
