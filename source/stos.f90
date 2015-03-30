SUBROUTINE stos(xi, eta, a, numgauss, valarray, m, qchemswitch)
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
REAL(KIND=DBL), PARAMETER :: sqrtpi2 = 1.2533141373155
COMPLEX(KIND=DBL), PARAMETER :: ci = (0.0d0,1.0d0)

! Input variables
INTEGER, INTENT(IN) :: numgauss                          ! Number of primative Gaussians
COMPLEX(KIND=DBL), INTENT(IN) :: xi                      ! Quadrature point in xi
COMPLEX(KIND=DBL), INTENT(IN) :: eta                     ! Quadrature point in eta
REAL(KIND=DBL), INTENT(IN) :: a                          ! Internuclear distance/2
INTEGER, INTENT(IN) :: m                                 ! Angular quantum number
INTEGER, INTENT(IN) :: qchemswitch                       ! What program generated the input?

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
CHARACTER(LEN=4) :: sym                       ! Read in the 2 character symmetry

! Local variables
INTEGER :: i                                  ! Loop variable
COMPLEX(KIND=DBL) :: gvalue                   ! Value of the Guassian part
REAL(KIND=DBL) :: norm                        ! Normalization
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
      gvalue = EXP(-expo*(rsq_2d+(z-zoff)*(z-zoff)))

      SELECT CASE(ABS(m))
         CASE(0)
            IF(sym .EQ. 'ss') THEN
               norm = (2.0d0*expo/pi)**(0.75)
               evalue = norm*gvalue*SQRT(2.0d0*pi)

            ELSEIF(sym .EQ. 'pz') THEN
               norm = 2.0d0*(2.0d0/pi)**(0.75)*expo**(1.25)
               evalue = norm*(z-zoff)*gvalue*SQRT(2.0d0*pi)

            ELSEIF(sym .EQ. 'xx') THEN
               norm = 4.0d0*(2.0d0/pi)**(0.75)*expo**(1.75)
               evalue = norm*rsq_2d*gvalue*sqrtpi2
            ELSEIF(sym .EQ. 'yy') THEN
               norm = 4.0d0*(2.0d0/pi)**(0.75)*expo**(1.75)
               evalue = norm*rsq_2d*gvalue*sqrtpi2
            ELSEIF(sym .EQ. 'zz') THEN
               norm = 4.0d0*(2.0d0/pi)**(0.75)*expo**(1.75)
               evalue = norm*(z-zoff)*(z-zoff)*gvalue*SQRT(2.0d0*pi)

            ELSEIF(sym .EQ. 'zzz') THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = norm*(z-zoff)*(z-zoff)*(z-zoff)*gvalue*SQRT(2.0d0*pi)
            ELSEIF(sym .EQ. 'xxz') THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = norm*rsq_2d*(z-zoff)*gvalue*sqrtpi2
            ELSEIF(sym .EQ. 'yyz') THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = norm*rsq_2d*(z-zoff)*gvalue*sqrtpi2

            ELSEIF(sym .EQ. 'xxxx') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*rsq_2d*gvalue*sqrtpi2*(0.75d0)
            ELSEIF(sym .EQ. 'yyyy') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*rsq_2d*gvalue*sqrtpi2*(0.75d0)
            ELSEIF(sym .EQ. 'xxyy') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*rsq_2d*gvalue*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'zzzz') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*(z-zoff)*(z-zoff)*(z-zoff)*(z-zoff)*gvalue*SQRT(2.0d0*pi)
            ELSEIF(sym .EQ. 'xxzz') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*(z-zoff)*(z-zoff)*gvalue*sqrtpi2
            ELSEIF(sym .EQ. 'yyzz') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*(z-zoff)*(z-zoff)*gvalue*sqrtpi2
               

            ELSE
                evalue = (0.0d0,0.0d0)
            ENDIF


         CASE(1)
            IF(sym .EQ. 'px') THEN
               norm = 2.0d0*(2.0d0/pi)**(0.75)*expo**(1.25)
               evalue = norm*xy*gvalue*sqrtpi2
            ELSEIF(sym .EQ. 'py' .AND. m .EQ. 1) THEN
               norm = 2.0d0*(2.0d0/pi)**(0.75)*expo**(1.25)
               evalue = norm*xy*gvalue*ci*sqrtpi2
            ELSEIF(sym .EQ. 'py' .AND. m .EQ. -1) THEN
               norm = 2.0d0*(2.0d0/pi)**(0.75)*expo**(1.25)
               evalue = -norm*xy*gvalue*ci*sqrtpi2

            ELSEIF(sym .EQ. 'xz') THEN
               norm = 4.0d0*(2.0d0/pi)**(0.75)*expo**(1.75)
               evalue = norm*xy*(z-zoff)*gvalue*sqrtpi2
            ELSEIF(sym .EQ. 'yz' .AND. m .EQ. 1) THEN
               norm = 4.0d0*(2.0d0/pi)**(0.75)*expo**(1.75)
               evalue = norm*xy*(z-zoff)*gvalue*ci*sqrtpi2
            ELSEIF(sym .EQ. 'yz' .AND. m .EQ. -1) THEN
               norm = 4.0d0*(2.0d0/pi)**(0.75)*expo**(1.75)
               evalue = -norm*xy*(z-zoff)*gvalue*ci*sqrtpi2

            ELSEIF(sym .EQ. 'xxx') THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = norm*rsq_2d*xy*gvalue*sqrtpi2*(0.75d0)
            ELSEIF(sym .EQ. 'yyy' .AND. m .EQ. 1) THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = norm*rsq_2d*xy*gvalue*ci*sqrtpi2*(0.75d0)
            ELSEIF(sym .EQ. 'yyy' .AND. m .EQ. -1) THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = -norm*rsq_2d*xy*gvalue*ci*sqrtpi2*(0.75d0)
            ELSEIF(sym .EQ. 'xxy' .AND. m .EQ. 1) THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = norm*rsq_2d*xy*gvalue*ci*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'xxy' .AND. m .EQ. -1) THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = -norm*rsq_2d*xy*gvalue*ci*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'xyy') THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = norm*rsq_2d*xy*gvalue*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'xzz') THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = norm*(z-zoff)*(z-zoff)*xy*gvalue*sqrtpi2
            ELSEIF(sym .EQ. 'yzz' .AND. m .EQ. 1) THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = norm*(z-zoff)*(z-zoff)*xy*gvalue*ci*sqrtpi2
            ELSEIF(sym .EQ. 'yzz' .AND. m .EQ. -1) THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = -norm*(z-zoff)*(z-zoff)*xy*gvalue*ci*sqrtpi2

            ELSEIF(sym .EQ. 'xxxz') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*xy*(z-zoff)*gvalue*sqrtpi2*(0.75d0)
            ELSEIF(sym .EQ. 'yyyz' .AND. m .EQ. 1) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*xy*(z-zoff)*gvalue*ci*sqrtpi2*(0.75d0)
            ELSEIF(sym .EQ. 'yyyz' .AND. m .EQ. -1) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = -norm*rsq_2d*xy*(z-zoff)*gvalue*ci*sqrtpi2*(0.75d0)
            ELSEIF(sym .EQ. 'xxyz' .AND. m .EQ. 1) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*xy*(z-zoff)*gvalue*ci*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'xxyz' .AND. m .EQ. -1) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = -norm*rsq_2d*xy*(z-zoff)*gvalue*ci*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'xyyz') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*xy*(z-zoff)*gvalue*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'xzzz') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*(z-zoff)*(z-zoff)*(z-zoff)*xy*gvalue*sqrtpi2
            ELSEIF(sym .EQ. 'yzzz' .AND. m .EQ. 1) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*(z-zoff)*(z-zoff)*(z-zoff)*xy*gvalue*ci*sqrtpi2
            ELSEIF(sym .EQ. 'yzzz' .AND. m .EQ. -1) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = -norm*(z-zoff)*(z-zoff)*(z-zoff)*xy*gvalue*ci*sqrtpi2

            ELSE
                evalue = (0.0d0,0.0d0)
            ENDIF


         CASE(2)
            IF(sym .EQ. 'xx') THEN
               norm = 4.0d0*(2.0d0/pi)**(0.75)*expo**(1.75)
               evalue = norm*rsq_2d*gvalue*sqrtpi2*(0.5d0)
            ELSEIF(sym .EQ. 'yy') THEN
               norm = 4.0d0*(2.0d0/pi)**(0.75)*expo**(1.75)
               evalue = -norm*rsq_2d*gvalue*sqrtpi2*(0.5d0)
            ELSEIF(sym .EQ. 'xy' .AND. m .EQ. 2) THEN
               norm = 4.0d0*(2.0d0/pi)**(0.75)*expo**(1.75)
               evalue = norm*rsq_2d*gvalue*ci*sqrtpi2*(0.5d0)
            ELSEIF(sym .EQ. 'xy' .AND. m .EQ. -2) THEN
               norm = 4.0d0*(2.0d0/pi)**(0.75)*expo**(1.75)
               evalue = -norm*rsq_2d*gvalue*ci*sqrtpi2*(0.5d0)

            ELSEIF(sym .EQ. 'xxz') THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = norm*rsq_2d*(z-zoff)*gvalue*sqrtpi2*(0.5d0)
            ELSEIF(sym .EQ. 'yyz') THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = -norm*rsq_2d*(z-zoff)*gvalue*sqrtpi2*(0.5d0)
            ELSEIF(sym .EQ. 'xyz' .AND. m .EQ. 2) THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = norm*rsq_2d*(z-zoff)*gvalue*ci*sqrtpi2*(0.5d0)            
            ELSEIF(sym .EQ. 'xyz' .AND. m .EQ. -2) THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = -norm*rsq_2d*(z-zoff)*gvalue*ci*sqrtpi2*(0.5d0)

            ELSEIF(sym .EQ. 'xxxx') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*rsq_2d*gvalue*sqrtpi2*(0.5d0)
            ELSEIF(sym .EQ. 'yyyy') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = -norm*rsq_2d*rsq_2d*gvalue*sqrtpi2*(0.5d0)         
            ELSEIF(sym .EQ. 'xxxy' .AND. m .EQ. 2) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*rsq_2d*gvalue*ci*sqrtpi2*(0.25d0)          
            ELSEIF(sym .EQ. 'xxxy' .AND. m .EQ. -2) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = -norm*rsq_2d*rsq_2d*gvalue*ci*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'xyyy' .AND. m .EQ. 2) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*rsq_2d*gvalue*ci*sqrtpi2*(0.25d0)          
            ELSEIF(sym .EQ. 'xyyy' .AND. m .EQ. -2) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = -norm*rsq_2d*rsq_2d*gvalue*ci*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'xxzz') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*(z-zoff)*(z-zoff)*gvalue*sqrtpi2*(0.5d0)
            ELSEIF(sym .EQ. 'yyzz') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = -norm*rsq_2d*(z-zoff)*(z-zoff)*gvalue*sqrtpi2*(0.5d0)
            ELSEIF(sym .EQ. 'xyzz' .AND. m .EQ. 2) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*(z-zoff)*(z-zoff)*gvalue*ci*sqrtpi2*(0.5d0)            
            ELSEIF(sym .EQ. 'xyzz' .AND. m .EQ. -2) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = -norm*rsq_2d*(z-zoff)*(z-zoff)*gvalue*ci*sqrtpi2*(0.5d0)

            ELSE
                evalue = (0.0d0,0.0d0)
            ENDIF


         CASE(3)
            IF(sym .EQ. 'xxx') THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = norm*rsq_2d*xy*gvalue*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'yyy' .AND. m .EQ. 3) THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = -norm*rsq_2d*xy*gvalue*ci*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'yyy' .AND. m .EQ. -3) THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = norm*rsq_2d*xy*gvalue*ci*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'xxy' .AND. m .EQ. 3) THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = norm*rsq_2d*xy*gvalue*ci*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'xxy' .AND. m .EQ. -3) THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = -norm*rsq_2d*xy*gvalue*ci*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'xyy') THEN
               norm = 8.0d0*(2.0d0/pi)**(0.75)*expo**(2.25)
               evalue = -norm*rsq_2d*xy*gvalue*sqrtpi2*(0.25d0)

            ELSEIF(sym .EQ. 'xxxz') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*xy*(z-zoff)*gvalue*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'yyyz' .AND. m .EQ. 3) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = -norm*rsq_2d*xy*(z-zoff)*gvalue*ci*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'yyyz' .AND. m .EQ. -3) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*xy*(z-zoff)*gvalue*ci*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'xxyz' .AND. m .EQ. 3) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*xy*(z-zoff)*gvalue*ci*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'xxyz' .AND. m .EQ. -3) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = -norm*rsq_2d*xy*(z-zoff)*gvalue*ci*sqrtpi2*(0.25d0)
            ELSEIF(sym .EQ. 'xyyz') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = -norm*rsq_2d*xy*(z-zoff)*gvalue*sqrtpi2*(0.25d0)
          
            ELSE
                evalue = (0.0d0,0.0d0)
            ENDIF

         CASE(4)
            IF(sym .EQ. 'xxxx') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*rsq_2d*gvalue*sqrtpi2*(0.125d0)
            ELSEIF(sym .EQ. 'yyyy') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*rsq_2d*gvalue*sqrtpi2*(0.125d0)
            ELSEIF(sym .EQ. 'xxxy' .AND. m .EQ. 4) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*rsq_2d*gvalue*ci*sqrtpi2*(0.125d0)          
            ELSEIF(sym .EQ. 'xxxy' .AND. m .EQ. -4) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = -norm*rsq_2d*rsq_2d*gvalue*ci*sqrtpi2*(0.125d0)
            ELSEIF(sym .EQ. 'xxyy') THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = -norm*rsq_2d*rsq_2d*gvalue*sqrtpi2*(0.125d0)
            ELSEIF(sym .EQ. 'xyyy' .AND. m .EQ. 4) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = -norm*rsq_2d*rsq_2d*gvalue*ci*sqrtpi2*(0.125d0)          
            ELSEIF(sym .EQ. 'xyyy' .AND. m .EQ. -4) THEN
               norm = 16.0d0*(2.0d0/pi)**(0.75)*expo**(2.75)
               evalue = norm*rsq_2d*rsq_2d*gvalue*ci*sqrtpi2*(0.125d0)

            ELSE
                evalue = (0.0d0,0.0d0)
            ENDIF

         CASE DEFAULT
             evalue = (0.0d0,0.0d0)
      END SELECT


      ! Normalization factors for Q.Chem Insertion
      IF(qchemswitch == 1) THEN
         IF(sym .EQ. 'xx' .OR. sym .EQ. 'yy' .OR. sym .EQ. 'zz') evalue = evalue/SQRT(3.0d0)

         IF(sym .EQ. 'xxx' .OR. sym .EQ. 'yyy' .OR. sym .EQ. 'zzz') evalue = evalue/SQRT(15.0d0)
         IF(sym .EQ. 'xxy' .OR. sym .EQ. 'xxz' .OR. sym .EQ. 'xyy') evalue = evalue/SQRT(3.0d0)
         IF(sym .EQ. 'yyz' .OR. sym .EQ. 'xzz' .OR. sym .EQ. 'yzz') evalue = evalue/SQRT(3.0d0)

         IF(sym .EQ. 'xxxx' .OR. sym .EQ. 'yyyy' .OR. sym .EQ. 'zzzz') evalue = evalue/SQRT(105.0d0)

         IF(sym .EQ. 'xxxy' .OR. sym .EQ. 'xxxz' .OR. sym .EQ. 'xyyy') evalue = evalue/SQRT(15.0d0)
         IF(sym .EQ. 'yyyz' .OR. sym .EQ. 'xzzz' .OR. sym .EQ. 'yzzz' ) evalue = evalue/SQRT(15.0d0)

         IF(sym .EQ. 'xxyy' .OR. sym .EQ. 'xxzz' .OR. sym .EQ. 'yyzz') evalue = evalue/3.0d0
         IF(sym .EQ. 'xxyz' .OR. sym .EQ. 'xyyz' .OR. sym .EQ. 'xyzz') evalue = evalue/SQRT(3.0d0)
      ENDIF




      valarray(i) = evalue
!      WRITE(1909,'(1x,I3,6ES16.7E2,1x,A5,5(ES15.6E3,1x))') m,xi,eta,zoff,a,sym,expo,evalue
   ENDDO readeval

! If there was a read error, let me know!
ELSE fileopen
   WRITE(*,*) " Error opening expos.dat for xi= ", REAL(xi), AIMAG(xi), " and eta= ", REAL(eta)
ENDIF fileopen
CLOSE(10)
ENDSUBROUTINE
