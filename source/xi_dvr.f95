SUBROUTINE xi_dvr(norder, nelements, ncscaled, scaling_angle, m_ang, xi_pts, xi_wts, ke_xi)
! This subroutine calculates the kinetic energy matrix in xi with proper normalization.
! We employ Gauss-Radau in the first element and Gauss-Lobatto in the subsequent elements
! Complex scaling is also preformed here
USE gaussquad
IMPLICIT NONE

! Parameters
INTEGER, PARAMETER :: DBL =  SELECTED_REAL_KIND(13,200)         ! Define Double Precision (precision,range)
REAL(KIND=DBL), PARAMETER :: pi = 3.14159265359

! Input variables
INTEGER, INTENT(IN) :: norder               ! order of the quadrature
INTEGER, INTENT(IN) :: nelements            ! number of finite elements
INTEGER, INTENT(IN) :: ncscaled             ! number of complex scaled elements
INTEGER, INTENT(IN) :: m_ang                ! m angular quantum number

REAL(KIND=DBL), INTENT(IN) :: scaling_angle ! Angle for complex scaling (deg)

! Output variables
! xi points and weights
COMPLEX(KIND=DBL), DIMENSION(1:nelements*(norder-1)+1), INTENT(OUT) :: xi_pts, xi_wts

! KE matrix in xi
COMPLEX(KIND=DBL), DIMENSION(1:nelements*(norder-1)+1,1:nelements*(norder-1)+1), INTENT(OUT) :: ke_xi

! Quadrature arrays
! Radau arrays
INTEGER :: ilegendre = 0, info
CHARACTER(LEN=1) :: radauchar = 'R'
REAL(KIND=DBL), DIMENSION(1:norder) :: rwork

! Lobotto arrays
CHARACTER(LEN=1) :: lobachar = 'B'
REAL(KIND=DBL), DIMENSION(1:norder) :: lwork

! Points and weights arrays
REAL(KIND=DBL), DIMENSION(1:norder, 1:2) :: quadpts, quadweights
REAL(KIND=DBL), DIMENSION(1:norder, 1:norder, 1:2) :: deriv
COMPLEX(KIND=DBL), DIMENSION(1:norder*nelements) :: ptstoobig, wtstoobig
COMPLEX(KIND=DBL), DIMENSION(1:nelements*(norder-1)+1,1:norder*nelements) :: xi_compder

! Read variables to find element boundaries
INTEGER :: ierror
REAL(KIND=DBL), DIMENSION(1:nelements) :: boundsizes

! Function to find the correct quadrature for each element
INTEGER :: el

! Local variables
INTEGER :: i, j, k
INTEGER :: ii, jj, kk
COMPLEX(KIND=DBL) :: integ, weightfac


! Read in the FEM boundaries
OPEN(UNIT=9, FILE='FEM.in', STATUS='OLD', ACTION='READ', IOSTAT=ierror)
openif: IF(ierror == 0) THEN
   DO i=1, nelements
      READ(9,*) boundsizes(i)
   ENDDO

   CLOSE(9)

ELSE openif 
   WRITE(6,*) "Couldn't open FEM.in"
   STOP
ENDIF openif

! Write out the FEM boundaries
WRITE(6,'(//," FEM DVR boundary sizes for xi are: ", 10ES15.3)') (boundsizes(i),i=1, nelements)

IF(ncscaled /= 0 .AND. ncscaled <= nelements) THEN
   WRITE(6,'(//," Complex scaling starts at boundary ", I3, " with the scaling at", 2ES15.5, " degrees.")') ncscaled, scaling_angle
ELSE
   WRITE(6,'(//," No Complex scaling")')
ENDIF

! Generate Radau points
CALL gauss_rule(ilegendre, norder, quadpts(:,1), quadweights(:,1), rwork, 0.0D0, 0.0D0, radauchar, info)

IF(info .NE. 0) THEN
     WRITE(6,'(//,"gauss_rule error in Radau step",I3)') info
     STOP
ENDIF

! Generate Lobotto points if needed
CALL gauss_rule(ilegendre, norder, quadpts(:,2), quadweights(:,2), lwork, 0.0D0, 0.0D0, lobachar, info)

IF(info .NE. 0) THEN
   WRITE(6,'(//,"gauss_rule error in Lobatto step",I3)') info
   STOP
ENDIF

! Generate the derivative matrices over the interpolating polynomials on [-1,1]
! for the ith function with the jth argument
! Uses real numbers
! ii=1 : Radau
! ii=2 : Lobatto
! Loop over quadratures
DO ii=1, 2
   ! Loop over functions
   DO i=1, norder
      ! Loop over points
      DO j=1, norder
         
         IF(j /= i) THEN
            deriv(i,j,ii) = 1.0d0/(quadpts(i,ii)-quadpts(j,ii))
            DO k=1, norder
               IF(k/=i .AND. k/=j) THEN
                  deriv(i,j,ii) = deriv(i,j,ii)*(quadpts(j,ii)-quadpts(k,ii))/(quadpts(i,ii)-quadpts(k,ii))
               ENDIF
            ENDDO

         ELSE
            deriv(i,i,ii) = 0.0d0
            DO k=1, norder
               IF(k/=i) THEN
                  deriv(i,i,ii) = deriv(i,i,ii) + 1.0d0/(quadpts(i,ii)-quadpts(k,ii))
               ENDIF
            ENDDO
         ENDIF
      ENDDO
   ENDDO
ENDDO

! Generate the grid of points and weights
! This is a Dan Haxton special
xi_wts = (0.0d0,0.0d0)
integ = (1.0d0,0.0d0)
j=0
jj=0
DO i=1, nelements
   weightfac = 1.0d0
   
   ! Are we on the ECS contour?
   IF(i >= ncscaled) THEN
      weightfac = EXP((0.0d0,1.0d0)*scaling_angle*pi/180.0d0)
   ENDIF

   DO k=1, norder
      jj=jj+1
      
      ! Don't count the first point in each element which overlaps with the previous element's last point
      IF((i==1) .OR. (k/=1)) j=j+1

      xi_pts(j) = weightfac*boundsizes(i)*0.5d0*(quadpts(k,el(i))+1.0d0)+integ
      xi_wts(j) = xi_wts(j)+quadweights(k,el(i))*boundsizes(i)*0.5d0*weightfac

      ptstoobig(jj) = xi_pts(j)
      wtstoobig(jj) = quadweights(k,el(i))*boundsizes(i)*0.5d0*weightfac
   ENDDO

   integ = integ + boundsizes(i)*weightfac
ENDDO


! Write to grid.out
OPEN(UNIT=20, FILE='grid.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
DO i=1,nelements*(norder-1)+1
      WRITE(20,'(1X,4ES17.8)') xi_pts(i), xi_wts(i)
ENDDO
CLOSE(20)

! Create the derivative matrix in xi at all points
xi_compder = (0.0d0,0.0d0)

DO i=1, nelements
   weightfac = 1.0d0
   
   ! Are we on the ECS contour?
   IF(i >= ncscaled) THEN
      weightfac = EXP((0.0d0,1.0d0)*scaling_angle*pi/180.0d0)
   ENDIF

   DO j=1, norder
      jj=(i-1)*(norder-1)+j
      DO k=1, norder
         kk = (i-1)*norder+k
         xi_compder(jj,kk) = xi_compder(jj,kk) + deriv(j,k,el(i))*2.0d0/weightfac/boundsizes(i)
      ENDDO
   ENDDO
ENDDO


! add in the appropriate factors if m is odd
IF(MOD(ABS(m_ang),2) == 1) THEN
   DO i=1, nelements*norder
      xi_compder(:,i) = xi_compder(:,i) * SQRT(ptstoobig(i)**2-1.0d0)
   ENDDO

   ! Evaulate the derivate of the ith function at the jth point
   ! Only extra factors comes from the derivative of the ith function at a point where it's defined
   j=0
   DO i=1, nelements*(norder-1)+1
      j=j+1
      xi_compder(i,j) = xi_compder(i,j) + ptstoobig(j)/SQRT(ptstoobig(j)**2-1.0d0)

      IF((MOD(i-1,norder-1) == 0) .AND. (i>1) .AND. (i < nelements*(norder-1)+1)) THEN
         j=j+1
         xi_compder(i,j) = xi_compder(i,j) + ptstoobig(j)/SQRT(ptstoobig(j)**2-1.0d0)
      ENDIF
   ENDDO
   
   ! Normalize
   DO i=1, nelements*(norder-1)+1
      xi_compder(i,:) = xi_compder(i,:)/SQRT(xi_wts(i))/SQRT(xi_pts(i)**2-1.0d0)
   ENDDO

! Else, normalize
ELSE
   DO i=1, nelements*(norder-1)+1
      xi_compder(i,:) = xi_compder(i,:)/SQRT(xi_wts(i))
   ENDDO
ENDIF

! Create the KE matrix
DO i=1, nelements*(norder-1)+1
   DO j=1, nelements*(norder-1)+1
      integ = (0.0d0,0.0d0)
      DO k=1, nelements*norder
         integ = integ - xi_compder(i,k)*xi_compder(j,k)*wtstoobig(k)*(ptstoobig(k)**2-1.0d0)
      ENDDO
      ke_xi(i,j) = integ

      IF(j==i) ke_xi(i,i) = ke_xi(i,i) - FLOAT(m_ang)**2/(xi_pts(i)**2-1.0d0)
   ENDDO
ENDDO

ke_xi = -0.5d0*ke_xi

END SUBROUTINE xi_dvr


FUNCTION el(index)
! This function returns which quadrature points to use based on which element one is on
IMPLICIT NONE
INTEGER :: el, index

IF(index == 1) THEN
   el=1
ELSE
   el=2
ENDIF
END FUNCTION el

