SUBROUTINE eta_dvr(norder, m_ang, eta_pts, eta_wts, ke_eta)
! This subroutine calculates the kinetic energy matrix in eta with proper normalization.
! We employ Gauss-Legendre quadrature in the interval (-1:1) since the volume element in 
! Prolate Spheroidal coordinates kills off the factors on the end points.
USE gaussquad
IMPLICIT NONE

! Parameter
INTEGER, PARAMETER :: DBL =  SELECTED_REAL_KIND(13,200)         ! Define Double Precision (precision,range)

! Input variables
INTEGER, INTENT(IN) :: norder               ! order of the quadrature
INTEGER, INTENT(IN) :: m_ang                ! m angular quantum number

! Output variables
! eta points and weights
COMPLEX(KIND=DBL), DIMENSION(1:norder), INTENT(OUT) :: eta_pts, eta_wts

! KE matrix in eta
COMPLEX(KIND=DBL), DIMENSION(1:norder, 1:norder), INTENT(OUT):: ke_eta

! Quadrature arrays
! Lobotto arrays
INTEGER :: ilegendre = 0, info
CHARACTER(LEN=1) :: legchar = 'N'
REAL(KIND=DBL), DIMENSION(1:norder) :: legpts, legweights, lwork
REAL(KIND=DBL), DIMENSION(1:norder, 1:norder) :: deriv

! Local variables
INTEGER :: i, j, k
COMPLEX(KIND=DBL) :: integ

! Generate Legendre points
CALL gauss_rule(ilegendre, norder, legpts, legweights, lwork, 0.0D0, 0.0D0, legchar, info)

IF(info .NE. 0) THEN
     WRITE(6,'(//,"gauss_rule error in Legendre step",I3)') info
     STOP
ENDIF

eta_pts = legpts
eta_wts = legweights

! Generate non-normalized derivative matrix
DO i=1, norder
   DO j=1, norder
      IF(j/=i) THEN
         deriv(i,j) = 1.0d0/(legpts(i)-legpts(j))
         DO k=1, norder
            IF((k/=i) .AND. (k/=j)) THEN
               deriv(i,j) = deriv(i,j) * (legpts(j)-legpts(k))/(legpts(i)-legpts(k))
            ENDIF
         ENDDO

      ELSE
         deriv(i,i) = 0.0d0
         DO k=1, norder
            IF(k/=i) THEN
               deriv(i,i) = deriv(i,i) + 1.0d0/(legpts(i)-legpts(k))
            ENDIF
         ENDDO
      ENDIF
   ENDDO
ENDDO

!If m is odd, multiply by the appropriate factor
IF(MOD(ABS(m_ang),2) == 1) THEN
   DO i=1, norder
      deriv(:,i) = deriv(:,i) * SQRT(1.0d0-legpts(i)**2)
   ENDDO
   
   DO i=1, norder
      deriv(i,i) = deriv(i,i) - legpts(i)/SQRT(1.0d0-legpts(i)**2)
   ENDDO

   DO i=1, norder
      deriv(i,:) = deriv(i,:) / SQRT(legweights(i))/SQRT(1.0d0-legpts(i)**2)
   ENDDO

ELSE
   DO i=1, norder
      deriv(i,:) = deriv(i,:) / SQRT(legweights(i))
   ENDDO
ENDIF


DO i=1, norder
   DO j=1, norder
      integ = (0.0d0,0.0d0)
      DO k=1, norder
         integ = integ - legweights(k)*(1.0d0-legpts(k)**2)*deriv(i,k)*deriv(j,k)
      ENDDO
      ke_eta(i,j) = integ

      IF(j==i) ke_eta(i,i) = ke_eta(i,i) - FLOAT(m_ang)**2/(1-legpts(i)**2)
   ENDDO
ENDDO

ke_eta = -0.5d0*ke_eta

END SUBROUTINE eta_dvr




