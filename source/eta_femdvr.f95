SUBROUTINE eta_femdvr(norder, m_ang, neta_divide, eta_pts, eta_wts, ke_eta)
! This subroutine calculates the kinetic energy matrix in eta using FEM DVR
! the subroutine assumes three elements [-1:-0.9,-0.9:0.9,0.9:1]
! Note, there is no complex scaling.
!
! Gauss-Lobotto quadrature is sandwiched between two Gauss-Radau quadratures
!
USE gaussquad
IMPLICIT NONE

! Parameters
INTEGER, PARAMETER :: DBL =  SELECTED_REAL_KIND(13,200)         ! Define Double Precision (precision,range)
REAL(KIND=DBL), PARAMETER :: pi = 3.14159265359

! Input variables
INTEGER, INTENT(IN) :: norder               ! order of the quadrature
INTEGER, INTENT(IN) :: m_ang                ! m angular quantum number
REAL(KIND=DBL), INTENT(IN) :: neta_divide   ! The size of the two end quadratures

! Output variables
! eta points and weights
COMPLEX(KIND=DBL), DIMENSION(1:3*(norder-1)+1), INTENT(OUT) :: eta_pts, eta_wts

! KE matrix in xi
COMPLEX(KIND=DBL), DIMENSION(1:3*(norder-1)+1,1:3*(norder-1)+1), INTENT(OUT) :: ke_eta

! Quadrature arrays
! Radau arrays
INTEGER :: ilegendre = 0, info
CHARACTER(LEN=1) :: radaucharr = 'R', radaucharl = 'L'
REAL(KIND=DBL), DIMENSION(1:norder) :: rwork, lwork

! Lobotto arrays
CHARACTER(LEN=1) :: lobachar = 'B'
REAL(KIND=DBL), DIMENSION(1:norder) :: mwork

! Points and weights arrays
REAL(KIND=DBL), DIMENSION(1:norder, 1:3) :: quadpts, quadweights
REAL(KIND=DBL), DIMENSION(1:norder, 1:norder, 1:2) :: deriv
COMPLEX(KIND=DBL), DIMENSION(1:norder*3) :: ptstoobig, wtstoobig
COMPLEX(KIND=DBL), DIMENSION(1:3*(norder-1)+1,1:norder*3) :: eta_compder

! Read variables to find element boundaries
INTEGER :: ierror
REAL(KIND=DBL), DIMENSION(1:3) :: boundsizes = (/ 0.1,1.8,0.1 /)

! Local variables
INTEGER :: i, j, k
INTEGER :: ii, jj, kk
REAL(KIND=DBL) :: integ

boundsizes(1) = neta_divide
boundsizes(2) = 2.0d0-2.0d0*neta_divide
boundsizes(3) = neta_divide

! Write out the FEM boundaries
WRITE(6,'(//," FEM DVR boundary sizes for eta are: ", 10ES15.3)') (boundsizes(i),i=1,3)

! Generate Radau points
CALL gauss_rule(ilegendre, norder, quadpts(:,1), quadweights(:,1), rwork, 0.0D0, 0.0D0, radaucharr, info)

IF(info .NE. 0) THEN
     WRITE(6,'(//,"gauss_rule error in Radau step",I3)') info
     STOP
ENDIF

! Generate Lobotto points if needed
CALL gauss_rule(ilegendre, norder, quadpts(:,2), quadweights(:,2), mwork, 0.0D0, 0.0D0, lobachar, info)

IF(info .NE. 0) THEN
   WRITE(6,'(//,"gauss_rule error in Lobatto step",I3)') info
   STOP
ENDIF

! Generate Lobotto points if needed
CALL gauss_rule(ilegendre, norder, quadpts(:,3), quadweights(:,3), lwork, 0.0D0, 0.0D0, radaucharl, info)

IF(info .NE. 0) THEN
   WRITE(6,'(//,"gauss_rule error in Radua step",I3)') info
   STOP
ENDIF

! Generate the derivative matrices over the interpolating polynomials on [-1,1]
! for the ith function with the jth argument
! Uses real numbers
! ii=1 : Radau
! ii=2 : Lobatto
! ii=3 : Radau
! Loop over quadratures
DO ii=1, 3
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
eta_wts = 0.0
integ = -1.0
j=0
jj=0
DO i=1, 3
   DO k=1, norder
      jj=jj+1
      
      ! Don't count the first point in each element which overlaps with the previous element's last point
      IF((i==1) .OR. (k/=1)) j=j+1

      eta_pts(j) = boundsizes(i)*0.5d0*(quadpts(k,i)+1.0d0)+integ
      eta_wts(j) = eta_wts(j)+quadweights(k,i)*boundsizes(i)*0.5d0

      ptstoobig(jj) = eta_pts(j)
      wtstoobig(jj) = quadweights(k,i)*boundsizes(i)*0.5d0
   ENDDO

   integ = integ + boundsizes(i)
ENDDO


! Write to etagrid.out
OPEN(UNIT=20, FILE='etagrid.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
DO i=1,3*(norder-1)+1
      WRITE(20,'(1X,2ES17.8)') REAL(eta_pts(i)), REAL(eta_wts(i))
ENDDO
CLOSE(20)

! Create the derivative matrix in xi at all points
eta_compder = 0.0d0

DO i=1, 3
   DO j=1, norder
      jj=(i-1)*(norder-1)+j
      DO k=1, norder
         kk = (i-1)*norder+k
         eta_compder(jj,kk) = eta_compder(jj,kk) + deriv(j,k,i)*2.0d0/boundsizes(i)
      ENDDO
   ENDDO
ENDDO


! add in the appropriate factors if m is odd
IF(MOD(ABS(m_ang),2) == 1) THEN
   DO i=1, 3*norder
      eta_compder(:,i) = eta_compder(:,i) * SQRT(1.0d0-ptstoobig(i)**2)
   ENDDO

   ! Evaulate the derivate of the ith function at the jth point
   ! Only extra factors comes from the derivative of the ith function at a point where it's defined
   j=0
   DO i=1, 3*(norder-1)+1
      j=j+1
      eta_compder(i,j) = eta_compder(i,j) + ptstoobig(j)/SQRT(1.0d0-ptstoobig(j)**2)

      IF((MOD(i-1,norder-1) == 0) .AND. (i>1) .AND. (i < 3*(norder-1)+1)) THEN
         j=j+1
         eta_compder(i,j) = eta_compder(i,j) + ptstoobig(j)/SQRT(1.0d0-ptstoobig(j)**2)
      ENDIF
   ENDDO
   
   ! Normalize
   DO i=1, 3*(norder-1)+1
      eta_compder(i,:) = eta_compder(i,:)/SQRT(eta_wts(i))/SQRT(1.0d0-eta_pts(i)**2)
   ENDDO

! Else, normalize
ELSE
   DO i=1, 3*(norder-1)+1
      eta_compder(i,:) = eta_compder(i,:)/SQRT(eta_wts(i))
   ENDDO
ENDIF

! Create the KE matrix
DO i=1, 3*(norder-1)+1
   DO j=1, 3*(norder-1)+1
      integ = 0.0d0
      DO k=1, 3*norder
         integ = integ - REAL(eta_compder(i,k)*eta_compder(j,k)*wtstoobig(k)*(1.0d0-ptstoobig(k)**2))
      ENDDO
      ke_eta(i,j) = integ

      IF(j==i) ke_eta(i,i) = ke_eta(i,i) - FLOAT(m_ang)**2/(1.0d0-eta_pts(i)**2)
   ENDDO
ENDDO

ke_eta = -0.5d0*ke_eta

END SUBROUTINE eta_femdvr
