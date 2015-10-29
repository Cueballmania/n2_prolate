! This function calculates the value of the partitioning function at a given x between [a,b]
FUNCTION beckepart(x, a, b)
IMPLICIT NONE

! Declare double precision
INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=13,r=200)
REAL(KIND=DBL) :: beckepart
REAL(KIND=DBL) :: becke

! Input variables: value and the boundaries a<b
REAL(KIND=DBL), INTENT(IN) :: x 
REAL(KIND=DBL), INTENT(IN) :: a, b

REAL(KIND=DBL) :: y, f

y = 2.0d0*(x-a)/(b-a)-1d0
f = becke(becke(becke(y)))
beckepart = 0.5d0*f+0.5d0

RETURN
END FUNCTION beckepart

FUNCTION becke(x)
IMPLICIT NONE

! Declare double precision
INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(p=13,r=200)

REAL(KIND=DBL) :: becke

! Input value
REAL(KIND=DBL), INTENT(IN) :: x

becke = 0.5*(3.0d0-x*x)*x

RETURN
END FUNCTION becke
