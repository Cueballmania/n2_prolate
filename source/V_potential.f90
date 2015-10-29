SUBROUTINE V_potential(xi,eta,R_over_two,vvalue)
!
!   potential for one electron problem specified completely here
!   needs to be cylindrically symmetric (otherwise phi dependence would couple 
!   m values
!
IMPLICIT NONE
INTEGER, PARAMETER :: DBL =  SELECTED_REAL_KIND(13,200)         ! Define Double Precision (precision,range)
COMPLEX(KIND=DBL), INTENT(IN) :: xi, eta
REAL(KIND=DBL), INTENT(IN) :: R_over_two
COMPLEX(KIND=DBL), INTENT(OUT) ::  vvalue

COMPLEX(KIND=DBL) :: r, z,xsquared_plus_ysquared

!  coordinates not dependent on phi
r = R_over_two*SQRT((xi**2-1.d0)*(1.d0-eta**2)+(xi**2)*(eta**2))
xsquared_plus_ysquared=(R_over_two**2)*(xi**2-1.d0)*(1.d0-eta**2)
z = R_over_two*xi*eta
!
! specify potential formula 
!
! potential from 1978 paper on complex basis function method
!  parameters of gaussian model potential from 1978 paper
!  note addition of the term proportional to sigma allows it to have bound states
!  REAL(KIND=DBL) :: lambda, gamma, A
!       lambda = 0.35d0
!       sigma = 1.25d0
!       gamma = 0.13d0
!       A = 2.d0
!       V_potential = lambda*(r**2)*(                                           
!     #   zexp(-gamma*(xsquared_plus_ysquared +(z-A)**2)) +
!     #   zexp(-gamma*(xsquared_plus_ysquared +(z+A)**2)) )
!       V_potential = (lambda*r**2 - sigma)*(
!     #   zexp(-gamma*(xsquared_plus_ysquared +(z-A)**2)) +
!     #   zexp(-gamma*(xsquared_plus_ysquared +(z+A)**2)) )
!
! old spherical model resonance problem
!V_potential =7.5d0*r**2*exp(-r)
!
! spherical harmonic oscillator
!        V_potential =.5d0*r**2
!
! Potential for n2
         vvalue = -2.0d0*7.0d0*xi/(xi**2-eta**2)/R_over_two 
!
! H atom: note you need to put the nucleus at one of the foci, pretty bad l=0 energies otherwise
!         V_potential = -1.d0/(xi+eta)/R_over_two
!V_potential = -1.0d0/r
!
! H2+ molecule  (nuclear repulsion is added in)
!         V_potential = -2.d0*xi/(xi**2-eta**2)/R_over_two + 1.d0/(2.d0*R_over_two)
END SUBROUTINE
