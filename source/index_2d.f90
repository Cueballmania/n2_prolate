INTEGER FUNCTION index_2d(m_xi,i_eta,nbas_xi)
IMPLICIT NONE
!
!  2D index for Hamiltonian and other arrays 
!      m_xi varies fastest
!      dimension of fastest varying array is nbas_xi
!      i_eta varies slowest
! 
INTEGER :: m_xi, i_eta, nbas_xi
!
index_2d = m_xi + (i_eta-1)*nbas_xi
END FUNCTION
