SUBROUTINE hamiltonian(xi_pts,ke_xi,eta_pts,ke_eta,eigbig,eta_wts,xi_wts)
!
!               Build 2D Hamiltonian
!
! The full Kinetic Energy matrix is passed in both variables.  Due to boundary conditions,
! one wants the DVR evaluation at the last endpoint in xi to be zero.  Hence, the grid and information
! will be left off and omitted from the full hamiltonian.
!
USE inputvars
IMPLICIT NONE

! Input variables
COMPLEX(KIND=DBL), INTENT(IN) :: xi_pts(1:n_xi), xi_wts(1:n_xi)     ! Arrays in Xi
COMPLEX(KIND=DBL), INTENT(IN) :: eta_pts(1:n_eta), eta_wts(1:n_eta) ! Arrays in Eta
COMPLEX(KIND=DBL), INTENT(IN) :: ke_xi(1:n_xi,1:n_xi)               ! KE from xi
COMPLEX(KIND=DBL), INTENT(IN) :: ke_eta(1:n_eta,1:n_eta)            ! KE from eta

! Output eigenvalue array
COMPLEX(KIND=DBL), INTENT(OUT) :: eigbig(1:(n_xi-1)*n_eta) ! Array of eigenvalues from the Hamiltonian
COMPLEX(KIND=DBL) :: eigvals(1:12)

! Functions
INTEGER :: index_2d                                   ! Quick index function for the hybrid DVR functions
COMPLEX(KIND=DBL) :: v_potential                      ! Array for DVR potentials

! Local Arrays
COMPLEX(KIND=DBL), ALLOCATABLE :: ham2d(:,:)                   ! 2-D hamiltonian
REAL(KIND=DBL), ALLOCATABLE :: rham2d(:,:), cham2d(:,:)        ! real and imaginary 2-D hamiltonian
COMPLEX(KIND=DBL), ALLOCATABLE :: insertpot(:,:)               ! Insertion potential

! Local variables
INTEGER :: i, j, k, n, m, ni, mi, nj
INTEGER :: nbas_xi

! Local normalization variables
COMPLEX(KIND=DBL) :: smetric_mj, smetric_mi
COMPLEX(KIND=DBL) :: smetric_ni, smetric_nj

!!$INTEGER(KIND=DBL),ALLOCATABLE :: sparsem(:,:)

! Scale the number of xi components to remove the endpoint
nbas_xi = n_xi-1

! Zero the hamiltonian

ALLOCATE(ham2d(1:ntot,1:ntot))

ham2d = (0.d0,0.d0)

!   Generate the DVR Kinetic energy.
!   As each element of the Hamiltonian is built, divide out the
!   appropriate combination of square roots of overlap matrix elements
!   denoted by smetric_ni etc.
!
! load the xi kinetic energy blocks, diagonal in eta
!
DO i=1,n_eta
   DO n=1,nbas_xi
      ni = index_2d(n,i,nbas_xi)
      smetric_ni = SQRT(xi_pts(n)**2 - eta_pts(i)**2)
      DO m=1,nbas_xi
         mi = index_2d(m,i,nbas_xi)
         smetric_mi = SQRT(xi_pts(m)**2 - eta_pts(i)**2)
         ham2d(ni,mi) = ham2d(ni,mi) + ke_xi(n,m)/(smetric_ni*smetric_mi)/R_over_two**2/amasse
      ENDDO
   ENDDO
ENDDO

!
!  load the eta kinetic energy blocks, diagonal in xi
!
DO n=1,nbas_xi
   DO i=1,n_eta
      ni = index_2d(n,i,nbas_xi)
      smetric_ni = SQRT(xi_pts(n)**2 - eta_pts(i)**2)
      DO j=1,n_eta
         nj =  index_2d(n,j,nbas_xi)
         smetric_nj = SQRT(xi_pts(n)**2 - eta_pts(j)**2)
         ham2d(ni,nj) = ham2d(ni,nj) + ke_eta(i,j)/(smetric_ni*smetric_nj)/R_over_two**2/amasse
         ENDDO
   ENDDO
ENDDO


! Check the switches to see what the Hamiltonian is made up from
IF(vswitch == 1) THEN
! DVR potential calculation
DO n=1,nbas_xi
   DO i=1,n_eta
      ni = index_2d(n,i,nbas_xi)
      ham2d(ni,ni) = ham2d(ni,ni) + V_potential(xi_pts(n),eta_pts(i),R_over_two)
   ENDDO
ENDDO

ELSE
   ALLOCATE(insertpot(1:ntot,1:ntot))
   CALL svd_insert(xi_pts,xi_wts,nbas_xi,eta_pts,eta_wts,n_eta,insertpot)

   IF(tswitch == 1) THEN
      DO n=1, nbas_xi*n_eta
         DO i=1, nbas_xi*n_eta
            ham2d(n,i) = ham2d(n,i) + insertpot(n,i)
         ENDDO
      ENDDO

   ELSE
      DO n=1, nbas_xi*n_eta
         DO i=1, nbas_xi*n_eta
            ham2d(n,i) = insertpot(n,i)
         ENDDO
      ENDDO
   ENDIF

   DEALLOCATE(insertpot)
ENDIF

WRITE(6,'(" finished building Hamiltonian")')


   OPEN(UNIT=10,FILE='rhamil.out',STATUS='UNKNOWN',ACTION='WRITE')
   DO i=1,ntot
       WRITE(10,'(1x,20000ES16.8)') (REAL(ham2d(i,j)),j=1, ntot)
  ENDDO
  CLOSE(10)
!!$
!!$
!!$   OPEN(UNIT=11,FILE='ihamil.out',STATUS='UNKNOWN',ACTION='WRITE')
!!$   DO i=1,ntot
!!$       WRITE(11,'(1x,20000ES16.8)') (AIMAG(ham2d(i,j)),j=1, ntot)
!!$   ENDDO
!!$   CLOSE(11)

!   Diagonalize ham2d, finding all eigenvalues, no eigenvectors
!   usiing LAPACK routine

!!$ALLOCATE(sparsem(1:ntot,1:ntot))
!!$DO i=1, ntot
!!$   DO j=1, ntot
!!$      IF(ABS(ham2d(i,j)) < 1E-6) THEN
!!$         sparsem(i,j) = 0
!!$      ELSE
!!$         sparsem(i,j) = 1
!!$      ENDIF
!!$      ENDDO
!!$ENDDO
!!$
!!$   OPEN(UNIT=10,FILE='lookhamil.out',STATUS='UNKNOWN',ACTION='WRITE')
!!$   DO i=1,ntot
!!$       WRITE(10,'(1x,200000I2)') (sparsem(i,j),j=1, ntot)
!!$   ENDDO
!!$   CLOSE(10)
!!$
!!$DEALLOCATE(sparsem)
!!$STOP

WRITE(6,'("begining diagonalization of hamiltonian ",i5," x",i5)') ntot,ntot


!!$ALLOCATE(rham2d(1:ntot,1:ntot), cham2d(1:ntot,1:ntot))
!!$
!!$rham2d = REAL(ham2d)
!!$cham2d = AIMAG(ham2d)
!!$
!!$OPEN(UNIT=32, FILE='rham.bin', FORM='unformatted', STATUS='UNKNOWN', ACTION='WRITE')
!!$WRITE(32) rham2d
!!$CLOSE(32)
!!$
!!$OPEN(UNIT=33, FILE='cham.bin', FORM='unformatted', STATUS='UNKNOWN', ACTION='WRITE')
!!$WRITE(33) cham2d
!!$CLOSE(33)
!!$
!!$
!!$DEALLOCATE(rham2d, cham2d)

WRITE(*,*) "The shift for the diagonalization is ", sigma

CALL arnoldi(ntot, ham2d, sigma, 12, eigvals)

!CALL diagwrap(ntot, ham2d, eigbig)

OPEN(UNIT=3737, FILE='spectrum.out', STATUS='UNKNOWN', ACTION='WRITE')
DO i=1, ntot
   WRITE(3737,'(1x,2ES19.9E4)') eigbig(i)
ENDDO
CLOSE(3737)

DEALLOCATE(ham2d)

END SUBROUTINE
