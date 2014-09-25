!Module to contain global variables
MODULE inputvars
IMPLICIT NONE
SAVE

!Global constants
INTEGER, PARAMETER :: DBL =  SELECTED_REAL_KIND(13,200)         ! Define Double Precision (precision,range)
REAL(KIND=DBL), PARAMETER :: amasse=1.0d0                       ! Mass of particle
REAL(KIND=DBL), PARAMETER :: pi=3.141592653589793238d0

! Input variables
! Grid
INTEGER :: nelements                                ! Number of elements in xi
INTEGER :: norder                                   ! Number of quadrature points per element
INTEGER :: ncscaled                                 ! Number of complex scaled elements
INTEGER :: n_xi                                     ! Number of pts in xi
INTEGER :: n_eta                                    ! Number of pts in eta
INTEGER :: ntot                                     ! Number of pair pts (n_xi-1)*n_eta
REAL(KIND=DBL) :: R_over_two                        ! +/- position of the nucleus 

! Switches
INTEGER :: tswitch                                  ! Type of calculation
                                                    !    tswitch = 1: DVR calc
                                                    !    tswitch = 2: Insertion

INTEGER :: vswitch                                  !    vswitch = 1: DVR calc
                                                    !    vswitch = 2: Insertion

INTEGER :: svdswitch                                ! SVDswitch
                                                    !    0: Insertion with pesudoinverse
                                                    !    1: SVD orthogonalization of orbitals

! Guassian
INTEGER :: numgauss                                 ! Number of contracted Gaussian basis
INTEGER :: numprimg                                 ! Number of primative Gaussians

! Angular
INTEGER :: j_ang                                    ! J angular momentum ***NOT REFERENCED***
INTEGER :: m_ang                                    ! M angular momentum

! Scaling
REAL(KIND=DBL) :: scaling_angle                     ! Scaling angle in degrees

! SVD
REAL(KIND=DBL) :: SVD_tol                           ! Tolerance of the SVD

END MODULE


PROGRAM prolate_dvr
!
!  This is the main program for one-electron calculations for any (cylindrially symmetric) potential 
!  in prolate spheroidal coordinate.
!
!  This program also allows for the insertion of a non-orthongonal basis to evaluate the
!  potential (and possibliy the kinetic energy) via an N^2 seperable potential.
!
!  It is meant to be a code for potential problems with a fixed distance between the foci 
!  (at R_over_two) of the prolate coordinate system.
!
!  The coordinates are 1 < xi < xi_max and -1 < eta < +1 where straight or exterior complex scaling
!  can be applied to xi.
!
!  This was modified from a working "improved adiabatic fixed nuclei" code in Jan 2013 by CWM
!  
!  This has been UPDATED to FORTRAN90/95 for array handeling and the insertion basis by BAA
!  This is also an update to include an FEM-DVR in xi
!
USE inputvars
IMPLICIT NONE

INTEGER :: ierror                                   ! Error flag
REAL(KIND=DBL) :: xi_max_real                       ! Max value of xi (unscaled)
COMPLEX(KIND=DBL) :: xi_min=(1.0d0,0.0d0)           ! Minimum of xi coordinates
COMPLEX(KIND=DBL) :: xi_max                         ! Complex maximum value of xi 
INTEGER :: i,j                                      ! Loop indices
REAL(KIND=DBL) :: rval                              ! R, the distance between focii

! arrays for 1D operators
INTEGER :: nbas_eta,nbas_xi                         ! Reduced size of the basis
COMPLEX(KIND=DBL) :: cmplx_scaling                  ! Complex scaling factor
COMPLEX(KIND=DBL), ALLOCATABLE :: ke_xi(:,:)        ! Kinetic energy matrix in xi
COMPLEX(KIND=DBL), ALLOCATABLE :: ke_eta(:,:)       ! Kinetic energy matrix in eta
COMPLEX(KIND=DBL), ALLOCATABLE :: der_xi(:,:)       ! Derivative matrix in xi
COMPLEX(KIND=DBL), ALLOCATABLE :: der_eta(:,:)      ! Derivative matrix in eta
COMPLEX(KIND=DBL), ALLOCATABLE :: xi_pts(:)         ! DVR pts in xi
COMPLEX(KIND=DBL), ALLOCATABLE :: eta_pts(:)        ! DVR pts in eta
COMPLEX(KIND=DBL), ALLOCATABLE :: xi_wts(:)         ! Weights in xi
COMPLEX(KIND=DBL), ALLOCATABLE :: eta_wts(:)        ! Weights in eta

! arrays for complex diagonalizations
COMPLEX(KIND=DBL), ALLOCATABLE :: eigbig(:)         ! Eigenvalues of the Hamiltonian

! Attempt to open the file indvr
OPEN(UNIT=9, FILE='indvr_complex', STATUS='OLD', ACTION='READ', IOSTAT=ierror)

openif: IF(ierror == 0) THEN

   ! Open was successful, read in inputs from file
   ! Angular values.  Note that j_ang is never referenced
   READ(9,*) j_ang, m_ang

   ! input information for grid in xi
   READ(9,*) nelements, norder

   ! input information for grid in eta
   READ(9,*) n_eta

   ! input information for r - the internuclear separation
   READ(9,*) rval

   ! input information for the complex scaling angle and number of complex elements
   READ(9,*) scaling_angle, ncscaled

   ! input information for the kinetic energy evaluation (DVR or insertion)
   READ(9,*) tswitch

   ! input information for the potential energy evaluation (DVR or insertion)
   READ(9,*) vswitch

   ! input information for how the insertion is evaluated and the SVD_tolerance
   !   svdswitch = 0 : Use pseudoinverse of the overlaps
   !   svdswitch = 1 : SVD orthogonalization
   READ(9,*) svdswitch, SVD_tol

   ! Read in number of Gaussians basis being inserted and number of primative gaussians
   READ(9,*) numgauss, numprimg
 
   CLOSE(9)

! If opening the file results in a error, print message and stop
ELSE openif
   WRITE(6,*) " Error opening input file: IOSTAT=", ierror
   STOP
ENDIF openif

! Print out input and preform some quick checks

WRITE(6,'(//,10x,"Prolate Spheroidal Coordinate DVR Code")')

! Check for even or odd 'm' and print out result
IF(m_ang == 2*(m_ang/2)) THEN
   WRITE(6,'(//," Even m case, using ordinary DVRs with m= ", I3)') m_ang
ELSE
   WRITE(6,'(//," Odd m case, using appropriate DVRs with m= ", I3)') m_ang
ENDIF

! Write out the type of scaling calculation
IF(nelements < ncscaled) THEN
   WRITE(6,'(//," Too many complex scaled elements (ncscaled > nelements).  Changing the calculation to straight scaling")')
   ncscaled = nelements
ELSE IF(nelements == ncscaled) THEN
   WRITE(6,'(//," Straight Complex Scaling calculation")')
ELSE IF(nelements > ncscaled .AND. ncscaled > 0) THEN
   WRITE(6,'(//," ECS calculation")')
ELSE
   WRITE(6,'(//," Real-valued coordinate calculation")')
ENDIF

! Write out the scaling angle
IF(ncscaled > 0) THEN
   WRITE(6,'(" Complex scaling by ",f10.5," degrees.",f10.5)') scaling_angle
   cmplx_scaling = EXP((0.d0,1.d0)*scaling_angle*pi/180.d0) 
   WRITE(6,'(" Complex_scaling factor = ",2e15.5) ') cmplx_scaling
ENDIF

! Print out the type of insertion calculation
IF(tswitch == 2 .AND. vswitch == 1) THEN
   WRITE(6,'(//," Input specified insertion for the T-matrix and DVR for the V-matrix.  Exiting")')
   STOP
ELSE IF(tswitch == 1 .AND. vswitch == 1) THEN
   WRITE(6,'(//," DVR calculation for both the T-matrix and the V-matrix.")')
ELSE IF(tswitch == 1 .AND. vswitch == 2) THEN
   WRITE(6,'(//," DVR calculation for the T-matrix and insertion for the V-matrix.")')
ELSE IF(tswitch == 2 .AND. vswitch == 2) THEN
   WRITE(6,'(//," Insertion calcuation for both the T-matrix and the V-matrix.")')
ENDIF

! Write out the SVD information if there is insertion
IF(vswitch == 2) THEN
   IF (svdswitch == 0) THEN
      WRITE(6, 400) numgauss, numprimg, svd_tol
400   FORMAT (1x, 'Insertion: the pseudoinverse of the overlap matrix for ', I3, ' contracted Gaussians and ',&
           &I3, ' primative Gaussians with a SVD_tolerance of', ES15.7)
   ELSEIF (svdswitch == 1) THEN
      WRITE(6, 401) numgauss, numprimg, svd_tol
401   FORMAT (1x, 'Insertion: Orthogonalizing ', I3, ' contracted Gaussians and ',&
           &I3, ' primative Gaussians using an SVD tolerance of ', ES15.7)
   ELSE
      WRITE(6,*) " an invalid SVD_switch = ", svdswitch
      STOP
   ENDIF
ENDIF

! Echo grid information
WRITE(6,'(//," R = ",2f12.6)') rval
WRITE(6,'(" Prolate spheroidal grid")')
WRITE(6,'(2x,i4," points in eta on -1, 1")') n_eta

n_xi = nelements*(norder-1)+1
WRITE(6,'(2x,i4," points in xi.")') n_xi
WRITE(6,'(//)')

! Now build the xi and eta kinetic energies without the (R/2)^(-2) overall factors
! We want to remove the last xi point to enforce the function to be zero at the boundary 
ntot = (n_xi-1)*n_eta 
R_over_two = rval/2.d0

! Build Xi grid with ECS and compute the Kinetic Energy
ALLOCATE(xi_pts(1:n_xi), xi_wts(1:n_xi), ke_xi(1:n_xi,1:n_xi))

CALL xi_dvr(norder, nelements, ncscaled, scaling_angle, m_ang, xi_pts, xi_wts, ke_xi)

WRITE(6,'(" Built grid and KE in xi")')

OPEN(UNIT=20, FILE='xi_1d.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
DO i=1, n_xi
   WRITE(20,'(1x,4ES16.6)') xi_pts(i), xi_wts(i)
ENDDO
CLOSE(20)

! Build eta grid and compute the ECS
ALLOCATE(eta_pts(1:n_eta), eta_wts(1:n_eta), ke_eta(1:n_eta,1:n_eta))

CALL eta_dvr(n_eta, m_ang, eta_pts, eta_wts, ke_eta)

WRITE(6,'(" Built grid and KE in eta")')

OPEN(UNIT=35, FILE='eta_1d.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
WRITE(*,*) ierror
DO i=1, n_eta
   WRITE(35,'(1x,4ES16.6)') eta_pts(i), eta_wts(i)
ENDDO
CLOSE(35)

! CAll the hamiltonian matrix
ALLOCATE(eigbig(1:(n_xi-1)*n_eta))

CALL hamiltonian(xi_pts,ke_xi,eta_pts,ke_eta,eigbig,eta_wts,xi_wts)

DO i=1, (n_xi-1)*n_eta
!   IF(REAL(eigbig(i)) < 0.24d0 .AND. REAL(eigbig(i)) > 0.05d0) WRITE(9101,'(2ES17.6)') REAL(eigbig(i)),AIMAG(eigbig(i))
   If(REAL(eigbig(i)) < -1.0D-10) WRITE(*,'(2ES17.6)') REAL(eigbig(i)), AIMAG(eigbig(i))
ENDDO

DEALLOCATE(xi_pts, xi_wts, ke_xi)
DEALLOCATE(eta_pts, eta_wts, ke_eta)
DEALLOCATE(eigbig)

END PROGRAM prolate_dvr
