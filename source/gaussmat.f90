!This subroutine evalues the integral, in the DVR approximation, of a DVR function
!multiplied by a primative gaussian function in prolate coordinates.
SUBROUTINE gaussmat(xi_pts,xi_wts,nbas_xi,eta_pts,eta_wts,nbas_eta,gaussdvrmat,mfac)
USE inputvars
IMPLICIT NONE

! Input variables
COMPLEX(KIND=DBL), INTENT(IN) :: xi_pts(1:n_xi)                ! DVR pts in xi
COMPLEX(KIND=DBL), INTENT(IN) :: xi_wts(1:n_xi)                ! weights in xi
COMPLEX(KIND=DBL), INTENT(IN) :: eta_pts(1:n_eta)              ! DVR pts in eta
COMPLEX(KIND=DBL), INTENT(IN) :: eta_wts(1:n_eta)              ! weights in eta
INTEGER, INTENT(IN) :: nbas_xi, nbas_eta                       ! number of basis functions in xi eta
INTEGER, INTENT(IN) :: mfac                                    ! m-factor for conjugation

! Output variable
COMPLEX(KIND=DBL), INTENT(OUT) :: gaussdvrmat(1:ntot,1:numprimg)  ! DVR matrix evaluation

! Local variables
INTEGER :: i, j, k                        ! loop variables
INTEGER :: ierror                         ! error variable
INTEGER :: ni, index_2d                   ! indexing variables
INTEGER :: meff                           ! Effective m
COMPLEX(KIND=DBL) :: var_jacob            ! Jacobian factor in prolate
COMPLEX(KIND=DBL) :: weights              ! Normalizing weight factor for DVR
COMPLEX(KIND=DBL) :: gaussarray(1:numprimg) ! Array holding the values of the gaussian at the DVR points
COMPLEX(KIND=DBL) :: local_xi             ! Calling variable
COMPLEX(KIND=DBL) :: local_eta            ! Calling variable

! Zero the array
gaussdvrmat = (0.0d0,0.0d0)

! Print out xi grid
OPEN(UNIT=3000, FILE='xigrid.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
IF (ierror == 0) THEN
   DO i=1, n_xi
      WRITE(3000,'(1x,80ES17.9)') xi_pts(i), xi_wts(i)
   ENDDO
ELSE
   WRITE(*,*) " Could not write xi grid to file"
ENDIF
CLOSE(3000)

! Print out eta grid
OPEN(UNIT=4000, FILE='etagrid.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
IF (ierror == 0) THEN
   DO i=1, n_eta
      WRITE(4000,'(1x,80ES17.9)') eta_pts(i), eta_wts(i)
   ENDDO
ELSE
   WRITE(*,*) " Could not write eta grid to file"
ENDIF
CLOSE(4000)

! Determine complex conjugation
meff=mfac*m_ang

! Main loop.  Loops over xi and eta, determines the index of the combination by using the
! user defined function index_2d.
! Then calculates the normalization factor and calls stos to get the primative gaussians
! evaluated at the particular xi and eta points.
! Finally, the matrix is made.
DO i=1,nbas_xi
   DO j=1,nbas_eta
      ni = index_2d(i,j,nbas_xi)
      local_xi = xi_pts(i)
      local_eta = eta_pts(j)
      weights = SQRT(xi_wts(i)*eta_wts(j))
      var_jacob = SQRT(R_over_two**3*(local_xi**2-local_eta**2))
      CALL stos(local_xi, local_eta, R_over_two, numprimg, gaussarray, meff)
      DO k=1, numprimg
         gaussdvrmat(ni,k) = weights*var_jacob*gaussarray(k)
      ENDDO
   ENDDO
ENDDO

! WRITE out the Gaussian evaluated DVR matrix
OPEN(UNIT=1000, FILE='gaussdvrmat.out', STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ierror)
IF (ierror == 0) THEN
   j=0
   DO i=1, nbas_xi
      ni = index_2d(i,j,nbas_xi)
      WRITE(1000,'(1x,296ES18.7E4)') xi_pts(i),(gaussdvrmat(ni,k),k=1,numprimg)
   ENDDO
ELSE
   WRITE(*,*) " Could not write gaussmat to file"
ENDIF

CLOSE(1000)

END SUBROUTINE
