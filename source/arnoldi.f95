SUBROUTINE arnoldi(nsize, mat, sigma, numeigen, eigenvs)
!
! This subroutine uses the arnoldi iteration for double-complex valued matrix diagonalization
! This particular scheme uses the shift and invert method to find eigenvalues near a particular value, sigma
! On input, the user must input:
!           int nsize: size of the matrix
!           complex*16 mat: nsize x nsize matrix of complex elements
!           complex*16 sigma: shifted value for eigenvalues near this point
!           int numeigen: number of eigenvalues sought
! Returns
!           complex*16 eigenvs: numeigen eigenvalues
!
! Note that ARPACK finds the eigenvalues less than sigma.  So, please use sigma as the upper bound.
!
! Also, no eigenvectors are given.  Please set rvec to true and look in v for the eigenvectors
!
!
IMPLICIT NONE 

INTEGER, PARAMETER :: DBL =  SELECTED_REAL_KIND(13,200)         ! Define Double Precision (precision,range)

! Input variables
INTEGER, INTENT(IN) :: nsize
INTEGER, INTENT(IN) :: numeigen
COMPLEX(KIND=DBL), INTENT(INOUT) :: mat(1:nsize,1:nsize)
COMPLEX(KIND=DBL), INTENT(IN) :: sigma                       ! Shift value -- referenced if iparam(7) = 3

! Output eigenvalues
COMPLEX(KIND=DBL) :: eigenvs(1:numeigen+1)

! ZNAUPD variables
INTEGER :: ido                    ! Reverse communication flag.  Must be equal to zero before we start the Arnoldi procedure
CHARACTER(LEN=1) :: bmat = 'I'    ! Specify that this is not a generalized eigenvalue equation, but rather A*x_i=E_i*x_i
INTEGER :: n                      ! Dimension of the input matrix
CHARACTER(LEN=2) :: which         ! Specify which eigenvalue(s) we're after
INTEGER :: nev                    ! Number of target eigenvalues (should be less than n or else should be using ZGEEV)
REAL(KIND=DBL) :: TOL=0.0d0       ! Tolerance needed for convergence.  Set higher if you're just testing things out.  

COMPLEX(KIND=DBL), ALLOCATABLE :: resid(:)   ! Contains residual vector of size n
INTEGER :: ncv                               ! Krylov Space size (number of arnoldi basis vectors)  Should be larger than nev. Recommended > 2*nev+1
COMPLEX(KIND=DBL), ALLOCATABLE :: v(:,:)     ! Contains the final Arnoldi vectors.  Size is n x ncv
INTEGER :: ldv                               ! Should be the same as n
INTEGER :: iparam(11)                        ! Parameters --  1, 3, 7 are the important ones as defined below
INTEGER :: ipntr(14)                         ! Pointer array to workd.  (1) points to the vector to be operated on.  (2) points to the result.
COMPLEX(KIND=DBL), ALLOCATABLE :: workd(:)   ! Work array used be ZNAUPD.  Pointer 1 is the vector created.  Pointer 2 points to where the resulting
                                             !        vector is located.  Size should be 3*n
COMPLEX(KIND=DBL), ALLOCATABLE :: workl(:)   ! Work array used by ZNAUPD and only referenced there.  Size should be lwork = 3* ncv**2 + 5* ncv
INTEGER :: lwork                             ! length of work array workl.  Size should be lwork = 3* ncv**2 + 5* ncv
REAL(KIND=DBL), ALLOCATABLE :: rwork(:)      ! Real valued work arrray of size ncv
INTEGER :: info                              ! Error flag

! iparam variables
INTEGER :: ishift=1     ! Implicit shifts?  Usually yes, thus ishift = 1
INTEGER :: maxitr=300   ! Number of Arnoldi iterations allowed
INTEGER :: mode=3       ! What kind of eigenproblem? Mode=3 works well finding things closest to the shift


! ZNEUPD variables
LOGICAL :: rvec                              ! Recover the eigenvectors
CHARACTER(LEN=1) :: HOWMNY = 'P'             ! Schur vectors, Ritz vectors do not work
INTEGER :: ldz                               ! Leading dimension of z, which is n
INTEGER :: ierr                              ! Error flag
LOGICAL, ALLOCATABLE :: selecte(:)           ! Work array to compute the Ritz vectors.  Size = ncv
COMPLEX(KIND=DBL), ALLOCATABLE :: d(:)       ! Array of eigenvaules.  Should be size = Nev+1
COMPLEX(KIND=DBL), ALLOCATABLE :: z(:,:)     ! Array of eigenvectors.  Size = n x nev
COMPLEX(KIND=DBL), ALLOCATABLE :: workev(:)  ! Work array of size 2*ncv

! LAPACK VARIABLES
INTEGER :: i
COMPLEX(KIND=DBL), ALLOCATABLE :: b(:)       ! zero array used for LU decomposition
INTEGER, ALLOCATABLE :: IPIV(:)              ! integer array not referenced but could hold pivots
INTEGER :: icant                             ! Error because it just can't
CHARACTER :: trans = 'N'                     ! Specifies the linear solver.  Might be a different letter if transposes/hermition congj are used.  See zgetrs


! Affix the input variables to ARPACK variables
n=nsize
nev = numeigen
ncv = nev*2


! Initialize ARPACK variables
! ZNAUPD variables
info = 0
ido = 0
which = 'SR'
ldv = n
iparam(1) = ishift
iparam(3) = maxitr
iparam(7) = mode
lwork = 3*ncv*ncv+5*ncv
ALLOCATE(resid(1:n))
ALLOCATE(v(1:n,1:ncv))
ALLOCATE(workd(1:3*n))
ALLOCATE(workl(1:lwork))
ALLOCATE(rwork(1:ncv))

! ZNEUPD variables
rvec = .FALSE.
ldz = n
ALLOCATE(selecte(1:ncv))
ALLOCATE(d(1:nev+1))
ALLOCATE(z(1:n,1:ncv))
ALLOCATE(workev(1:2*ncv))

! Create the pivot point for the diagonalization
! That is B = (A - sigma * I)
! where A is the input matrix and sigma is the value that we want eigenvalues near
! (I is the identity matrix)
DO i=1, n
   mat(i,i) = mat(i,i) - sigma
ENDDO


! Factor the matrix using LU decomposition calling the LAPACK routine
ALLOCATE(IPIV(1:n))

CALL ZGETRF(n, n, mat, n, IPIV, icant)

IF( icant .NE. 0) THEN
   WRITE(*,*) "Error in ZGETRF, info =", icant
   STOP
ENDIF

! ARPACK Arnoldi implicit restart algorithm

DO WHILE(ido == 0 .OR. ido == 1 .OR. ido == -1)

     CALL znaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, &
          ipntr, workd, workl, lwork, rwork, info)

     IF (ido .eq. -1 .or. ido .eq. 1) THEN

        ! Copy the array over
        CALL zcopy (n, workd(ipntr(1)), 1, workd(ipntr(2)), 1)

        CALL ZGETRS('N', n, 1, mat, n, IPIV, workd(ipntr(2)), n, icant)
        IF( icant .NE. 0) THEN
           WRITE(*,*) "Error in ZGETRS, info =", icant
           STOP
        ENDIF

     END IF
END DO

IF ( info < 0 ) THEN

    PRINT *, ' '
    PRINT *, ' Error with _naupd, info = ', info
    PRINT *, ' Check the documentation of _naupd'
    PRINT *, ' '

ELSE

   ! Cannot get HOWMNY='A' to work.  Instead, use the trick where we send v instead of z in parameter 5
    CALL zneupd (rvec, HOWMNY, selecte, d, v, ldv, sigma, workev, bmat, n, which, &
         nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lwork, rwork, ierr)

    IF ( ierr .NE. 0) THEN
        PRINT *, ' '
        PRINT *, ' Error with _neupd, info = ', ierr
        PRINT *, ' Check the documentation of _neupd. '
        PRINT *, ' '

    END IF
!!$c
!!$c        %-------------------------------------------%
!!$c        | Print additional convergence information. |
!!$c        %-------------------------------------------%
!!$c
   IF (info == 1) THEN
      PRINT *, ' '
      PRINT *, ' Maximum number of iterations reached.'
      PRINT *, ' '
   ELSEIF ( info == 3) THEN
      PRINT *, ' ' 
      PRINT *, ' No shifts could be applied during implicit',' Arnoldi update, try increasing NCV.'
      PRINT *, ' '
   ENDIF

   PRINT *, ' '
   PRINT *, '_NDRV2 '
   PRINT *, '====== '
   PRINT *, ' '
   PRINT *, ' Size of the matrix is ', n
   PRINT *, ' The number of Ritz values requested is ', nev
   PRINT *, ' The number of Arnoldi vectors generated (NCV) is ', ncv
   PRINT *, ' What portion of the spectrum: ', which
   PRINT *, ' The number of Implicit Arnoldi update iterations taken is ', iparam(3)
   PRINT *, ' The number of OP*x is ', iparam(9)
   PRINT *, ' The convergence criterion is ', tol
   PRINT *, ' '
ENDIF

DO i=1, nev
   eigenvs(i) = d(i)
   PRINT *, d(i)
ENDDO

! DEALLOCATE
DEALLOCATE(resid)
DEALLOCATE(v)
DEALLOCATE(workd)
DEALLOCATE(workl)
DEALLOCATE(rwork)
DEALLOCATE(selecte)
DEALLOCATE(d)
DEALLOCATE(z)
DEALLOCATE(workev)
DEALLOCATE(IPIV)

END SUBROUTINE arnoldi
