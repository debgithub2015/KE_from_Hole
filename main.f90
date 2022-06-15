
!A simple program to test the *.wfn reader and the grid.  It also checks the
!kinetic energy functional program.

PROGRAM testKE

USE kinds; USE variables_wfn, ONLY: wfn_filename
USE inout_module; USE KE_module; USE KE_vars; USE normalize_module
USE grid_vars, ONLY: radscale,radfac,Wbecke

IMPLICIT NONE

INTEGER :: i,cnt_failures,cnt_grid    
REAL(dbl) :: normalization_error,normalization(1:2)

!cnt_failures counts the number of times the program has failed.
!cnt_expand counts the number of times the program has increased the grid size.
!normalization_error estimates the error in normalization.
!normalization(1:2) -- the normalization of the spin-densities.

!The first item in the command line is the wfn file.
call getarg(1,wfn_filename)

!Put in wfn filename.
wfn_filename = TRIM(wfn_filename)//'.wfn'
WRITE(*,*) 'Reading wfn file: ',wfn_filename
WRITE(*,*) ' '

normalization_tolerance = 1.0D-3
max_iter = 100
cnt_failures = 0
! Jvectors -- the number of vectors kept to store the Jacobian in limited-memory        !
!             Broyden's bad method.  This doesn't affect the cost of the iteration much,!
!             but round-off error makes this converge more slowly if chosen too large.  !
!             If you use a negative number, this defaults to 5.  Zero is diagonal-only  !
!             update.  
! Jdiag_update -- the diagonal is updated every this-many iterations.  Doubles the cost !
!                 of the iteration but increases accuracy.                              !
! Jtrust_radius -- controls when and how the local trust radii are updated.             !
!                             0  --  no local trust radius used.                        !
!                         other   --  local trust radii used to construct the next step.!
! A previous version of the program allowed one to reject certain portions of a previous!
! step using the the trust radius; this never seems to help.                            !



!This reads the input file containing the grid parameters and also parses the
!*.wfn file.  It evaluates the density and the density gradient.
call startup()


DO i=0,5

cnt_failures = 0
!We loop through the rest of the program until we converge or give up.
DO
   IF (cnt_failures > 0) THEN
      !Change this line to decide how many "failures" you allow.
      WRITE(*,*) "We weren't able to find a grid that allowed us to converge"
      WRITE(*,*) "the calculation."
      WRITE(*,*) "*************TOTAL AND IRREVERSIBLE CONVERGENCE FAILURE*******************"
      EXIT
   ELSE IF (cnt_failures < 0) THEN
      !Successful convergence.  Call output subroutine.
      EXIT
   ELSE IF (cnt_failures > 0) THEN
      !Scale the radfac and radscale and generate the new grid/density/gradient
      radscale = radscale*.80_dbl   !Scale the radius down by 20%
      !Choose a grid that gives sufficient accuracy.
      cnt_grid = 0
      DO
         !Compute normalization.
         cnt_grid = cnt_grid + 1
         call re_start() 
         normalization(1) = DOT_PRODUCT(Wbecke,rho(:,1))
         normalization(2) = DOT_PRODUCT(Wbecke,rho(:,2))
         normalization_error = MAXVAL(ABS(NINT(normalization(1:2))-normalization(1:2)))
         WRITE(*,*) ' '
         WRITE(*,*) 'Normalization of alpha & beta spin ',normalization(1:2)
         WRITE(*,*) 'Estimated Normalization Error      ', normalization_error        
         IF (normalization_error < normalization_tolerance) THEN
            !The grid is accurate enough.
            WRITE(*,*) 'Compressing grid to obtain more compact representation of points.'
            WRITE(*,*) 'atomic radii are scaled by: ',radscale
            WRITE(*,*) 'number of radial points relative to default Becke grid: ',radfac
            WRITE(*,*) ' '
            EXIT         
         ELSE IF (cnt_grid > 5) THEN
            !We couldn't find a good enough grid.  Proceed anyway.
            WRITE(*,*) 'Compressing grid to obtain more compact representation of points.'
            WRITE(*,*) 'atomic radii are scaled by: ',radscale
            WRITE(*,*) 'number of radial points relative to default Becke grid: ',radfac
            WRITE(*,*) ' '
            WRITE(*,*) '*****WARNING: This grid may not be accurate enough to meet your******'
            WRITE(*,*) '*****         accuracy specifications.                         ******'
            WRITE(*,*) ' '
            EXIT    
         ELSE
            !Increase the number of radial grid points.
            radfac = radfac/SQRT(.80_dbl) !This scaling means that the most distant shell of 
                                          !grid points is roughly the same distance away as  
                                          !in the initial grid.  
         ENDIF
      ENDDO
   ELSE
      !This is the first iteration; no need to do anything
      CONTINUE
   ENDIF
      
   !Compute conventional kinetic energy formulae
   call compute_ke_conventional()

   !Compute 1-point and LDA-style normalization.  These methods give the same
   !result for any p-value.
   WRITE(*,*) '-------------------------------------------------------------------'
   WRITE(*,*) 'Uniform Electron Gas normalization and 1-pt. normalization gives   '
   WRITE(*,*) 'kinetic energies and effective Fermi momenta that are invariant to '
   WRITE(*,*) 'the value of p that is used in the generalize p-mean.              '
   WRITE(*,*) ' '

   !Normalize the holes.  
   call hole_normalize(0,cnt_failures)      !the number here is the number of points used in the
                               !normalization.  0 = LDA; 1 = average density approx.
                               !2 = exact normalization.
   call hole_normalize(1,cnt_failures)

   Jvectors = 7
   Jdiag_update = 1
   Jtrust_radius = 0
   p_mean = 2.0_dbl**i    !The p_mean is 1,2,4,8,16,32

   WRITE(*,*) '-------------------------------------------------------------------'
   WRITE(*,*) 'Bad-Broyden with 7 vectors and Diagonal updates; No Local Trust-Radius Pruning.'
   WRITE(*,*) 'p value for generalized p-mean: ',p_mean
   WRITE(*,*) ' '
   !Compute improved "normalized hole" kinetic energy formulae
   call hole_normalize(2,cnt_failures)
   WRITE(*,*) ' '
   WRITE(*,*) '-------------------------------------------------------------------'
   WRITE(*,*) 'p value for generalized p-mean: ',p_mean
   WRITE(*,*) ' '
   call compute_ke_improved()
   WRITE(*,*) ' '
   WRITE(*,*) 'Bad-Broyden with 7 vectors and Diagonal updates; No Local Trust-Radius Pruning.'
   WRITE(*,*) 'p value for generalized p-mean: ',p_mean
   WRITE(*,*) '-------------------------------------------------------------------'
ENDDO

ENDDO

END PROGRAM testKE