!-------------------------------------------------------------------------------!
!                               NORMALIZE_MODULE                                !
!                                                                               !
! This module contains the subroutines for computing the normalization factors  !
! for the exchange-correlation hole.                                            !
!-------------------------------------------------------------------------------!

MODULE NORMALIZE_MODULE

CONTAINS

!-------------------------------------------------------------------------------!
!                       hole_normalize                                          !
!                                                                               !
! This subroutine performs normalizes the exchange-correlation hole.            !
!                                                                               !
! First it chooses the same normalization, kf, that would be used for the       !
! uniform electron gas.                                                         !
!                                                                               !
! Second, it refines this by writing the normalization condition,               !
!         -1 = INT (rho(r) hxc(kf;r,r') dr) for kf(r').                         !
! which on the grid gives:                                                      !
!        SUM(i) rho(i) hxc(i,j,kf) w(i) + 1 = 0 solved for kf(j)                !
! This equation can be solved by a one-dimensional Newton's method with a trust !
! radius that prevents kf from increase/decreasing more than a set amount. The  !
! uniform electron gas normalization is used as an initial guess.               !
!                                                                               !
! The problem with method two is that it uses only one point to define the      !
! width of the hole, so it is not symmetric.  It is better to use two points.   !
!         -1 = INT(rho(r) hxc(kf(r,r');r,r') dr                                 !
! where                                                                         !
!         kf(r,r') = {(kf(r)^(1/p) + kf(r')^(1/p))/2}^p                         !
! is the generalized p-mean of the kf values.  This equation then is solved for !
! kf(r).                                                                        !
!                                                                               !
! The equations are ill-conditioned when the density is very small.  The energy !
! is also not very sensitive to the behavior to the behavior when the energy    !
! is small.  For this reason, the *ACTUAL* equations that are implemented are   !
! weighted by the density at the reference point.  I.e., we actually solve:     !
!     -rho(r') = rho(r') * INT(rho(r) hxc(kf(r,r');r,r') dr                     !
!-------------------------------------------------------------------------------!
!                           DICTIONARY OF LOCAL VARIABLES                       !
!                                                                               !
! refpoints -- How many reference points do we use to define the hole.          !
!              This can be 1 (method 2) or 2 (method 3).                        !
! spin -- counter for spin variables.                                           !
! debug -- =.true. prints debug print info.                                     !
! localKE -- local kinetic energy for debug print.                              !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM OTHER MODULES                            !
! ----The definitive statement of these definitions is found in            ---- !
! ----another module!!!                                                    ---- !
! n_grid -- number of points on the grid, including the points used to compute  !
!          the derivatives.                                                     !
! rho(1:n_grid,1:2) -- the spin-density.                                        !
! a_LDA(1:ngrid,1:3) -- the alpha values used to control the curvature of the   !
!                       hole at same-spin electron coalescence.                 !
!                   1.  From uniform electron gas.                              !
!                   2.  Based on one-point normalization.                       !
!                   3.  Based on exact normalization.                           !
! Jacobian_type -- The type of Jacobian being used.                             !
!                      0  -- diagonal.                                          !
!                     -1  -- Broyden's bad method.                              !
!                      1  -- Broyden's good method.                             !
!                     -2  -- Broyden's bad method with diagonal-only update     !
!                      2  -- Broyden's good method with diagonal-only update.   !
!                     -3  -- Broyden's bad method with diagonal update and      !
!                            perturbative off-diagonal update.                  !
!                      3  -- Broyden's good method with diagonal-only update.   ! 
!                            perturbative off-diagonal update.                  !
!                     -4  -- Broyden's bad method with perturbative update to   !
!                            capture changes on the diagonal.                   !
!                      4  -- Broyden's good method with perturbative update to  !
!                            capture changes on the diagonal.                   !  
!-------------------------------------------------------------------------------!

subroutine hole_normalize(refpoints)

USE kinds
USE constants
USE grid_vars, ONLY: n_grid,Wbecke
USE KE_vars, ONLY: a_LDA,rho,KE_conv,D2g_LDA,normalization_tolerance,Jacobian_type
USE variables_wfn, ONLY: multiplicity

IMPLICIT NONE

INTEGER(istd) :: spin,refpoints,i
REAL(dbl), ALLOCATABLE :: localKE(:)
LOGICAL       :: debug

debug = .false.

IF (debug) THEN
   ALLOCATE(localKE(1:n_grid))
ENDIF

!Initialize *EVERYTHING* to LDA.
IF (refpoints == 0) THEN
   FORALL(spin=1:2)
         a_LDA(:,1,spin) = (6*Pi**2*rho(:,spin))**(1.0_dbl/3)
   ENDFORALL
   !Initialize the other cases to LDA.
   a_LDA(:,2,:) = a_LDA(:,1,:)
   a_LDA(:,3,:) = a_LDA(:,1,:) 
ELSE IF (refpoints == 1) THEN
   !Call Newton's method to do the normalization,
   !First make sure the effective kf values are not so close to zero that we
   !will have convergence problems.
   WHERE (a_LDA(:,2,:) < SQRT(EPSILON(normalization_tolerance)))
         a_LDA(:,2,:) = SQRT(EPSILON(normalization_tolerance))
   ENDWHERE
   !alpha-spin
   call Newton1pt(a_LDA(:,2,1),rho(:,1),normalization_tolerance)
   IF (multiplicity == 1) THEN
      a_LDA(:,2,2) = a_LDA(:,2,1)
   ELSE
      !Solve for beta-spin.
      call Newton1pt(a_LDA(:,2,2),rho(:,2),normalization_tolerance)
   ENDIF
   !Initialize the 2-point formula to the one-point value.
   !This is sometimes so small that it causes problems.  We increase these
   !elements when they are much smaller than the accuracy we can expect.
   a_LDA(:,3,:) = a_LDA(:,2,:)
ELSE
   !Call Newton's method to do the normalization.  The initial guess is the
   !previous 2-point optimization that you performed *unless* you have 
   !previously done a 2-point optimization, in which case you reuse that value. 
   !alpha-spin
   !First make sure the effective kf values are not so close to zero that we
   !will have convergence problems.
   a_LDA(:,3,:) = a_LDA(:,2,:)
   WHERE (a_LDA(:,3,:) < EPSILON(normalization_tolerance))
         a_LDA(:,3,:) = EPSILON(normalization_tolerance)
   ENDWHERE
   call Newton2pt(a_LDA(:,3,1),rho(:,1),normalization_tolerance,Jacobian_type)
   IF (multiplicity == 1) THEN
      a_LDA(:,3,2) = a_LDA(:,3,1)
   ELSE
      !Solve for beta-spin.
      call Newton2pt(a_LDA(:,3,2),rho(:,2),normalization_tolerance,Jacobian_type)
   ENDIF
ENDIF

IF (debug) THEN
   localKE(:) = -3/10*rho(:,1)*D2g_LDA*a_LDA(:,refpoints,1)**2                   &
                    - 3/10*rho(:,2)*D2g_LDA*a_LDA(:,refpoints,2)**2
   WRITE(*,*) ' '
   WRITE(*,*) 'KE with ',refpoints,'normalization:  ',DOT_PRODUCT(Wbecke,localKE) + KE_conv(2)
ENDIF

IF (debug) THEN
   DEALLOCATE(localKE)
ENDIF

end subroutine hole_normalize

!-------------------------------------------------------------------------------!
!                           Newton1pt                                           !
!                                                                               !
! This solver uses Newton's method, with a trust radius, to solve for the       !
! normalization condition in the "exact" case.  This is a "fullgrid" Newtons    !
! method because all of the coefficients a_???(:,3) are solved for at           !
! the same time.  Since refpoints = 1, the Jacobian is diagonal and             !
! this is equivalent to n_grid 1-D Newton solves.                               !
!-------------------------------------------------------------------------------!
! anow -- the current vector of "effective fermi momenta."                      !
! fnow -- the current normalization errors.                                     !
! trust_radii -- the vector of trust radii for the steps.                       !
! fpredicted -- the vector of predicted values for the normalization errors.    !
! fnew -- the vector of new function values.                                    !
! anew -- the vector of new Fermi momenta (the guess).                          !
! trust_ratio -- the ratio that is used to control the trust radius.  If the    !
!                ratio between the "expected improvement" and the "actual       !
!                improvement is close to zero, we decrease the trust radius.    !
!                If it is close to one, we increase the trust radius.  Otherwise!
!                we leave the trust radius unchanged.                           !
! value -- the initial guess on input; the solution upon output.                !
! refpoints -- the number of reference points used.                             !
! Jnew -- the diagonal element of the Jacobian in the next step.                !
! Jnow -- the diagonal elements of the current Jacobian.                        !
! updated -- keeps track of which steps were accepted and which were rejected.  !
! rho -- the spin-density of the channel currently under consideration.         !
! debug -- a debug print flag.                                                  !
! max_error -- maximum error in nonlinear equations.                            !
! avg_error -- average error in normalization.                                  !
! L_lingood -- .true. if the linear model is trustworthy.                       !
! L_linbad -- .true. if the linear model is bad.                                !
! max_iter -- the maximum number of iterations.                                 !
! norm_error -- the error in the normalization.                                 !
!-------------------------------------------------------------------------------!

subroutine Newton1pt(value,rho,normalization_tolerance)

USE kinds
USE grid_vars, ONLY: n_grid,Wbecke

IMPLICIT NONE

REAL(dbl) :: normalization_tolerance,value(1:n_grid)
REAL(dbl) :: max_error,avg_error,norm_error
INTEGER(istd) :: iter,i,max_iter
REAL(dbl) :: anow(1:n_grid),fnow(1:n_grid),anew(1:n_grid),fnew(1:n_grid)
REAL(dbl) :: trust_radii(1:n_grid),fpredicted(1:n_grid),temp(1:n_grid)
REAL(dbl) :: trust_ratio(1:n_grid),Jnew(1:n_grid),Jnow(1:n_grid),rho(1:n_grid)
LOGICAL :: updated(1:n_grid),debug,L_linbad(1:n_grid),L_lingood(1:n_grid)

!Initialize to guessed "effective kf" value.  Choose predicted function value
!to be zero, as this is the extreme case and is unlikely to be attained.
anew(:) = value(:)
fpredicted(:) = 0.0_dbl 
fnow(:) = 100000.0_dbl
Jnow(:) = 1000000.0_dbl
trust_radii(:) = MAX(ABS(anew(:))/3,normalization_tolerance/10)
iter = 1
max_iter = 50

!Set debug=.true. for debug printing.
debug = .true.

IF (debug) THEN
   WRITE(*,*) ' '
   WRITE(*,*) 'For 1 point normalization.'
   WRITE(*,*) 'Number of electrons of this spin, ', DOT_PRODUCT(Wbecke,rho)
   WRITE(*,*) 'Iteration info:'
ENDIF

DO   
   !Default is to update all points
   updated(:) = .true.
   L_lingood(:) = .false.
   L_linbad(:) = .false.

   !Call subroutine to evaluate the error in the normalization of the hole.
   call normalizationUEG1pt(rho,anew,fnew,Jnew)
   !End if converged.
   !Compute the maximum error in normalization and the integrated error.
   max_error = MAXVAL(abs(fnew))
   avg_error = SUM(abs(fnew)*Wbecke)
   norm_error = ABS(SUM(fnew*Wbecke))
   
   IF (debug) THEN   !Print debug info if not converged.
      WRITE(*,*) iter,' the maximum error: ',max_error
      WRITE(*,*) '                   total error: ',avg_error
      WRITE(*,*) '           normalization error: ',norm_error
   ENDIF
   
   IF (max_error < normalization_tolerance                                      &
       .and. avg_error < normalization_tolerance) THEN
      WRITE(*,*) '1-point hole normalization converged.'
      WRITE(*,*) ' '
      EXIT
   ENDIF
   
   !If not converged, then we should evaluate whether the "guessed" improvement
   !was at all close to the "actual" improvement.
   FORALL(i=1:n_grid,fpredicted/=fnow)
         trust_ratio(i) = min(abs(fpredicted(i)-fnow(i)), abs(fnew(i)-fnow(i)))   &
                           /max(abs(fpredicted(i)-fnow(i)), abs(fnew(i)-fnow(i)))
   ENDFORALL
   
   WHERE (abs(fnew) >= abs(fnow))
         !Reject move if further from solution than before
         !do *NOT* update guess; try again with smaller step.
         !update trust radius to force smaller step (eventually)
         trust_radii = trust_radii/2.0_dbl
         updated = .false.
   ELSEWHERE (trust_ratio > .8_dbl)
         !The Newton step gave an accurate prediction of the step.  Take the step
         !and increase the trust radius.  Compute the Jacobian.  Because every iteration
         !is expensive, it is very bad to reject an iteration.  So we expand the trust radius
         !conservatively but reduce the trust radius aggressively.
         trust_radii =  trust_radii*1.5_dbl
         anow = anew; fnow = fnew; Jnow = Jnew
         L_lingood = .true.         
   ELSEWHERE (trust_ratio > .2_dbl) 
         !The Newton step was not that accurate but it didn't mess us up too badly.  Keep
         !the old trust radius and make the Newton step.
         anow = anew; fnow = fnew; Jnow = Jnew
   ELSEWHERE 
         !The Newton step was quite inaccurate; reduce the trust radius but take the step
         !since at least it got us closer (ever-so-slightly) to the solution.
         trust_radii = trust_radii/1.5_dbl
         anow = anew; fnow = fnew; Jnow = Jnew
         L_linbad = .true.
   END WHERE
   
   IF (debug) THEN
      WRITE(*,*) '             ',COUNT(updated), ' of ',n_grid,'points improved.'
      WRITE(*,*) '             ',COUNT(L_lingood), ' of ',n_grid,'points had an good linear update.'
      WRITE(*,*) '             ',COUNT(L_linbad), ' of ',n_grid,'points had an bad linear update.'
   ENDIF
                              
   !Compute the predicted step.  The output of this is a new choice for a (anew)
   !and a predicted value for the normalization error after updating a (fpredicted)
   call diagNewton_Step_UEG(rho,anow,fnow,Jnow,1,trust_radii,anew,fpredicted,   &
                            normalization_tolerance)
   
   iter = iter + 1
   IF (iter > max_iter) THEN
      WRITE(*,*) 'WARNING:  Failure to converge the normalization constraint.'
      WRITE(*,*) 'This is usually because your integration grid is not accurate'
      WRITE(*,*) 'enough. The program will continue but you should be careful '
      WRITE(*,*) 'about interpreting the results.'
      WRITE(*,*) 'Program failed in solving full-grid Newton method.'
      WRITE(*,*) 'Consider printing debug info.'
      WRITE(101,*) 'Program filed in solving full-grid Newton method for the so-called'
      WRITE(101,*) 'exact normalization scheme.  Consider printing'
      WRITE(101,*) 'debug info.'
      EXIT
   ENDIF
   !Every now and then we should reinitialize the trust radii because they can 
   !get "trapped" at a value that is too small.  So
   IF (mod(iter,20) == 0) THEN
      trust_radii(:) = ABS(anow(:))/2
      WRITE(*,*) '......resetting trust radii.......'
   ENDIF
ENDDO

!After we exit the loop, we set the latest a-value to return:
value = anew

end subroutine Newton1pt

!-------------------------------------------------------------------------------!
!                           normalizationUEG1pt                                 !
!                                                                               !
! This subroutine evaluates the normalization error and the Jacobian of the     !
! nonlinear equations to be solved for the 1-point normalization case.          !
!-------------------------------------------------------------------------------!
!                               DICTIONARY OF VARIABLES                         !
!                                                                               !
! a -- effective kf value.                                                      !
! f -- error in normalization of the hole, = <rho(r)*h(r,r')> - (-1).           !
! Jdiag -- the diagonal element of the Jacobian of the nonlinear equations,     !
!          f = 0.                                                               !
! gUEG(1:n_grid,1:n_grid) -- the square root of the exchange hole in the        !
!                            uniform electron gas.                              !
! dgUEG(1:n_grid,1:n_grid) -- the derivative of gUEG with respect to the a      !
!                             (the effective kf value).  Used to compute the    !
!                             Jacobian.                                         !
! dij(1:n_grid,1:n_grid) -- the distance between two grid points.               !
!-------------------------------------------------------------------------------!
! n_grid -- the number of grid points.                                          !
! Wbecke -- the integration weights.                                            !
!-------------------------------------------------------------------------------!

subroutine normalizationUEG1pt(rho,a,f,Jdiag)

USE kinds
USE grid_vars, ONLY: n_grid,Wbecke,XYZbecke

IMPLICIT NONE

REAL(dbl) :: rho(1:n_grid),a(1:n_grid),f(1:n_grid),Jdiag(1:n_grid)
REAL(dbl) :: gUEG(1:n_grid,1:n_grid),dgUEG(1:n_grid,1:n_grid)
REAL(dbl) :: dij(1:n_grid,1:n_grid)
INTEGER(istd) :: k,l

LOGICAL :: debug

debug =.false.

!Compute the distance between the grid points
FORALL(k=1:n_grid,l=1:n_grid)
      dij(k,l) = SQRT(DOT_PRODUCT(XYZbecke(1:3,k)-XYZbecke(1:3,l),                  &
                             XYZbecke(1:3,k)-XYZbecke(1:3,l)))
ENDFORALL

!We do a simple normalization integral with a one-point
!"width" of the hole.  This gives us, as equations to solve,
! 0 = [ INT(rho(r)h(a|r-r'|)dr) - (-1) ] * rho(r')

!The default value of the g function is 1.0, which occurs at coalescence (dij= 0)
FORALL(k=1:n_grid,l=1:n_grid, (a(k)*dij(k,l)<=.0000001_dbl))
      gUEG(k,l) = 1.0_dbl
      dgUEG(k,l) = 0.0_dbl
ENDFORALL
!If the points are not too close together, then compute the term when necessary.
FORALL(k=1:n_grid,l=1:n_grid, (a(k)*dij(k,l)>.0000001_dbl))
      gUEG(k,l) = 3*(sin(a(k)*dij(k,l)) - a(k)*dij(k,l)*cos(a(k)*dij(k,l))) &
                    /(a(k)*dij(k,l))**3 
ENDFORALL
FORALL(k=1:n_grid,l=1:n_grid, (a(k)*dij(k,l)>.0000001_dbl))
      dgUEG(k,l) = 3*(sin(a(k)*dij(k,l)) - a(k)*dij(k,l)*gUEG(k,l))         &
                    /(a(k)*dij(k,l))**2
ENDFORALL
FORALL(k=1:n_grid)
      f(k) = (DOT_PRODUCT(Wbecke(:),rho(:)*(-1)*gUEG(k,:)**2) + 1)*rho(k)
      Jdiag(k) = (2*DOT_PRODUCT(Wbecke(:),rho(:)*(-1)*gUEG(k,:)*dij(k,:)*dgUEG(k,:)))*rho(k)
ENDFORALL

end subroutine normalizationUEG1pt

!-------------------------------------------------------------------------------!
!                           normalizationUEG2pt                                 !
!                                                                               !
! This subroutine evaluates the normalization error for the 2-pt hole           !
! normalization with the p-mean.                                                !
!                                                                               !
! This particular subroutine is for the uniform electron gas hole, but the      !
! modifications to consider other holes are minor.                              !
!-------------------------------------------------------------------------------!
!                               DICTIONARY OF VARIABLES                         !
!                                                                               !
! a -- effective kf value.                                                      !
! f -- error in normalization of the hole, = <rho(r)*h(r,r')> - (-1).           !
! gUEG(1:n_grid,1:n_grid) -- the square root of the exchange hole in the        !
!                            uniform electron gas.                              !
! dij(1:n_grid,1:n_grid) -- the distance between two grid points.               !
! pdist(1:n_grid,1:n_grid) -- the "effective" kF value times the distance       !
!                             between points in the 2-point model.              !
!-------------------------------------------------------------------------------!
! n_grid -- the number of grid points.                                          !
! Wbecke -- the integration weights.                                            !
! p_mean -- the p-value that specifies the appropriate "generalized mean."      !
!-------------------------------------------------------------------------------!

subroutine normalizationUEG2pt(rho,a,f)

USE kinds
USE grid_vars, ONLY: n_grid,Wbecke,XYZbecke
USE KE_vars, ONLY: p_mean

IMPLICIT NONE

REAL(dbl) :: rho(1:n_grid),a(1:n_grid),f(1:n_grid)
REAL(dbl) :: gUEG(1:n_grid,1:n_grid)
REAL(dbl) :: dij(1:n_grid,1:n_grid),pdist(1:n_grid,1:n_grid)
INTEGER(istd) :: refpoints,k,l

LOGICAL :: debug

debug =.false.

!Compute the distance between the grid points
FORALL(k=1:n_grid,l=1:n_grid)
      dij(k,l) = SQRT(DOT_PRODUCT(XYZbecke(1:3,k)-XYZbecke(1:3,l),                  &
                             XYZbecke(1:3,k)-XYZbecke(1:3,l)))
ENDFORALL

!Two reference points will be used but we are using the quasi-Newton update and
!don't need the diagonal element of the Jacobian.
FORALL(k=1:n_grid,l=1:n_grid)
      pdist(k,l) = ((a(k)**p_mean + a(l)**p_mean)/2.0_dbl)**(1.0_dbl/p_mean)*dij(k,l)
ENDFORALL
!The default value of the g function is 1.0, which occurs at coalescence (dij= 0)
FORALL(k=1:n_grid,l=1:n_grid, (pdist(k,l)<=.0000001_dbl))
         gUEG(k,l) = 1.0_dbl
ENDFORALL   
!If the points are not too close together, use both points.
FORALL(k=1:n_grid,l=1:n_grid, (pdist(k,l)>.0000001_dbl))
         gUEG(k,l) = 3.0_dbl * (sin(pdist(k,l))                                 &
                           - pdist(k,l)*cos(pdist(k,l))) &
                         /(pdist(k,l))**3 
ENDFORALL
FORALL(k=1:n_grid)
         f(k) = (DOT_PRODUCT(Wbecke(:),rho(:)*(-1.0_dbl)*gUEG(k,:)**2) + 1.0_dbl)*rho(k)
ENDFORALL

end subroutine normalizationUEG2pt



!-------------------------------------------------------------------------------!
!                           normalizationUEG2ptJacobian                         !
!                                                                               !
! This subroutine evaluates the normalization error and the diagonal element    !
! of the Jacobian of the nonlinear equations to be solved for the 2-point       !
! normalization.                                                                !
!                                                                               !
! This particular subroutine is for the uniform electron gas hole, but the      !
! modifications to consider other holes are minor.                              !
!-------------------------------------------------------------------------------!
!                               DICTIONARY OF VARIABLES                         !
!                                                                               !
! a -- effective kf value.                                                      !
! f -- error in normalization of the hole, = <rho(r)*h(r,r')> - (-1).           !
! Jdiag -- the diagonal element of the Jacobian of the nonlinear equations,     !
!          f = 0.                                                               !
! gUEG(1:n_grid,1:n_grid) -- the square root of the exchange hole in the        !
!                            uniform electron gas.                              !
! dgUEG(1:n_grid,1:n_grid) -- the derivative of gUEG with respect to the a      !
!                             (the effective kf value).  Used to compute the    !
!                             Jacobian.                                         !
! dij(1:n_grid,1:n_grid) -- the distance between two grid points.               !
! pdist(1:n_grid,1:n_grid) -- the "effective" kF value times the distance       !
!                             between points in the 2-point model.              !
!-------------------------------------------------------------------------------!
! n_grid -- the number of grid points.                                          !
! Wbecke -- the integration weights.                                            !
! p_mean -- the p-value that specifies the appropriate "generalized mean."      !
!-------------------------------------------------------------------------------!

subroutine normalizationUEG2ptJacobian(rho,a,f,Jdiag)

USE kinds
USE grid_vars, ONLY: n_grid,Wbecke,XYZbecke
USE KE_vars, ONLY: p_mean

IMPLICIT NONE

REAL(dbl) :: rho(1:n_grid),a(1:n_grid),f(1:n_grid),Jdiag(1:n_grid)
REAL(dbl) :: gUEG(1:n_grid,1:n_grid),dgUEG(1:n_grid,1:n_grid)
REAL(dbl) :: dij(1:n_grid,1:n_grid),pdist(1:n_grid,1:n_grid)
INTEGER(istd) :: k,l

LOGICAL :: debug

debug =.false.

!Compute the distance between the grid points
FORALL(k=1:n_grid,l=1:n_grid)
      dij(k,l) = SQRT(DOT_PRODUCT(XYZbecke(1:3,k)-XYZbecke(1:3,l),              &
                             XYZbecke(1:3,k)-XYZbecke(1:3,l)))
ENDFORALL

!Two reference points will be used and the diagonal of the Jacobian will be
!generated.
FORALL(k=1:n_grid,l=1:n_grid)
      pdist(k,l) = ((a(k)**p_mean + a(l)**p_mean)/2.0_dbl)**(1.0_dbl/p_mean)*dij(k,l)
ENDFORALL

!The default value of the g function is 1.0, which occurs at coalescence (dij= 0)
FORALL(k=1:n_grid,l=1:n_grid, ((pdist(k,l)<=.0000001_dbl)                       &
                              .or. (max(a(k),a(l)) <= EPSILON(a(k))**2)))
      gUEG(k,l) = 1.0_dbl
      dgUEG(k,l) = 0.0_dbl
ENDFORALL   
!If the points are not too close together, use both points.
FORALL(k=1:n_grid,l=1:n_grid, ((pdist(k,l) >.0000001_dbl)                       &
                              .and. (max(a(k),a(l)) > EPSILON(a(k))**2)))
      gUEG(k,l) = 3.0_dbl * (sin(pdist(k,l)) - pdist(k,l)*cos(pdist(k,l)))      &
                      /(pdist(k,l))**3 
ENDFORALL
FORALL(k=1:n_grid,l=1:n_grid, ((pdist(k,l) >.0000001_dbl)                       &
                              .and. (max(a(k),a(l)) > EPSILON(a(k))**2)))
      dgUEG(k,l) = 3.0_dbl*(sin(pdist(k,l)) - pdist(k,l)*gUEG(k,l))             &
                     /(pdist(k,l))**2
ENDFORALL
FORALL(k=1:n_grid)
      f(k) = (DOT_PRODUCT(Wbecke(:),rho(:)*(-1.0_dbl)*gUEG(k,:)**2) + 1.0_dbl)*rho(k)
ENDFORALL
FORALL(k=1:n_grid,(a(k) > EPSILON(a(k))**2))
      Jdiag(k) = (DOT_PRODUCT(Wbecke(:),rho(:)*(-1.0_dbl)*gUEG(k,:)*dgUEG(k,:)  &
                                       *a(k)**(p_mean-1)*pdist(k,:)             &
                                       /((a(k)**p_mean+a(:)**p_mean)/2)))*rho(k)
ENDFORALL
FORALL(k=1:n_grid,(a(k) <= EPSILON(a(k))**2))
      Jdiag(k) = EPSILON(a(k))**2
ENDFORALL

IF (debug) THEN
   DO k=1,n_grid,100
      WRITE(*,*) rho(k),a(k),f(k),Jdiag(k)
   ENDDO
ENDIF

end subroutine normalizationUEG2ptJacobian



!-------------------------------------------------------------------------------!
!                           diagNewton_Step_UEG                                 !
!                                                                              !
! This subroutine determines the Newton step when only the diagonal is used.    !
! This is the case for the 1-pt. formula and also for the "diagonal only" case  !
! of the two-point normalization.                                               !
!-------------------------------------------------------------------------------!
! anow -- the current vector of "effective fermi momenta."                      !
! fnow -- the current normalization errors.                                     !
! trust_radii -- the vector of trust radii for the steps.                       !
! fpredicted -- the vector of predicted values for the normalization errors.    !
! anew -- the vector of new Fermi momenta (the guess).                          !
! trust_ratio -- the ratio that is used to control the trust radius.  If the    !
!                ratio between the "expected improvement" and the "actual       !
!                improvement is close to zero, we decrease the trust radius.    !
!                If it is close to one, we increase the trust radius.  Otherwise!
!                we leave the trust radius unchanged.                           !
! refpoints -- the number of reference points used.                             !
! Jii -- the diagonal element of the Jacobian                                   !
! rho -- the spin-density of the channel currently under consideration.         !
! step -- the step that will be taken.  At the beginning it is the full step,   !
!         then it is scaled back.                                               !
! n_grid -- the number of grid points.                                          !
! fullstep -- the full Newton step.                                             !
!-------------------------------------------------------------------------------!

subroutine diagNewton_Step_UEG(rho,anow,fnow,Jii,refpoints,trust_radii,         &
                               anew,fpredicted,normalization_tolerance)

USE kinds
USE grid_vars, ONLY: n_grid,Wbecke

IMPLICIT NONE

REAL(dbl) :: rho(1:n_grid),anow(1:n_grid),fnow(1:n_grid),Jii(1:n_grid)
REAL(dbl) :: trust_radii(1:n_grid),anew(1:n_grid),fpredicted(1:n_grid)
REAL(dbl) :: step(1:n_grid),fullstep(1:n_grid),normalization_tolerance
INTEGER(istd) :: refpoints,j
LOGICAL :: debug

debug = .true.

!Don't allow extremely small trust radii
FORALL(j=1:n_grid)
      trust_radii(j) = MAX(trust_radii(j),normalization_tolerance**2)
ENDFORALL

!First compute the Newton step.
IF (refpoints == 2) THEN
   !We need to update the Jacobian
   call normalizationUEG2ptJacobian(rho,anow,fnow,Jii)   
ENDIF

WHERE(Jii /= 0)
     !Compute Newton step.
     fullstep = -1*fnow/Jii
ELSEWHERE
     !If the Jacobian is zero, then this is a strange situation unless the
     !point has almost no density.  We can choose the step in the right direction
     !(decrease a if fnow > 0; increase otherwise)
     fullstep = -1*fnow
ENDWHERE

!Start by trying the full step
step = fullstep

IF (debug) THEN
   WRITE(*,*) SUM(Wbecke*ABS(fnow+Jii*step))
ENDIF

!Now we need to scale back the step to agree with the trust radius
WHERE (ABS(fullstep) > trust_radii)
      step = SIGN(trust_radii,step)
END WHERE

IF (debug) THEN
   WRITE(*,*) SUM(Wbecke*ABS(fnow+Jii*step))
ENDIF

!This forall loop ensures that we never take a step that makes the
!effective Fermi momentum negative.  This is important because a<=0 is
!absurd.
FORALL (j=1:n_grid, step(j) < -9*anow(j)/10)
       step(j) = -9*anow(j)/10
       trust_radii(j) = anow(j)/2   
ENDFORALL
    
anew = anow + step
!Make sure that anew is positive
anew = MAX(anew,0.0_dbl)

fpredicted = fnow + Jii*step

end subroutine diagNewton_Step_UEG

!-------------------------------------------------------------------------------!
!                           Newton2pt                                           !
!                                                                               !
! This solver uses Newton's method, with a trust radius, to solve for the       !
! normalization condition in the "exact" case.  This is a "fullgrid" Newtons    !
! method because all of the coefficients a_???(:,3) are solved for at           !
! the same time.                                                                !
!-------------------------------------------------------------------------------!
! anow -- the current vector of "effective fermi momenta."                      !
! da -- the change in the a values.                                             !
! fnow -- the current normalization errors.                                     !
! trust_radii -- the vector of trust radii for the steps.                       !
! fpredicted -- the vector of predicted values for the normalization errors.    !
! fnew -- the vector of new function values.                                    !
! anew -- the vector of new Fermi momenta (the guess).                          !
! aold -- the "effective Fermi momenta" from the previous step.                 !
! fold -- the vector of function values from the previous step.                 !
! trust_ratio -- the ratio that is used to control the trust radius.  If the    !
!                ratio between the "expected improvement" and the "actual       !
!                improvement is close to zero, we decrease the trust radius.    !
!                If it is close to one, we increase the trust radius.  Otherwise!
!                we leave the trust radius unchanged.                           !
! value -- the initial guess on input; the solution upon output.                !
! updated -- keeps track of which steps were accepted and which were rejected.  !
! rho -- the spin-density of the channel currently under consideration.         !
! debug -- a debug print flag.                                                  !
! max_error -- maximum error in nonlinear equations.                            !
! avg_error_old -- average absolute error from previous step.                   !
! avg_error_new -- average absolute error from present step.                    !
! norm_error -- the error in the normalization.                                 !
! pred_error -- the predicted average absolute error.                           !
! list_error(1:5) -- the list of the errors in the most recent iterations.      !
! L_lingood -- .true. if the linear model is trustworthy.                       !
! L_linbad -- .true. if the linear model is bad.                                !
! max_iter -- the maximum number of iterations.                                 !
! iter_backtrace -- number of iterations in the backtrace.                      !
! Jdiag_old(1:n_grid) -- the "old" diagonal from the previous iteration.        !
! Jacobian_type -- The type of Jacobian being used.                             !
!                      0  -- diagonal.                                          !
!                     -1  -- Broyden's bad method.                              !
!                      1  -- Broyden's good method.                             !
!                     -2  -- Broyden's bad method with diagonal-only update     !
!                      2  -- Broyden's good method with diagonal-only update.   !
!                     -3  -- Broyden's bad method with diagonal update and      !
!                            perturbative off-diagonal update.                  !
!                      3  -- Broyden's good method with diagonal-only update.   ! 
!                            perturbative off-diagonal update.                  !
!                     -4  -- Broyden's bad method with perturbative update to   !
!                            capture changes on the diagonal.                   !
!                      4  -- Broyden's good method with perturbative update to  !
!                            capture changes on the diagonal.                   !
!-------------------------------------------------------------------------------!

subroutine Newton2pt(value,rho,normalization_tolerance,Jacobian_type)

USE kinds
USE grid_vars, ONLY: n_grid,Wbecke

IMPLICIT NONE

REAL(dbl) :: normalization_tolerance,value(1:n_grid)
REAL(dbl) :: max_error,norm_error,list_error(1:5),avg_error_new,avg_error_old,pred_error
INTEGER(istd) :: iter, refpoints,i,max_iter,Jacobian_type,iter_backtrace
REAL(dbl) :: anow(1:n_grid),fnow(1:n_grid),anew(1:n_grid),fnew(1:n_grid),da(1:n_grid)
REAL(dbl) :: aold(1:n_grid),fold(1:n_grid)
REAL(dbl) :: trust_radii(1:n_grid),fpredicted(1:n_grid),temp(1:n_grid)
REAL(dbl) :: trust_ratio(1:n_grid),Jnew(1:n_grid),Jnow(1:n_grid),rho(1:n_grid)
REAL(dbl) :: JBroydenInv(1:n_grid,1:n_grid),Jdiag_old(1:n_grid)
LOGICAL :: updated(1:n_grid),debug,L_linbad(1:n_grid),L_lingood(1:n_grid)

!Initialize to guessed "effective kf" value.  Choose predicted function value
!to be zero, as this is the extreme case and is unlikely to be attained.
anew(:) = value(:)
aold(:) = value(:)
fpredicted(:) = 0.0_dbl 
fnow(:) = 100000.0_dbl
list_error(:) = 100000.0_dbl
avg_error_old = 100000.0_dbl
avg_error_new = 100000.0_dbl
FORALL(i=1:n_grid)
      trust_radii(i) = MAX(anow(i)/3,normalization_tolerance/10)
ENDFORALL
iter = 1
max_iter = 75

!Set debug=.true. for debug printing.
debug = .true.

IF (debug) THEN
   WRITE(*,*) ' '
   WRITE(*,*) 'For 2 point normalization.'
   WRITE(*,*) 'Number of electrons of this spin, ', DOT_PRODUCT(Wbecke,rho)
   WRITE(*,*) ' ' 
   WRITE(*,*) 'Iteration info:'
ENDIF

DO   
   !Default is to update all points
   updated(:) = .true.
   L_lingood(:) = .false.
   L_linbad(:) = .false.

   !Call subroutine to evaluate the error in the normalization of the hole.  
   !This is used only to accept/reject moves, and the Jacobian is not needed.
   call normalizationUEG2pt(rho,anew,fnew)
   !End if converged.
   !Compute the maximum error in normalization and the integrated error. 
   max_error = MAXVAL(abs(fnew))
   avg_error_new = SUM(Wbecke*abs(fnew))
   norm_error = ABS(SUM(Wbecke*fnew))
   pred_error = SUM(Wbecke*abs(fpredicted))
   !Update array of errors.
   DO i=1,4
         list_error(i) = list_error(i+1)
   ENDDO
   list_error(5) = avg_error_new
   
   IF (debug) THEN   !Print debug info.
      WRITE(*,*) '        maximum absolute error: ',max_error
      WRITE(*,*) '           normalization error: ',norm_error
      WRITE(*,*) '    predicted total abs. error: ',pred_error
      WRITE(*,*) '          total absolute error: ',avg_error_new
      WRITE(*,*) '                     step size: ',SQRT(SUM((anew-aold)*(anew-aold)))

   ENDIF
   
   IF (max_error < normalization_tolerance                                      &
       .and. avg_error_new < normalization_tolerance) THEN
      WRITE(*,*) 'SUCCESS!!!  2-point hole normalization converged.'
      WRITE(*,*) ' '
      EXIT
   ENDIF
   
   
   !If the rms error increase, reject the step and reduce the trust radii.
   IF (avg_error_new > avg_error_old) THEN
      WRITE(*,'(I3,A)') iter,'  Avg. Absolute Error increased; will backtrace to find good step.'   
      iter_backtrace = 1
      da = anew - aold
      DO
         trust_radii = trust_radii/2
         anew = aold + da*2**(-1*iter_backtrace)
         call normalizationUEG2pt(rho,anew,fnew)  
         max_error = MAXVAL(abs(fnew))
         avg_error_new = SUM(Wbecke*abs(fnew))
         norm_error = ABS(SUM(Wbecke*fnew))
         IF (avg_error_new <= avg_error_old) THEN
            WRITE(*,*) iter_backtrace,'          max absolute error: ',max_error
            WRITE(*,*) '                     normalization error: ',norm_error
            WRITE(*,*) '                    total absolute error: ',avg_error_new
            avg_error_old = avg_error_new
            list_error(5) = avg_error_new
            anow = anew
            EXIT
         ELSE IF (iter_backtrace > 10) THEN
            anow = anew
            WRITE(*,*) 'Backtrace did not converge.'
            WRITE(*,*) iter_backtrace,'          max absolute error: ',max_error
            WRITE(*,*) '                     normalization error: ',norm_error
            WRITE(*,*) '                    total absolute error: ',avg_error_new
            EXIT
         ELSE
            !WRITE(*,*) '      ....halving step size....'
            iter_backtrace = iter_backtrace + 1            
         ENDIF         
      ENDDO 
   ELSE
      avg_error_old = avg_error_new  
      
      !Global adjustments of trust radii
      IF ((MAXVAL(list_error)-MINVAL(list_error)) < normalization_tolerance) THEN
         !We are having trouble converging.  Try shorter steps but accept the move.
         WRITE(*,'(I3,A)') iter,'  convergence is slow.  Try reducing trust radius.'     
         trust_radii = trust_radii/1.5
      ELSE
         !Since the error went down, slightly increase all the trust radii.
         trust_radii = MIN(1.5,(avg_error_old/avg_error_new))*trust_radii
      ENDIF
      
      !Now Adjust the trust radii individually.
      WHERE (abs(fnew) >= abs(fold))
            !Reject move if further from solution than before
            !do *NOT* update guess; try again with smaller step.
            !update trust radius to force smaller step (eventually)
            trust_radii = trust_radii/2
            anow = aold
            updated = .false.
      ELSEWHERE((abs(fnew) > abs(fold)/2) .and. (fold*fnew < 0))
            !The step was OK, but it was too far and it changed the sign of f.
            trust_radii = trust_radii/1.5
            anow = anew
            L_linbad = .true.   
      ELSEWHERE((abs(fnew) > abs(fold)/2) .and. (fold*fnew > 0))
            !The step was OK, but it wasn't far enough.  
            trust_radii =  trust_radii*1.5_dbl
            anow = anew
            L_linbad = .true. 
      ELSEWHERE((abs(fnew) < abs(fold)/2) .and. (fold*fnew > 0))
            !The step was quite good, but it wasn't far enough.  
            trust_radii =  trust_radii*2
            anow = anew
            L_lingood = .true.
      ELSEWHERE      
            !The step was quite good, but we slightly overshot the solution.  Keep the
            !trust radius the same.
            anow = anew
            L_lingood = .true.
      END WHERE
      IF (debug) THEN
         WRITE(*,'(I3,A,3(2X,F3.0))') iter,'  (a) % points that improved (b) w/ good update (c) w/ bad update:',  &
                              100.0_dbl*COUNT(updated)/n_grid,                                 &
                              100.0_dbl*COUNT(L_lingood)/n_grid,                                 &
                              100.0_dbl*COUNT(L_linbad)/n_grid                                 
      ENDIF
   ENDIF                           
   !Compute the predicted step.  The output of this is a new choice for a (anew)
   !and a predicted value for the normalization error after updating a (fpredicted)
   IF (Jacobian_type == 0) THEN
      !We use only diagonal element of the Jacobian.  fnow and Jdiag_old are only
      !dummy arguments, and immediately overwritten inside the subroutine.
      call diagNewton_Step_UEG(rho,anow,fnow,Jdiag_old,2,trust_radii,           &
                               anew,fpredicted,normalization_tolerance)
   ELSE
      !We are going to be using a quasi-Newton method *IF* it seems likely to
      !help.  The basic idea is to use the diagonal approximation until it stagnates
      !or for 10 iterations, whichever comes first.  
      !At that stage, the diagonal element of the Jacobian should be about right,
      !and it is probably sensible to start working on the off-diagonal.
      !We also periodically reset the Jacobian to its diagonal.  
      IF (((MAXVAL(list_error)-MINVAL(list_error)) > normalization_tolerance    &
           .and. iter < 10)                                                     &
           .or. (MOD(iter,40) == 0)) THEN
         PRINT*, '====Diagonal Update====='
         !In the last 5 iterations we improved the calculation significantly or
         !it is time to reset the diagonal.  So
         call diagNewton_Step_UEG(rho,anow,fnow,Jdiag_old,2,trust_radii,        &
                                  anew,fpredicted,normalization_tolerance)  
         !Update the Jacobian in case you want to do real quasi-Newton next time.
         JBroydenInv = 0.0_dbl
         !First get rid of NaN's and 0's in the diagonal of the Jacobian:
         WHERE(ABS(Jdiag_old) > (EPSILON(normalization_tolerance))**(-1))
              Jdiag_old = (SIGN(EPSILON(normalization_tolerance),Jdiag_old))**(-1)
         ELSEWHERE (ABS(Jdiag_old) < EPSILON(normalization_tolerance))
              Jdiag_old = SIGN(EPSILON(normalization_tolerance),Jdiag_old)
         ENDWHERE
         !Then update the inverse Jacobian.
         FORALL(i=1:n_grid)
               JBroydenInv(i,i) = Jdiag_old(i)**(-1)
         ENDFORALL
      ELSE 
         PRINT*, '====Quasi Newton Update Type',Jacobian_type,'===='
         !Do quasi-Newton step.  We need to give the step (anow - anew)
         !and then we will (in the subroutine) compute the update to the Jacobian
         call Broyden_Step_UEG(rho,anow,fnow,aold,fold,anew,fpredicted,         &
                               Jdiag_old,JBroydenInv,Jacobian_type,             &
                               trust_radii,normalization_tolerance)
      ENDIF
   ENDIF
   !Store the values from the previous step
   aold = anow
   fold = fnow   
   iter = iter + 1
   IF (iter > max_iter) THEN
      WRITE(*,*) 'WARNING:  Failure to converge the normalization constraint.'
      WRITE(*,*) 'This is usually because your integration grid is not accurate'
      WRITE(*,*) 'enough. The program will continue but you should be careful '
      WRITE(*,*) 'about interpreting the results.'
      WRITE(*,*) 'Program failed in solving full-grid Newton method.'
      WRITE(*,*) 'Consider printing debug info.'
      WRITE(101,*) 'Program filed in solving full-grid Newton method for the so-called'
      WRITE(101,*) 'exact normalization scheme.  Consider printing'
      WRITE(101,*) 'debug info.'
      EXIT
   ENDIF
   
   !Every now and then we should reinitialize the trust radii because they can 
   !get "trapped" at a value that is too small.  It is also useful to reset the
   !trust radius if we have gone five iterations without much improvement.
   IF (mod(iter,20) == 1)  THEN
      FORALL(i=1:n_grid)
            trust_radii(i) = MAX(anow(i)/2,normalization_tolerance)
      ENDFORALL
      WRITE(*,*) '......resetting trust radii.......'
   ENDIF
   
ENDDO

!Return the best value you found:
value = anew

end subroutine Newton2pt


!-------------------------------------------------------------------------------!
!                           BroydenStep_UEG                                     !
!                                                                               !
! This subroutine determines the quasi-Newton step from the Broyden-based       !
! Jacobian.                                                                     !
!-------------------------------------------------------------------------------!
! anow -- the current vector of "effective fermi momenta."                      !
! fnow -- the current normalization errors.                                     !
! trust_radii -- the vector of trust radii for the steps.                       !
! fpredicted -- the vector of predicted values for the normalization errors.    !
! fpredicteddiag -- the vector of predicted values based on the diagonal approx.!
!                   to the Jacobian.                                            !
! anew -- the vector of new Fermi momenta (the guess).                          !
! aold -- the "effective Fermi momenta" from the previous step.                 !
! fold -- the vector of function values from the previous step.                 !
! rho -- the spin-density of the channel currently under consideration.         !
! debug -- a debug print flag.                                                  !
! Jdiag_old(1:n_grid) -- the "old" diagonal from the previous iteration.        !
! Jdiag_new(1:n_grid) -- the "new" diagonal from the present iteration.         !
! Japprox(1:n_grid,1:n_grid) -- the approximate Jacobian used to predict the    !
!                               improvement in the solution.                    !
! dx(1:n_grid) -- an intermediate quantity used to construct the Jacobian       !
!                       update.                                                 !
! tmp(1:n_grid) -- a temporary vector used in Broyden's "good" update.          !
! denominator -- the denominator in the Hessian updates.                        !
! dJdiag(1:n_grid) -- the change in the diagonal of the Jacobian.               !
! Jacobian_type -- The type of Jacobian being used.                             !
!                      0  -- diagonal.                                          !
!                     -1  -- Broyden's bad method.                              !
!                      1  -- Broyden's good method.                             !
!                     -2  -- Broyden's bad method with diagonal-only update     !
!                      2  -- Broyden's good method with diagonal-only update.   !
!                     -3  -- Broyden's bad method with diagonal update and      !
!                            perturbative off-diagonal update.                  !
!                      3  -- Broyden's good method with diagonal-only update.   ! 
!                            perturbative off-diagonal update.                  !
!                     -4  -- Broyden's bad method with perturbative update to   !
!                            capture changes on the diagonal.                   !
!                      4  -- Broyden's good method with perturbative update to  !
!                            capture changes on the diagonal.                   !
! step -- the step that will be taken.  At the beginning it is the full step,   !
!         then it is scaled back.                                               !
! n_grid -- the number of grid points.                                          !
! fullstep -- the full Newton step.                                             !
! df(1:n_grid) -- the change in the normalization function.                     !
! da(1:n_grid) -- the change in the effective Fermi momenta.                    !
!-------------------------------------------------------------------------------!

subroutine Broyden_Step_UEG(rho,anow,fnow,aold,fold,anew,fpredicted,            &
                               Jdiag_old,Jinv,Jacobian_type,                    &
                               trust_radii,normalization_tolerance)

USE kinds
USE grid_vars, ONLY: n_grid

IMPLICIT NONE

REAL(dbl) :: rho(1:n_grid)
REAL(dbl) :: anow(1:n_grid),aold(1:n_grid),anew(1:n_grid)
REAL(dbl) :: fnow(1:n_grid),fold(1:n_grid)
REAL(dbl) :: fpredicted(1:n_grid),fpredicteddiag(1:n_grid)
REAL(dbl) :: df(1:n_grid),da(1:n_grid),dx(1:n_grid),denominator,tmp(1:n_grid)
REAL(dbl) :: Jdiag_old(1:n_grid),Jdiag_new(1:n_grid),dJdiag(1:n_grid)
REAL(dbl) :: Jinv(1:n_grid,1:n_grid),Japprox(1:n_grid,1:n_grid)
REAL(dbl) :: trust_radii(1:n_grid),normalization_tolerance
REAL(dbl) :: step(1:n_grid),fullstep(1:n_grid)
INTEGER(istd) :: i,j,k,l,Jacobian_type
LOGICAL :: debug

debug = .true.

!We first need to compute the function value for the "current" value the
!effective Fermi momenta, anow(:).  
IF (ABS(Jacobian_type) == 1) THEN
   !We don't need the diagonal element of the Jacobian.  So
   call normalizationUEG2pt(rho,anow,fnow)
ELSE
   !We do need the diagonal element of the Jacobian.
   call normalizationUEG2ptJacobian(rho,anow,fnow,Jdiag_new)
   !Get rid of NaN's and 0's in the diagonal of the Jacobian.
   WHERE(ABS(Jdiag_new) > (EPSILON(normalization_tolerance))**(-1))
        Jdiag_new = (SIGN(EPSILON(normalization_tolerance),Jdiag_new))**(-1)
   ELSEWHERE (ABS(Jdiag_new) < EPSILON(normalization_tolerance))
        Jdiag_new = SIGN(EPSILON(normalization_tolerance),Jdiag_new)
   ENDWHERE
ENDIF

!Now we update the Jacobian.  We start by updating the diagonal.  We then
!update the rest of the matrix.  This ensures that the matrix is consistent
!with the secant approximation, even though it does (slightly) mess up the
!diagonal.
IF (ABS(Jacobian_type) == 2 .or. ABS(Jacobian_type) == 3) THEN
   !Update the diagonal element.  So
   FORALL(i=1:n_grid)
         JInv(i,i) = (Jdiag_new(i))**(-1)
   ENDFORALL
   IF (ABS(Jacobian_type) == 3) THEN
      !Also update the off-diagonal elements.
      FORALL(k=1:n_grid,l=1:n_grid,(k /= l))
            Jinv(k,l) = (Jinv(k,k)*Jinv(l,l)*Jdiag_new(k)*Jdiag_new(l))**(-1)   &
                         * Jinv(k,l)
      ENDFORALL
   ENDIF
   Jdiag_old = Jdiag_new
ELSE IF (ABS(Jacobian_type) == 4) THEN
   !We will use perturbation theory to update the diagonal.  Compute the change
   !in the Jacobian
   dJdiag = Jdiag_new - Jdiag_old
   Jdiag_old = Jdiag_new      
ENDIF

!Compute the updated Jacobian
da = anow - aold
df = fnow - fold
dx = MATMUL(Jinv,df)
IF (Jacobian_type > 0) THEN
   !Use Broyden's Good Method to do the update.
   denominator = SUM(df*dx)
   !Use Japprox to store the low-rank update matrix.
   tmp(:) = MATMUL(df,Jinv)/denominator
   FORALL(k=1:n_grid,l=1:n_grid)
         Jinv(k,l) = Jinv(k,l) + (da(k)-dx(k))*tmp(l)
   ENDFORALL   
ELSE 
   !Use Broyden's Bad Method to do the update.
   denominator = SUM(df*df)
   FORALL(k=1:n_grid,l=1:n_grid)
         Jinv(k,l) = Jinv(k,l) + (da(k)-dx(k))*df(l)/denominator
   ENDFORALL
ENDIF

!Compute the Newton Step:
fullstep = -1*MATMUL(Jinv,fnow)

!If we are using the Diagonal Perturbation update, there is a second term.  There
!are third and higher order terms, but we leave them out because it gets very
!expensive if we include them.
IF (ABS(Jacobian_type) == 4) THEN
   fullstep = fullstep - MATMUL(Jinv,dJdiag*fullstep)
ENDIF

!Start by trying the full step
step = fullstep

!Now we need to scale back the step to agree with the trust radius
WHERE (ABS(fullstep) > trust_radii)
      step = SIGN(trust_radii,step)
END WHERE

!This forall loop ensures that we never take a step that reduces the value
!of a by more than a factor of 1/2.  This is important because a<=0 is
!absurd, and we want to approach a .approx. zero very gingerly. . . ..
!The trust radius must be reset b/c otherwise one "falsely accepts" a short
!step as very reliable and increases the trust radius.
FORALL (j=1:n_grid, step(j) < -1*anow(j)/2)
       step(j) = -1*anow(j)/2
       trust_radii(j) = anow(j)/3  
ENDFORALL

!Make the step.    
anew = anow + step

!Now we have a problem.  We really need a prediction of the step, but we
!only have the *INVERSE* Jacobian and not the real Jacobian.  We approximate the
!Jacobian by assuming that the inverse is diagonally dominant.  We know the
!exact diagonal element, so we only need to approximate the off-diagonal.
FORALL(k=1:n_grid)
      !Approx. diagonal element.
      Japprox(k,k) = Jdiag_new(k)
ENDFORALL
FORALL(k=1:n_grid,l=1:n_grid,(k /= l))
      Japprox(k,l) = -1 * (Jinv(k,k)*Jinv(l,l))**(-1) * Jinv(k,l)
ENDFORALL

fpredicted = fnow + MATMUL(Japprox,step)
fpredicteddiag = fnow + Jdiag_new*step

!Choose the smaller of the two predictions; this if anything demands "too much"
!rather than "too little" from the calculation.
WHERE (ABS(fpredicted) > ABS(fpredicteddiag))
      fpredicted = fpredicteddiag
ENDWHERE

end subroutine Broyden_Step_UEG

end module NORMALIZE_MODULE