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
!-------------------------------------------------------------------------------!

subroutine hole_normalize(refpoints)

USE kinds
USE constants
USE grid_vars, ONLY: n_grid,Wbecke
USE KE_vars, ONLY: a_LDA,rho,KE_conv,D2g_LDA

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
   !alpha-spin
   call Newton1pt(a_LDA(:,2,1),rho(:,1))
   IF (multiplicity == 1) THEN
      a_LDA(:,2,2) = a_LDA(:,2,1)
   ELSE
      !Solve for beta-spin.
      call Newton1pt(a_LDA(:,2,2),rho(:,2))
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
   call Newton2pt(a_LDA(:,3,1),rho(:,1))
   IF (multiplicity == 1) THEN
      a_LDA(:,3,2) = a_LDA(:,3,1)
   ELSE
      !Solve for beta-spin.
      call Newton2pt(a_LDA(:,3,2),rho(:,2))
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

subroutine Newton1pt(value,rho)

USE kinds
USE KE_vars, ONLY: normalization_tolerance
USE grid_vars, ONLY: n_grid,Wbecke
USE holesubroutines

IMPLICIT NONE

REAL(dbl) :: value(1:n_grid)
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
   call normalization1pt(rho,anew,fnew,Jnew)
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
   call diagNewton_Step_UEG(anow,fnow,Jnow,trust_radii,anew,fpredicted)
   
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
      trust_radii(:) = MAX(ABS(anow(:))/2,normalization_tolerance)
      WRITE(*,*) '......resetting trust radii.......'
   ENDIF
ENDDO

!After we exit the loop, we set the latest a-value to return:
value = anew

end subroutine Newton1pt


!-------------------------------------------------------------------------------!
!                           diagNewton_Step_UEG                                 !
!                                                                               !
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

subroutine diagNewton_Step_UEG(anow,fnow,Jii,trust_radii,         &
                               anew,fpredicted)

USE kinds
USE grid_vars, ONLY: n_grid,Wbecke
USE holesubroutines
USE KE_vars, ONLY: normalization_tolerance

IMPLICIT NONE

REAL(dbl) :: rho(1:n_grid),anow(1:n_grid),fnow(1:n_grid),Jii(1:n_grid)
REAL(dbl) :: trust_radii(1:n_grid),anew(1:n_grid),fpredicted(1:n_grid)
REAL(dbl) :: step(1:n_grid),fullstep(1:n_grid)
INTEGER(istd) :: refpoints,j
LOGICAL :: debug

debug = .false.

WHERE(Jii /= 0)
     !Compute Newton step.
     fullstep = -1*fnow/Jii
ELSEWHERE
     !If the Jacobian is zero, then this is a strange situation unless the density
     !is also zero.  Try choosing the "a" value for this point as zero.
     fullstep = -1*anow
ENDWHERE

!Start by trying the full step
step = fullstep

IF (debug) THEN
   WRITE(*,*) SUM(Wbecke*ABS(fnow+Jii*step))
ENDIF

!Now we need to scale back the step to agree with the trust radius
FORALL(j=1:n_grid,(ABS(fullstep(j)) > trust_radii(j)))
      step(j) = SIGN(trust_radii(j),step(j))
ENDFORALL

IF (debug) THEN
   WRITE(*,*) SUM(Wbecke*ABS(fnow+Jii*step))
ENDIF

!This forall loop ensures that we never take a step that makes the
!effective Fermi momentum negative.  This is important because a<=0 is
!absurd.
FORALL (j=1:n_grid, step(j) < -1*anow(j)/2)
       step(j) = -1*anow(j)/2
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
! value -- the initial guess on input; the solution upon output.                !
! rho -- the spin-density of the channel currently under consideration.         !
!                                                                               !
! anew(:) -- the newest value of the effective Fermi momentum.                  !
! fnew(:) -- the newest value of the function.                                  !
! aold(:) -- the previous value of the function.                                !
! fold(:) -- the previous value of the function.                                !
! Jiinew(:) -- the newest value of the Jacobian's diagonal.                     !
! abest(:) -- the effective fermi momentum for the best point found thus far.   !
! fbest(:) -- the normalization error for the best point found thus far.        !
! Jiibest(:) -- the diagonal of the Jacobian for the best point found thus far. !
! a(:,1:Jvectors) -- the effective Fermi momenta from the previous Jvectors steps.!
! f(:,1:Jvectors) -- the effective Fermi momenta from the previous Jvectors steps.!
! fguess(:) -- the guess for the next step based on the diagonal approx. of the !
!              Jacobian. This is just for debugging.                            !
!                                                                               !
! trust_radii(:) -- local trust radii to control individual components of the   !
!                   step.                                                       !
! Global_trust -- the global trust radius.                                      !
! step_size -- stores the step size in cases where it is annoying to compute it !
!              in situ.                                                         !
! step -- the step when we use the diagonal approximation.                      !
!                                                                               !
! fratio(:) -- the ratio between the "latest" f value and the "best previous" f.!
! max_error -- maximum error in nonlinear equations.                            !
! avg_error_old -- average absolute error from previous step.                   !
! avg_error_new -- average absolute error from present step.                    !
! avg_error_best -- the lowest average absolute error found so far.             !
! norm_error -- the error in the normalization.                                 !
! pred_error -- the predicted average absolute error.                           !
! list_error(1:5) -- the list of the errors in the most recent iterations.      !
! integration_error -- an estimate of the error in the numerical integration.   !
!                                                                               !
! Lbwardstep(:) -- a logical variable that says whether the step was backwards, !
!                  i.e., if the step was in the opposite direction that         !
!                  intuition would suggest.                                     !
! Lcutstep(:) -- a logical variable that says whether the step was cut in order !
!                to agree with the trust radius.                                !
! Lincrtrust(:) -- a logical variable that says whether the trust radius was    !
!                  increased.                                                   !
! Ldecrtrust(:) -- a logical variable that says whether the trust radius was    !
!                  decreased.                                                   !   
! Limprove(:) -- keeps track of whether the objective function improved.        !
! Lshort(:) -- steps that are too short.                                        !
! Llong(:) -- steps that are too long.                                          !
! Lgoodbward(:) -- steps that were "backwards" but still reduced the error.     !
! Lscalestep -- .true. if the global trust radius was used.                     !
! debug -- a debug print flag.                                                  !
!                                                                               !
! iter -- current iteration number.                                             !
! max_iter -- the maximum number of iterations.                                 !
! iter_bad -- this is the iteration when a step is rejected.                    !
! iter_slow -- this is the iteration when the program is slow.                  !
!                                                                               !
! Nelectrons -- the number of electrons.                                        !
!                                                                               !
! Jdiag_old(1:n_grid) -- the "old" diagonal from the previous iteration.        !
! Jvectors -- the number of vectors kept to store the Jacobian in limited-memory!
!             Broyden's bad method.                                             !
! Jdiag_update -- the diagonal is updated every this-many iterations.           !
! Jtrust_radius -- controls when and how the local trust radii are updated.     !
!                    0  --  no local trust radius used.                         !
!                   >0   --  local trust radii used to construct the next step. !
!                   <0  --  local trust radii used to reject certain portions   !
!                                    of the current step.                       !
! normalization_tolerance -- the calculations converges when the average        !
!                            absolute error and the maximum error are less than !
!                            this.                                              !
!                                                                               !
! i,j,k,l -- counters.                                                          !
!-------------------------------------------------------------------------------!

subroutine Newton2pt(value,rho)

USE kinds
USE grid_vars, ONLY: n_grid,Wbecke
USE holesubroutines
USE KE_vars, ONLY: Jtrust_radius,Jdiag_update,Jvectors,normalization_tolerance,  &
                   max_iter

IMPLICIT NONE

REAL(dbl) :: value(1:n_grid),rho(1:n_grid)

REAL(dbl) :: abest(1:n_grid),fbest(1:n_grid),Jiibest(1:n_grid)
REAL(dbl) :: anew(1:n_grid),fnew(1:n_grid),Jiinew(1:n_grid)
REAL(dbl) :: aold(1:n_grid),fold(1:n_grid)
REAL(dbl) :: a(1:n_grid,1:Jvectors),f(1:n_grid,1:Jvectors)
REAL(dbl) :: fguess(1:n_grid),pred_error,fratio(1:n_grid)

REAL(dbl) :: trust_radii(1:n_grid),Global_trust,step_size,step(1:n_grid)
REAL(dbl) :: max_error,norm_error,list_error(1:5),Nelectrons
REAL(dbl) :: avg_error_new,avg_error_old,avg_error_best,integration_error

INTEGER(istd) :: iter,i,j,k,l,iter_bad,iter_slow

LOGICAL :: debug,Lbwardstep(1:n_grid),Lcutstep(1:n_grid),Lincrtrust(1:n_grid)
LOGICAL :: Ldecrtrust(1:n_grid),Limprove(1:n_grid),Lshort(1:n_grid),LLong(1:n_grid)
LOGICAL :: Lgoodbward(1:n_grid),Lscalestep

!Initialize to guessed "effective kf" value.  Choose predicted function value
!to be zero, as this is the extreme case and is unlikely to be attained.
a(:,:) = -100.0_dbl
f(:,:) = -100.0_dbl
anew(:) = value(:)
aold(:) = -100.0_dbl
abest(:) = -100.0_dbl
fbest(:) = SQRT(HUGE(list_error))
fold(:) = SQRT(HUGE(list_error))
fguess(:) = 0.0_dbl
list_error(:) = SQRT(HUGE(list_error))
avg_error_old = SQRT(HUGE(avg_error_old))
avg_error_new = SQRT(HUGE(avg_error_old))
avg_error_best = SQRT(HUGE(avg_error_old))
trust_radii(:) = MAX(anew(:)/3,normalization_tolerance**2)
Global_trust = SUM(ABS(trust_radii))
iter_bad = 1
iter_slow = 1

!Certain parameters can be used as flags.
IF (Jdiag_update < 1) THEN
   !By default, update the diagonal of the Jacobian in every iteration.
   Jdiag_update = 1
ENDIF
IF (Jvectors < 0) THEN
   !By default, keep 5 vectors.
   Jvectors = 5
ENDIF

iter = 1

!Set debug=.true. for debug printing.
debug = .true.

!Estimate the integration error
Nelectrons = DOT_PRODUCT(Wbecke,rho)
integration_error = ABS(NINT(Nelectrons)-Nelectrons)

IF (debug) THEN
   WRITE(*,*) ' '
   WRITE(*,*) 'For 2 point normalization.'
   WRITE(*,*) 'Number of electrons of this spin, ', Nelectrons
   WRITE(*,*) 'Iteration info:'
ENDIF

1110 FORMAT(1X,I3,3X,A,ES12.5)
1111 FORMAT(10X,A,ES12.5)
1112 FORMAT(10X,8(2X,A8))
1113 FORMAT(10X,8(7X,I3))
1115 FORMAT(11X,I3,3X,A,ES12.5)
1116 FORMAT(20X,A,ES12.5)

DO   
   !Default is to update all points
   Lincrtrust(:) = .false.
   Ldecrtrust(:) = .false.
   Lbwardstep(:) = .false.
   Lcutstep(:) = .false.
   Limprove(:) = .true.
   Lshort(:) = .false.
   Llong(:) = .false.
   Lgoodbward(:) = .false.
   Lscalestep = .false.
   
   WHERE((anew-abest)*fbest > 0)
        !If the previous value of the normalization (fbest(i)) is positive(negative), 
        !then the hole is normalized to a number that is too small(big) in absolute
        !magnitude and da = anew - abest should decrease(increase).  When this is 
        !not true it is strange, so we flag this case. 
        Lbwardstep = .true.
   ENDWHERE
   
   !Experience shows that backwards steps are almost *ALWAYS* a bad idea.  Occassionally
   !they don't hurt you in the first few iterations, however.
   IF (iter > 3) THEN
      WHERE(Lbwardstep)
           !If the step is backwards, then we make a very small step in the 
           !direction that causes f to decrease in size.  
           anew = abest - SIGN( MIN(abest*ABS(fbest), abest/50), fbest)
      ENDWHERE
      anew = MIN(abest,0.0_dbl)
   ENDIF
   
   !Call subroutine to evaluate the error in the normalization of the hole.
   IF (iter == 1                                                                &
       .or. (MOD(iter,Jdiag_update) == 0 .and. .not. (Jtrust_radius > 0))) THEN
      !We need to compute the Jacobian.  If Jtrust_radius > 0, then we wait until
      !after we compute the "scaled back" step to compute the Jacobian.
      WRITE(*,1111) "   Updated Jacobian's Diagonal"
      call normalization2ptJacobian(rho,anew,fnew,Jiinew)
   ELSE 
      !We only need to compute the function value.
      WRITE(*,1111) '   No Jacobian Update This Iteration'
      call normalization2pt(rho,anew,fnew)
   ENDIF
   
   fratio(:) = fnew(:)/fbest(:)

   WHERE (ABS(fratio) > 1.0_dbl)
         Limprove = .false.
   ENDWHERE

   !If the result isn't 10 times better than it was before, it is either too
   !short or too long.
   WHERE(fratio > .1_dbl .and. Limprove)
        Lshort = .true.
   ENDWHERE
   WHERE(fratio < -.1_dbl)
        Llong = .true.
   ENDWHERE

   WHERE (Limprove .and. Lbwardstep)
         Lgoodbward = .true.
   ENDWHERE
   
   !Compute average absolute error.
   avg_error_new = SUM(Wbecke*abs(fnew))   
   max_error = MAXVAL(ABS(fnew))
   step_size = SUM(ABS(anew-abest))
      
   !Check for convergence.
   !EXITLOOPEXITLOOPEXITLOOPEXITLOOPEXITLOOPEXITLOOPEXITLOOPEXITLOOPEXITLOOPEXITLOOP
   IF (max_error < normalization_tolerance                                   &
             .and. avg_error_new < normalization_tolerance) THEN
      WRITE(*,*) 'SUCCESS!!!  2-point hole normalization converged.'
      norm_error = ABS(SUM(Wbecke*fnew))
      pred_error = SUM(Wbecke*abs(fguess))  
      WRITE(*,1110) iter, ' In this step the normalization improved by ',     &
                             avg_error_old-avg_error_best
      WRITE(*,1111) '    maximum absolute error: ',max_error
      WRITE(*,1111) '      total absolute error: ',avg_error_new
      WRITE(*,1111) 'predicted total abs. error: ',pred_error
      WRITE(*,1111) '       normalization error: ',norm_error
      WRITE(*,1111) '                 step size: ',step_size
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      EXIT
   ELSE IF (iter > max_iter) THEN
      WRITE(*,*) '***FAILURE***  2-point hole normalization did not converge due to'
      WRITE(*,*) 'too many iterations.'  
      norm_error = ABS(SUM(Wbecke*fnew))
      pred_error = SUM(Wbecke*abs(fguess))  
      step_size = SUM(ABS(anew-abest))
      WRITE(*,1110) iter, ' In this step the normalization improved by ',     &
                             avg_error_old-avg_error_best
      WRITE(*,1111) '    maximum absolute error: ',max_error
      WRITE(*,1111) '      total absolute error: ',avg_error_new
      WRITE(*,1111) 'predicted total abs. error: ',pred_error
      WRITE(*,1111) '       normalization error: ',norm_error
      WRITE(*,1111) '                 step size: ',step_size
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      EXIT
      EXIT
   ELSE
      CONTINUE  !Try again.
   ENDIF     
   !EXITLOOPEXITLOOPEXITLOOPEXITLOOPEXITLOOPEXITLOOPEXITLOOPEXITLOOPEXITLOOPEXITLOOP


   !Decide whether the new step is better or worse than the old step.
   IF ((avg_error_new - avg_error_best) > integration_error) THEN
      !BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP
      !There was no improvement in this step.  Try again with smaller trust radius.
      step_size = SUM(ABS(anew-abest))
      Global_trust = MIN(Global_trust,step_size)
      Global_trust = MAX(Global_trust,normalization_tolerance**2)   
                                                   !Reset trust radius to current step in
                                                   !cases where this was not already true.
      Global_trust = Global_trust/2
      IF (avg_error_new < avg_error_old) THEN
         avg_error_old = avg_error_new
         !It is the second-best value, at least; maybe don't reduce the trust radius
         !**TOO** aggressively.
         Global_trust = Global_trust*1.5
      ENDIF

      WRITE(*,1115) iter_bad, ' The error went ***UP*** by  ',avg_error_new-avg_error_best
      WRITE(*,1116) '               step size: ',step_size
      
      IF (Jvectors > 0) THEN
         IF (iter_bad == 1) THEN
            !We need to store the data from previous iterations to use in the 
            !limited memory quasi-Newton.  We only do this in step 1; subsequent
            !steps are collinear with step 1; we should not update discard data to
            !use those points
            DO i=1,Jvectors-1
               a(:,i) = a(:,i+1)
               f(:,i) = f(:,i+1)
            ENDDO
         ENDIF
         !We always update the last point.  As the step gets smaller and smaller
         !the "secant" approximation we are making gets better and better.
         a(:,Jvectors) = anew(:)   !Keep the "best" point but do use the new
         f(:,Jvectors) = fnew(:)   !data to improve the Quasi-Newton.
      ENDIF  
      
      IF (iter_slow > 1) THEN
         !Try just scaling back the step.
         step = anew - abest
         step = step/1.5
         anew = abest + step
         anew = MAX(0.0_dbl,anew)
         WRITE(*,1116) '           new step size: ',SUM(ABS(anew-abest))         
      ELSE IF (iter_bad < 11) THEN
         !Try Broyden step.
         call Broyden_Step(abest,fbest,Jiibest,a,f,anew,fguess,                 &
                           trust_radii,Global_trust,Lcutstep,Lscalestep)
      ELSE IF (iter_bad < 21) THEN
         WRITE(*,*) '           ***Failed to converge.  Could not find a good trust radius.***'
         WRITE(*,*) '           ***Try steepest descent step.***'
         IF (iter_bad == 11) THEN
            !We're searching in a new direction; reset the trust radius.
            Global_trust = SUM(ABS(abest(:)))/100
         ENDIF
         step = -1*fbest
         step = MAX(step,-1*abest)
         step_size = SUM(ABS(step))
         anew = abest + (Global_trust/step_size)*step
         anew = MAX(anew,0.0_dbl)
         WRITE(*,1116) '           new step size: ',SUM(ABS(anew-abest))
      ELSE IF (iter_bad < 31) THEN
         WRITE(*,*) '           ***Failed to converge.  Could not find a good trust radius.***'
         WRITE(*,*) '           ***Try diagonal Jacobian step.***'
         IF (iter_bad == 21) THEN
            !We're searching in a new direction; reset the trust radius.
            Global_trust = SUM(ABS(abest(:)))/100
         ENDIF
         step = 0.0_dbl
         WHERE (ABS(Jiibest) > EPSILON(normalization_tolerance)**2)
               step = -1*fbest/Jiibest
         ELSEWHERE
               step = -1*abest     !assume that Jiibest being very small means that a should
                                   !be nearly zero.
         ENDWHERE
         step = MAX(step,-1*abest)
         step_size = SUM(ABS(step))
         anew = abest + (Global_trust/step_size)*step
         anew = MAX(anew,0.0_dbl)
         WRITE(*,1116) '           new step size: ',SUM(ABS(anew-abest))
      ELSE 
         !This seems hopeless :-( 
         WRITE(*,*) '*** TOTAL FAILURE TO CONVERGE TO DESIRED ACCURACY ***.'
         WRITE(*,*) '*** :-( :-( :-( :-( :-( :-( :-( :-( :-( :-( :-(   ***.'
         EXIT
      ENDIF

      iter_bad = iter_bad + 1
      
      !BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP BADSTEP
   ELSE
      !The error may have increased, yes, but the error is close enough to the same to be within
      !"rounding error" for this grid.
      !GOODSTEP GOODSTEP GOODSTEP GOODSTEP GOODSTEP GOODSTEP GOODSTEP GOODSTEP GOODSTEP GOODSTEP GOODSTEP 
      !Reset the counter for inner iterations.
      iter_bad = 1

      !This is the best step thus far.  Compute errors.
      avg_error_old = avg_error_best
      avg_error_best = avg_error_new
      norm_error = ABS(SUM(Wbecke*fnew))
      pred_error = SUM(Wbecke*abs(fguess)) 
      step_size = SUM(ABS(anew-abest)) 
      DO i=1,4
         list_error(i) = list_error(i+1)
      ENDDO
      list_error(5) = avg_error_best
        
      WRITE(*,1110) iter, ' In this step the normalization improved by ',     &
                             avg_error_old-avg_error_best
      WRITE(*,1111) '    maximum absolute error: ',max_error
      WRITE(*,1111) '      total absolute error: ',avg_error_best
      WRITE(*,1111) 'predicted total abs. error: ',pred_error
      WRITE(*,1111) '       normalization error: ',norm_error
      WRITE(*,1111) '                 step size: ',step_size
      WRITE(*,1111) '       global trust radius: ',Global_trust
      WRITE(*,*) ' '
      
      IF (Jvectors > 0) THEN
         !We need to store the data from previous iterations to use in the 
         !limited memory quasi-Newton:
         DO i=1,Jvectors-1
            a(:,i) = a(:,i+1)
            f(:,i) = f(:,i+1)
         ENDDO      
         a(:,Jvectors) = abest(:)
         f(:,Jvectors) = fbest(:)
      ENDIF
      aold = abest
      fold = fbest
      abest(:) = anew(:)
      fbest(:) = fnew(:)
      Jiibest(:) = Jiinew(:)
             
      !Check for failure to converge
      IF ( ((MAXVAL(list_error) - MINVAL(list_error))/5*(max_iter-iter+1)        &
           < avg_error_best) .and. (MOD(iter,10) == 0)) THEN
         !Very slow convergence.  Too slow to converge even if linear from here on
         !out.  Try to see if we can find a new direction.
         iter_slow = iter_slow + 1
         WRITE(*,*) '===Very Slow Convergence==='
         IF (iter_slow == 2) THEN
            !Try brute strength step
            WRITE(*,*) '===Try steepest descent step.==='
            step = -1*fbest
            step = MAX(10*step,-1*abest)
            anew = abest + 10*step
            anew = MAX(anew,0.0_dbl)
            WRITE(*,1116) '           new step size: ',SUM(ABS(anew-abest))
            CYCLE
         ELSE IF (iter_slow == 3) THEN
            !Try diagonal step.
            WRITE(*,*) '===Try diagonal Jacobian step.==='
            step = 0.0_dbl
            WHERE (ABS(Jiibest) > EPSILON(normalization_tolerance)**2)
                  step = -1*fbest/Jiibest
            ELSEWHERE
                  step = -1*abest     !assume that Jiibest being very small means that a should
                                      !be nearly zero.
            ENDWHERE
            step = MAX(10*step,-1*abest)
            anew = abest + step
            anew = MAX(anew,0.0_dbl)
            WRITE(*,1116) '           new step size: ',SUM(ABS(anew-abest))
            CYCLE
         ELSE
            CONTINUE
         ENDIF
      ENDIF
      
      iter_slow = 1

      !UPDATE TRUSTRADII --  UPDATE TRUSTRADII --  UPDATE TRUSTRADII --  UPDATE TRUSTRADII --  UPDATE TRUSTRADII
      !Update the trust radii.
      !For the Global Trust Radius.
      IF (COUNT(Lshort) > 2*COUNT(Llong) .and. Lscalestep) THEN
         !There are many more "too short" steps than "too long" steps, and we
         !scaled back the step.  Try increasing it.
         Global_Trust = Global_Trust*MAX(COUNT(Lshort)/(2*COUNT(Llong)+1),2)
      ELSE IF (COUNT(LLong) > 1.5*COUNT(Lshort)) THEN
         !There are many more "too long" steps than "too short" steps.  Reduce
         !the trust radius.
         IF (Lscalestep) THEN
            !We did scale the step back before; we just need to truncate the step
            !more agressively.
            Global_Trust = Global_Trust/MAX( COUNT(Llong)/(1.5*COUNT(Lshort)+1), 2.0_dbl)
         ELSE
            !We need to scale back the trust radius quite a bit; we aren't even
            !scaling back the step yet.  Make the next step a fraction of the length
            !of the current step.
            step_size = SUM(abs(abest-aold))
            Global_Trust = step_size/MAX( COUNT(Llong)/(1.5*COUNT(Lshort)+1), 2.0_dbl)
         ENDIF
      ELSE
         !Keep old trust radius
         CONTINUE
      ENDIF
      WRITE(*,1111) 'updated global trust radius ',Global_trust
    
      WHERE (.not. Limprove) 
            !These are cases where the error became bigger
            !Reduce the trust radius.
            trust_radii = trust_radii/MIN(3.0_dbl, ABS(fratio)+1)
            Ldecrtrust = .true.
      ELSEWHERE (fratio > .5 .and. Lcutstep)
           !These are cases where (a) the error became smaller :-)
           !                      (b) the step was controlled by the trust radius
           !                      (c) the function did not change sign, so probably
           !                          a bigger step would be better.
           trust_radii = trust_radii*1.5
           Lincrtrust = .true.
      ELSEWHERE (fratio < -.5)
           !These are cases where (a) the error became smaller :-)
           !                      (b) the step was too big, and the sign changed.
           trust_radii = trust_radii/1.5
           Ldecrtrust = .true.
      ENDWHERE
      !In other cases you can argue that the method did as it should.  It either
      !(a) reduced the error by 1/2 *because* of the trust radius (good trust radius!)
      !(b) or it didn't use the trust radius at all (no problem with trust radius!)
      !UPDATE TRUSTRADII --  UPDATE TRUSTRADII --  UPDATE TRUSTRADII --  UPDATE TRUSTRADII --  UPDATE TRUSTRADII
      
      !REJECT MOVES?? ??  REJECT MOVES?? ??  REJECT MOVES?? ??  REJECT MOVES?? ??  REJECT MOVES??  ?? REJECT MOVES??
      !**IF** we are screening the points, then do so now.
      IF ((Jtrust_radius > 0) .and. (COUNT(Limprove) < n_grid)) THEN
         anew = abest
         WHERE (.not. Limprove)
               !Whereever the value did not improve, go back to the previous 
               !choice of vectors.
               anew = aold
         ENDWHERE
         
         !We have changed the a-vector, so now we need to *RECOMPUTE* the function
         !and possibly also the Jacobian.  So:
         IF (MOD(iter,Jdiag_update) == 0) THEN
            !We need to compute the Jacobian.  If Jtrust_radius > 0, then we wait until
            !after we compute the "scaled back" step to compute the Jacobian.
            call normalization2ptJacobian(rho,anew,fnew,Jiinew)
            Jiibest(:) = Jiinew(:)
            WRITE(*,1111) ' Updated Jacobian'
         ELSE 
            !We only need to compute the function value.
            call normalization2pt(rho,anew,fnew)
         ENDIF          
         !Compute error and update vectors
         !Compute average absolute error.
         avg_error_new = SUM(Wbecke*abs(fnew))   
         
         IF (avg_error_new >= avg_error_best) THEN
            WRITE(*,1111) 'After pruning step with local trust radius, the error went ***UP*** by ', &
                                avg_error_new-avg_error_best
            !Update the vectors:
            IF (Jvectors > 0) THEN
               !We need to store the data from previous iterations to use in the 
               !limited memory quasi-Newton:
               DO i=1,Jvectors-1
                  a(:,i) = a(:,i+1)
                  f(:,i) = f(:,i+1)
               ENDDO      
               a(:,Jvectors) = anew(:)
               f(:,Jvectors) = fnew(:)
            ENDIF
         ELSE
            WRITE(*,1111) 'Pruning the step with local trust radius decreased the error by ', &
                                avg_error_best-avg_error_new
            !The error improved, as expected.
            avg_error_best = avg_error_new
            list_error(5) = avg_error_best
            IF (Jvectors > 0) THEN
               !We need to store the data from previous iterations to use in the 
               !limited memory quasi-Newton:
               DO i=1,Jvectors-1
                  a(:,i) = a(:,i+1)
                  f(:,i) = f(:,i+1)
               ENDDO      
               a(:,Jvectors) = abest(:)
               f(:,Jvectors) = fbest(:)
            ENDIF
            aold(:) = abest(:)
            fold(:) = fbest(:)
            abest(:) = anew(:)
            fbest(:) = fnew(:)
         ENDIF
      ENDIF
      
      !REJECT MOVES?? ??  REJECT MOVES?? ??  REJECT MOVES?? ??  REJECT MOVES?? ??  REJECT MOVES??  ?? REJECT MOVES??
      
      !COMPUTE NEWMOVE !!  COMPUTE NEWMOVE !!  COMPUTE NEWMOVE !!  COMPUTE NEWMOVE !!  COMPUTE NEWMOVE !!  COMPUTE NEWMOVE
      !Write out some statistics
      IF (debug) THEN
         WRITE(*,*) '          Percentage of',n_grid, ' points that: ' 
         WRITE(*,1112) 'improved','in-trust','de-trust','backward','truncate','tooshort',' toolong','goodbwrd'
         WRITE(*,1113) NINT(COUNT(Limprove)*100.0_dbl/n_grid),                 &
                       NINT(COUNT(Lincrtrust)*100.0_dbl/n_grid),               &
                       NINT(COUNT(Ldecrtrust)*100.0_dbl/n_grid),               &
                       NINT(COUNT(Lbwardstep)*100.0_dbl/n_grid),               &
                       NINT(COUNT(Lcutstep)*100.0_dbl/n_grid),                 &
                       NINT(COUNT(Lshort)*100.0_dbl/n_grid),                   &
                       NINT(COUNT(Llong)*100.0_dbl/n_grid),                    &
                       NINT(COUNT(Lgoodbward)*100.0_dbl/COUNT(Lbwardstep))
         WRITE(*,*) ' '
         WRITE(*,*) ' '
      ENDIF
       
      !Make next step.
      call Broyden_Step(abest,fbest,Jiibest,a,f,anew,fguess,                   &
                       trust_radii,Global_trust,Lcutstep,Lscalestep)
      !call diagNewton_Step_UEG(abest,fbest,Jiibest,trust_radii,                 &
      !                        anew,fguess)

      !COMPUTE NEWMOVE !!  COMPUTE NEWMOVE !!  COMPUTE NEWMOVE !!  COMPUTE NEWMOVE !!  COMPUTE NEWMOVE !!  COMPUTE NEWMOVE
       
      !Reset the trust radii every now and then
      IF (MOD(iter,MAX(5*Jdiag_update,max_iter/3-5)) == 0) THEN
         !Reset the trust radii.
         trust_radii(:) = MAX(ABS(anew(:))/2,normalization_tolerance)
         Global_trust = SUM(ABS(trust_radii))          
         WRITE(*,*) '.........resetting trust radii......'
      ENDIF

      !Increment iteration.
      iter = iter + 1
      !GOODSTEP GOODSTEP GOODSTEP GOODSTEP GOODSTEP GOODSTEP GOODSTEP GOODSTEP GOODSTEP GOODSTEP GOODSTEP 
   ENDIF
   
   !Back to the top to check the new step!

ENDDO

!Return the best value you found:
value = anew

end subroutine Newton2pt


!-------------------------------------------------------------------------------!
!                           Broyden_Step                                        !
!                                                                               !
! This subroutine determines the limited-memory Broyden step for the bad        !
! Broyden method with Jvectors vectors stored.                                  !
!                                                                               !
! Define:                                                                       !
!   f(:) = current function value.                                              !
!   da(:,1:Jvectors) = vector of change in variable values.                     !
!   df(:,1:Jvectors) = vector of change in function values.                     !
!   dfsq(1:Jvectors) = df {dot} df                                              !
!   dy(1:Jvec) = df(:,j) {dot} f                                                !
!                                                                               !
! The step is computed recursively using:                                       !
!  Jinv(k) = (da(k) X dF(k))/dfsq(k) + Jinv(k-1)(Id - (df X df)/dfsq            !
! Where a X b is the exterior product of two vectors.                           !
!                                                                               !
!-------------------------------------------------------------------------------!
! anow -- the current vector of "effective fermi momenta."                      !
! fnow -- the current normalization errors.                                     !
! a -- previous effective fermi momenta.                                        !
! f -- previous normalization errors.                                           !
! Jii -- diagonal of the Jacobian.                                              !
! da = anow - a.                                                                !
! df = fnow - f                                                                 !
! dfsq -- df {dot} df                                                           !
! Nrecur -- the number of vectors available in the recursion.                   !
!                                                                               !
! anew -- the output vector of fermi momenta.                                   !
! fguess -- the guessed error in the step; this is based on the diagonal.       !
! fullstep -- the full Newton step.                                             !
! step -- the Newton step after using the trust radius.                         !
! part1 -- the first part of the Broyden step.  (sum part, above.)              !
! part2 -- the second part of the Broyden step. (product part, above.)          !
! trust_radii -- the vector of trust radii for the steps.                       !
! Global_trust -- the global trust radius.                                      !
! trust_compare -- a number that is compared to the trust radius.               !
!                                                                               !
! debug -- a debug print flag.                                                  !
!                                                                               !
! n_grid -- the number of grid points.                                          !
! Jvectors -- the number of vectors kept to store the Jacobian in limited-memory!
!             Broyden's bad method.                                             !
! Jdiag_update -- the diagonal is updated every this-many iterations.           !
! Jtrust_radius -- controls when and how the local trust radii are updated.     !
!                    0  --  no local trust radius used.                         !
!                   >0   --  local trust radii used to construct the next step. !
!                   <0  --  local trust radii used to reject certain portions   !
!                                    of the current step.                       !
! normalization_tolerance -- the calculations converges when the average        !
!                            absolute error and the maximum error are less than !
!                            this.                                              !
!                                                                               !
! i,j,k,l -- counters.                                                          !
!-------------------------------------------------------------------------------!

subroutine Broyden_Step(anow,fnow,Jii,a,f,anew,fguess,                          &
                        trust_radii,Global_trust,Lcutstep,Lscalestep)

USE kinds
USE grid_vars, ONLY: n_grid
USE KE_vars, ONLY: Jtrust_radius,Jvectors

IMPLICIT NONE

INTEGER :: i,j,k,l,Nrecur

REAL(dbl) :: anow(1:n_grid),anew(1:n_grid),a(1:n_grid,1:Jvectors)
REAL(dbl) :: fnow(1:n_grid),fguess(1:n_grid),f(1:n_grid,1:Jvectors)
REAL(dbl) :: Jii(1:n_grid)
REAL(dbl) :: da(1:n_grid,1:Jvectors),df(1:n_grid,Jvectors),dfsq(1:Jvectors)
REAL(dbl) :: part1(1:n_grid),part2(1:n_grid),step(1:n_grid),fullstep(1:n_grid)

REAL(dbl) :: trust_radii(1:n_grid),Global_trust,trust_compare

LOGICAL :: debug,Lcutstep(1:n_grid),Lscalestep

1114 FORMAT(10X,A,G12.5)
debug = .true.
Lcutstep = .false.

!Compute the differences of the vectors.
FORALL(j=1:Jvectors)
      da(:,j) = anow(:) - a(:,j)
      df(:,j) = fnow(:) - f(:,j)
ENDFORALL
FORALL(j=1:Jvectors)
      dfsq(j) = DOT_PRODUCT(df(:,j),df(:,j))
ENDFORALL

Nrecur = Jvectors + 1

DO j=1,Jvectors
   IF (a(1,j) >= 0) THEN
      !There is real data here. Notice that this depends on the way we
      !initialized a!!.
      Nrecur = j 
      EXIT
   ENDIF
ENDDO

IF (Nrecur > Jvectors) THEN
   WRITE(*,*) '          Diagonal step.'
ENDIF

!For part 2, let's first compute the terms as if the Jacobian was an identity matrix.
part1 = 0.0_dbl
part2(:) = fnow(:)
DO j=Jvectors,Nrecur,-1
   part1(:) = part1(:) + da(:,j)*DOT_PRODUCT(df(:,j),part2(:))/dfsq(j)
   part2(:) = part2(:) - df(:,j)*DOT_PRODUCT(df(:,j),part2(:))/dfsq(j)
ENDDO

!Now use the diagonal element.
WHERE(Jii > EPSILON(Jii))
     part2 = part2/Jii
ELSEWHERE
     !In other cases, we assume that Jii is the identity matrix, but with the sign
     !given to be that of the function.  
     part2 = SIGN(part2,fnow)
ENDWHERE

!We are implicitly assuming that Jii is the identity matrix for the other elements.

!We have construct Jinverse * f; the step is -1 times this:
fullstep = -1*(part1 + part2)

!Now we need to prune according to the trust radius.
trust_compare = SUM(ABS(fullstep))

!Sometimes you wind up with a very short step.  Thatt's problematic.  Don't allow
!steps less than about SQRT(epsilon)
trust_compare = MAX(trust_compare,SQRT(EPSILON(trust_compare)))

IF (trust_compare <= Global_trust) THEN
   !We will accept the step.
   step = fullstep
   Lscalestep = .false.
ELSE
   Lscalestep = .true.
   step = fullstep*Global_trust/trust_compare
   WRITE(*,1114) 'Scaling back whole step by factor of',trust_compare/Global_trust
ENDIF   

!We will prune things back even further if we allow local trust radii to be
!considered.
IF (Jtrust_radius /= 0) THEN
   !Use local trust radii to prune the step.
   WHERE (ABS(step) > trust_radii)
         step = SIGN(trust_radii,step)
         Lcutstep = .true.
   END WHERE
ENDIF

!Finally we want to be absolutely sure that we never let the value of the
!variable become negative.  (a = 0 is absurd!!)
WHERE (step < -1*anow/2)
      step = -1*anow/2
      Lcutstep = .true.
      trust_radii = anow/2
ENDWHERE

!Make the step.    
anew = anow + step
anew = MAX(anew,0.0_dbl)   !Use this just to be absolutely positively sure you never
                           !get a negative number.

!Now we have a problem.  We really need a prediction of the step, but we
!only have the *INVERSE* Jacobian and not the real Jacobian.  We approximate the
!Jacobian by assuming that the inverse is diagonally dominant.  
fguess = fnow + Jii*step

end subroutine Broyden_Step

end module NORMALIZE_MODULE

