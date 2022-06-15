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

subroutine hole_normalize()

USE kinds
USE constants
USE grid_vars, ONLY: n_grid,Wbecke
USE KE_vars, ONLY: a_LDA,rho,KE_conv,D2g_LDA,normalization_tolerance
USE variables_wfn, ONLY: multiplicity

IMPLICIT NONE

INTEGER(istd) :: spin,refpoints
REAL(dbl), ALLOCATABLE :: localKE(:)
LOGICAL       :: debug

debug = .false.

IF (debug) THEN
   ALLOCATE(localKE(1:n_grid))
ENDIF

FORALL(spin=1:2)
      a_LDA(:,1,spin) = (6*Pi**2*rho(:,spin))**(1.0_dbl/3)
ENDFORALL

!If debug print, then print LDA k.e. here.
IF (debug) THEN
   localKE(:) = rho(:,1)*D2g_LDA*a_LDA(:,1,1)**2 + rho(:,2)*D2g_LDA*a_LDA(:,1,2)**2
   WRITE(*,*) 'TF the other way  ', DOT_PRODUCT(Wbecke,localKE)
   WRITE(*,*) 'UEG normalization ', DOT_PRODUCT(Wbecke,localKE) + KE_conv(2)
   WRITE(*,*) ' '
ENDIF

!Do one-reference-point solution.
refpoints = 1

!Initialize to UEG values.
a_LDA(:,2,:) = a_LDA(:,1,:)

!Solve for 1-reference-point hole
!alpha-spin
call Newton(a_LDA(:,2,1),rho(:,1),normalization_tolerance,refpoints)
IF (multiplicity == 1) THEN
   a_LDA(:,2,2) = a_LDA(:,2,1)
ELSE
   !Solve for beta-spin.
   call Newton(a_LDA(:,2,2),rho(:,2),normalization_tolerance,refpoints)
ENDIF

IF (debug) THEN
   localKE(:) = rho(:,1)*D2g_LDA*a_LDA(:,2,1)**2 + rho(:,2)*D2g_LDA*a_LDA(:,2,2)**2
   WRITE(*,*) ' '
   WRITE(*,*) '1-pt. normalization ',DOT_PRODUCT(Wbecke,localKE) + KE_conv(2)
ENDIF

!Do two reference-point solution.
refpoints = 2
!Initialize to 1-reference-point hole
a_LDA(:,3,:) = a_LDA(:,2,:)

!Solve for alpha-spin
call Newton(a_LDA(:,3,1),rho(:,1),normalization_tolerance,refpoints)
IF (multiplicity == 1) THEN
   a_LDA(:,3,2) = a_LDA(:,3,1)
ELSE
   !Solve for beta-spin.
   call Newton(a_LDA(:,3,2),rho(:,2),normalization_tolerance,refpoints)
ENDIF

IF (debug) THEN
   localKE(:) = rho(:,1)*D2g_LDA*a_LDA(:,3,1)**2 + rho(:,2)*D2g_LDA*a_LDA(:,3,2)**2
   WRITE(*,*) ' '
   WRITE(*,*) '2-pt. normalization ', DOT_PRODUCT(Wbecke,localKE) + KE_conv(2)
ENDIF

IF (debug) THEN
   DEALLOCATE(localKE)
ENDIF

end subroutine hole_normalize

!-------------------------------------------------------------------------------!
!                           Newton                                              !
!                                                                               !
! This solver uses Newton's method, with a trust radius, to solve for the       !
! normalization condition in the "exact" case.  This is a "fullgrid" Newtons    !
! method because all of the coefficients a_???(:,3) are solved for at           !
! the same time.  If refpoints = 1, however, then the Jacobian is diagonal and  !
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

subroutine Newton(value,rho,normalization_tolerance,refpoints)

USE kinds
USE grid_vars, ONLY: n_grid,Wbecke

IMPLICIT NONE

REAL(dbl) :: normalization_tolerance,value(1:n_grid)
REAL(dbl) :: max_error,avg_error,norm_error
INTEGER(istd) :: iter, refpoints,i,max_iter
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
trust_radii(:) = ABS(anew(:))/3
iter = 1
max_iter = 100

!Set debug=.true. for debug printing.
debug = .true.

IF (debug) THEN
   WRITE(*,*) ' '
   WRITE(*,*) 'For ', refpoints, ' point normalization.'
   WRITE(*,*) 'Number of electrons, ', DOT_PRODUCT(Wbecke,rho)
   WRITE(*,*) 'Iteration info:'
ENDIF

DO   
   !Default is to update all points
   updated(:) = .true.
   L_lingood(:) = .false.
   L_linbad(:) = .false.

   !Call subroutine to evaluate the error in the normalization of the hole.
   call normalizationUEG(rho,anew,fnew,Jnew,refpoints,updated)
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
      value = anew
      WRITE(*,*) refpoints,'-point hole normalization converged.'
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
         trust_radii = trust_radii/1.4_dbl
         updated = .false.
   ELSEWHERE (trust_ratio > .8_dbl)
         !The Newton step gave an accurate prediction of the step.  Take the step
         !and increase the trust radius.  Compute the Jacobian.  Because every iteration
         !is expensive, it is very bad to reject an iteration.  So we expand the trust radius
         !conservatively but reduce the trust radius aggressively.
         trust_radii =  trust_radii*1.1_dbl
         anow = anew; fnow = fnew; Jnow = Jnew
         L_lingood = .true.         
   ELSEWHERE (trust_ratio > .2_dbl) 
         !The Newton step was not that accurate but it didn't mess us up too badly.  Keep
         !the old trust radius and make the Newton step.
         anow = anew; fnow = fnew; Jnow = Jnew
   ELSEWHERE 
         !The Newton step was quite inaccurate; reduce the trust radius but take the step
         !since at least it got us closer (ever-so-slightly) to the solution.
         trust_radii = trust_radii/1.2_dbl
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
   call diagNewton_Step_UEG(rho,anow,fnow,Jnow,refpoints,updated,trust_radii,  &
                            anew,fpredicted,normalization_tolerance)
   
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
      trust_radii(:) = ABS(anow(:))/10
      WRITE(*,*) '......resetting trust radii.......'
   ENDIF
ENDDO

end subroutine Newton

!-------------------------------------------------------------------------------!
!                           normalizationUEG                                    !
!                                                                               !
! This subroutine evaluates the normalization error and the Jacobian of the     !
! nonlinear equations to be solved, but only the diagonal element of the        !
! Jacobian is used.                                                             !
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
! refpoints -- the number of "reference points" used to evaluate the hole.      !
! gUEG(1:n_grid,1:n_grid) -- the square root of the exchange hole in the        !
!                            uniform electron gas.                              !
! dgUEG(1:n_grid,1:n_grid) -- the derivative of gUEG with respect to the a      !
!                             (the effective kf value).  Used to compute the    !
!                             Jacobian.                                         !
! dij(1:n_grid,1:n_grid) -- the distance between two grid points.               !
! pdist(1:n_grid,1:n_grid) -- the "effective" kF value times the distance       !
!                             between points in the 2-point model.              !
! computeflag -- .true. if this value is computed. .false. otherwise.           !
!-------------------------------------------------------------------------------!
! n_grid -- the number of grid points.                                          !
! Wbecke -- the integration weights.                                            !
! p_mean -- the p-value that specifies the appropriate "generalized mean."      !
!-------------------------------------------------------------------------------!

subroutine normalizationUEG(rho,a,f,Jdiag,refpoints,computeflag)

USE kinds
USE grid_vars, ONLY: n_grid,Wbecke,XYZbecke
USE KE_vars, ONLY: p_mean

IMPLICIT NONE

REAL(dbl) :: rho(1:n_grid),a(1:n_grid),f(1:n_grid),Jdiag(1:n_grid)
REAL(dbl) :: gUEG(1:n_grid,1:n_grid),dgUEG(1:n_grid,1:n_grid)
REAL(dbl) :: dij(1:n_grid,1:n_grid),pdist(1:n_grid,1:n_grid)
INTEGER(istd) :: refpoints,k,l

LOGICAL :: computeflag(1:n_grid),debug

debug =.true.

!Compute the distance between the grid points
FORALL(k=1:n_grid,l=1:n_grid, computeflag(k))
      dij(k,l) = SQRT(DOT_PRODUCT(XYZbecke(1:3,k)-XYZbecke(1:3,l),                  &
                             XYZbecke(1:3,k)-XYZbecke(1:3,l)))
ENDFORALL

!When refpoints = 1, then we do a simple normalization integral with a one-point
!"width" of the hole.  This gives us, as equations to solve,
! 0 = [ INT(rho(r)h(a|r-r'|)dr) - (-1) ] * rho(r')

IF (refpoints == 1) THEN

   !The default value of the g function is 1.0, which occurs at coalescence (dij= 0)
   FORALL(k=1:n_grid,l=1:n_grid, (a(k)*dij(k,l)<=.0000001_dbl .AND. computeflag(k)))
         gUEG(k,l) = 1.0_dbl
         dgUEG(k,l) = 0.0_dbl
   ENDFORALL
   !If the points are not too close together, then compute the term when necessary.
   FORALL(k=1:n_grid,l=1:n_grid, (a(k)*dij(k,l)>.0000001_dbl .AND. computeflag(k)))
         gUEG(k,l) = 3*(sin(a(k)*dij(k,l)) - a(k)*dij(k,l)*cos(a(k)*dij(k,l))) &
                         /(a(k)*dij(k,l))**3 
   ENDFORALL
   FORALL(k=1:n_grid,l=1:n_grid, (a(k)*dij(k,l)>.0000001_dbl .AND. computeflag(k)))
         dgUEG(k,l) = 3*(sin(a(k)*dij(k,l)) - a(k)*dij(k,l)*gUEG(k,l))         &
                        /(a(k)*dij(k,l))**2
   ENDFORALL
   FORALL(k=1:n_grid, computeflag(k))
         f(k) = (DOT_PRODUCT(Wbecke(:),rho(:)*(-1)*gUEG(k,:)**2) + 1)*rho(k)
         Jdiag(k) = (2*DOT_PRODUCT(Wbecke(:),rho(:)*(-1)*gUEG(k,:)*dij(k,:)*dgUEG(k,:)))*rho(k)
   ENDFORALL
   
ELSE
   !Assume that two reference points will be used.  First make the appropriate
   !weighting value.
   FORALL(k=1:n_grid,l=1:n_grid, computeflag(k))
         pdist(k,l) = ((a(k)**p_mean + a(l)**p_mean)/2)**(1.0_dbl/p_mean)*dij(k,l)
   ENDFORALL

   !The default value of the g function is 1.0, which occurs at coalescence (dij= 0)
   FORALL(k=1:n_grid,l=1:n_grid, (pdist(k,l)<=.0000001_dbl .AND. computeflag(k)))
         gUEG(k,l) = 1.0_dbl
         dgUEG(k,l) = 0.0_dbl
   ENDFORALL   
   !If the points are not too close together, use both points.
   FORALL(k=1:n_grid,l=1:n_grid, (pdist(k,l)>.0000001_dbl .AND. computeflag(k)))
         gUEG(k,l) = 3 * (sin(pdist(k,l))                                 &
                           - pdist(k,l)*cos(pdist(k,l))) &
                         /(pdist(k,l))**3 
   ENDFORALL
   FORALL(k=1:n_grid,l=1:n_grid, (pdist(k,l)>.0000001_dbl .AND. computeflag(k)))
         dgUEG(k,l) = 3*(sin(pdist(k,l)) - pdist(k,l)*gUEG(k,l))          &
                        /(pdist(k,l))**2
   ENDFORALL
   FORALL(k=1:n_grid, computeflag(k))
         f(k) = (DOT_PRODUCT(Wbecke(:),rho(:)*(-1)*gUEG(k,:)**2) + 1)*rho(k)
         Jdiag(k) = (DOT_PRODUCT(Wbecke(:),rho(:)*(-1)*gUEG(k,:)*dgUEG(k,:)      &
                                          *a(k)**(p_mean-1)*pdist(k,:)        &
                                          /((a(k)**p_mean+a(:)**p_mean)/2)))*rho(k)
    ENDFORALL
ENDIF

end subroutine normalizationUEG

!-------------------------------------------------------------------------------!
!                           diagNewton_Step_UEG                                 !
!                                                                               !
! This subroutine determines the Newton step.                                   !
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
! updated -- keeps track of which steps were accepted and which were rejected.  !
! rho -- the spin-density of the channel currently under consideration.         !

! step -- the step that will be taken.  At the beginning it is the full step,   !
!         then it is scaled back.                                               !
! compute_flag(1:ngrid) -- the elements of the Jacobian and the matrix to be    !
!                          updated.                                             !
! n_grid -- the number of grid points.                                          !
! fullstep -- the full Newton step.                                             !
!-------------------------------------------------------------------------------!

subroutine diagNewton_Step_UEG(rho,anow,fnow,Jii,refpoints,updated,trust_radii, &
                           anew,fpredicted,normalization_tolerance)

USE kinds
USE grid_vars, ONLY: n_grid

IMPLICIT NONE

REAL(dbl) :: rho(1:n_grid),anow(1:n_grid),fnow(1:n_grid),Jii(1:n_grid)
REAL(dbl) :: trust_radii(1:n_grid),anew(1:n_grid),fpredicted(1:n_grid)
REAL(dbl) :: step(1:n_grid),fullstep(1:n_grid),normalization_tolerance
INTEGER(istd) :: refpoints,j
LOGICAL :: updated(1:n_grid),compute_flag(1:n_grid),debug

debug = .true.

!First compute the Newton step.
IF (refpoints == 1) THEN
   !In this case the Newton step is very straightforward, because the different
   !"a" values are uncoupled
   WHERE(Jii > 0)
        !Compute Newton step.
        fullstep = -1*fnow/Jii
   ELSEWHERE
        !Either this point has converged or the Jacobian is zero.
        fullstep = 0.0_dbl
   ENDWHERE
       
ELSE
   !The problem is that some (hopefully most) a values will have been updated,
   !while others will revert back to their previous values.  If almost all of the
   !a values were updated, then it is reasonable to *only* compute the values 
   !that are *not* updated.  So 
   IF (COUNT(updated) > .95_dbl*n_grid) THEN
      !More than 95% of the points were updated. Update the Jacobian only for
      !the points that were *not* updated.
      compute_flag(:) = .not. updated(:)  
      call normalizationUEG(rho,anow,fnow,Jii,refpoints,compute_flag)
   ELSE
      !A lot of points were not updated.  Update all points.
      compute_flag(:) = .true.
      call normalizationUEG(rho,anow,fnow,Jii,refpoints,compute_flag)
   ENDIF
   
   WHERE(Jii > 0)
        !Compute Newton step.
        fullstep = -1*fnow/Jii
   ELSEWHERE
        !Either this point has converged or the Jacobian is zero.
        fullstep = 0.0_dbl
   ENDWHERE
   
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
!The trust radius must be reduced b/c otherwise one "falsely accepts" a short
!step as very reliable and increases the trust radius.
FORALL (j=1:n_grid, step(j) < -1*anow(j)/2)
       step(j) = -1*anow(j)/2
       trust_radii(j) = anow(j)/3    !Can't use anow(j)/2 without risking a limit cycle.
ENDFORALL
    
anew = anow + step
fpredicted = fnow + Jii*step

end subroutine diagNewton_Step_UEG

end module NORMALIZE_MODULE
