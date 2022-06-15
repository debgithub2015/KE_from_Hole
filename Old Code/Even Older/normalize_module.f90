!-------------------------------------------------------------------------------!
!                               NORMALIZE_MODULE                                !
!                                                                               !
! This module contains the subroutines for computing the normalization factors  !
! for the exchange-correlation hole.                                            !
!-------------------------------------------------------------------------------!

MODULE NORMALIZE_MODULE

USE kinds; USE constants; USE vars_density; USE vars_grid

CONTAINS

!-------------------------------------------------------------------------------!
! hole_normalize -- driver routine.                                             !
! UEG_normalization -- computes the normalization based on the uniform electron !
!                      gas model.                                               !
! Simple_normalization -- a 1D-Newton's method solver for the normalization     !
!                         constant.                                             !
! Exact_normalization -- a 3D-Newton's method solver for the normalization      !
!                        constant.                                              !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                           hole_normalize                                      !
! driver routine for hole normalization.                                        !
!-------------------------------------------------------------------------------!

subroutine hole_normalize

call UEG_normalization

call Simple_normalization

call Exact_normalization

end subroutine hole_normalize

!-------------------------------------------------------------------------------!
!                               UEG_normalization                               !
!                                                                               !
! This subroutine computes the "effective kF" value based on the uniform        !
! electron gas model for the energy.                                            !
!-------------------------------------------------------------------------------!
!                     DICTIONARY OF LOCAL VARIABLES                             !
!                                                                               !
! rho_one3rd(1:ngrid) -- the density to the one-third power.                    !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_GRID                                !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_grid module!!!                                                  ---- !
! ngrid -- number of points on the grid, including the points used to compute   !
!          the derivatives.                                                     !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_density                             !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_density module!!!                                               ---- !
! rho(:) -- the density at a point                                              !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_KE                                  !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_KE module!!!                                                    ---- !
! a_LDA(1:ngrid,1:3) -- the alpha values used to control the curvature of the   !
!                       hole at same-spin electron coalescence.                 !
!                   1.  From uniform electron gas.                              !
!                   2.  Based on one-point normalization.                       !
!                   3.  Based on exact normalization.                           !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!

subroutine UEG_normalization

IMPLICIT NONE

REAL(dbl) :: rho_one3rd(1:ngrid)

rho_one3rd = rho(:)**(1_dbl/3)

a_LDA(:,1) = (6*Pi)**(1_dbl/3)*rho_one3rd(:)

end subroutine UEG_normalization

!-------------------------------------------------------------------------------!
!                       simple_normalization                                    !
!                                                                               !
! This subroutine performs the "simple" approximate normalization of the        !
! of the exchange hole.  This is based on solving the equation:                 !
!         -1 = INT (rho(r) hxc(kf;r,r') dr) for kf(r').                         !
! which on the grid gives:                                                      !
!        SUM(i) rho(i) hxc(i,j,kf) w(i) + 1 = 0 solved for kf(j)                !
! This equation can be solved by a one-dimensional Newton's method with a trust !
! radius that prevents kf from increase/decreasing more than a set amount. The  !
! uniform electron gas normalization is used as an initial guess.               !
!-------------------------------------------------------------------------------!
!                           DICTIONARY OF LOCAL VARIABLES                       !
!                                                                               !
! j -- counter                                                                  !
! guess -- the guessed value for the "effective kf" of the exchange hole.       !
! answer -- the value for the "effective kf" obtained using Newton's method.    !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_GRID                                !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_grid module!!!                                                  ---- !
! ngrid -- number of points on the grid, including the points used to compute   !
!          the derivatives.                                                     !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_KE                                  !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_KE module!!!                                                    ---- !
! a_LDA(1:ngrid,1:3) -- the alpha values used to control the curvature of the   !
!                       hole at same-spin electron coalescence.                 !
!                   1.  From uniform electron gas.                              !
!                   2.  Based on one-point normalization.                       !
!                   3.  Based on exact normalization.                           !
!-------------------------------------------------------------------------------!

subroutine simple_normalization

IMPLICIT NONE

INTEGER :: j
REAL(dbl) :: normalization_tolerance,guess,answer
CHARACTER(len=4) :: holetype

!The normalization tolerance could be sensibly updated.
normalization_tolerance = 1.0D-5

!Solve the equation for each grid point.  There is probably a way to do every 
!grid point at the same time using some fancy FORALL construct, but since the
!"good normalization" routine will be so much slower, I will not worry about 
!that
DO j = 1,ngrid
   !Solve for LDA hole
   guess= a_LDA(j,1)
   holetype = 'LDAp'
   call Newton1D(guess,answer,j,normalization_tolerance)
   a_LDA(j,2) = answer
ENDDO

end subroutine simple_normalization

!-------------------------------------------------------------------------------!
!                           exact_normalization                                 !
!                                                                               !
! This subroutine uses Newton's method to compute the *exact* two-point         !
! normalization factor for the exchange-correlation hole.  The "simple"         !
! one-point normalization condition is used as an initial guess.                !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                           DICTIONARY OF LOCAL VARIABLES                       !
!                                                                               !
! j -- counter                                                                  !
! guess -- the guessed value for the "effective kf" of the exchange hole.       !
! answer -- the value for the "effective kf" obtained using Newton's method.    !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_GRID                                !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_grid module!!!                                                  ---- !
! ngrid -- number of points on the grid, including the points used to compute   !
!          the derivatives.                                                     !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_KE                                  !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_KE module!!!                                                    ---- !
! a_LDA(1:ngrid,1:3) -- the alpha values used to control the curvature of the   !
!                       hole at same-spin electron coalescence.                 !
!                   1.  From uniform electron gas.                              !
!                   2.  Based on one-point normalization.                       !
!                   3.  Based on exact normalization.                           !
!-------------------------------------------------------------------------------!

subroutine exact_normalization()

IMPLICIT NONE

REAL(dbl) :: normalization_tolerance,guess(1:ngrid),answer(1:ngrid)
CHARACTER(len=4) :: holetype

!The normalization tolerance could be sensibly updated.
normalization_tolerance = 1.0D-4

!Solve for LDA hole
guess= a_LDA(:,2)
holetype = 'LDAp'
call Newton_fullgrid(guess,answer,normalization_tolerance)
a_LDA(:,3) = answer
   
end subroutine exact_normalization

!-------------------------------------------------------------------------------!
!                               Newton1D                                        !
!                                                                               !
! This is a one-dimensional Newton solver using a trust radius to scale back    !
! questionable steps and reject steps that increase the value of the objective  !
! function.                                                                     !
!-------------------------------------------------------------------------------!
! j_pt -- the point at which the hole normalization is being determined.        !
! iter -- number of iterations.                                                 !
! anow -- the current value for the alpha constant.                             !
! fnow -- the current deviation of the hole's normalization from -1             !
!       f = SUM(i) rho(i) hxc(i,j,kf) w(i) + 1                                  !
! dfnow -- the current gradient value.                                          !
! trust_radius -- the trust radius of the current iterate.                      !
! fullstep -- the full step that would be taken according to Newton's method.   !
! fpredicted -- the predicted value of f according to Newton's method.          !
! anew -- newton value for the next point.                                      !
! fnew -- newton value for the function evaluation.                             !
! dfnew -- newton value for the gradient.                                       !
! normalization_tolerance -- the tolerance factor for the normalization.  This  !
!                            MIGHT be adjusted by the user to make the code     !
!                            faster but here I will choose something rather     !
!                            tight for testing purposes.                        !
! trust_ratio -- the ratio that is used to control the trust radius.  If the    !
!                ratio between the "expected improvement" and the "actual       !
!                improvement is close to zero, we decrease the trust radius.    !
!                If it is close to one, we increase the trust radius.  Otherwise!
!                we leave the trust radius unchanged.                           !
!-------------------------------------------------------------------------------!

subroutine Newton1D(guess,answer,j_pt,normalization_tolerance)

IMPLICIT NONE

INTEGER :: j_pt
REAL(dbl) :: guess, answer, normalization_tolerance, trust_ratio

INTEGER :: iter
REAL(dbl) :: anow,fnow,dfnow,trust_radius,fullstep,fpredicted
REAL(dbl) :: fnew,anew,dfnew

!Initialize to guessed "effective kf" value.  Choose predicted function value
!to be zero, as this is the extreme case and is unlikely to be attained.
anew = guess
fpredicted = 0.0_dbl 
fnow = 100_dbl         !essentially ensures that we accept the first step!
trust_radius = anew/2  !this ensures that the first step does not change the "a"
                       !value too much.

DO   
   !Call subroutine to evaluate the function and the gradient.
   call normalization_evaluate(anew,fnew,dfnew,j_pt)
     
   !End if converged.
   IF (abs(fnew) < normalization_tolerance) THEN
      answer = anew
      EXIT
   ENDIF
   
   trust_ratio = min(abs(fpredicted-fnow),abs(fnew-fnow))     &
                 /max(abs(fpredicted-fnow),abs(fnew-fnow))
     
   !If not converged....
   !Reject move if further from solution than before
   IF (abs(fnew) > abs(fnow)) THEN
      !do *NOT* update guess; try again with smaller step.
      !update trust radius to force smaller step (eventually)
      trust_radius = abs(anow-anew)/2
   ELSE IF (trust_ratio > .8_dbl) THEN
      !The Newton step gave an accurate prediction of the step.  Take the step
      !and increase the trust radius.  However, never let the trust radius exceed .9
      !times the value of a, so that you can never get absurd negative values of a.
      trust_radius = min(2*trust_radius,9*anew/10)
      anow = anew; fnow = fnew; dfnow = dfnew
   ELSE IF (trust_ratio > .2_dbl) THEN
      !The Newton step was not that accurate but it didn't mess us up too badly.  Keep
      !the old trust radius and make the Newton step.
      anow = anew; fnow = fnew; dfnow = dfnew
   ELSE 
      !The Newton step was quite inaccurate; reduce the trust radius but take the step
      !since at least it got us closer (ever-so-slightly) to the solution.
      trust_radius = MIN(trust_radius/2,9*anew/10)
      anow = anew; fnow = fnew; dfnow = dfnew
   ENDIF
        
   !Compute the full Newton step.
   fullstep = -1_dbl*fnow/dfnow        
        
   !Scale back the Newton step to agree with the trust radius if needed.
   IF (abs(fullstep) > trust_radius) THEN
      fullstep = SIGN(trust_radius,fullstep)
      fpredicted = fnow + dfnow*fullstep
   ELSE
      fpredicted = 0_dbl
   ENDIF
 
   !Compute Newton guess:
   anew = anow + fullstep

   iter = iter + 1
   IF (iter > 500) THEN
     WRITE(*,*) 'Program filed in solving 1-D Newton method for the so-called'
     WRITE(*,*) 'simple approximate normalization scheme.  Consider printing'
     WRITE(*,*) 'debug info.'
     WRITE(101,*) 'Program filed in solving 1-D Newton method for the so-called'
     WRITE(101,*) 'simple approximate normalization scheme.  Consider printing'
     WRITE(101,*) 'debug info.'
   STOP
   ENDIF
ENDDO

end subroutine Newton1D

!-------------------------------------------------------------------------------!
!                           Newton_fullgrid                                     !
!                                                                               !
! This solver uses Newton's method, with a trust radius, to solve for the       !
! normalization condition in the "exact" case.  This is a "fullgrid" Newtons    !
! method because all of the coefficients a_????(:,3) have to be solved for at   !
! the same time, because they are coupled together in this formulation.         !
!                                                                               !
!-------------------------------------------------------------------------------!
! anow -- the current vector of "effective fermi momenta."                      !
! fnow -- the current normalization errors.                                     !
! df -- the Jacobian matrix, dfnow/danow.                                       !
! trust_radius -- the trust radius for the steps.                               !
! fullstep -- the vector of "steps" for the next trial vector of effective      !
!             Fermi momenta.                                                    !
! fpredicted -- the vector of predicted values for the normalization errors.    !
! fnew -- the vector of new function values.                                    !
! anew -- the vector of new Fermi momenta (the guess).                          !
! trust_ratio -- the ratio that is used to control the trust radius.  If the    !
!                ratio between the "expected improvement" and the "actual       !
!                improvement is close to zero, we decrease the trust radius.    !
!                If it is close to one, we increase the trust radius.  Otherwise!
!                we leave the trust radius unchanged.                           !
!-------------------------------------------------------------------------------!

subroutine Newton_fullgrid(guess,answer,normalization_tolerance)

IMPLICIT NONE

REAL(dbl) :: normalization_tolerance,guess(1:ngrid),answer(1:ngrid)
INTEGER :: iter
REAL(dbl) :: anow(1:ngrid),fnow(1:ngrid),df(1:ngrid,1:ngrid)
REAL(dbl) :: trust_radius,fullstep(1:ngrid),fpredicted(1:ngrid)
REAL(dbl) :: fnew(1:ngrid),anew(1:ngrid),trust_ratio

!Initialize to guessed "effective kf" value.  Choose predicted function value
!to be zero, as this is the extreme case and is unlikely to be attained.
anew = guess
fpredicted = 0.0_dbl 
fnow = 100_dbl
trust_radius = SUM(ABS(anew))/(2*ngrid)

DO   
   !Call subroutine to evaluate the function and the gradient.
   call normalization_fullgrid(anew,fnew)
     
   !End if converged.
   IF (MAXVAL(abs(fnew)) < normalization_tolerance) THEN
      answer = anew
      EXIT
   ENDIF
   
   trust_ratio = min(SUM(abs(fpredicted-fnow)),SUM(abs(fnew-fnow)))          &
                 /max(SUM(abs(fpredicted-fnow)),SUM(abs(fnew-fnow)))
     
   !If not converged....
   !Reject move if further from solution than before
   IF (SUM(abs(fnew)) > SUM(abs(fnow))) THEN
      !do *NOT* update guess; try again with smaller step.
      !update trust radius to force smaller step (eventually)
      trust_radius = trust_radius/2
   ELSE IF (trust_ratio > .8_dbl) THEN
      !The Newton step gave an accurate prediction of the step.  Take the step
      !and increase the trust radius.  Compute the Jacobian.  Because every iteration
      !is expensive, it is very bad to reject an iteration.  So we expand the trust radius
      !conservatively but reduce the trust radius aggressively.
      call normalization_Jacobian(df)
      trust_radius = 1.5*trust_radius
      anow = anew; fnow = fnew
   ELSE IF (trust_ratio > .2_dbl) THEN
      !The Newton step was not that accurate but it didn't mess us up too badly.  Keep
      !the old trust radius and make the Newton step.
      call normalization_Jacobian(df)
      anow = anew; fnow = fnew
   ELSE 
      !The Newton step was quite inaccurate; reduce the trust radius but take the step
      !since at least it got us closer (ever-so-slightly) to the solution.
      trust_radius = trust_radius/2
      call normalization_Jacobian(df)
      anow = anew; fnow = fnew
   ENDIF
        
   !Compute the full Newton step.
   call Newton_LinEq_Solve(df,fnow,fullstep)
   
   !Scale back the Newton step to agree with the trust radius if needed.
   FORALL (j=1:ngrid, ABS(fullstep(j)) > trust_radius)
       fullstep(j) = SIGN(trust_radius,fullstep(j))
   ENDFORALL
   
   !This forall loop ensures that we never take a step that reduces the value
   !of a by more than a factor of 1/2.  This is important because a<=0 is
   !absurd.
   FORALL (j=1:ngrid, fullstep(j) < 0)
       fullstep(j) = MAX(-1*anow(j)/2,fullstep(j))
   ENDFORALL
    
   !Compute Newton guess.  I am not 100% sure this matrix multiplication is in the
   !right order....
   anew = anow + fullstep
   fpredicted = fnow + MATMUL(df,fullstep)

   iter = iter + 1
   IF (iter > 1000) THEN
      WRITE(*,*) 'Program filed in solving full-grid Newton method for the so-called'
      WRITE(*,*) 'exact normalization scheme.  Consider printing'
      WRITE(*,*) 'debug info.'
      WRITE(101,*) 'Program filed in solving full-grid Newton method for the so-called'
      WRITE(101,*) 'exact normalization scheme.  Consider printing'
      WRITE(101,*) 'debug info.'
      STOP
   ENDIF
ENDDO

!-------------------------------------------------------------------------------!
!                           normalization_evaluate                              !
!                                                                               !
! This subroutine evaluates the normalization error and the gradient of the     !
! normalization error for a given model of the exchange-correlation hole.       !
! It uses a *simple* one-point normalization, rather than the more sophisticated!
! symmetric two-point normalization factor.                                     !
!-------------------------------------------------------------------------------!

subroutine normalization_evaluate(anew,fnew,dfnew,j_pt)

??????

end subroutine normalization_evaluate

!-------------------------------------------------------------------------------!
!                           normalization_fullgrid                              !
!                                                                               !
! This subroutine evaluates the normalization error in the two-point symmetric  !
! hole normalization convention.                                                !
!-------------------------------------------------------------------------------!

subroutine normalization_fullgrid(anew,fnew)

??????

end subroutine full_normalization_evaluate


!-------------------------------------------------------------------------------!
!                           normalization_Jacobian                              !
!                                                                               !
! This subroutine evaluates the Jacobian for the two-point normalization        !
! convention.                                                                   !
!-------------------------------------------------------------------------------!

subroutine normalization_Jacobian(df)

??????

end subroutine normalization_Jacobian


subroutine Newton_LinEq_Solve(df,fnow,fullstep)

???????

end subroutine Newton_LinEq_Solve

!-------------------------------------------------------------------------------!
!                               DM_finder                                       !
!                                                                               !
! This subroutine takes the input exchange hole and uses it to find the density !
! matrix for a given spin channel.  All of the eigenvectors of the matrix       !
!   SQRT(-rho(r)rho(r')h(r,r'))                                                 !
! are computed and the N largest ones are set to one, while the rest are set to !
! zero.  This gives a highly approximate density matrix that can be numerically !
! differentiated to obtain the kinetic energy.                                  !
!    It would be faster to write this program so that only the N largest        !
! eigenvalues are computed but this is a little too involved for us right now.  !
! Alternatively, one could project onto a basis set and proceed analytically,   !
! in the normal quantum-chemistry way.                                          !
!-------------------------------------------------------------------------------!






