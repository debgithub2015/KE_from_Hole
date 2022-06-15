!-------------------------------------------------------------------------------!
!                               holesubroutines (unscaled)                      !
!                                                                               !
! This module contains the subroutines for computing the normalization factors  !
! for the exchange-correlation hole and their Jacobians.                        !
!                                                                               !
! This is the only module that needs to be replaced if you want to do the       !
! normalization for a different case.                                           !
!                                                                               !
! normalization1pt -- computes the normalization error and the Jacobian for the !
!                     1-point normalization case.                               !
! normalization2pt -- computes the normalization error for the 2-point          !
!                     normalization case.                                       !
! normalization2ptJacobian -- computes the normalization error and the          !
!                     Jacobian for the 2-pt normalization case.                 !
!                                                                               !
! The function value that we use to test the normalization is:                  !
! f(r') =  (<rho(r)*h(r,r')> - (-1))  * rho(r')**(1/3)                          !
!       =  (<rho(r)*h(r,r')> + 1) * rho(r')**(1/3)                              !
!                                                                               !
! The Jacobian is df(r')/da(r') where a(r') is the effective Fermi momentum at  !
! r'.                                                                           !
!                                                                               !
! When f(r') > 0, this means that the hole is normalized to a number that is    !
!                 greater than -1.  In this case we need the hole to be bigger, !
!                 and we expect that the change in a should be negative (da<0). !
! When f(r') < 0, this means hat the hole is normalized to a number that is     !
!                 less than -1, and we need the hole to be smaller, so the      !
!                 change in a should be positive. (da > 0)                      !
! Since the Newton step, in the diagonal case, is                               !
!         f(i) + J(i,i)*da(i) = 0  ==>  da(i) = -f(i)/Jii(i)                    !
!                                  ==> J(i,i) = -f(i)/da(i)                     !
! we expect that J(i,i) is positive.                                            !
!-------------------------------------------------------------------------------!

MODULE holesubroutines

CONTAINS

!-------------------------------------------------------------------------------!
!                           normalization1pt                                    !
!                                                                               !
! This subroutine evaluates the normalization error and the Jacobian of the     !
! nonlinear equations to be solved for the 1-point normalization case.          !
!-------------------------------------------------------------------------------!
!                               DICTIONARY OF VARIABLES                         !
!                                                                               !
! a -- effective kf value.                                                      !
! f -- error in normalization of the hole, = [<rho(r)*h(r,r')> - (-1)]rho(r')   !
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

subroutine normalization1pt(rho,a,f,Jdiag)

USE kinds
USE grid_vars, ONLY: n_grid,Wbecke,XYZbecke

IMPLICIT NONE

REAL(dbl) :: rho(1:n_grid),a(1:n_grid),f(1:n_grid),Jdiag(1:n_grid)
REAL(dbl) :: gUEG(1:n_grid,1:n_grid),dgUEG(1:n_grid,1:n_grid)
REAL(dbl) :: dij(1:n_grid,1:n_grid),pdist(1:n_grid,1:n_grid)
INTEGER(istd) :: k,l

LOGICAL :: debug

debug = .false.

!Compute the distance between the grid points
FORALL(k=1:n_grid,l=1:n_grid)
      dij(k,l) = SQRT(DOT_PRODUCT(XYZbecke(1:3,k)-XYZbecke(1:3,l),                  &
                             XYZbecke(1:3,k)-XYZbecke(1:3,l)))
ENDFORALL

!only one reference function is used.
FORALL(k=1:n_grid,l=1:n_grid)
      pdist(k,l) = a(k)*dij(k,l)
ENDFORALL

!We do a simple normalization integral with a one-point
!"width" of the hole.  This gives us, as equations to solve,
! 0 = [ INT(rho(r)h(a|r-r'|)dr) - (-1) ] * rho(r')

!The default value of the g function is 1.0, which occurs at coalescence (dij= 0)
!This cutoff comes from my experimentation.  It seems that round-off is *killing*
!us here.  Once pdist < epsilon**(1/4) then we only have that g = 1 to within 
!SQRT(epsilon) and dg = 0 to within (epsilon)**(1/4).  The following Taylor 
!Series seem more accurate up to about .03.
FORALL(k=1:n_grid,l=1:n_grid, (pdist(k,l) <= .03_dbl))
      gUEG(k,l) = 1.0_dbl - pdist(k,l)**2/10 + pdist(k,l)**4/280 - pdist(k,l)**6/15120
      dgUEG(k,l) = 0.0_dbl - pdist(k,l)/5 + pdist(k,l)**3/70 - pdist(k,l)**5/2520      &
                   + pdist(k,l)**7/166320
ENDFORALL   

!If the points are not too close together, use both points.
FORALL(k=1:n_grid,l=1:n_grid, (pdist(k,l) > .03_dbl))
      gUEG(k,l) = 3 * (sin(pdist(k,l))/pdist(k,l)**3 - cos(pdist(k,l))/pdist(k,l)**2) 
ENDFORALL
FORALL(k=1:n_grid,l=1:n_grid, (pdist(k,l) > .03_dbl))
      dgUEG(k,l) = 3*(sin(pdist(k,l))/pdist(k,l)**2 - gUEG(k,l)/pdist(k,l))
ENDFORALL

FORALL(k=1:n_grid)
      f(k) = (DOT_PRODUCT(Wbecke(:),rho(:)*(-1)*gUEG(k,:)**2) + 1)*rho(k)*Wbecke(k)
      Jdiag(k) = (2*DOT_PRODUCT(Wbecke(:),rho(:)*(-1)*gUEG(k,:)*dij(k,:)*dgUEG(k,:)))*rho(k)*Wbecke(k)
ENDFORALL

IF (debug) THEN
   DO k=1,n_grid
      IF (Jdiag(k) < 0.0_dbl) THEN
         WRITE(*,*) 'Found negative value of Jacobian.  This is unexpected.'
         WRITE(*,*) rho(k),a(k),Jdiag(k)
      ENDIF
   ENDDO
ENDIF

end subroutine normalization1pt

!-------------------------------------------------------------------------------!
!                           normalization2pt                                    !
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
! f -- error in normalization of the hole, = [<rho(r)*h(r,r')> - (-1)]rho(r')   !
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

subroutine normalization2pt(rho,a,f)

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
!This cutoff comes from my experimentation.  It seems that round-off is *killing*
!us here.  Once pdist < epsilon**(1/4) then we only have that g = 1 to within 
!SQRT(epsilon) and dg = 0 to within (epsilon)**(1/4).  The following Taylor 
!Series seem more accurate up to about .03.
FORALL(k=1:n_grid,l=1:n_grid, (pdist(k,l) <= .03_dbl))
      gUEG(k,l) = 1.0_dbl - pdist(k,l)**2/10 + pdist(k,l)**4/280 - pdist(k,l)**6/15120
ENDFORALL   

!If the points are not too close together, use both points.
FORALL(k=1:n_grid,l=1:n_grid, (pdist(k,l) > .03_dbl))
      gUEG(k,l) = 3 * (sin(pdist(k,l))/pdist(k,l)**3 - cos(pdist(k,l))/pdist(k,l)**2) 
ENDFORALL

FORALL(k=1:n_grid)
         f(k) = (DOT_PRODUCT(Wbecke(:),rho(:)*(-1.0_dbl)*gUEG(k,:)**2) + 1.0_dbl)*rho(k)*Wbecke(k)
ENDFORALL

end subroutine normalization2pt



!-------------------------------------------------------------------------------!
!                           normalization2ptJacobian                            !
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
! f -- error in normalization of the hole, = [<rho(r)*h(r,r')> - (-1)]rho(r')   !
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

subroutine normalization2ptJacobian(rho,a,f,Jdiag)

USE kinds
USE grid_vars, ONLY: n_grid,Wbecke,XYZbecke
USE KE_vars, ONLY: p_mean

IMPLICIT NONE

REAL(dbl) :: rho(1:n_grid),a(1:n_grid),f(1:n_grid),Jdiag(1:n_grid)
REAL(dbl) :: gUEG(1:n_grid,1:n_grid),dgUEG(1:n_grid,1:n_grid)
REAL(dbl) :: dij(1:n_grid,1:n_grid),pdist(1:n_grid,1:n_grid)
INTEGER(istd) :: k,l

LOGICAL :: debug

debug = .false.

!Compute the distance between the grid points
FORALL(k=1:n_grid,l=1:n_grid)
      dij(k,l) = SQRT(DOT_PRODUCT(XYZbecke(1:3,k)-XYZbecke(1:3,l),              &
                             XYZbecke(1:3,k)-XYZbecke(1:3,l)))
ENDFORALL

!Two reference points will be used and the diagonal of the Jacobian will be
!generated.
FORALL(k=1:n_grid,l=1:n_grid)
      pdist(k,l) = ((a(k)**p_mean + a(l)**p_mean)/2)**(1.0_dbl/p_mean)*dij(k,l)
ENDFORALL

DO k=1,n_grid
   DO l=1,n_grid
      IF (pdist(k,l) < 0) THEN
         WRITE(*,*) k,l,dij(k,l),a(k),a(l),p_mean
      ENDIF
   ENDDO
ENDDO

!The default value of the g function is 1.0, which occurs at coalescence (dij= 0)
!This cutoff comes from my experimentation.  It seems that round-off is *killing*
!us here.  Once pdist < epsilon**(1/4) then we only have that g = 1 to within 
!SQRT(epsilon) and dg = 0 to within (epsilon)**(1/4).  The following Taylor 
!Series seem more accurate up to about .03.
FORALL(k=1:n_grid,l=1:n_grid, (pdist(k,l) <= .03_dbl))
      gUEG(k,l) = 1.0_dbl - pdist(k,l)**2/10 + pdist(k,l)**4/280 - pdist(k,l)**6/15120
      dgUEG(k,l) = 0.0_dbl - pdist(k,l)/5 + pdist(k,l)**3/70 - pdist(k,l)**5/2520  &
                    + pdist(k,l)**7/166320
ENDFORALL   

!If the points are not too close together, use both points.
FORALL(k=1:n_grid,l=1:n_grid, (pdist(k,l) > .03_dbl))
      gUEG(k,l) = 3 * (sin(pdist(k,l))/pdist(k,l)**3 - cos(pdist(k,l))/pdist(k,l)**2) 
ENDFORALL

FORALL(k=1:n_grid,l=1:n_grid, (pdist(k,l) > .03_dbl))
      dgUEG(k,l) = 3*(sin(pdist(k,l))/pdist(k,l)**2 - gUEG(k,l)/pdist(k,l))
ENDFORALL

f(:) = 0.0_dbl
FORALL(k=1:n_grid)
      f(k) = (DOT_PRODUCT(Wbecke(:),rho(:)*(-1)*gUEG(k,:)**2) + 1)*rho(k)*Wbecke(k)
ENDFORALL

!Initialize Jacobian to zero.
Jdiag(:) = 0.0_dbl  
FORALL(k=1:n_grid)
      !The +EPSILON(a(k))**2 terms are a way to avoid dividing by zero, and should cause
      !negligible error, since we won't ever converge the a's that accurately anyway.
         Jdiag(k) = -4*DOT_PRODUCT(Wbecke(:),rho(:)*gUEG(k,:)*dgUEG(k,:)                     &
                                        *(a(k)+EPSILON(a(k))**2)**(p_mean-1)                        &
                                        *pdist(k,:)/(a(k)**p_mean+a(:)**p_mean+EPSILON(a(k))**2))*rho(k)*Wbecke(k)
ENDFORALL

IF (debug) THEN
   DO k=1,n_grid
      IF (.not. (Jdiag(k) > 0)) THEN
         WRITE(10,*) 'density, eff. Fermi momentum, Jacobian diagonal, fvalue',rho(k),a(k),Jdiag(k),f(k)
      ELSE IF (Jdiag(k) > 1/EPSILON(Jdiag(k)) .or. f(k) > 1/EPSILON(f(k))) THEN
         WRITE(10,*) 'density, eff. Fermi momentum, Jacobian diagonal, fvalue',rho(k),a(k),Jdiag(k),f(k)
      ELSE IF (Jdiag(k) < 0) THEN   
         WRITE(10,*) 'density, eff. Fermi momentum, Jacobian diagonal, fvalue',rho(k),a(k),Jdiag(k),f(k)
      ELSE
         CONTINUE
      ENDIF
   ENDDO
ENDIF

!The Jacobian should never be negative.  Get rid of that problem.  Zero values are
!also problematic because they can lead to gigantic steps.
!FORALL(k=1:n_grid, (.not. (Jdiag(k) >= 0)))
!      Jdiag(k) = 0.0_dbl
!ENDFORALL

end subroutine normalization2ptJacobian




END MODULE holesubroutines
