!-------------------------------------------------------------------------------!
!                               holesubroutines                                 !
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
! f(r') = rho(r') * (<rho(r)*h(r,r')> - (-1))                                   !
!       = rho(r') * (<rho(r)*h(r,r')> + 1)                                      !
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
! gUEG(1:n_grid) -- the square root of the exchange hole in the                 !
!                            uniform electron gas.                              !
! dgUEG(1:n_grid) -- the derivative of gUEG with respect to the a               !
!                             (the effective kf value).  Used to compute the    !
!                             Jacobian.                                         !
! pdist(1:n_grid) -- the distance in the weighted density approximation.        !
! dR(1:n_grid) -- the distance between two grid points.                         !
! dij(1:n_grid) -- the distance between two points.                             !
!-------------------------------------------------------------------------------!
! n_grid -- the number of grid points.                                          !
! Wbecke -- the integration weights.                                            !
!-------------------------------------------------------------------------------!

subroutine normalization1pt(rho,a,f,Jdiag)

USE kinds
USE grid_vars, ONLY: n_grid,Wbecke,XYZbecke

IMPLICIT NONE

REAL(dbl) :: rho(1:n_grid),a(1:n_grid),f(1:n_grid),Jdiag(1:n_grid)
REAL(dbl) :: gUEG(1:n_grid),dgUEG(1:n_grid)
REAL(dbl) :: dij(1:n_grid),pdist(1:n_grid),dR(1:3,1:n_grid)
INTEGER(istd) :: k,l

LOGICAL :: debug

debug =.false.


DO k=1,n_grid
   !Compute the distance between the grid points
   dR(1,:) = XYZbecke(1,:)-XYZbecke(1,k)
   dR(2,:) = XYZbecke(2,:)-XYZbecke(2,k)
   dR(3,:) = XYZbecke(3,:)-XYZbecke(3,k)
   FORALL(l=1:n_grid)
         dij(l) = SQRT(DOT_PRODUCT(dR(1:3,l),dR(1:3,l)))
   ENDFORALL
   pdist(:) = a(k)*dij(:)

   !We do a simple normalization integral with a one-point
   !"width" of the hole.  This gives us, as equations to solve,
   ! 0 = [ INT(rho(r)h(a|r-r'|)dr) - (-1) ] * rho(r')

   !The default value of the g function is 1.0, which occurs at coalescence (dij= 0)
   !This cutoff comes from my experimentation.  It seems that round-off is *killing*
   !us here.  Once pdist < epsilon**(1/4) then we only have that g = 1 to within 
   !SQRT(epsilon) and dg = 0 to within (epsilon)**(1/4).  The following Taylor 
   !Series seem more accurate up to about .03.
   WHERE(pdist <= .03_dbl) 
        gUEG = 1.0_dbl - pdist**2/10 + pdist**4/280 - pdist**6/15120
   ELSEWHERE
        gUEG = 3 * (sin(pdist)/pdist**3 - cos(pdist)/pdist**2) 
   ENDWHERE
   WHERE(pdist <= .03_dbl)
        dgUEG = 0.0_dbl - pdist/5 + pdist**3/70 - pdist**5/2520 + pdist**7/166320
   ELSEWHERE
        dgUEG = 3*(sin(pdist)/pdist**2 - gUEG/pdist)
   ENDWHERE

   f(k) = Wbecke(k)*rho(k)*(DOT_PRODUCT(Wbecke(:),rho(:)*(-1)*gUEG(:)**2) + 1)
   Jdiag(k) = -2*Wbecke(k)*rho(k)*DOT_PRODUCT(Wbecke(:),rho(:)*gUEG(:)*dij(:)*dgUEG(:))
ENDDO   

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
! gUEG(1:n_grid) -- the square root of the exchange hole in the                 !
!                            uniform electron gas.                              !
! pdist(1:n_grid) -- the "effective" kF value times the distance                !
!                             between points in the 2-point model.              !
! dR(1:n_grid) -- the distance vector between two grid points.                  !
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
REAL(dbl) :: gUEG(1:n_grid)
REAL(dbl) :: pdist(1:n_grid),dR(1:3,1:n_grid)
INTEGER(istd) :: k,l

LOGICAL :: debug

debug =.false.

DO k=1,n_grid
   !Compute the distance between the grid points
   dR(1,:) = XYZbecke(1,:)-XYZbecke(1,k)
   dR(2,:) = XYZbecke(2,:)-XYZbecke(2,k)
   dR(3,:) = XYZbecke(3,:)-XYZbecke(3,k)
   FORALL(l=1:n_grid)
          pdist(l) = ((a(k)**p_mean + a(l)**p_mean)/2)**(1.0_dbl/p_mean)        &
                      *SQRT(DOT_PRODUCT(dR(1:3,l),dR(1:3,l)))
   ENDFORALL

   !The default value of the g function is 1.0, which occurs at coalescence (dij= 0)
   !This cutoff comes from my experimentation.  It seems that round-off is *killing*
   !us here.  Once pdist < epsilon**(1/4) then we only have that g = 1 to within 
   !SQRT(epsilon) and dg = 0 to within (epsilon)**(1/4).  The following Taylor 
   !Series seem more accurate up to about .03.
   WHERE(pdist <= .03_dbl) 
        gUEG = 1.0_dbl - pdist**2/10 + pdist**4/280 - pdist**6/15120
   ELSEWHERE
        gUEG = 3 * (sin(pdist)/pdist**3 - cos(pdist)/pdist**2) 
   ENDWHERE

   f(k) = Wbecke(k)*rho(k)*(DOT_PRODUCT(Wbecke(:),rho(:)*(-1)*gUEG(:)**2) + 1)
   
ENDDO
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
! gUEG(1:n_grid) -- the square root of the exchange hole in the                 !
!                            uniform electron gas.                              !
! dgUEG(1:n_grid) -- the derivative of gUEG with respect to the a               !
!                             (the effective kf value).  Used to compute the    !
!                             Jacobian.                                         !
! pdist(1:n_grid) -- the "effective" kF value times the distance                !
!                             between points in the 2-point model.              !
! dR(1:n_grid) -- distance between grid points.                                 !
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
REAL(dbl) :: gUEG(1:n_grid),dgUEG(1:n_grid)
REAL(dbl) :: pdist(1:n_grid),dR(1:3,1:n_grid)
INTEGER(istd) :: k,l

LOGICAL :: debug

debug = .false.

DO k=1,n_grid
   !Compute the distance between the grid points
   dR(1,:) = XYZbecke(1,:)-XYZbecke(1,k)
   dR(2,:) = XYZbecke(2,:)-XYZbecke(2,k)
   dR(3,:) = XYZbecke(3,:)-XYZbecke(3,k)
   FORALL(l=1:n_grid)
          pdist(l) = ((a(k)**p_mean + a(l)**p_mean)/2)**(1.0_dbl/p_mean)        &
                      *SQRT(DOT_PRODUCT(dR(1:3,l),dR(1:3,l)))
   ENDFORALL

   !The default value of the g function is 1.0, which occurs at coalescence (dij= 0)
   !This cutoff comes from my experimentation.  It seems that round-off is *killing*
   !us here.  Once pdist < epsilon**(1/4) then we only have that g = 1 to within 
   !SQRT(epsilon) and dg = 0 to within (epsilon)**(1/4).  The following Taylor 
   !Series seem more accurate up to about .03.
   WHERE(pdist <= .03_dbl) 
        gUEG = 1.0_dbl - pdist**2/10 + pdist**4/280 - pdist**6/15120
   ELSEWHERE
        gUEG = 3 * (sin(pdist)/pdist**3 - cos(pdist)/pdist**2) 
   ENDWHERE

   WHERE(pdist <= .03_dbl)
        dgUEG = 0.0_dbl - pdist/5 + pdist**3/70 - pdist**5/2520 + pdist**7/166320
   ELSEWHERE
        dgUEG = 3*(sin(pdist)/pdist**2 - gUEG/pdist)
   ENDWHERE

   f(k) = Wbecke(k)*rho(k)*(DOT_PRODUCT(Wbecke(:),rho(:)*(-1)*gUEG(:)**2) + 1)
   
   Jdiag(k) = -2*Wbecke(k)*rho(k)*DOT_PRODUCT(Wbecke(:),rho(:)*gUEG(:)*dgUEG(:)                 &
                                    *(a(k)+EPSILON(a(k))**2)**(p_mean-1)                        &
                                    *pdist(:)/(a(k)**p_mean+a(:)**p_mean+EPSILON(a(k))**2))
ENDDO


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

end subroutine normalization2ptJacobian

END MODULE holesubroutines
