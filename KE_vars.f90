!---------------------------------------------------------------------------------------!
!                               KE_vars                                                 !
!                                                                                       !
! This module contains the kinetic energy from various functionals computed with various!
! normalization conventions.                                                            !
!---------------------------------------------------------------------------------------!
!                           DICTIONARY                                                  !
!                                                                                       !
! KE_LDA(1:3) -- the LDA kinetic energy with                                            !
!                   1.  The hole normalized as for a uniform electron gas.              !
!                   2.  The hole normalized approximately, point-by-point.              !
!                   3.  The hole normalized exactly, using pairs of points.             !
! KE_conv(1:10) -- **CONVENTIONAL** kinetic energy functionals.                         !
!                   1.  Thomas-Fermi.                                                   !
!                   2.  Weizsacker.                                                     !
!                   3.  TF + 1/9 W                                                      !
!                   4.  TF + 1/5 W                                                      !
!                       other functionals.                                              !
! a_LDA(1:ngrid,1:3,1:2) -- the alpha values used to control the curvature of the hole  !
!                           at same-spin electron coalescence.  1 (alpha) and 2 (beta)  !
!                   1.  From uniform electron gas.                                      !
!                   2.  Based on one-point normalization.                               !
!                   3.  Based on exact normalization.                                   !
! KE_exact -- the exact kinetic energy                                                  !
! p_mean -- the p value used to define the generalized mean.                            !
! rho(1:n_grid,1:2) -- the spin densities.                                              !
! drho(1:3,1:n_grid,1:2) -- the gradient of the spin-densities.                         !
! D2G_LDA -- the second derivative of the square root of the exchange-correlation hole  !
!            of the UEG at the origin.                                                  !
! normalization_tolerance -- the normalization tolerance that is used for the hole.  It !
!                            seems that the error in the k.e. is about equal to this    !
!                            tolerance (to within a factor of 2 or 3, anyway).          !
! Jacobian_type -- The type of Jacobian being used.                                     !
! Jacobian_type -- The type of Jacobian being used.                                     !
!                      0  -- diagonal.                                                  !
!                     -1  -- Broyden's bad method.                                      !
!                      1  -- Broyden's good method.                                     !
!                     -2  -- Broyden's bad method with diagonal-only update             !
!                      2  -- Broyden's good method with diagonal-only update.           !
!                     -3  -- Broyden's bad method with diagonal update and              !
!                            perturbative off-diagonal update.                          !  
!                      3  -- Broyden's good method with diagonal-only update.           ! 
!                            perturbative off-diagonal update.                          !
!                     -4  -- Broyden's bad method with perturbative update to           !
!                            capture changes on the diagonal.                           !
!                      4  -- Broyden's good method with perturbative update to          !
!                            capture changes on the diagonal.                           ! 
! Jvectors -- the number of vectors kept to store the Jacobian in limited-memory        !
!             Broyden's bad method.                                                     !
! Jdiag_update -- the diagonal is updated every this-many iterations.                   !
! Jtrust_radius -- controls when and how the local trust radii are updated.             !
!                             0  --  no local trust radius used.                        !
!                         other   --  local trust radii used to construct the next step.!
! A previous version of the program allowed one to reject certain portions of a previous!
! step using the the trust radius; this never seems to help.                            !
! max_iter -- the maximum number of iterations used to converge.                        !
!---------------------------------------------------------------------------------------!

MODULE KE_vars

USE kinds

IMPLICIT NONE

SAVE

INTEGER :: Jacobian_type,Jtrust_radius,Jdiag_update,Jvectors,max_iter

REAL(dbl) :: p_mean

REAL(dbl) :: KE_LDA(1:3),KE_conv(1:10),KE_exact,D2g_LDA,normalization_tolerance

REAL(dbl), ALLOCATABLE :: a_LDA(:,:,:)
REAL(dbl), ALLOCATABLE :: rho(:,:),drho(:,:,:)

END MODULE KE_vars
