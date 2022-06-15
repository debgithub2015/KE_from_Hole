
!-------------------------------------------------------------------------------!
!                           KE_MODULE                                           !
!                                                                               !
! Contains subroutines for computing the kinetic energy.                        !
!-------------------------------------------------------------------------------!

MODULE KE_MODULE

CONTAINS

!-------------------------------------------------------------------------------!
! compute_ke_conventional -- computes classical kinetic energy functionals.     !
! compute_ke_improved -- computes kinetic energy functionals with improved      !
!                        normalization constants.                               !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!*******************************************************************************!
!-------------------------------------------------------------------------------!
!                               compute_ke_conventional                         !
!                                                                               !
! This subroutine computes classical kinetic energy functionals for comparison  !
! purposes.  The emphasis is on very simple functionals:  TF, Weizsacker, and   !
! the gradient expansion.                                                       !
!-------------------------------------------------------------------------------!
!                          INTERNAL DICTIONARY                                  !
!                                                                               !
! local_KE(:) -- the local kinetic energy.                                      !
! rho(1:ngrid,1:2) -- the density at a point                                    !
! Drho(1:3,1:ngrid,1:2) -- the gradient of the density at a grid point.         !
! debug -- .true. gives extra printing.                                         !
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_GRID                                !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_grid module!!!                                                  ---- !
! n_grid -- number of points on the grid.                                       !
! Wbecke(1:ngrid) -- the integration weights.                                   !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_density                             !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_density module!!!                                               ---- !

!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_KE                                  !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_KE module!!!                                                    ---- !
! KE_conv(1:10) -- **CONVENTIONAL** kinetic energy functionals.                 !
!                   1.  Thomas-Fermi.                                           !
!                   2.  Weizsacker.                                             !
!                   3.  Gradient expansion (TF + 1/9 Weizsacker)                !
!                   4.  Thomas-Fermi+(1/5)Weizsacker.                           !
!                   5.  Other?                                                  !
!                   6.  Other?                                                  !
!                       ....                                                    !
!-------------------------------------------------------------------------------!

SUBROUTINE compute_ke_conventional()

USE kinds
USE grid_vars, ONLY: Wbecke,n_grid
USE constants
USE KE_vars, ONLY: KE_conv,rho,drho
USE variables_wfn, ONLY: KE


IMPLICIT NONE


REAL(dbl)     :: localKE(1:n_grid)
INTEGER(istd) :: i  !just counters
LOGICAL :: debug

debug = .true.

!Compute the Thomas-Fermi kinetic energy.  The two lines are for alpha and beta
!spin, respectively.
!If problems in this statement use BLAS and DDOT.

localKE(:) = rho(:,1)**(5.0_dbl/3)*(3.0_dbl/10)*(6*pi**2)**(2.0_dbl/3)            &
             +rho(:,2)**(5.0_dbl/3)*(3.0_dbl/10)*(6*pi**2)**(2.0_dbl/3)

KE_conv(1) = DOT_PRODUCT(Wbecke,localKE)

!Compute the Weizsacker k.e..  If problems, use DDOT from BLAS.
!The two lines are for the two spin types.
FORALL(i=1:n_grid,(rho(i,1)>0 .and. rho(i,2)>0))
      localKE(i) = DOT_PRODUCT(Drho(1:3,i,1),Drho(1:3,i,1))/(8*rho(i,1))     &
                   + DOT_PRODUCT(Drho(1:3,i,2),Drho(1:3,i,2))/(8*rho(i,2))
ENDFORALL

KE_conv(2) = DOT_PRODUCT(Wbecke,localKE)

KE_conv(3) = KE_conv(1) + KE_conv(2)/9
KE_conv(4) = KE_conv(1) + KE_conv(2)/5

IF (debug) THEN
   WRITE(*,*) 'Thomas-Fermi k.e.             ',KE_conv(1)
   WRITE(*,*) 'von Weizsacker k.e.           ',KE_conv(2)
   WRITE(*,*) 'grad. expansion (TF + 1/9 vW) ',KE_conv(3)
   WRITE(*,*) 'TF + 1/5 vW                   ',KE_conv(4)
   WRITE(*,*) 'exact k.e. from wfn file      ',KE
   WRITE(*,*) ' '
ENDIF

!Other k.e. functionals can be evaluated here.
!*******************
!*******************

end subroutine compute_ke_conventional

!-------------------------------------------------------------------------------!
!                               compute_ke_improved                             !
!                                                                               !
! This subroutine computes the "improved" kinetic energy functionals with       !
! different choices for the normalization.                                      !
!-------------------------------------------------------------------------------!
!                          INTERNAL DICTIONARY                                  !
!                                                                               !
! local_KE(1:n_grid) -- the local kinetic energy.                               !
! D2g_LDA -- g"(0) from my notes--a prefactor for LDA.                          !
! rho(1:ngrid,1:2) -- the density at a point                                    !
! Drho(1:3,1:ngrid,1:2) -- the gradient of the density at a grid point.         !
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_GRID                                !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_grid module!!!                                                  ---- !
! n_grid -- number of points on the grid.                                       !
! Wbecke(1:n_grid) -- the integration weights.                                  !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_KE                                  !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_KE module!!!                                                    ---- !
! KE_LDA(1:3) -- the LDA kinetic energy with                                    !
!                   1.  The hole normalized as for a uniform electron gas.      !
!                   2.  The hole normalized approximately, point-by-point.      !
!                   3.  The hole normalized exactly, using pairs of points.     !
! a_LDA(1:ngrid,1:3,1:2) -- the alpha values used to control the curvature of   !
!                           the hole at same-spin electron coalescence.         !
!                   1.  From uniform electron gas.                              !
!                   2.  Based on one-point normalization.                       !
!                   3.  Based on exact normalization.                           !
! KE_conv(2) -- the Weizsacker kinetic energy.                                  !
!-------------------------------------------------------------------------------!

subroutine compute_ke_improved()

USE kinds
USE grid_vars, ONLY: Wbecke,n_grid
USE KE_vars, ONLY: KE_LDA,KE_conv,a_LDA,rho,D2g_LDA
USE variables_wfn, ONLY: KE
USE normalize_module

IMPLICIT NONE

REAL(dbl)     :: localKE(1:n_grid)
INTEGER(istd) :: i      !just a counter
LOGICAL :: debug

D2g_LDA = -3.0_dbl/5
debug = .true.

!Normalize the holes.
call hole_normalize()

!Compute LDA K.E. values:
!Uniform electron gas normalization.  The two lines are for the two spin types.
localKE(:) = -.5_dbl*(rho(:,1)*D2g_LDA*a_LDA(:,1,1)**2 + rho(:,2)*D2g_LDA*a_LDA(:,1,2)**2)

KE_LDA(1) = DOT_PRODUCT(Wbecke,localKE) 

!One-point hole normalization.  The two lines are for the two spin types.
localKE(:) = -.5_dbl*(rho(:,1)*D2g_LDA*a_LDA(:,2,1)**2 + rho(:,2)*D2g_LDA*a_LDA(:,2,2)**2)
            
KE_LDA(2) = DOT_PRODUCT(Wbecke,localKE)

!Exact hole normalization.  The two lines are for the two spin types.
localKE(:) = -.5_dbl*(rho(:,1)*D2g_LDA*a_LDA(:,3,1)**2 + rho(:,2)*D2g_LDA*a_LDA(:,3,2)**2)

KE_LDA(3) = DOT_PRODUCT(Wbecke,localKE)

IF (debug) THEN
   WRITE(*,*) 'for g"(0) = -3/5'
   WRITE(*,*) 'Thomas-Fermi-ish k.e.       ',KE_LDA(1)
ENDIF

!Add on Weizsacker contribution
KE_LDA(:) = KE_LDA(:) + KE_conv(2)

IF (debug) THEN
   WRITE(*,*) 'LDA-based k.e. from hole    ',KE_LDA(1)
   WRITE(*,*) 'one-point normalization KE  ',KE_LDA(2)
   WRITE(*,*) 'two-point normalization KE  ',KE_LDA(3)
   WRITE(*,*) 'exact k.e. from wfn file    ',KE
   WRITE(*,*) ' '
ENDIF

end subroutine compute_ke_improved

end module KE_module


