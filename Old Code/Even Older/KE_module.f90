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
! compute_ke_idempotent -- computes kinetic energy functionals for the          !
!                          idempotent density matrix.                           !
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
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_GRID                                !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_grid module!!!                                                  ---- !
! ngrid -- number of points on the grid, including the points used to compute   !
!          the derivatives.                                                     !
! Wint(1:ngrid) -- the integration weights.                                     !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_density                             !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_density module!!!                                               ---- !
! rho(:) -- the density at a point                                              !
! Drho(1:3,:) -- the gradient of the density at a grid point.                   !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_KE                                  !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_KE module!!!                                                    ---- !
! KE_conv(1:4) -- **CONVENTIONAL** kinetic energy functionals.                  !
!                   1.  Thomas-Fermi.                                           !
!                   2.  Weizsacker.                                             !
!                   3.  Gradient expansion (TF + 1/9 Weizsacker)                !
!                   4.  Thomas-Fermi-(1/5)Weizsacker.                           !
!-------------------------------------------------------------------------------!

SUBROUTINE compute_ke_conventional

IMPLICIT NONE

USE kinds
USE vars_grid; USE vars_density; USE constants

REAL(dbl) :: localKE(1:ngrid)
INTEGER :: i  !just a counter

!Compute the Thomas-Fermi kinetic energy.  If problems in this statement use BLAS and DDOT.
localKE(:) = rho(:)**(5_dbl/3)*(3_dbl/10)*(6*pi**2)**(2_dbl/3)
KE_conv(1) = DOT_PRODUCT(Wint,localKE)

!Compute the Weizsacker k.e..  If problems, use DDOT from BLAS.
FORALL(i=1,ngrid)
      localKE(i) = DOT_PRODUCT(Drho(1:3,i),Drho(1:3,i))/(8*rho(i))
ENDFORALL
KE_conv(2) = DOT_PRODUCT(Wint,localKE)

KE_conv(3) = KE_conv(1) + KE_conv(2)/9
KE_conv(4) = KE_conv(1) + KE_conv(2)/5

end subroutine compute_ke_conventional

!-------------------------------------------------------------------------------!
!                               compute_ke_improved                             !
!                                                                               !
! This subroutine computes the "improved" kinetic energy functionals with       !
! different choices for the normalization.                                      !
!-------------------------------------------------------------------------------!
!                          INTERNAL DICTIONARY                                  !
!                                                                               !
! local_KE(:) -- the local kinetic energy.                                      !
! D2g_LDA -- g"(0) from my notes--a prefactor for LDA.                          !
! D2g_LeeP -- g"(0) from my notes--a prefactor for the Gaussian hole.           !
! D2g_Lrnz -- g"(0) from my notes--a prefactor for the squard Lorentzian hole.  !
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_GRID                                !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_grid module!!!                                                  ---- !
! ngrid -- number of points on the grid, including the points used to compute   !
!          the derivatives.                                                     !
! Wint(1:ngrid) -- the integration weights.                                     !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_density                             !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_density module!!!                                               ---- !
! rho(:) -- the density at a point                                              !
! Drho(1:3,:) -- the gradient of the density at a grid point.                   !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       VARIABLES FROM VARS_KE                                  !
! ----The definitive statement of these definitions is found in the        ---- !
! ----vars_KE module!!!                                                    ---- !
! KE_LDA(1:4) -- the LDA kinetic energy with                                    !
!                   1.  The hole normalized as for a uniform electron gas.      !
!                   2.  The hole normalized approximately, point-by-point.      !
!                   3.  The hole normalized exactly, using pairs of points.     !
!                   4.  The hole using idempotency to force N-representability. !
! KE_LeeP(1:4) -- the Lee-Parr "Gaussian" hole model. Indices as in KE_LDA.     !
! KE_Lrnz(1:4) -- the Lorentzian hole model.  Indices as in KE_LDA.             !
! a_LDA(1:ngrid,1:3) -- the alpha values used to control the curvature of the   !
!                       hole at same-spin electron coalescence.                 !
!                   1.  From uniform electron gas.                              !
!                   2.  Based on one-point normalization.                       !
!                   3.  Based on exact normalization.                           !
! a_LeeP(1:ngrid,1:3) -- same as a_LDA but for the Lee-Parr Gaussian hole model.!
! a_Lrnz(1:ngrid,1:3) -- same as a_LDA but for the squared Lorentzian hole model!
!-------------------------------------------------------------------------------!

subroutine compute_ke_improved()

IMPLICIT NONE

USE kinds
USE vars_grid; USE vars_density; USE constants

REAL(dbl) :: localKE(1:ngrid),D2g_LDA,D2g_LeeP,D2g_Lrnz
INTEGER :: i  !just a counter

D2g_LDA = ???
D2g_LeeP = ???
D2g_Lrnz = ???

!Compute LDA K.E. values:
!Uniform electron gas normalization
local_KE(:) = rho(:)*D2g_LDA*a_LDA(:,1)**2
KE_LDA(1) = DOT_PRODUCT(Wint,localKE)

!One-point hole normalization
local_KE(:) = rho(:)*D2g_LDA*a_LDA(:,2)**2
KE_LDA(2) = DOT_PRODUCT(Wint,localKE)

!Exact hole normalization
local_KE(:) = rho(:)*D2g_LDA*a_LDA(:,3)**2
KE_LDA(3) = DOT_PRODUCT(Wint,localKE)

!Compute Gaussian hole K.E. values:
!Uniform electron gas normalization
local_KE(:) = rho(:)*D2g_LeeP*a_LeeP(:,1)**2
KE_LeeP(1) = DOT_PRODUCT(Wint,localKE)

!One-point hole normalization
local_KE(:) = rho(:)*D2g_LeeP*a_LeeP(:,2)**2
KE_LeeP(2) = DOT_PRODUCT(Wint,localKE)

!Exact hole normalization
local_KE(:) = rho(:)*D2g_LeeP*a_LeeP(:,3)**2
KE_LeeP(3) = DOT_PRODUCT(Wint,localKE)

!Compute squared Lorentzian hole K.E. values:
!Uniform electron gas normalization
local_KE(:) = rho(:)*D2g_Lrnz*a_Lrnz(:,1)**2
KE_Lrnz(1) = DOT_PRODUCT(Wint,localKE)

!One-point hole normalization
local_KE(:) = rho(:)*D2g_Lrnz*a_Lrnz(:,2)**2
KE_Lrnz(2) = DOT_PRODUCT(Wint,localKE)

!Exact hole normalization
local_KE(:) = rho(:)*D2g_Lrnz*a_Lrnz(:,3)**2
KE_Lrnz(3) = DOT_PRODUCT(Wint,localKE)

!Add on Weizsacker contribution
KE_LDA(:) = KE_LDA(:) + T_conv(2)
KE_LeeP(:) = KE_LeeP(:) + T_conv(2)
KE_Lrnz(:) = KE_Lrnz(:) + T_conv(2)

end subroutine compute_ke_improved

!-------------------------------------------------------------------------------!
!                           compute_ke_idempotent                               !
!                                                                               !
! This subroutine computes the idempotent kinetic energies.  The density matrix !
! is determined and differentiated "on the fly" in some sense, which ensures    !
! that the memory demands are not too large.                                    !
!-------------------------------------------------------------------------------!
! Sequence of events:                                                           !
! 1.  Starting from the exactly normalized exchange hole, compute the density   !
!     matrix.                                                                   !
! 2.  Use finite differencing to evaluate the local kinetic energy and then the !
!     kinetic energy.                                                           !
!-------------------------------------------------------------------------------!

subroutine compute_ke_idempotent

IMPLICIT NONE

USE kinds
USE vars_grid; USE vars_density; USE constants

REAL(dbl) :: localKE(1:ngrid),DM1(1:ngrid,1:ngrid)


!Compute idempotent density matrix.

!Compute local kinetic energy for every grid point
localKE = 0.0_dbl
j = 1
DO 
   localKE(j) = (SUM(DM1(j,j+1:j+6))-6*DM1(j,j))/h**2
   !The next point is seven points later in the grid.
   j = j + 7
ENDDO








