!---------------------------------------------------------------------------------------!
!                               VARS_KE                                                 !
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
! occ_nums(1:ngrid,1:2) -- the occupation numbers of the 1-electron reduced density     !
!                          matrix.                                                      !
!---------------------------------------------------------------------------------------!

MODULE VARS_KE

IMPLICIT NONE

SAVE

USE kinds

INTEGER(istd) :: p_mean

REAL(dbl) :: KE_LDA(1:3),KE_conv(1:10),KE_exact

REAL(dbl), ALLOCATABLE :: a_LDA(:,:,:)
REAL(dbl), ALLOCATABLE :: occ_nums(:)

END MODULE vars_KE