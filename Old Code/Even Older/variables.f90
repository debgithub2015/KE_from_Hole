!***************************************************************************************!
!									  MODULE Kinds                                      !
!                                                                                       !
!                                       Paul Ayers                                      !
!                                                                                       !
!                                       July 4, 2000                                    !
!                                                                                       !
!     Provenance; last modified July 6, 2000. (streamlined module)                      !
!                 renamed and revised, May 12, 2002.                                    !
! This module assigns various integer and real variable types to meaningful expressions.!
!***************************************************************************************!

MODULE kinds

IMPLICIT NONE

! First we define three levels of integer precision
  INTEGER, PARAMETER ::  ihalf  = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER ::  isngl = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER ::  idbl = SELECTED_INT_KIND(8) 
  INTEGER, PARAMETER ::  itpl = SELECTED_INT_KIND(12)

! Now we define 3 levels of real precision.  
  INTEGER, PARAMETER ::  sngl = SELECTED_REAL_KIND(p=5)
  INTEGER, PARAMETER ::  dbl  = SELECTED_REAL_KIND(p=12)
! We assume that requiring 18 digits of precision is enough to force the program 
! into more the double precision; when that exists.  If it does not exist,
! then SELECTED_REAL_KIND < 0, so trpl <= dbl.
  INTEGER, PARAMETER ::  trpl = MAX(SELECTED_REAL_KIND(p=18),dbl)
! We try for quadruple precision:
  INTEGER, PARAMETER ::  qdpl = MAX(SELECTED_REAL_KIND(p=25),trpl)

END MODULE kinds

!-----------------------------------------------------------------------!
!                               CONSTANTS                               !
!                                                                       !
! This module contains the values of key constants like pi.             !
!                                                                       !
!***********************************************************************!
! Dependencies:                                                         !
!    kinds -- module containing variable types.                         !
!                                                                       !
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
!                                                                       !
! Contains:                                                             !
!     (none)                                                            !
!                                                                       !
!     OBITER DICTUM: Variable-storage modules should not contain any    !
!                    subroutines.                                       !
!                                                                       !
!-----------------------------------------------------------------------!
! Authors:                                                              !
!     Juan I. Rodriguez (rodrigji@mcmaster.ca)                          !
!     Paul W. Ayers (ayers@mcmaster.ca)                                 !
!                                                                       !
! Provenance:                                                           !
!     May 1, 2007 (PWA) documentation of existing module.               !
!                                                                       !
!-----------------------------------------------------------------------!

MODULE constants

USE kinds

IMPLICIT NONE

SAVE

REAL(dbl) :: pi = 3.14159265358979323846_dbl

END MODULE constants



!---------------------------------------------------------------------------------------!
!                               VARS_GRID                                               !
!                                                                                       !
! This module contains the key variables for defining the grid.  The values of these    !
! parameters are read in from the file InputFile.inp                                    !
!---------------------------------------------------------------------------------------!
!                           DICTIONARY                                                  !
!                                                                                       !
! Nr -- minimum number of radial points.                                                !
! Rm -- Rm parameter for constructing the radial grid.  Similar to Becke.               !
! Lebedev_order -- Number of Lebedev points.                                            !
! ngrid -- number of points on the grid, including the points used to compute the       !
!          derivatives.                                                                 !
! X(1:3,1:ngrid) -- the grid points.                                                    !
! Wint(1:ngrid) -- the integration weights.                                             !
! exp_step -- the step size will be EPSILON**exp_step.                                  !
! h -- the step size for evaluating the derivatives, EPSILON**exp_step.                 !
!---------------------------------------------------------------------------------------!

MODULE vars_grid

IMPLICIT NONE

SAVE

USE kinds

INTEGER   :: Nr,Lebedev_order,PointsGrid,ngrid
REAL(dbl) :: Rm,exp_step,h
REAL(dbl), ALLOCATABLE :: X(:,:),Wint(:)

END MODULE vars_grid

!---------------------------------------------------------------------------------------!
!                               VARS_density                                            !
!                                                                                       !
! This module contains the density and the density gradient.                            !
!---------------------------------------------------------------------------------------!
!                           DICTIONARY                                                  !
!                                                                                       !
! rho(:,1:2) -- the spin-density at a point                                             !
! Drho(1:3,1:ngrid,1:2) -- the gradient of the spin-density at a grid point.            !
! Z_Atomic -- atomic number                                                             !
!---------------------------------------------------------------------------------------!

MODULE vars_density

IMPLICIT NONE

SAVE

USE kinds

INTEGER :: Z_Atomic

REAL(dbl), ALLOCATABLE :: Drho(:,:,:),rho(:,:)

END MODULE vars_density

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
! KE_conv(1:3) -- **CONVENTIONAL** kinetic energy functionals.                          !
!                   1.  Thomas-Fermi.                                                   !
!                   2.  Weizsacker.                                                     !
!                   3.  A better k.e. functional.                                       !
! a_LDA(1:ngrid,1:3) -- the alpha values used to control the curvature of the hole at   !
!                       same-spin electron coalescence.                                 !
!                   1.  From uniform electron gas.                                      !
!                   2.  Based on one-point normalization.                               !
!                   3.  Based on exact normalization.                                   !
! KE_exact -- the exact kinetic energy                                                  !
! g(1:ngrid,1:ngrid) -- the model for the square root of the hole.                      !
! occ_nums(:) -- the occupation numbers of the electron density.                        !
! DM1(1:ngrid,1:ngrid) -- the model density matrix.                                     !
!---------------------------------------------------------------------------------------!

MODULE VARS_KE

IMPLICIT NONE

SAVE

USE kinds

REAL(dbl) :: KE_LDA(1:3),KE_conv(1:3),KE_exact

REAL(dbl), ALLOCATABLE :: a_LDA(:,:)
REAL(dbl), ALLOCATABLE :: g(:,:),occ_nums(:)

END MODULE vars_KE