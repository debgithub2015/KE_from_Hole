!-------------------------------------------------------------------------------!
!8888888888888888888888888888888888888888888888888888888888888888888888888888888!
!8888888888888888888888888888888888888888888888888888888888888888888888888888888!
!8888888888888888888888888888888888888888888888888888888888888888888888888888888!
!8888888888888888888888888888888888888888888888888888888888888888888888888888888!
!8888888888888888888888888888888888888888888888888888888888888888888888888888888!
!-------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------!
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
!-------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------!
!                       grid_variables                                          !
!                                                                               !
! A module containing variables needed to define the grid.                      !
! Originally by Patrick Bultinck; streamlined by Paul Ayers.                    !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                       DICTIONARY OF VARIABLES                                 !
!                                                                               !
! radfac -- controls the number of radial grid points.  1.0 is Becke's          !
!           recommendation.                                                     !
! n_leb -- the number of Lebedev points.  The choices are: 6,14,26,38,50,74,86, !
!          110,146,170,194,230,266,302,350,434,590,770,974,1202,1454,1730,2030, !
!          2354,2702,3074,3470,3890,4334,4802,5294,5810                         !
! stepf_iter -- controls the stepfunction smoothing between atomic regions.     !
!               Becke recommends 3.0.                                           !
! filename -- filename.wfn is the wavefunction file; this contains the geometry !
!             info we need to make the Becke grid.                              !
! n_grid -- the number of grid points.                                          !
! XYZbecke(1:3,1:n_grid) -- the grid points.                                    !
! Wbecke(1:n_grid) -- the grid weights.                                         !
! radscale -- a real number that lets you scale the Bragg-Slater radii to put   !
!             more (or fewer) points close to the nucleus.  A number smaller    !
!             than 1.0 causes the grid to be more compact; a number greater than!
!             1.0 causes the grid to be more diffuse.                           !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!

module grid_vars

USE kinds

implicit none

SAVE

REAL(dbl), ALLOCATABLE :: Wbecke(:),XYZbecke(:,:)

INTEGER(istd) :: stepf_iter,n_leb,n_grid
REAL(dbl) :: radfac,radscale

end module grid_vars


