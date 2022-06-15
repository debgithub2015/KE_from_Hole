!-------------------------------------------------------------------------------!
!                           INOUT_MODULE                                        !
!                                                                               !
! Contains subroutines for input and output.                                    !
!-------------------------------------------------------------------------------!

MODULE INOUT_MODULE

CONTAINS

!-------------------------------------------------------------------------------!
! read_input -- driver routine for input.                                       !
! input_file -- reads in the parameters that define the grid and the system.    !
! getgrid -- constructs the grid.                                               !
! get_density -- constructs the atomic density and its gradient on the grid.    !
! check_normalization -- checks the normalization of the density.               !
! write_output -- summarizes the output.                                        !
!-------------------------------------------------------------------------------!
!*******************************************************************************!
!                             STARTUP                                           !
!                                                                               !
! A driver subroutine that reads the input file, constructs grids for           !
! integration and differentiation, etc..                                        !
!                                                                               !
! n_atoms -- number of atoms.                                                   !
! Ratom(1:3,1:n_atoms) -- positions of the atoms.                               !
! Zatom(1:n_atoms) -- atomic number of the atoms.                               !
! n_grid -- number of grid points.                                              !
! XYZbecke(1:n_grid) -- grid points on the Becke grid.                          !
! rho(1:n_grid,1:2) -- the spin-density.                                        !
! drho(1:n_grid,1:2) -- the gradient of the spin-density.                       !
! a_LDA(1:n_grid,1:3,1:2) -- the "effective kF"                                 !
! istat -- status for allocating/deallocating.                                  !
!-------------------------------------------------------------------------------!

subroutine startup()

USE kinds
USE grid_vars, ONLY: n_grid,XYZbecke
USE gridmaker
USE process_wfn
USE variables_wfn, ONLY: n_atoms,Zatom,atm_position,wfn_filename
USE KE_vars, ONLY: a_LDA,rho,drho

IMPLICIT NONE

INTEGER :: istat

DEALLOCATE(a_LDA,rho,drho,stat=istat)

call input_grid()   !Reads the input file containing the grid-parameters

call read_wfn()     !Reads the wfn file.

call Becke(n_atoms,atm_position,Zatom)    !Generates the grid.

!ALLOCATE ARRAYS
ALLOCATE(a_LDA(1:n_grid,1:3,1:2),rho(1:n_grid,1:2))
ALLOCATE(drho(1:3,1:n_grid,1:2))

call densitygradient_wfn(rho,drho,XYZbecke,n_grid)  !Evaluates the density and density gradient.

end subroutine startup

!-------------------------------------------------------------------------------!
!*******************************************************************************!
!                             re_start                                          !
!                                                                               !
! A driver subroutine that generates a new grid and evaluates the density and   !
! gradient on the new grid.                                                     !
!                                                                               !
! n_atoms -- number of atoms.                                                   !
! Ratom(1:3,1:n_atoms) -- positions of the atoms.                               !
! Zatom(1:n_atoms) -- atomic number of the atoms.                               !
! n_grid -- number of grid points.                                              !
! XYZbecke(1:n_grid) -- grid points on the Becke grid.                          !
! rho(1:n_grid,1:2) -- the spin-density.                                        !
! drho(1:n_grid,1:2) -- the gradient of the spin-density.                       !
! a_LDA(1:n_grid,1:3,1:2) -- the "effective kF"                                 !
! istat -- status for allocating/deallocating.                                  !
!-------------------------------------------------------------------------------!

subroutine re_start()

USE kinds
USE grid_vars, ONLY: n_grid,XYZbecke
USE gridmaker
USE process_wfn
USE variables_wfn, ONLY: n_atoms,Zatom,atm_position,wfn_filename
USE KE_vars, ONLY: a_LDA,rho,drho

IMPLICIT NONE

INTEGER :: istat

DEALLOCATE(a_LDA,rho,drho,stat=istat)

call Becke(n_atoms,atm_position,Zatom)    !Generates the grid.

!ALLOCATE ARRAYS
ALLOCATE(a_LDA(1:n_grid,1:3,1:2),rho(1:n_grid,1:2))
ALLOCATE(drho(1:3,1:n_grid,1:2))

call densitygradient_wfn(rho,drho,XYZbecke,n_grid)  !Evaluates the density and density gradient.

end subroutine re_start







!-------------------------------------------------------------------------------!
!                           write_output                                        !
!                                                                               !
! This subroutine writes an output file containing the evaluated kinetic        !
! energies.  It also investigates the eigenvalue spectrum of the model density  !
! matrix.                                                                       !
!-------------------------------------------------------------------------------!








!USE CalcMatrix			! Module that contain subroutines to operate with matrices
						! particularly the computation of eigenvectors and eigenvalues

!IMPLICIT NONE


!REAL(dp), ALLOCATABLE :: AMatrix(:,:) ! Array for the xc-hole matrix
!REAL(dp), ALLOCATABLE :: EigenVectorsAMatrix(:,:) ! Array for the eigenvectors of the xc-hole matrix
!REAL(dp), ALLOCATABLE :: EigenValuesAMatrix(:) ! Array for the eigenvalues of the xc-hole matrix
!INTEGER :: i, j, k !! Just counters
!REAL(dp) :: Value3DIntegral !! Value of the 3D integral
!REAL(dp) :: Value6DIntegral ! Value of the 6D integral



! The matrix to be diagonalizes is allocated as well as their eigenvalues and eigenvectors 
! arrays.
! In this particular example the dimension of the matrix used is the same as the number of
! grid points.
!ALLOCATE(AMatrix(1:PointsGrid,1:PointsGrid),STAT=istat)
!ALLOCATE(EigenVectorsAMatrix(1:PointsGrid,1:PointsGrid),STAT=istat)
!ALLOCATE(EigenValuesAMatrix(1:PointsGrid),STAT=istat)

!DO i = 1 , PointsGrid
!	DO j = 1 , PointsGrid
!		AMatrix(i,j) = i+j
!	END DO
!END DO


! This subroutine performs the computation of the eigenvalues and eigenvectors of the chosen
! matrix to be diagionalized.
!CALL EigenMatrix(PointsGrid,AMatrix,EigenValuesAMatrix,EigenVectorsAMatrix)


end module INOUT_MODULE