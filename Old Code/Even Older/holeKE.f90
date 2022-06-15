!-------------------------------------------------------------------------------!
!                               holeKE                                          !
!                                                                               !
! This program is designed to compute the kinetic energy from an approximate    !
! exchange-correlation hole.  As an initial approach, it will consider the      !
! density to be a sum of spherical Gaussians, and it will only work for atoms.  !
!                                                                               !
! Given the density, the program first computes several standard k.e. formulae, !
! including Thomas-Fermi and Weizsacker, and some formulas we derived that      !
! are related to these.  Next, it normalizes the exchange-correlation hole      !
! either (1) approximately or (2) exactly.  These give new k.e. fucntionals.    !
! Finally, it forces idempotency of the density matrix by a brute strength      !
! procedure and thereby obtains yet a third family of k.e. functionals.         !
!-------------------------------------------------------------------------------!
!                               FILES                                           !
!                                                                               !
! general.out -- an echo of the "program status" from the screen.               !
! table.out -- a table in easy-to-import column format.                         !
!-------------------------------------------------------------------------------!

PROGRAM hole KE

IMPLICIT NONE

USE inout_module
USE normalize_module
USE grid_module

!Open the output file
OPEN(UNIT=101,STATUS='UNKNOWN',FILE='general.out')
OPEN(UNIT=102,STATUS='UNKNOWN',FILE='table.out')

!read in the parameters that define the calculation.
call read_input()

!Construct the grid.
call GetGrid()

!Construct the density and the density gradient.
call get_density()  

!compute "standard" k.e. formulae.
call compute_ke_conventional()

!determine approximate and exact hole normalization.
call hole_normalize()

!compute improved k.e. formulae
call compute_ke_improved()

!force density matrix to be idempotent and recompute k.e.
call compute_ke_idempotent()

!print output
call write_output()

CLOSE(101);CLOSE(102)

END PROGRAM hole KE

USE nrtype				! Module where the type of variables is defined
						! Gotten from numerical receipes
USE InputFile			! Module that contain the subroutine that reads the
						! input file containing the parameters of the grid
USE GridSubroutines		! Module that contains the soubroutines related with the
						! constructed grid. Among others: creating the grid points file.
USE CheckNormalization  ! Module that contain a subroutine to check the normalization
						! of the chosen atomic density.
USE IntegralsOnGrid		! Module that contain the subroutines that perform the
						! 3D and 6D integrals.
USE DensityFunction		! Module containing some of the atomic densities as well
						! as the test functions to perform the 3D and 6D integrals
USE CalcMatrix			! Module that contain subroutines to operate with matrices
						! particularly the computation of eigenvectors and eigenvalues

IMPLICIT NONE

!===========												============
!===========  Set of variables for reading the input file ==============
!=====														============

INTEGER :: Nr		        ! Minimum number of radial points
REAL(dp):: Rm				! Rm parameter when constructing the radial grid
INTEGER :: Lebedev_order    ! Number of Lebedev points
INTEGER :: Z_Atomic		    ! Atomic number
INTEGER :: PointsGrid		! Number of Points in the whole grid (radial times angular)
INTEGER :: istat
REAL(dp), ALLOCATABLE :: AMatrix(:,:) ! Array for the xc-hole matrix
REAL(dp), ALLOCATABLE :: EigenVectorsAMatrix(:,:) ! Array for the eigenvectors of the xc-hole matrix
REAL(dp), ALLOCATABLE :: EigenValuesAMatrix(:) ! Array for the eigenvalues of the xc-hole matrix
INTEGER :: i, j, k !! Just counters
REAL(dp) :: Value3DIntegral !! Value of the 3D integral
REAL(dp) :: Value6DIntegral ! Value of the 6D integral

!===========												============
!===========  ============================================= =============
!=====														=============



write(*,*) 'Auxiliary code for the KEF project'
write(*,*) 'February 11th, 2010'
write(*,*) '						'



! The input file with all the parameters is read
CALL Input_File(Nr,Rm,Lebedev_order,Z_Atomic)


! The file containing the grid to perform the calculation is created
! The name of such a file is grid.txt
CALL GetGrid_File(Nr,Lebedev_order,Rm,Z_Atomic,PointsGrid)


! The subroutine to check the normalization with the constructed grid is invoked
! The value of the normalization is printed out on screen.
CALL Check_Normalization(Nr,Lebedev_order,PointsGrid,Rm,Z_Atomic)


! This subroutine determines the value of the 3D integral for the function FuncTest
! IMPORTANT:The used function should be in spherical coordinates since the grid points
! are in that coordinate system/
! The function FuncTest is defined in the MODULE DensityFunction.
! This function depends on three real variables (r,theta,phi) 
! The result of the integral is stogared in Value3DIntegral
CALL IntegrateFunction3D(FuncTest3D,Nr,Rm,Lebedev_order,Value3DIntegral)


write(*,*) Value3DIntegral


! This subroutine determines the value of the 3D integral for the function FuncTest2
! IMPORTANT:The used function should be in spherical coordinates since the grid points
! are in that coordinate system
! The function FuncTest2 is defined in the MODULE DensityFunction 
! This functions depends on six real varaibles (r1,theta1,phi1,r2,theta2,phi2)
! The result of the integral is stogared in Value6DIntegral
CALL IntegrateFunction6D(FuncTest6D,Nr,Rm,Lebedev_order,Value6DIntegral)


write(*,*) Value6DIntegral



! The matrix to be diagonalizes is allocated as well as their eigenvalues and eigenvectors 
! arrays.
! In this particular example the dimension of the matrix used is the same as the number of
! grid points.
ALLOCATE(AMatrix(1:PointsGrid,1:PointsGrid),STAT=istat)
ALLOCATE(EigenVectorsAMatrix(1:PointsGrid,1:PointsGrid),STAT=istat)
ALLOCATE(EigenValuesAMatrix(1:PointsGrid),STAT=istat)

DO i = 1 , PointsGrid
	DO j = 1 , PointsGrid
		AMatrix(i,j) = i+j
	END DO
END DO


! This subroutine performs the computation of the eigenvalues and eigenvectors of the chosen
! matrix to be diagionalized.
CALL EigenMatrix(PointsGrid,AMatrix,EigenValuesAMatrix,EigenVectorsAMatrix)


