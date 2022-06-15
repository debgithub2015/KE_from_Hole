! This MODULE contains the following SUBROUTINES:
!
! SUBROTUINE IntegrateFunction3D(func3D,Nr,Rm,Lebedev_order,Value3DIntegral)
! This subroutine computes the integral of a 3D function
!
! SUBROUTINE IntegrateFunction6D(func6D,Nr,Rm,Lebedev_order,Value6DIntegral)
! This subroutine computes the integral of a 6D function
!

MODULE IntegralsOnGrid

CONTAINS


!=================================================================================================
!=================================================================================================

! SUBROTUINE IntegrateFunction3D(func3D,Nr,Rm,Lebedev_order,Value3DIntegral)
! This subroutine computes the integral of a 3D function

SUBROUTINE IntegrateFunction3D(func3D,Nr,Rm,Lebedev_order,Value3DIntegral)

USE GridSubroutines		! Module that contains the suboutine to get the grid points.
USE nrtype				! Module that contains the type of used variables

IMPLICIT NONE

REAL(dp), EXTERNAL :: func3D		 !(input,function) The previously defined function that is going to be 
									 ! integrated.
									 ! This function is defined in the MODULE DensityFunction
INTEGER, INTENT(IN) :: Nr			 !(input) Number of radial points
INTEGER, INTENT(IN) :: Lebedev_order !(input) Number of angular points
REAL(dp), INTENT(IN) :: Rm			 !(input) Rm parameter to construct the radial grid
REAL(dp), INTENT(OUT) :: Value3DIntegral ! (output) Value of the 3D integral


REAL(dp), ALLOCATABLE :: Wtot(:), Xtot(:), Ytot(:), Ztot(:) ! Allocatable arrays
															! Parameters of the grid
															! weights and (r,theta,phi) coordinates
															! respectively
REAL(dp), ALLOCATABLE :: Integrand(:) ! Auxiliary array used to construct the "integrand"
INTEGER :: i						  ! Just a counter
INTEGER :: NoPoints					  ! Number of points in the grid


! Since we are deling with the tensor product of two quadratures, the total number of
! points is just the number of radial points times the number of angular ones.
NoPoints = Nr*Lebedev_order

! The arrays that will contain the weights, the coordinates of the grid points and the
! auxiliary Integrand array are allocated to have dimension NoPoints.
ALLOCATE(Wtot(1:NoPoints),Xtot(1:NoPoints),Ytot(1:NoPoints),Ztot(1:NoPoints),Integrand(1:NoPoints))

! The gird points and weigth are called by invoking this subroutine
! which is contained in the module GridSubroutines
CALL Grid_Points(NoPoints,Wtot,Xtot,Ytot,Ztot)

!For each of the grid points the function is evaluated and multiplied by the weights
DO i = 1 , NoPoints
	Integrand(i) = Wtot(i)*func3D(Xtot(i),Ytot(i),Ztot(i))
END DO

! The value of the integral is computed by just adding all the elements in the 
! array Integrand. The factor of 4*Pi is, again, just a consequence of how the Lebedev
! quadrature is constructed
Value3DIntegral = 4.d0*Pi*Sum(Integrand(1:NoPoints))


! All the previously allocated allocatable arrays are deallocated
! (this sentence it was intented to be a "joke". If you do not find it funny, get rid of it and/or
! change it to a real/better one)
DEALLOCATE(Wtot,Xtot,Ytot,Ztot,Integrand)


END SUBROUTINE IntegrateFunction3D


!=================================================================================================
!=================================================================================================

! SUBROUTINE IntegrateFunction6D(func6D,Nr,Rm,Lebedev_order,Value6DIntegral)
! This subroutine computes the integral of a 6D function
!

SUBROUTINE IntegrateFunction6D(func6D,Nr,Rm,Lebedev_order,Value6DIntegral)

USE GridSubroutines		! This module contains the subroutine to get the grid used
USE nrtype				! This module contains the type of used variables

IMPLICIT NONE

REAL(dp), EXTERNAL :: func6D			 !(input,function) The previously 6D defined function 
										 ! to be integrated
INTEGER, INTENT(IN) :: Nr				 !(input) Number of radial points
INTEGER, INTENT(IN) :: Lebedev_order	 !(input) Number of angular points
REAL(dp), INTENT(IN) :: Rm				 !(input) The Rm parameter for constructing the radial grid
REAL(dp), INTENT(OUT) :: Value6DIntegral !(output) Value of the 6D integral


REAL(dp), ALLOCATABLE :: Wtot(:), Xtot(:), Ytot(:), Ztot(:) !Allocatble arrays that contain the weigths
															! and the points of the grid.
REAL(dp), ALLOCATABLE :: Integrand(:,:)		! Auxiliary allocatable array to perform the integral
INTEGER :: i, j			! Just counters
INTEGER :: NoPoints		! Number of grid points

! Since we are dealing with the tensor product of two quadratures, the number of total points is
! just the product of the radial and angular number of points
NoPoints = Nr*Lebedev_order

! The arrays that will contain the weights, the coordinates of the grid points and the
! auxiliary Integrand array are allocated to have dimension NoPoints.
ALLOCATE(Wtot(1:NoPoints),Xtot(1:NoPoints),Ytot(1:NoPoints),Ztot(1:NoPoints))
ALLOCATE(Integrand(1:NoPoints,1:NoPoints))

! The gird points and weigth are called by invoking this subroutine
! which is contained in the module GridSubroutines
CALL Grid_Points(NoPoints,Wtot,Xtot,Ytot,Ztot)

! The 6D function is evaluated at all the grid points
DO i = 1 , NoPoints
	DO j = 1, NoPoints
		Integrand(i,j) = Wtot(i)*Wtot(j)*func6D(Xtot(i),Ytot(i),Ztot(i),Xtot(j),Ytot(j),Ztot(j))
	END DO
END DO


! The value of the integral. Remember there is a factor of 4*Pifor each 3D integral
Value6DIntegral = ((4.d0*Pi)**(2.d0))*Sum(Integrand(1:NoPoints,1:NoPoints))

! All the previously allocated allocatable arrays are deallocated
! (this sentence it was intented to be a "joke". If you do not find it funny, get rid of it and/or
! change it to a real/better one)
DEALLOCATE(Wtot,Xtot,Ytot,Ztot,Integrand)

END SUBROUTINE IntegrateFunction6D


!=================================================================================================
!=================================================================================================


END MODULE IntegralsOnGrid