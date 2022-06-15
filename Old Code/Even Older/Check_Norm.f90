! This MODULE contains the following SUBROUTINE:
!
! SUBROUTINE Check_Normalization(Nr,Lebedev_order,GPoints,Rm,Z_Atomic)
! This subroutine computes the normalization of the density as a check.
!

MODULE CheckNormalization

CONTAINS

!==========================================================================
!=========================================================================

! SUBROUTINE Check_Normalization(Nr,Lebedev_order,GPoints,Rm,Z_Atomic)
! This subroutine computes the normalization of the density as a check.

SUBROUTINE Check_Normalization(Nr,Lebedev_order,GPoints,Rm,Z_Atomic)

USE nrtype			! it uses the module nrtype where the kind of used variables is decleared
USE DensityFunction ! It uses the module DensityFunction where the functional form of the
					! used densities is computed
USE GridSubroutines	! It uses the module GridSubroutines where the construction and reading
					! of the grid is performed


IMPLICIT NONE


INTEGER, INTENT(IN) :: Nr			 !(input) Number of radial points
INTEGER, INTENT(IN) :: Lebedev_order !(input) Number of Lebedev points
REAL(dp), INTENT(IN) :: Rm			 !(input) Rm parameter for the construction of the radial grid
INTEGER, INTENT(IN) :: Z_Atomic		 !(input) Atomic number
INTEGER, INTENT(IN) :: GPoints		 !(input) The number of points in the grid.

REAL(dp) :: NumPoint(1:GPoints) ! Array containing the labels for the grid points
REAL(dp) :: W(1:GPoints)		! Array containing the weights of the grid points
REAL(dp) :: X(1:GPoints)		! Array containing the R coordiante of the grid
REAL(dp) :: Y(1:GPoints)		! Array containing the THETA coodinate of the grid
REAL(dp) :: Z(1:GPoints)		! Array containing the PHI coordinate of the grid
REAL(dp) :: Norm_density(1:GPoints) ! Auxiliary array. Integrand

INTEGER :: i, j  !! counters
INTEGER :: ierror !! status integer varaible

! A helpful message
write(*,*) 'You are in the checking Normalization  part of the program!'

! The subroutine to get the grid points is invoked.
! This subroutine is contained in the MODULE GridSubroutines
CALL Grid_Points(GPoints,W,X,Y,Z)

! For all the grid points the "integrand" is evaluated
! That is, the function at the grid points times the corresponding wiegth
DO i = 1, GPoints
	Norm_density(i) = W(i)*Density(Z_Atomic,X(i),Y(i),Z(i)) 
END DO

! The file "Check_Normalization.out" is open and the result of the
! density integration is printed out. This file is labeld as UNIT 1000
OPEN(UNIT=1000,STATUS='UNKNOWN',FILE='Check_Normalization.out', &
ACTION='write')

! It prints useful information on the file
WRITE(1000,*) 'Atom = ' , Z_Atomic 
WRITE(1000,*) 'Nr = ' , Nr
WRITE(1000,*) 'NLeb = ' , Lebedev_order
WRITE(1000,*) 'Ntot = ' , GPoints

! The result of the integral is printed out
! The extra factor 4*Pi is a consquence of the way how the Lebedev quadrature is constructed
WRITE(1000,*) 'Normalization value = ' , SUM(Norm_density(1:GPoints)) * (4.d0*Pi) 
  
CLOSE(1000)

! The value of the integral is printed out on the screen
write(*,*) 'Norm value = ' , SUM(Norm_density(1:GPoints)) * (4.d0*Pi)  

! A useful message
write(*,*) 'The checking Normalization part of the program! is finished!'
write(*,*) '						'


END SUBROUTINE Check_Normalization



END MODULE CheckNormalization
