! This module contains the following SUBROUTINES:
!
! SUBROUTINE GetGrid_File(Nr,Lebedev_order,Rm,Z_Atomic,Ntot) which reads the file
! UNIT=15 'Lebedev_RECgrid.txt' This file contains the angular points in cartesian 
! coordinates (x,y,z). Once the points are read they are changed to points in spherical
! coordinates (r,theta,phi).
! Finally the tensor product of the radial and angular grids is performed and the
! complete grid is printed out in the file UNIT=20 'grid.txt'
!
! SUBROUTINE Grid_Points(NoPoints,Wtot,Xtot,Ytot,Ztot)
! This subroutine reads the grid points of the complete grid from the file
! "grid_points.txt" that has been identifed as UNIT 20
!

MODULE GridSubroutines

USE nrtype	! Use the module nrtype where the type of variables used is decleared

CONTAINS

!========================================================================
!========================================================================


! This SUBROUTINE reads the file
! UNIT=15 'Lebedev_RECgrid.txt' This file contains the angular points in cartesian 
! coordinates (x,y,z). Once the points are read they are changed to points in spherical
! coordinates (r,theta,phi).
! Finally the tensor product of the radial and angular grids is performed and the
! complete grid is printed out in the file UNIT=20 'grid.txt'

SUBROUTINE GetGrid_File(Nr,Lebedev_order,Rm,Z_Atomic,Ntot)

IMPLICIT NONE


INTEGER, INTENT(IN) :: Nr			 !(input) Number of points in the radial grid
INTEGER, INTENT(IN) :: Lebedev_order !(input) Number of Lebedev points in the angular grid
REAL(dp), INTENT(IN) :: Rm			 !(input) The value of the parameter Rm to generate the
									 ! radial grid. This value is atom dependent.
INTEGER, INTENT(IN) :: Z_Atomic		 !(input) Atomic number
INTEGER, INTENT(OUT) :: Ntot		 !(output) Number of points of the complete 
									 !(radial times angular) grid

REAL(dp), ALLOCATABLE :: Wang(:)		!Allocatable array where the weights of the
										!angular grid are storaged 
REAL(dp), ALLOCATABLE :: Xang(:)		!Allocatable array. 
										!X coordiante of the Lebedev grid in cartesian coordiantes
REAL(dp), ALLOCATABLE :: Yang(:)		!Allocatable array
										!Y coordiante of the Lebedev grid in cartesian coordiantes
REAL(dp), ALLOCATABLE :: Zang(:)		!Allocatable array
										!Z coordiante of the Lebedev grid in cartesian coordiantes
REAL(dp), ALLOCATABLE :: NangPoint(:)	!Allocatable array
										!Number of points in the angular grid
REAL(dp), ALLOCATABLE :: PHI(:)			!Allocatable array
										!Phi coordiante of each grid point (spherical coordinates)
REAL(dp), ALLOCATABLE :: THETA(:)		!Allocatable array
										!Theta coordinte of each grid point (spherical coordiantes)
REAL(dp), ALLOCATABLE :: R(:)			!Allocatable array
										!Radial coordiante of each grid point
REAL(dp), ALLOCATABLE :: Wr(:)			!Allocatable array
										!Weights of the radial grid
REAL(dp), ALLOCATABLE :: Wtot(:)		!Allocatable array	
										!Weigth of the complete grid
REAL(dp), ALLOCATABLE :: Xtot(:)		!Allocatable array
										!Radial coordinate of the complete grid (r coodiante)
REAL(dp), ALLOCATABLE :: Ytot(:)		!Allocatable array
										!Theta coordiante of the complete grid (theta coordiante)
REAL(dp), ALLOCATABLE :: Ztot(:)		!Allocatable array
										!Phi coordinate of the complete grid (phi coordiante)
REAL(dp), ALLOCATABLE :: NpointTot(:)	!Allocatable array
										!An auxiliary array to label the number of grid point 
										!we are dealing with
INTEGER :: i, j, k						! Just counters
INTEGER :: ierror						! ierror status variable
REAL(dp) :: Rmax						!Auxiliary real value to construct the radial grid
CHARACTER(len=20) :: Lebedev_RECgrid	!Charcater type variable. Determines the Lebedev file to be opened.
INTEGER :: Nang							!Number of angular points


! The following IF condition is just to establish what is the Lebedev file that is going to be
! opened.

IF (Lebedev_order == 38) THEN

	Lebedev_RECgrid = 'Lebedev_38_REC.txt'

 ELSE IF (Lebedev_order == 50) THEN

	Lebedev_RECgrid = 'Lebedev_50_REC.txt'

 ELSE IF (Lebedev_order == 86) THEN

	Lebedev_RECgrid = 'Lebedev_86_REC.txt'

 ELSE IF (Lebedev_order == 110) THEN

	Lebedev_RECgrid = 'Lebedev_110_REC.txt'

 ELSE IF (Lebedev_order == 170) THEN

	Lebedev_RECgrid = 'Lebedev_170_REC.txt'

 ELSE IF (Lebedev_order == 350) THEN

	Lebedev_RECgrid = 'Lebedev_350_REC.txt'

 ELSE IF (Lebedev_order == 1202) THEN

	Lebedev_RECgrid = 'Lebedev_1202_REC.txt'

END IF

! End of the IF 

! Since we will be counting the number of angular points, such a number is set to zero. After that
! the counting sequence will be started.
Nang = 0

!! Open the file where the Lebedev grid points (RECTANGULAR) are stored and counts the number of points in the grid.
!! The name of the file is Lebdev_RECgrid and is declared as UNIT 15
OPEN(UNIT=15,STATUS='old',FILE=Lebedev_RECgrid,ACTION='read',IOSTAT=ierror)
 DO 
   READ(15,*,IOSTAT=ierror)
   IF(ierror /= 0) EXIT     !If the end of file has been reached, then it exists the file
   Nang = Nang + 1			!If a data is contained in that line, then the counter is increased.
 END DO
CLOSE(15)
! The files is closed

! The arrays involving the angular grid are allocated to have dimension Nang, where Nang is the number
! of angular points read in the Lebedev file.
ALLOCATE(NangPoint(1:Nang),Wang(1:Nang),Xang(1:Nang),Yang(1:Nang),Zang(1:Nang),PHI(1:Nang),THETA(1:Nang))


!! Open the file where the Lebedev grid points (RECTANGULAR) are stored and counts the number of points in the grid.
!! The name of the file is Lebdev_RECgrid and is declared as UNIT 15
OPEN(UNIT=15,STATUS='unknown',FILE=Lebedev_RECgrid,ACTION='read',IOSTAT=ierror)
 DO i = 1, Nang
   ! The data contained in the file are read, this are j, Wj, Xj, Yj, Zj. These being
   ! the label of the point, the weight, the x coordianted, y coordinate and z coodinate
   ! respectively
   READ(15,*) NangPoint(i), Wang(i), Xang(i), Yang(i), Zang(i)
 END DO
CLOSE(15)
! the file is closed.


! The following DO is to convert the angular points from cartesian coordiantes (x,y,z) to
! spherical coordiantes (1,theta,phi). Notice that the Lebedev grids are constructed on a unit
! sphere and then in all of the following sequence the value of r is r=1

DO i = 1, Nang

 IF (Xang(i) == 0.d0 .AND. Yang(i) == 0.d0 .AND. Zang(i) /= 0.d0) THEN

    PHI(i) = 0.d0

	IF (Zang(i) > 0.d0) THEN

     THETA(i) = 0.d0

	ELSE
  
     THETA(i) = Pi

	END IF

 ELSE IF (Xang(i) == 0.d0 .AND. Yang(i) /= 0.d0 .AND. Zang(i) == 0.d0) THEN

    THETA(i) = Pi/2.d0

	IF (Yang(i) > 0.d0) THEN

     PHI(i) = Pi/2.d0

	ELSE 

	 PHI(i) = (3.d0/2.d0)*Pi

	END IF

 ELSE IF (Xang(i) == 0.d0 .AND. Yang(i) /= 0.d0 .AND. Zang(i) /= 0.d0) THEN	

    IF (Yang(i) > 0.d0) THEN

	 PHI(i) = Pi/2.d0

	ELSE

	 PHI(i) = (3.d0/2.d0)*Pi

	END IF


	IF (Zang(i) > 0.d0) THEN

	 THETA(i) = ACos(Abs(Zang(i)))

	ELSE

     THETA(i) = ACos(Abs(Zang(i))) + Pi
	
	END IF


 ELSE IF (Xang(i) /= 0.d0 .AND. Yang(i) == 0.d0 .AND. Zang(i) == 0.d0) THEN	

   THETA(i) = Pi/2.d0

   IF (Xang(i) > 0.d0) THEN

    PHI(i) = 0.d0

   ELSE 

    PHI(i) = Pi

   END IF


 ELSE IF (Xang(i) /= 0.d0 .AND. Yang(i) == 0.d0 .AND. Zang(i) /= 0.d0) THEN	


   IF (Xang(i) > 0.d0) THEN

    PHI(i) = 0.d0

   ELSE 

    PHI(i) = Pi

   END IF

   IF (Zang(i) > 0.d0) THEN

    THETA(i) = ACos(Abs(Zang(i)))

   ELSE 

    THETA(i) = ACos(Abs(Zang(i))) + Pi

   END IF


 ELSE IF (Xang(i) /= 0.d0 .AND. Yang(i) /= 0.d0 .AND. Zang(i) == 0.d0) THEN

   THETA(i) = Pi/2.d0

   IF ( (Xang(i) > 0.d0) .AND. (Yang(i) > 0.d0) ) THEN

    PHI(i) = ATan(Yang(i)/Xang(i))

   ELSE IF ( (Xang(i) < 0.d0) .AND. (Yang(i) > 0.d0) ) THEN

    PHI(i) = ATan(Yang(i)/Xang(i)) + (Pi/4.d0)


   ELSE IF ( (Xang(i) < 0.d0) .AND. (Yang(i) < 0.d0) ) THEN

    PHI(i) = ATan(Yang(i)/Xang(i)) + (Pi/2.d0)


   ELSE IF ( (Xang(i) > 0.d0) .AND. (Yang(i) < 0.d0) ) THEN

    PHI(i) = ATan(Yang(i)/Xang(i)) + (3.d0*Pi/2.d0)

   END IF


 ELSE IF (Xang(i) /= 0.d0 .AND. Yang(i) /= 0.d0 .AND. Zang(i) /= 0.d0) THEN

   IF ( (Xang(i) > 0.d0) .AND. (Yang(i) > 0.d0) ) THEN

    PHI(i) = ATan(Yang(i)/Xang(i))

   ELSE IF ( (Xang(i) < 0.d0) .AND. (Yang(i) > 0.d0) ) THEN

    PHI(i) = ATan(Yang(i)/Xang(i)) + (Pi/4.d0)


   ELSE IF ( (Xang(i) < 0.d0) .AND. (Yang(i) < 0.d0) ) THEN

    PHI(i) = ATan(Yang(i)/Xang(i)) + (Pi/2.d0)


   ELSE IF ( (Xang(i) > 0.d0) .AND. (Yang(i) < 0.d0) ) THEN

    PHI(i) = ATan(Yang(i)/Xang(i)) + (3.d0*Pi/2.d0)

   END IF

   IF (Zang(i) > 0.d0) THEN

    THETA(i) = ACos(Abs(Zang(i)))

   ELSE 

    THETA(i) = ACos(Abs(Zang(i))) + Pi

   END IF


 END IF

END DO


! This parameter is defined based on Rm. It is set to avoid very small values of the density
! for the DCF project. However I do not think this is relevant for the KEF project and then
! it can be modified at will.
Rmax = (4.d0*Rm)/(Nr**2.d0)

! The arrays R and Wr are allocated to have dimension Nr, the number of radial points.
ALLOCATE(R(1:Nr),Wr(1:Nr))


! The radial grid is generated
! This grid can be find in the following reference:
!
!
DO i = 1, Nr
 R(i) = Rmax*(i**2.d0)*(((Nr+1.d0) + 1.d0 - i)**(-2.d0))
 Wr(i) = 2.d0*(Rmax**3.d0)*((Nr+1.d0) + 1.d0)*(i**5.d0)*(((Nr+1.d0) + 1.d0 - i)**(-7.d0))
END DO

! Since we are doing tensor product of radial and angular quadratures, the total number of points
! is just the product of the angular and radial ones.
Ntot = Nr*Nang

! The arrays Wtot, Xtot, Ytot, Ztot and NPointTot are allocated to have dimension Ntot, the total number
! of points
! Those arrays are, respectively, the weight, r,theta,phi coordiantes of the whole grid and the labels
! of the grid points
ALLOCATE(Wtot(1:Ntot),Xtot(1:Ntot),Ytot(1:Ntot),Ztot(1:Ntot),NPointTot(1:Ntot))

!Auxiliary counter to keep track of the total number of grid points.
!The counter is set to 1 at first
k = 1
! While the counter is less than or equal to the total number of points, the tensor product
! rule is done
DO WHILE (k <= Ntot)

 DO i = 1, Nr
  !The following DO sequence makes the tensor product of one radial point with all the angular ones
  DO j = 1, Nang
	Wtot(k) = Wr(i)*Wang(j) ! The weigths of the complete gird is just the product of the radial
							! and angular ones.
	Xtot(k) = R(i)			! The r coodiante is the one obtained from the radial grid (remember that
							! the Lebedev grid is constructed in the unit sphere
	Ytot(k) = THETA(j)		
	Ztot(k) = PHI(j)
	!By increasing this counter the next radial point is obtained to perform the previous
	!tensor product of a radial point with all the angular ones.
	k = k + 1
  END DO
 END DO

END DO


!! Open the file where the complete grid is going to be printed out on. 
!! The name of the file is "grid_points.txt" and is declared as UNIT 20
OPEN(UNIT=20,STATUS='unknown',FILE='grid_points.txt',ACTION='write',IOSTAT=ierror)
 ! For all the grid points do
 DO i = 1, Ntot
   !Writes on the file the label of the point, the weight, the r , theta and phi coordinate. This is
   ! j, Wj, Rj, THETAj, PHIj
   WRITE(20,5) i, Wtot(i), Xtot(i), Ytot(i), Ztot(i)
   ! Format in which the grid points are verted in the file.
   5 FORMAT (I7,2E27.17,2E27.17,2E27.17,2E27.17)
 END DO
CLOSE(20)
!The file is closed.


! A helpful message
write(*,*) 'The grid file has been created!'
write(*,*) '												  '


! All the previously allocated arrays are deallocated
DEALLOCATE(R,Wr)
DEALLOCATE(Wtot,Xtot,Ytot,Ztot,NPointTot)



END SUBROUTINE GetGrid_File



!========================================================================
!========================================================================

! SUBROUTINE Grid_Points(NoPoints,Wtot,Xtot,Ytot,Ztot)
! This subroutine reads the grid points of the complete grid from the file
! "grid_points.txt" that has been identifed as UNIT 20
!

SUBROUTINE Grid_Points(NoPoints,Wtot,Xtot,Ytot,Ztot)


IMPLICIT NONE


INTEGER, INTENT(IN) :: NoPoints				!(input) Number of grid points in the grid
REAL(dp), INTENT(OUT) :: Wtot(1:NoPoints)	!(output) Array where the weights of the grid are storaged
REAL(dp), INTENT(OUT) :: Xtot(1:NoPoints)	!(output) Array where the r coordiante of the grid are storaged	
REAL(dp), INTENT(OUT) :: Ytot(1:NoPoints)	!(output) Array where the theta coord of the grid are storaged
REAL(dp), INTENT(OUT) :: Ztot(1:NoPoints)	!(output) Array where the phi coord of the grid are storaged

REAL(dp) :: NumPoint(1:NoPoints)	! Array that contains the labels for the grid points
INTEGER :: i		! Just a counter
INTEGER ::  ierror	! integer of status of files


!! Open the file where the complete grid has been previously printed out.
!! The name of the file is "grid_points.txt" and is declared as UNIT 20
OPEN(UNIT=20,STATUS='old',FILE='grid_points.txt',ACTION='read',IOSTAT=ierror)
 DO i = 1, NoPoints
   READ(20,*) NumPoint(i), Wtot(i), Xtot(i), Ytot(i), Ztot(i)
 END DO
CLOSE(20)
!The file is closed.


END SUBROUTINE Grid_Points



!==============================================================================
!==============================================================================



END MODULE GridSubroutines