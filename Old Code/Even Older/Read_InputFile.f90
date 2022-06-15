! This module contains the SUBROUTINE Input_File(Nr,Rm,Lebedev_order,Z_Atomic)
! It gets the parameters to construct the grid of interest.
! The parameters are read from the UNIT=5 file InputFile.inp


MODULE InputFile

USE nrtype


CONTAINS


SUBROUTINE Input_File(Nr,Rm,Lebedev_order,Z_Atomic)


IMPLICIT NONE

INTEGER, INTENT(OUT) :: Nr				! (output) Number of radial points in the grid
REAL(dp), INTENT(OUT) :: Rm				! (output) Parameter to generate the radial grid
										! This parameter is atom dependent
INTEGER, INTENT(OUT) :: Lebedev_order	! (output) Number of angual points in the Lebedev grid
INTEGER, INTENT(OUT) :: Z_Atomic		! (output) Atomic number 
CHARACTER(len=1) :: comment				! (internal) Just a character variable to identify comments in the file


! The file is open
OPEN(UNIT=5, STATUS='OLD',FILE='InputFile.inp',ACTION='READ')

DO
   READ(5,*) comment

   IF (comment /= '!')THEN !This line contains data, not a comment!

      backspace(5)
      READ(5,*) Nr		   ! The number of radial points is read

      EXIT

   ENDIF
ENDDO




DO
   READ(5,*) comment

   IF (comment /= '!')THEN !This line contains data, not a comment!

      backspace(5)
      READ(5,*) Rm			! The parameter for creating the radial grid is read			
							! Recommended values for each of the corresponding atoms are
							! included in the input file itself
      EXIT

   ENDIF
ENDDO


DO
   READ(5,*) comment

   IF (comment /= '!')THEN !This line contains data, not a comment!

      backspace(5)
      READ(5,*) Lebedev_order  ! Number of angular grids in the Lebedev grid
							   ! The code permits to choose Lebedev grids with 38, 50, 86, 110,
							   ! 170, 350 and 1202 angular points

      EXIT

   ENDIF
ENDDO



DO
   READ(5,*) comment

   IF (comment /= '!')THEN !This line contains data, not a comment!

      backspace(5)
      READ(5,*) Z_Atomic		! The atomic number is read

      EXIT

   ENDIF
ENDDO

! The file is closed.
CLOSE(5)


END SUBROUTINE Input_File


END MODULE InputFile
