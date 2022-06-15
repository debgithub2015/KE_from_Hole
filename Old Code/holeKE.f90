!-------------------------------------------------------------------------------!
!                               holeKE                                          !
!                                                                               !
! This program is designed to compute the kinetic energy from an approximate    !
! density matrix and exchange hole.                                             !
!                                                                               !
! Given the density, the program first computes several standard k.e. formulae, !
! including Thomas-Fermi and Weizsacker, and some formulas we derived that      !
! are related to these.  Next, it normalizes the exchange-correlation hole      !
! either (1) approximately or (2) exactly.  These give new k.e. fucntionals.    !
! Finally, it computes the eigenvalues of the density matrix.                   !
!                                                                               !
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
USE KE_module
USE variables_wfn, ONLY: wfn_filename

!Open the output file
OPEN(UNIT=101,STATUS='UNKNOWN',FILE='general.out')
OPEN(UNIT=102,STATUS='UNKNOWN',FILE='table.out')

!Read the name of the molecule
call getarg(1,wfn_filename)
wfn_filename = TRIM(ADJUSTL(wfn_filename))//'.wfn'

!Start up the calculation.  This:
!   1.  Reads in the parameters that define the grid.
!   2.  Reads the *.wfn file.  This gives the molecular geometry and 
!       exact noninteracting kinetic energy.
!   3.  Generates the Becke-Lebedev grid
!   4.  Evaluates the molecular electron density and density gradient.
call startup()

!compute "standard" k.e. formulae.
call compute_ke_conventional()

!determine approximate and exact hole normalization.
call hole_normalize()

!compute improved k.e. formulae
call compute_ke_improved()

!print output.  The eigenvalues of the effective density matrix can be 
!computed here, if desired.
call write_output()

CLOSE(101); CLOSE(102)

END PROGRAM hole KE

