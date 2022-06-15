!****************************************************************************
!  	Program written by Dieter Ghillemijn and Patrick Bultinck
!     	in cooperation with (mainly) Dimitri Van Neck and Paul W. Ayers
!     	www.quantum.ugent.be
!****************************************************************************

!-------------------------------------------------------------------------------!
!                                   variables_wfn                               !
!                                                                               !
! This module contains the main variables that are contained in a .wfn file.    !
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                           DICTIONARY OF VARIABLES                             !
!-------------------------------------------------------------------------------!
! atm_charge(1:n_atoms) -- the atomic numbers of the atoms.                     !
! atm_position(1:3,1:n_atoms) -- the position of each atom.                     !
! basis_center(1:n_prim) -- the center to which the gaussian primitive is       !
!                           assigned.                                           !
! basis_exp(1:n_prim) -- exponents of the gaussian primitives.                  !
! basis_type(1:n_prim) -- the type of gaussian primitive.                       !
! basis_ijk(1:3,1:nprim) -- assigns exponents to each Gaussian primitive        !
!                              using type_exponents.                            !
! multiplicity -- the spin-multiplicity.                                        !
! mo_energy(1:n_orbs,1:2) -- molecular orbital energies for 1. alpha and 2.     !
!                            beta-spin orbitals                                 !
! mo_coeff(1:n_prim,1:n_orb,1:2) -- molecular orbital coefficients of the spin  !
!                                   molecular orbitals (1 = alpha; 2 = beta)    !
! mo_occ(1:n_orbs,1:2) -- molecular orbital occupation numbers for orbitals of  !
!                        each spin.                                             !
! n_atoms -- the number of atoms.                                               !
! n_el -- the number of electrons                                               !
! n_orbs -- the number of molecular orbitals read in.                           !
! n_prim -- the number of Gaussian primitives.                                  !
! Energy -- the energy of the system.                                           !
! virial -- the virial ration, -V/T.                                            !
! KE -- the kinetic energy.                                                     !
! spinorbs -- .true. if we have spin orbitals. .false. if only spatial orbitals.!
!-------------------------------------------------------------------------------!

MODULE variables_wfn

USE kinds

IMPLICIT NONE

SAVE

INTEGER(istd)   :: n_orbs,n_prim,n_atoms,n_el,multiplicity
INTEGER(istd), ALLOCATABLE :: basis_center(:), basis_type(:)
INTEGER(istd), ALLOCATABLE :: basis_ijk(:,:)
INTEGER(istd), ALLOCATABLE :: atm_prim(:,:),Zatom(:)

REAL(dbl) :: Energy,virial,KE
REAL(dbl), ALLOCATABLE :: mo_occ(:,:),mo_coeff(:,:,:),mo_energy(:,:)
REAL(dbl), ALLOCATABLE :: atm_position(:,:),atm_charge(:),basis_exp(:)

LOGICAL :: spinorbs
 
CHARACTER(len=80) :: wfn_filename

END MODULE variables_wfn