!-------------------------------------------------------------------------------!
!                                   process_wfn                                 !
!                                                                               !
! This module contains the main subroutines that are used to access and use a   !
! wfn file.                                                                     !
!                                                                               !
! INCLUDED SUBROUTINES:                                                         !
!       read_wfn -- reads a wfn file.                                           !
!       density_wfn -- gives the density at a set of grid points.               !
!       densitygradient_wfn -- gives the gradient of the electron density and   !
!                              the electron density.                            !
!       densityders_wfn -- gives the derivatives of the density at a set of     !
!                          grid points.                                         !
!       KEdensity_wfn -- gives the kinetic energy density at a set of grid pts  !
!       stress_wfn -- gives the stress tensor and tension force at a set of     !
!                     grid points.                                              !
!       orbs_wfn -- gives the value of a set of MOs at a set of grid points.    !
!       wavefunction_wfn -- gives the value of the wavefunction at a set of     !
!                           grid points.  Only for a Slater determinant.        !
!       DM1_wfn -- gives the values of the density matrix at a set of grid      !
!                  points.                                                      !
!       dDM1_wfn -- gives the gradient of the density matrix at a set of grid   !
!                   points.                                                     !
!       d2DM1_wfn -- gives the 2nd derivative information about the 1-electron  !
!                   reduced density matrix at a set of grid points.             !
!       d3DM1_wfn -- gives the 3rd derivative information about the 1-electron  !
!                   reduced density matrix at a set of grid points.             !
!       DM2_wfn -- gives the 2-RDM at a set of grid points.  Only for a Slater  !
!                  determinant.                                                 !
!       hx_wfn -- the exchange hole.                                            !
!       localIP_wfn -- the local ionization potential.                          !
!       KSresponse_wfn -- the Kohn-Sham response function.                      !
!-------------------------------------------------------------------------------!
!*******************************************************************************!
! There are several important IOPs to be used in Gaussian.                      !
!   1.  Make a wfn file!!                                                       !
!   2.  If you are doing an ab initio calculation and you want the ab initio    !
!       density instead of Hartree-Fock, remember to specify the keyword        !
!               DENSITY=CURRENT                                                 !
!   3.  Use IOp(99/6=100) to make a wfn file.  (Can also be done with keyword)  !
!   4.  Use IOp(99/6=1000) to use natural orbitals in the wfn file.             !
!   5.  OR use IOP(99/10=256) to print natural orbitals in the wfn file.        !
!   6.  Use IOp(99/18=-1) to print all virtual orbitals in the wfn file.        !
!                                                                               !
!*******************************************************************************!
!-------------------------------------------------------------------------------!
!                            DICTIONARY OF VARIABLES                            !
!                      used repeatedly*in this module.                          !
!                                                                               !
! *some other variables are "local" to individual subroutines.                  !
!                                                                               !
! alpha,beta -- terms that define the stress tensor.                            !
!               The stress tensor definition has a characteristic element:      !
!                 -1/2*(alpha*(d/dx d/dy' + d/dx' d/dy)                         !
!                       - (1-alpha)(d/dx d/dy + d/dx' d/dy'))*D(r,r')           !
!                 - 1/2*delta(i,j)*beta*Laplacian(rho(r))                       !
! i -- counter.                                                                 !
! n_pts -- the number of points at which the density should be evaluated.       !
! XYZ(1:3,1:npts) -- the location of the point in Cartesian coordinates.        !
! R6d(1:6,1:npts) -- the 6-dimensional diagonal coordinate of the 1-electron    !
!                    reduced density matrix.                                    !
! R12d(1:6,1:npts) -- the 12-dimensional diagonal coordinate of the 2-electron  !
!                     reduced density matrix.                                    !
! rho(1:npts,1:2) -- the spin-density at the points.                            !
! drho(1:3,1:npts,1:2) -- the gradient of the spin-density at the points.       !
! d2rho(1:6,1:npts,1:2) -- the second derivative matrix of the spin-density at  !
!                          the input points.                                    !
!                         1 -- d^2/dx^2 rho(r)                                  !
!                         2 -- d^2/dy^2 rho(r)                                  !
!                         3 -- d^2/dz^2 rho(r)                                  !
!                         4 -- d/dx d/dy rho(r)                                 !
!                         5 -- d/dx d/dz rho(r)                                 !
!                         6 -- d/dy d/dz rho(r)                                 !
! d3rho(1:10,1:npts,1:2) -- the third derivative matrix of the spin-density at  !
!                           the points                                          !
!                         1 -- d^3/dx^3 rho(r)                                  !
!                         2 -- d^3/dy^3 rho(r)                                  !
!                         3 -- d^3/dz^3 rho(r)                                  !
!                         4 -- d^2/dx^2 d/dy rho(r)                             !
!                         5 -- d^2/dx^2 d/dz rho(r)                             !
!                         6 -- d^2/dy^2 d/dx rho(r)                             !
!                         7 -- d^2/dy^2 d/dz rho(r)                             !
!                         8 -- d^2/dz^2 d/dx rho(r)                             !
!                         9 -- d^2/dz^2 d/dy rho(r)                             !
!                        10 -- d/dx d/dy d/dz rho(r)                            !
! Lap(1:npts,1:2) -- the Laplacian of the electron density at the points.  The  !
!                    trace of the second derivative matrix.                     !
! stress(1:6,1:npts,1:2) -- the electronic ("kinetic energy") stress tensor at  !
!                           the points.  The elements are ordered in the same   !
!                           as the second derivative matrix, d2rho(1:6,:,:)     !
! stressForce(1:3,1:npts,1:2) -- the spin-resolved tension (or Ehrenfest) force !
!                                at the points. This is the divergence of the   !
!                                stress tensor.                                 !
! KEdensity -- the kinetic energy density that corresponds to this stress       !
!              tensor. Equal to -1/2*Tr[stress]                                 !
! tpos(1:npts,1:2) -- the positive-definite kinetic energy density at the       !
!                     points.                                                   !
! DM(1:npts,1:2) -- the spin-resolved one-electron reduced density matrix.      !
! DM2nd(1:npts,1:2) -- the spin-resolved two-electron reduced density matrix.   !
! hx(1:npts,1:2) -- the spin-resolved exchange hole.                            !
! localIP(1:npts,1:2) -- the spin-resolved local ionization potential.          !
! PHIa(1:nlist,1:npts) -- alpha-spin molecular orbitals at the grid points.     !
! PHIb(1:nlist,1:npts) -- beta-spin MOs at the grid points.                     !
! KSresponse(1:npts,1:2) -- the Kohn-Sham response at a set of grid points.     !
! dDM(1:3,1:npts,1:2) -- the spin (1:2) components of the gradient of the first-!
!                        order density matrix at the given grid points (1:npts).!
!                        The components are:                                    !
!                         1 -- d/dx D(r',r)                                     !
!                         2 -- d/dy D(r',r)                                     !
!                         3 -- d/dz D(r',r)                                     !
! d2DM(1:12,1:npts,1:2) -- the spin (1:2) components of the second derivatives  !
!                         of the first-order density matrix at the given        !
!                         grid points (1:npts).  The components are:            !
!                         1 -- d^2/dx^2 D(r',r)                                 !
!                         2 -- d^2/dy^2 D(r',r)                                 !
!                         3 -- d^2/dz^2 D(r',r)                                 !
!                         4 -- d/dx d/dy D(r',r)                                !
!                         5 -- d/dx d/dz D(r',r)                                !
!                         6 -- d/dy d/dz D(r',r)                                !
!                         7 -- d/dx d/dx' D(r',r)                               !
!                         8 -- d/dy d/dy' D(r',r)                               !
!                         9 -- d/dz d/dz' D(r',r)                               !
!                        10 -- d/dx d/dy' D(r',r)                               !
!                        11 -- d/dx d/dz' D(r',r)                               !
!                        12 -- d/dy d/dz' D(r',r)                               !
! d3DM(1:28,1:npts,1:2) -- the spin (1:2) components of the third derivatives   !
!                         of the first-order density matrix at the given        !
!                         grid points (1:npts).  The components are:            !
!                         1 -- d^3/dx^3 D(r',r)                                 !
!                         2 -- d^3/dy^3 D(r',r)                                 !
!                         3 -- d^3/dz^3 D(r',r)                                 !
!                         4 -- d^2/dx^2 d/dy D(r',r)                            !
!                         5 -- d^2/dx^2 d/dz D(r',r)                            !
!                         6 -- d^2/dy^2 d/dx D(r',r)                            !
!                         7 -- d^2/dy^2 d/dz D(r',r)                            !
!                         8 -- d^2/dz^2 d/dx D(r',r)                            !
!                         9 -- d^2/dz^2 d/dy D(r',r)                            !
!                        10 -- d/dx d/dy d/dz D(r',r)                           !
!                        11 -- d^2/dx^2 d/dx' D(r',r)                           !
!                        12 -- d^2/dy^2 d/dy' D(r',r)                           !
!                        13 -- d^2/dz^2 d/dz' D(r',r)                           !
!                        14 -- d^2/dx^2 d/dy' D(r',r)                           !
!                        15 -- d^2/dx^2 d/dz' D(r',r)                           !
!                        16 -- d^2/dy^2 d/dx' D(r',r)                           !
!                        17 -- d^2/dy^2 d/dz' D(r',r)                           !
!                        18 -- d^2/dz^2 d/dx' D(r',r)                           !
!                        19 -- d^2/dz^2 d/dy' D(r',r)                           !
!                        20 -- d/dx d/dy d/dx' D(r',r)                          !
!                        21 -- d/dx d/dz d/dx' D(r',r)                          !
!                        22 -- d/dy d/dx d/dy' D(r',r)                          !
!                        23 -- d/dy d/dx d/dy' D(r',r)                          !
!                        24 -- d/dz d/dx d/dz' D(r',r)                          !
!                        25 -- d/dz d/dy d/dz' D(r',r)                          !
!                        26 -- d/dx' d/dy d/dz D(r',r)                          !
!                        27 -- d/dx d/dy' d/dz D(r',r)                          !
!                        28 -- d/dx d/dy d/dz' D(r',r)                          !
!-------------------------------------------------------------------------------!
!                           ITEMS READ FROM *.WFN FILE                          !
!                                                                               !
! atm_charge(1:n_atoms) -- the atomic numbers of the atoms.                     !
! atm_position(1:3,1:n_atoms) -- the position of each atom.                     !
! basis_center(1:n_prim) -- the center to which the gaussian primitive is       !
!                           assigned.                                           !
! basis_exp(1:n_prim) -- exponents of the gaussian primitives.                  !
! basis_type(1:n_prim) -- the type of gaussian primitive.                       !
! basis_ijk(1:3,1:nprim) -- assigns exponents to each Gaussian primitive        !
!                              using type_exponents.                            !
! multiplicity -- the spin-multiplicity.                                        !
! mo_energy(1:n_orbs,1:2) -- molecular orbital energies for 1. alpha and 2. beta-!
!                           spin orbitals                                       !
! mo_coeff(1:n_prim,1:n_orbs,1:2) -- molecular orbital coefficients of the spin  !
!                                   molecular orbitals (1 = alpha; 2 = beta)    !
! mo_occ(1:n_orbs,1:2) -- molecular orbital occupation numbers for orbitals of   !
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


MODULE process_wfn

USE kinds
USE variables_wfn

CONTAINS

!-------------------------------------------------------------------------------!
!                           read_wfn                                            !
!                                                                               !
! This subroutines reads a *.wfn file and stores the data in the appropriate    !
! variables.                                                                    !
!-------------------------------------------------------------------------------!
!                       DICTIONARY OF LOCAL VARIABLE                            !
!                                                                               !
! line -- this reads in the lines from the *.wfn file.  The line is then        !
!         parsed for content.                                                   !
! j -- counters                                                                 !
! ifile -- to find end of file.                                                 !
! iorb -- counts the orbitals.                                                  !
! ibss -- counts the primitives.                                                !
! nlines -- counts the number of lines to be read.                              !
! ndata -- the number of data items to be read on the line.                     !
! label -- the label for the line.  Used to debug and head off errors.          !
! orb_occ -- a temporary variable holding the orbital occupation # for the 1st  !
!            MO.                                                                !
! orb_energy -- a temporary variable holding the orbital energy for the 1st MO. !
! type_exponents(1:3,1:20) -- the values of (i,j,k) for the Cartesian           !
!                         prefactor, x**i * y**j * z**k, of the Gausian         !
!                         function in the primitive basis functions.            !
!-------------------------------------------------------------------------------!
!                       VARIABLES from variables_wfn MODULE                     !
!                                                                               !
! atm_charge(1:n_atoms) -- the atomic numbers of the atoms.                     !
! atm_position(1:3,1:n_atoms) -- the position of each atom.                     !
! basis_center -- the center to which the gaussian primitive is assigned.       !
! basis_exp -- exponents of the gaussian primitives.                            !
! basis_ijk(1:3,1:nprim) -- assigns exponents to each Gaussian primitive        !
!                              using type_exponents.                            !
! basis_type -- the type of gaussian primitive.                                 !
! type_exponents(1:3,1:20) -- the values of (i,j,k) for the Cartesian           !
!                         prefactor, x**i * y**j * z**k, of the Gausian         !
!                         function in the primitive basis functions.            !
! multiplicity -- the spin-multiplicity.                                        !
! mo_energy -- molecular orbital energy                                         !
! mo_coeff -- molecular orbital coefficients                                    !
! mo_occ -- molecular orbital occupation numbers                                !
! n_atoms -- the number of atoms.                                               !
! n_el -- the number of electrons                                               !
! n_orbs -- the number of molecular orbitals read in.                           !
! n_prim -- the number of Gaussian primitives.                                  !
! Energy -- the energy of the system.                                           !
! virial -- the virial ration, -V/T.                                            !
! KE -- the kinetic energy.                                                     !
! spinorbs -- logical variable that tells us whether we have spin orbitals or   !
!             only spatial orbitals.                                            !
!-------------------------------------------------------------------------------!

subroutine read_wfn

USE utilities_wfn

IMPLICIT NONE

CHARACTER(len=80) :: line
INTEGER(istd) :: j,ifile,iorb,ibss,iatm,ndata,nlines
CHARACTER(len=3) :: label
REAL(dbl) :: orb_occ,orb_energy,type_exponents(1:3,1:20)

!Open the wfn file.
OPEN(unit=10,FILE=wfn_filename,STATUS='old',IOSTAT=ifile)

!The next do loop just reads finds the first "content" line in the *.wfn file.
!The rest of the file is not read with the loop; it is read as fixed-format
!instead.
DO
   !Read the next line in the *.wfn file.
   READ(10,'(A80)',IOSTAT=ifile) line

   !The first line of the file with "content" starts with 'GAUSSIAN'
   IF (line(1:8) == 'GAUSSIAN') THEN
      READ (line,'(18X,I5)') n_orbs
      READ(line,'(38X,I5)') n_prim
      READ(line,'(59X,I4)') n_atoms
      !Now we read the rest of the wfn file.
      EXIT
   ENDIF
   
   !Control Statement to end the loop if not GAUSSIAN line is found.
   IF (ifile /= 0) THEN
      WRITE(*,*) "Unexpected format for WFN file."
      WRITE(*,*) "You probably typed in the name of your wfn file incorrectly."
      STOP
   ENDIF
   
ENDDO
   
!The rest of this file reads the "other data" in the *.wfn file.
   
   !Allocate arrays for defining the basis functions.
   ALLOCATE(basis_exp(1:n_prim),basis_center(1:n_prim),basis_type(1:n_prim))
   ALLOCATE(basis_ijk(1:3,1:n_prim))
   ALLOCATE(atm_charge(1:n_atoms),atm_position(1:3,1:n_atoms),Zatom(1:n_atoms))
   
   !The next items in the file are the positions and types of the atoms.
   DO iatm=1,n_atoms
      READ(10,'(24X,3(F12.8),10X,F5.1)') (atm_position(j,iatm),j=1,3), atm_charge(iatm)
   ENDDO
   
   !Gaussian uses real numbers for atomic numbers.  Let's convert to integers
   Zatom(:) = NINT(atm_charge(:))
   
   !The next thing should be the atomic centers.  These are printed out 20
   !per line.
   nlines = n_prim/20 
   ndata = mod(n_prim,20)
   IF (ndata == 0) THEN
      !Actually we have a full line of data to read.
      ndata = 20
      nlines = nlines - 1
   ENDIF
      
   DO ibss = 0,nlines-1
      !This loop reads in the "complete" lines.
      READ(10,'(A3,17X,20(I3))') label,(basis_center(20*ibss+j),j=1,20)
      IF (label /= 'CEN') THEN
         WRITE(*,*) 'Error reading atomic center assigned to basis functions.'
         WRITE(*,*) 'Expected CEN but got a different label instead ',label
         STOP
      ENDIF
   ENDDO
   READ(10,'(A3,17X,20(I3))') label,(basis_center(20*nlines+j),j=1,ndata)
   IF (label /= 'CEN') THEN
      WRITE(*,*) 'Error reading atomic center assigned to basis functions.'
      WRITE(*,*) 'Expected CEN but got a different label instead ',label
      STOP
   ENDIF   
   
   !Next comes the "type assignments."  This tells us the angular momentum factor
   !for each primitive.
   DO ibss = 0,nlines-1
      !This loop reads in the "complete" lines.
      READ(10,'(A3,17X,20(I3))') label,(basis_type(20*ibss+j),j=1,20)
      IF (label /= 'TYP') THEN
         WRITE(*,*) 'Error reading types of basis functions.'
         WRITE(*,*) 'Expected TYP but got a different label instead ',label
         STOP
      ENDIF
   ENDDO
   READ(10,'(A3,17X,20(I3))') label,(basis_type(20*nlines+j),j=1,ndata)  
   IF (label /= 'TYP') THEN
      WRITE(*,*) 'Error reading types of basis functions.'
      WRITE(*,*) 'Expected TYP but got a different label instead ',label
      STOP
   ENDIF   

   !The exponents of the Gaussians in the basis set are printed out 5 per line.
   nlines = n_prim/5
   ndata = mod(n_prim,5)
   IF (ndata == 0) THEN
      !Actually we have a full line of data to read.
      ndata = 5
      nlines = nlines - 1
   ENDIF
   DO ibss = 0,nlines-1
      !This loop reads in the "complete" lines.
      READ(10,'(A3,7X,5(D14.7))') label,(basis_exp(5*ibss+j),j=1,5)
      IF (label /= 'EXP') THEN
         WRITE(*,*) 'Error reading exponents of basis functions.'
         WRITE(*,*) 'Expected EXP but got a different label instead ',label
         STOP
      ENDIF
   ENDDO
   READ(10,'(A3,7X,5(D14.7))') label,(basis_exp(5*nlines+j),j=1,ndata)
   IF (label /= 'EXP') THEN
      WRITE(*,*) 'Error reading exponents of basis functions.'
      WRITE(*,*) 'Expected EXP but got a different label instead ',label
      STOP
   ENDIF   

   !Next come the molecular orbitals.  This is tricky because we might
   !have spin-orbitals or spatial orbitals, and we want to be able to tell the
   !difference.  
   spinorbs = .true.  !default is to have spin orbitals.
   IF (MOD(n_orbs,2) == 1) THEN
      !We have an odd number of orbitals, so obviously we aren't printing out
      !the spin orbitals.
      spinorbs = .false.
   ELSE
      !It is still possible that we just have an even number of spatial orbitals.
      !We read the first line of the MO orbitals.
      READ(10,'(A3,31X,D13.7,15X,D12.6)') label,orb_occ,orb_energy
      
      !If the occupation number is greater than 1, then it is safe to say that
      !we have spatial orbitals.
      IF (orb_occ > 1.01) THEN
         spinorbs = .false.
      ENDIF
      BACKSPACE(10)  !We will read this line again!
   ENDIF
   
   !Now we can allocate the arrays for the MOs.
   IF (spinorbs) THEN
      !The number of spatial orbitals is 1/2 the total number of orbitals.
      n_orbs = n_orbs/2
   ENDIF

   ALLOCATE(MO_coeff(1:n_prim,1:n_orbs,1:2),MO_energy(1:n_orbs,1:2),           &
            MO_occ(1:n_orbs,1:2))

   !MO coefficients are printed out 5 per line.
   nlines = n_prim/5
   ndata = mod(n_prim,5)
   IF (ndata == 0) THEN
      !Actually we have a full line of data to read.
      ndata = 5
      nlines = nlines - 1
   ENDIF

   !Now we read the MOs.
   DO iorb = 1,n_orbs
      !Read orbital occupation number and orbital energy
      READ(10,'(A3,31X,D13.7,15X,D12.6)') label,MO_occ(iorb,1),MO_energy(iorb,1)
      IF (label /= 'MO ') THEN
         WRITE(*,*) 'Error reading molecular orbitals from wfn file.'
         WRITE(*,*) 'Expected MO[space] but got a different label instead ',label
         STOP
      ENDIF
      
      !Read MO coefficients, 5 per line.
      DO ibss = 0,nlines-1
         !This loop reads in the "complete" lines.
         READ(10,'(5(D16.8))') (MO_coeff(5*ibss+j,iorb,1),j=1,5)
      ENDDO
      READ(10,'(5(D16.8))') (MO_coeff(5*nlines+j,iorb,1),j=1,ndata)
   
   ENDDO      
   
   !If we have only spatial MOs, then the other spin-MOs are the same and
   !the occupation numbers should all be halved.
   IF (.not. spinorbs) THEN
      !halve the occupation numbers
      MO_occ(:,2) = MO_occ(:,1)/2
      MO_occ(:,1) = MO_occ(:,2)
      !copy the orbitals and orbital energies over
      MO_coeff(:,:,2) = MO_coeff(:,:,1)
      MO_energy(:,2) = MO_energy(:,1)
   ELSE
      !We need to read in the second set of spin orbitals.
      !Now we read the beta-spin MOs.
      DO iorb = 1,n_orbs
         !Read orbital occupation number and orbital energy
         READ(10,'(A3,31X,D13.7,15X,D12.6)') label,MO_occ(iorb,2),MO_energy(iorb,2)
         IF (label /= 'MO ') THEN
            WRITE(*,*) 'Error reading molecular orbitals from wfn file.'
            WRITE(*,*) 'Expected MO[space] but got a different label instead ',label
            STOP
         ENDIF
         
         !Read MO coefficients, 5 per line.
         DO ibss = 0,nlines-1
            !This loop reads in the "complete" lines.
            READ(10,'(5(D16.8))') (MO_coeff(5*ibss+j,iorb,2),j=1,5)
         ENDDO
         READ(10,'(5(D16.8))') (MO_coeff(5*nlines+j,iorb,2),j=1,ndata)
      ENDDO 
   ENDIF
   
   !Now we read in the energy and the virial ratio.
   READ(10,'(A3)') label
   IF (label /= 'END') THEN
      WRITE(*,*) 'There is something fishy here.  We should be done reading'
      WRITE(*,*) 'the MO coefficients but it appears that we are not.'
      WRITE(*,*) 'Expected END but got a different label instead ',label
   ENDIF
   READ(10,'(A3)') label
   BACKSPACE(10)
   IF (label == ' TH') THEN
      READ(10,'(17X,G20.12,18X,G13.8)') energy,virial
   ELSE IF (label == ' TO') THEN
      READ(10,'(17X,G18.12,18X,G13.8)') energy,virial
   ELSE
      WRITE(*,*) 'Error reading the energy and virial.'
      WRITE(*,*) 'Expected [space]TH but got a different label instead',label
      STOP 'error in wfn read'
   ENDIF
   !Compute kinetic energy.  -V/T = virial. T + V = energy
   ! T = energy/(1-virial)
   KE = energy/(1-virial)

!This is the last line in the .wfn file
close(10)
!The *.wfn file is now complete finished.


!Decipher the right "types" for each basis function
call primitive_type_assignments(type_exponents)

FORALL(ibss=1:n_prim)
      basis_ijk(1:3,ibss) = type_exponents(1:3,basis_type(ibss))
ENDFORALL

!Determine the multiplicity and number of electrons
n_el = NINT(SUM(MO_occ))

! multiplicity = 1 + |Nalpha - Nbeta|
multiplicity = 1 + NINT( ABS( SUM(MO_occ(:,1)) - SUM(MO_occ(:,2)) ) )

end subroutine read_wfn

!-------------------------------------------------------------------------------!
!                           densitygradient_wfn                                 !
!                                                                               !
! This subroutine outputs the density and the gradient of the density at a      !
! select group of points.                                                       !
!-------------------------------------------------------------------------------!


subroutine densitygradient_wfn(rho,drho,XYZ,npts)

IMPLICIT NONE

INTEGER(istd) :: npts
REAL(dbl)     :: XYZ(1:3,1:npts),rho(1:npts,1:2),drho(1:3,1:npts,1:2)
REAL(dbl)     :: R6d(1:6,1:npts)

R6d(1:3,:) = XYZ(1:3,:)
R6d(4:6,:) = XYZ(1:3,:)

call density_wfn(rho,XYZ,npts)
call dDM1_wfn(drho,R6d,npts)

drho = 2*drho

end subroutine densitygradient_wfn

!-------------------------------------------------------------------------------!
!                           densityders_wfn                                     !
!                                                                               !
! This subroutine outputs the density and the first three derivatives of the    !
! density at a set of grid points.                                              !
!-------------------------------------------------------------------------------!

subroutine densityders_wfn(rho,drho,d2rho,d3rho,XYZ,npts)

IMPLICIT NONE

INTEGER(istd) :: npts,i
REAL(dbl)     :: XYZ(1:3,1:npts),rho(1:npts,1:2),drho(1:3,1:npts,1:2)
REAL(dbl)     :: d2rho(1:6,1:npts,1:2),d3rho(1:10,1:npts,1:2)
REAL(dbl)     :: d2DM(1:12,1:npts,1:2),d3DM(1:28,1:npts,1:2)
REAL(dbl)     :: R6d(1:6,1:npts)

R6d(1:3,:) = XYZ(1:3,:)
R6d(4:6,:) = XYZ(1:3,:)

!Evaluate density
call density_wfn(rho,XYZ,npts)

!Evaluate first derivative
call dDM1_wfn(drho,R6d,npts)
drho = 2*drho

!Evaluate second derivative
call d2DM1_wfn(d2DM,R6d,npts)

FORALL(i=1:6)
      d2rho(i,:,:) = 2*(d2DM(i,:,:)+d2DM(i+6,:,:))
ENDFORALL

!Evaluate third derivative
call d3DM1_wfn(d3DM,R6d,npts)

FORALL(i=1:3)
      d3rho(i,:,:) = 2*(d3DM(i,:,:)+3*d3DM(10+i,:,:))
ENDFORALL
FORALL(i=4:9)
      d3rho(i,:,:) = 2*(d3DM(i,:,:)+d3DM(10+i,:,:)+2*d3DM(16+i,:,:))
ENDFORALL

d3rho(10,:,:) = 2*(d3DM(10,:,:) + d3DM(26,:,:) + d3DM(27,:,:) + d3DM(28,:,:))

end subroutine densityders_wfn

!-------------------------------------------------------------------------------!
!                           KEdensity_wfn                                       !
!                                                                               !
! Computes the density, gradient of the density, Laplacian of the density, and  !
! the positive-definite kinetic energy density                                  !
!-------------------------------------------------------------------------------!

subroutine KEdensity_wfn(rho,drho,Lap,tpos,XYZ,npts)

IMPLICIT NONE

INTEGER(istd) :: npts,i
REAL(dbl)     :: XYZ(1:3,1:npts),rho(1:npts,1:2),drho(1:3,1:npts,1:2)
REAL(dbl)     :: d2rho(1:6,1:npts,1:2),d3rho(1:10,1:npts,1:2)
REAL(dbl)     :: d2DM(1:12,1:npts,1:2),Lap(1:npts,1:2),tpos(1:npts,1:2)
REAL(dbl)     :: R6d(1:6,1:npts)

R6d(1:3,:) = XYZ(1:3,:)
R6d(4:6,:) = XYZ(1:3,:)

!Evaluate density
call density_wfn(rho,XYZ,npts)

!Evaluate first derivative
call dDM1_wfn(drho,R6d,npts)
drho = 2*drho

!Evaluate second derivative
call d2DM1_wfn(d2DM,R6d,npts)

!Get the Laplacian
FORALL(i=1:3)
      d2rho(i,:,:) = 2*(d2DM(i,:,:)+d2DM(i+6,:,:))
ENDFORALL
Lap(:,:) = d2rho(1,:,:) + d2rho(2,:,:) + d2rho(3,:,:)

tpos(:,:) = (d2DM(7,:,:)+d2DM(8,:,:)+d2DM(9,:,:))/2

end subroutine KEdensity_wfn

!-------------------------------------------------------------------------------!
!                           stress_wfn                                          !
!                                                                               !
! This subroutine outputs the density, gradient of the density, second          !
! derivative matrix of the density, and the KE density, stress tensor, and      !
! tension (Ehrenfest) force at a set of grid points.                            !
!    The stress tensor definition has a characteristic element:                 !
!   -1/2*(alpha*(d/dx d/dy' + d/dx' d/dy)                                       !
!        - (1-alpha)(d/dx d/dy + d/dx' d/dy'))*D(r,r')                          !
!   -1/2*delta(i,j)*beta*Laplacian(rho(r))                                      !
!                                                                               !
! See                                                                           !
!   J.S.M. Anderson, P.W.Ayers, J.I. Rodriguez Hernandez; JPCA (2010) for       !
! details.                                                                      !
!   The tension force is the divergence of the stress tensor.                   !
!-------------------------------------------------------------------------------!

subroutine stress_wfn(rho,drho,d2rho,stress,stressForce,KEdensity,Lap,    &
                           alpha,beta,XYZ,npts)

IMPLICIT NONE

INTEGER(istd) :: npts,i
REAL(dbl)     :: XYZ(1:3,1:npts),rho(1:npts,1:2),drho(1:3,1:npts,1:2)
REAL(dbl)     :: d2rho(1:6,1:npts,1:2),KEdensity(1:npts,1:2),Lap(1:npts,1:2)
REAL(dbl)     :: stress(1:6,1:npts,1:2),stressForce(1:3,1:npts,1:2)
REAL(dbl)     :: d3rho(1:10,1:npts,1:2)
REAL(dbl)     :: d2DM(1:12,1:npts,1:2),d3DM(1:28,1:npts,1:2)
REAL(dbl)     :: R6d(1:6,1:npts),alpha,beta

R6d(1:3,:) = XYZ(1:3,:)
R6d(4:6,:) = XYZ(1:3,:)

!Evaluate density
call density_wfn(rho,XYZ,npts)

!Evaluate first derivative
call dDM1_wfn(drho,R6d,npts)
drho = 2*drho

!Evaluate second derivative
call d2DM1_wfn(d2DM,R6d,npts)

FORALL(i=1:6)
      d2rho(i,:,:) = 2*(d2DM(i,:,:)+d2DM(i+6,:,:))
ENDFORALL
Lap = d2rho(1,:,:) + d2rho(2,:,:) + d2rho(3,:,:)

!Evaluate the stress tensor
FORALL(i=1:6)
      stress(i,:,:) = (1-alpha)*d2DM(i,:,:)-alpha*d2DM(i+6,:,:)
ENDFORALL
!Now we need to add the Laplacian contribution
FORALL(i=1:3)
      stress(i,:,:) = stress(i,:,:) - beta*Lap(:,:)/2
ENDFORALL

!The KE density is the trace of the stress tensor, times -1/2.  So
KEdensity(:,:) = (stress(1,:,:)+stress(2,:,:)+stress(3,:,:))/(-2)

!Evaluate third derivative
call d3DM1_wfn(d3DM,R6d,npts)

FORALL(i=1:3)
      d3rho(i,:,:) = 2*(d3DM(i,:,:)+3*d3DM(10+i,:,:))
ENDFORALL
FORALL(i=4:9)
      d3rho(i,:,:) = 2*(d3DM(i,:,:)+d3DM(10+i,:,:)+2*d3DM(16+i,:,:))
ENDFORALL

d3rho(10,:,:) = 2*(d3DM(10,:,:) + d3DM(26,:,:) + d3DM(27,:,:) + d3DM(28,:,:))

!Now we need to evaluate the tension (Ehrenfest) force.

FORALL(i=1:3)
      d3rho(i,:,:) = 2*(d3DM(i,:,:)+3*d3DM(10+i,:,:))
ENDFORALL
FORALL(i=4:9)
      d3rho(i,:,:) = 2*(d3DM(i,:,:)+d3DM(10+i,:,:)+2*d3DM(16+i,:,:))
ENDFORALL

d3rho(10,:,:) = 2*(d3DM(10,:,:) + d3DM(26,:,:) + d3DM(27,:,:) + d3DM(28,:,:))

!The third-derivative density is not printed out here, but it is a useful
!intermediate quantity.  We now compute the "Laplacian" piece of the stress force
!as
stressForce(1,:,:) = beta/(-2)*(d3rho(1,:,:)+d3rho(6,:,:)+d3rho(8,:,:))
stressForce(2,:,:) = beta/(-2)*(d3rho(2,:,:)+d3rho(4,:,:)+d3rho(9,:,:))
stressForce(3,:,:) = beta/(-2)*(d3rho(3,:,:)+d3rho(5,:,:)+d3rho(7,:,:))

!The "hard part" of the tension force is the stress tensor piece.  It is
stressForce(1,:,:) = stressForce(1,:,:)                                         &
                   +(1-alpha)/2*(d3DM(1,:,:)+d3DM(11,:,:)+d3DM(6,:,:)           &
                                +d3DM(8,:,:)+d3DM(22,:,:)+d3DM(24,:,:))         &
                   +alpha/2*(2*d3DM(11,:,:)+d3DM(16,:,:)+d3DM(18,:,:)           &
                            +d3DM(22,:,:)+d3DM(24,:,:))
stressForce(2,:,:) = stressForce(2,:,:)                                         &
                   +(1-alpha)/2*(d3DM(2,:,:)+d3DM(12,:,:)+d3DM(4,:,:)           &
                                +d3DM(9,:,:)+d3DM(20,:,:)+d3DM(25,:,:))         &
                   +alpha/2*(2*d3DM(12,:,:)+d3DM(14,:,:)+d3DM(19,:,:)           &
                            +d3DM(20,:,:)+d3DM(25,:,:))
stressForce(3,:,:) = stressForce(3,:,:)                                         &
                   +(1-alpha)/2*(d3DM(3,:,:)+d3DM(13,:,:)+d3DM(5,:,:)           &
                                +d3DM(7,:,:)+d3DM(21,:,:)+d3DM(23,:,:))         &
                   +alpha/2*(2*d3DM(13,:,:)+d3DM(15,:,:)+d3DM(17,:,:)           &
                            +d3DM(21,:,:)+d3DM(23,:,:))

end subroutine stress_wfn

!-------------------------------------------------------------------------------!
!                               density_wfn                                     !
!                                                                               !
! This subroutine computes the electron density at a set of points.             !
!      rho(r) = SUM(i=1:n_orbs) mo_occ(i) * PSI(r)**2                           !
!      PSI(r) = SUM(j=1:n_prim) MO_coeff(j)*primitive(j,r)                      !
!-------------------------------------------------------------------------------!
! mo_coeff -- molecular orbital coefficients                                    !
! mo_occ -- molecular orbital occupation numbers                                !
! n_orbs -- the number of molecular orbitals read in.                           !
! n_prim -- the number of Gaussian primitives.                                  !
! spinorbs -- logical variable that tells us whether we have spin orbitals or   !
!             only spatial orbitals.                                            !
!-------------------------------------------------------------------------------!
! prim_r(1:n_prim) -- the value of the primitives at r.                         !
! orb_contribution(j) -- the contribution of the jth MO to the density matrix.  !
!-------------------------------------------------------------------------------!

subroutine density_wfn(rho,XYZ,npts)

USE utilities_wfn

IMPLICIT NONE

INTEGER(istd) :: ipt,npts,j
REAL(dbl)     :: rho(1:npts,1:2),XYZ(1:3,1:npts)
REAL(dbl)     :: prim_r(1:n_prim)
REAL(dbl)     :: orb_contribution(1:n_orbs)

rho(:,:) = 0.0_dbl

print*, 'in density subroutine'

DO ipt=1,npts
   !Evaluate the values of the Gaussian primitives at these points
   call Gau_prim(prim_r,XYZ(1:3,ipt))
   
   !We will start by computing the alpha-spin density.
   orb_contribution(:) = 0.0_dbl
   DO j=1,n_orbs
      IF(MO_occ(j,1) > EPSILON(MO_occ(j,1)))THEN
         orb_contribution(j) = MO_occ(j,1)                                      &
                                  *DOT_PRODUCT(MO_coeff(:,j,1),prim_r(:))**2 
      ENDIF
   ENDDO
   rho(ipt,1) = SUM(orb_contribution)
   
   IF (.not. spinorbs) THEN
      !The second component of the density is the same as the first.
      rho(ipt,2) = rho(ipt,1)
   ELSE
      !Evaluate the beta-spin part of the density in the same way.
      orb_contribution(:) = 0.0_dbl
      DO j=1,n_orbs
         IF(MO_occ(j,2) > EPSILON(MO_occ(j,2)))THEN
            orb_contribution(j) = MO_occ(j,2)                                      &
                                     *DOT_PRODUCT(MO_coeff(:,j,2),prim_r(:))**2 
         ENDIF
      ENDDO
      rho(ipt,2) = SUM(orb_contribution)  
   ENDIF
ENDDO

print*, 'end density subroutine.'

end subroutine density_wfn
   
!-------------------------------------------------------------------------------!
!                               localIP_wfn                                     !
!                                                                               !
! This subroutine computes the local ionization potential of Politzer at a set  !
! of points.                                                                    !
!      localIP = SUM(i=1:n_orbs) mo_occ(i)*mo_energy(i) * PSI(r)*PSI(r)         !
!      PSI(r) = SUM(j=1:n_prim) MO_coeff(j)*primitive(j,r)                      !
!-------------------------------------------------------------------------------!
! mo_energy -- molecular orbital energy                                         !
! mo_coeff -- molecular orbital coefficients                                    !
! mo_occ -- molecular orbital occupation numbers                                !
! n_orbs -- the number of molecular orbitals read in.                           !
! n_prim -- the number of Gaussian primitives.                                  !
! spinorbs -- logical variable that tells us whether we have spin orbitals or   !
!             only spatial orbitals.                                            !
!-------------------------------------------------------------------------------!
! prim_r1(1:n_prim) -- the value of the primitives at r1.                       !
! prim_r2(1:n_prim) -- the value of the primitve basis functions at r2.         !
! orb_contribution(j) -- the contribution of the jth MO to the density matrix.  !
! localIP(1:npts) -- the local ionization potential.                            !
!-------------------------------------------------------------------------------!

subroutine localIP_wfn(localIP,XYZ,npts)

USE utilities_wfn

IMPLICIT NONE

INTEGER(istd) :: ipt,npts,j
REAL(dbl)     :: localIP(1:npts,1:2),XYZ(1:3,1:npts)
REAL(dbl)     :: prim_r(1:n_prim)
REAL(dbl)     :: orb_contribution(1:n_orbs)

DO ipt=1,npts
   !Evaluate the values of the Gaussian primitives at these points
   call Gau_prim(prim_r,XYZ(1:3,ipt))
   
   !We will start by computing the alpha-spin local IP.
   orb_contribution(:) = 0.0_dbl
   DO j=1,n_orbs
      IF(MO_occ(j,1) > EPSILON(MO_occ(j,1)))THEN
         orb_contribution(j) = MO_occ(j,1)*MO_energy(j,1)                       &
                                  *DOT_PRODUCT(MO_coeff(:,j,1),prim_r(:))**2 
      ENDIF
   ENDDO
   localIP(ipt,1) = SUM(orb_contribution)
   
   IF (.not. spinorbs) THEN
      !The second component of the localIP is the same as the first.
      localIP(ipt,2) = localIP(ipt,1)
   ELSE
      !Evaluate the beta-spin part of the localIP in the same way.
      orb_contribution(:) = 0.0_dbl
      DO j=1,n_orbs
         IF(MO_occ(j,2) > EPSILON(MO_occ(j,2)))THEN
            orb_contribution(j) = MO_occ(j,2)*MO_energy(j,2)                    &
                                     *DOT_PRODUCT(MO_coeff(:,j,2),prim_r(:))**2 
         ENDIF
      ENDDO
      localIP(ipt,2) = SUM(orb_contribution)  
   ENDIF
ENDDO

end subroutine localIP_wfn

!-------------------------------------------------------------------------------!
!                               DM1_wfn                                         !
!                                                                               !
! This subroutine computes the 1-electron reduced density matrix at a set of    !
! points.                                                                       !
!      DM(r1,r2) = SUM(i=1:n_orbs) mo_occ(i) * PSI(r1)*PSI(r2)                  !
!      PSI(r) = SUM(j=1:n_prim) MO_coeff(j)*primitive(j,r)                      !
!-------------------------------------------------------------------------------!
! mo_coeff -- molecular orbital coefficients                                    !
! mo_occ -- molecular orbital occupation numbers                                !
! n_orbs -- the number of molecular orbitals read in.                           !
! n_prim -- the number of Gaussian primitives.                                  !
! spinorbs -- logical variable that tells us whether we have spin orbitals or   !
!             only spatial orbitals.                                            !
!-------------------------------------------------------------------------------!
! prim_r1(1:n_prim) -- the value of the primitives at r1.                       !
! prim_r2(1:n_prim) -- the value of the primitives at r2.                       !
! orb_contribution(j) -- the contribution of the jth MO to the density matrix.  !
! DM(1:npts,1:2) -- the value of the density matrix at each point.              !
!-------------------------------------------------------------------------------!

subroutine DM1_wfn(DM,R6d,npts)

USE utilities_wfn

IMPLICIT NONE

INTEGER(istd) :: ipt,npts,j
REAL(dbl)     :: DM(1:npts,1:2),R6d(1:6,1:npts)
REAL(dbl)     :: prim_r1(1:n_prim),prim_r2(1:n_prim)
REAL(dbl)     :: orb_contribution(1:n_orbs)

DO ipt=1,npts
   !Evaluate the values of the Gaussian primitives at these points
   call Gau_prim(prim_r1,R6d(1:3,ipt))
   
   !Don't evaluate the primitives at r2 if r2 is very close to r1.
   IF (MAXVAL(ABS(R6d(1:3,ipt)-R6d(4:6,ipt))) < EPSILON(R6d(1,ipt))**(.75)) THEN
      prim_r2 = prim_r1
   ELSE
      call Gau_prim(prim_r2,R6d(4:6,ipt))
   ENDIF
      
   !We will start by computing the alpha-spin density matrix.
   orb_contribution(:) = 0.0_dbl
   DO j=1,n_orbs
      IF(MO_occ(j,1) > EPSILON(MO_occ(j,1)))THEN
         orb_contribution(j) = MO_occ(j,1)                                      &
                                  *DOT_PRODUCT(MO_coeff(:,j,1),prim_r1(:))      &
                                  *DOT_PRODUCT(MO_coeff(:,j,1),prim_r2(:)) 
      ENDIF
   ENDDO
   DM(ipt,1) = SUM(orb_contribution)
   
   IF (.not. spinorbs) THEN
      !The second component of the density matrix is the same as the first.
      DM(ipt,2) = DM(ipt,1)
   ELSE
      !Evaluate the beta-spin part of the density matrix in the same way.
      orb_contribution(:) = 0.0_dbl
      DO j=1,n_orbs
         IF(MO_occ(j,2) > EPSILON(MO_occ(j,2)))THEN
            orb_contribution(j) = MO_occ(j,2)                                   &
                                  *DOT_PRODUCT(MO_coeff(:,j,2),prim_r1(:))      &
                                  *DOT_PRODUCT(MO_coeff(:,j,2),prim_r2(:)) 
         ENDIF
      ENDDO
      DM(ipt,2) = SUM(orb_contribution)  
   ENDIF
ENDDO

end subroutine DM1_wfn
   
!-------------------------------------------------------------------------------!
!                               DM2_wfn                                         !
!                                                                               !
! Evaluates the two-electron reduced density matrix at a set of grid points.    !
! The formula that is used is only valid for Slater determinants.  It is        !
!   DM2(r1,r2,r1',r2',alpha,alpha) = 1/2*{DM1(r1,r1',alpha)*DM1(r2,r2',alpha)   !
!                                  - DM1(r2,r1',alpha)*DM1(r1,r2',alpha)}       !
!   DM2(r1,r2,r1',r2',alpha,beta) = 1/2*{DM1(r1,r1',alpha)*DM1(r2,r2',beta)     !
!-------------------------------------------------------------------------------!
! R12d -- the input points in the order (r1,r2,r1',r2')                         !
! R6a -- the first collection of points for evaluating the 1RDM.               !
! R6b -- the second collection of points for evaluating the 1RDM.              !
! DM1a -- the first density matrix.                                             !
! DM2a -- the second density matrix.                                            !
!-------------------------------------------------------------------------------!

subroutine DM2_wfn(DM2,R12d,npts)

IMPLICIT NONE

INTEGER(istd) :: npts
REAL(dbl)     :: DM2(1:npts,1:2,1:2),DM1a(1:npts,1:2),DM1b(1:npts,1:2)
REAL(dbl)     :: R6a(1:6,1:npts),R6b(1:6,1:npts),R12d(1:12,1:npts)

!Start by computing the positive componetnts.

R6a(1:3,:) = R12d(1:3,:); R6a(4:6,:) = R12d(7:9,:)
R6b(1:3,:) = R12d(4:6,:); R6b(4:6,:) = R12d(10:12,:)

call DM1_wfn(DM1a,R6a,npts)
call DM1_wfn(DM1b,R6b,npts)

!Opposite spin components
DM2(1:npts,1,2) = DM1a(:,1)*DM1b(:,2)
DM2(1:npts,2,1) = DM1a(:,2)*DM1b(:,1)

!The independent-electron piece of the same-spin components is
DM2(1:npts,1,1) = DM1a(:,1)*DM1b(:,1)
DM2(1:npts,2,2) = DM1a(:,2)*DM1b(:,2)

!The exchange correction is now computed.
R6a(1:3,:) = R12d(4:6,:)
R6b(1:3,:) = R12d(1:3,:)

call DM1_wfn(DM1a,R6a,npts)
call DM1_wfn(DM1b,R6b,npts)
DM2(1:npts,1,1) = DM2(1:npts,1,1) - DM1a(:,1)*DM1b(:,1)
DM2(1:npts,2,2) = DM2(1:npts,2,2) - DM1a(:,2)*DM1b(:,2)
 
end subroutine DM2_wfn

!-------------------------------------------------------------------------------!
!                               hx_wfn                                          !
!                                                                               !
! This subroutine computes the spin-resolved exchange hole.                     !
!   hx(r,r',spin) = -1*DM1(r,r',spin)**2/(rho(r,spin)*rho(r',spin))             !
!-------------------------------------------------------------------------------!
! mo_coeff -- molecular orbital coefficients                                    !
! mo_occ -- molecular orbital occupation numbers                                !
! n_orbs -- the number of molecular orbitals read in.                           !
! n_prim -- the number of Gaussian primitives.                                  !
! spinorbs -- logical variable that tells us whether we have spin orbitals or   !
!             only spatial orbitals.                                            !
!-------------------------------------------------------------------------------!
! prim_r1(1:n_prim) -- the value of the primitives at r1.                       !
! prim_r2(1:n_prim) -- the value of the primitives at r2.                       !
! orb_contDM(j) -- the contribution of the jth MO to the density matrix.        !
! orb_contdensity1(j,1:2) -- the contributions of the jth MO to the electron    !
!                           density at r1                                       !
! orb_contdensity2(j,1:2) -- the contributions of the jth MO to the electron    !
!                           density at r2                                       !
! hx(1:npts,1:2) -- the spin components of the exchange hole.                   !
! R6d -- the a point (r,r') in 6 dimensions.                                    !
! MO_r1 -- the value of the MO at r1.                                           !
! MO_r2 -- the value of the MO at r2.                                           !
!-------------------------------------------------------------------------------!

subroutine hx_wfn(hx,R6d,npts)

use utilities_wfn

IMPLICIT NONE
INTEGER(istd) :: ipt,npts,j
REAL(dbl)     :: DM(1:npts,1:2),R6d(1:6,1:npts)
REAL(dbl)     :: prim_r1(1:n_prim),prim_r2(1:n_prim)
REAL(dbl)     :: orb_contDM(1:n_orbs)
REAL(dbl)     :: orb_contdensity1(1:n_orbs),orb_contdensity2(1:n_orbs)
REAL(dbl)     :: MO_r1,MO_r2,hx(1:npts,1:2)

DO ipt=1,npts
   !Evaluate the values of the Gaussian primitives at these points
   call Gau_prim(prim_r1,R6d(1:3,ipt))
   call Gau_prim(prim_r2,R6d(4:6,ipt))
      
   !We will start by computing the alpha-spin density matrix and the alpha-spin
   !densities
   orb_contDM(:) = 0.0_dbl
   orb_contdensity1(:) = 0.0_dbl
   orb_contdensity2(:) = 0.0_dbl
   DO j=1,n_orbs
      IF(MO_occ(j,1) > EPSILON(MO_occ(j,1)))THEN
         MO_r1 = DOT_PRODUCT(MO_coeff(:,j,1),prim_r1(:))
         MO_r2 = DOT_PRODUCT(MO_coeff(:,j,1),prim_r2(:))
         orb_contDM(j) = MO_occ(j,1)*MO_r1*MO_r2
         orb_contdensity1(j) = MO_occ(j,1)*MO_r1**2
         orb_contdensity2(j) = MO_occ(j,1)*MO_r2**2       
      ENDIF
   ENDDO
   
   !construct the alpha-spin exchange hole
   hx(ipt,1) = -1*SUM(orb_contDM)**2                                 &
                     /(SUM(orb_contdensity1)*SUM(orb_contdensity2))
   
   IF (.not. spinorbs) THEN
      !The second component of the exchange hole is the same as the first.
      hx(ipt,2) = hx(ipt,1)
   ELSE
      !Evaluate the beta-spin part of the exchange hold in the same way.
      orb_contDM(:) = 0.0_dbl
      orb_contdensity1(:) = 0.0_dbl
      orb_contdensity2(:) = 0.0_dbl
      DO j=1,n_orbs
         IF(MO_occ(j,2) > EPSILON(MO_occ(j,2)))THEN
            MO_r1 = DOT_PRODUCT(MO_coeff(:,j,2),prim_r1(:))
            MO_r2 = DOT_PRODUCT(MO_coeff(:,j,2),prim_r2(:))
            orb_contDM(j) = MO_occ(j,2)*MO_r1*MO_r2
            orb_contdensity1(j) = MO_occ(j,2)*MO_r1**2
            orb_contdensity2(j) = MO_occ(j,2)*MO_r2**2       
         ENDIF
      ENDDO
      hx(ipt,2) = -1*SUM(orb_contDM)**2                                 &
                        /(SUM(orb_contdensity1)*SUM(orb_contdensity2))
   ENDIF
ENDDO

end subroutine hx_wfn

!-------------------------------------------------------------------------------!
!                              KSresponse_wfn                                   !
!                                                                               !
! This subroutine computes the Kohn-Sham response.                              !
! SUM(i) SUM(j/=i) (n(j)-n(i))/(e(j)-e(i)) PSI(i,r)PSI(j,r)PSI(i,r')PSI(j,r')   !
!-------------------------------------------------------------------------------!
! mo_energy -- molecular orbital energy.                                        !
! mo_coeff -- molecular orbital coefficients                                    !
! mo_occ -- molecular orbital occupation numbers                                !
! n_orbs -- the number of molecular orbitals read in.                           !
! n_prim -- the number of Gaussian primitives.                                  !
! spinorbs -- logical variable that tells us whether we have spin orbitals or   !
!             only spatial orbitals.                                            !
!-------------------------------------------------------------------------------!
! prim_r1(1:n_prim) -- the value of the primitives at r1.                       !
! prim_r2(1:n_prim) -- the value of the primitives at r2.                       !
! orb_contribution(1:n_orbs) -- the contribution of the jth MO to the KSresponse!
! KSresponse(1:npts,1:2) -- the spin components of the KS response function.    !
! R6d -- a point (r,r') in 6 dimensions.                                        !
! MO_r1 -- the value of the MOs at r1.                                          !
! MO_r2 -- the value of the MOs  at r2.                                         !
!-------------------------------------------------------------------------------!

subroutine KSresponse_wfn(KSresponse,R6d,npts)

use utilities_wfn

IMPLICIT NONE
INTEGER(istd) :: ipt,npts,j,k
REAL(dbl)     :: DM(1:npts,1:2),R6d(1:6,1:npts)
REAL(dbl)     :: prim_r1(1:n_prim),prim_r2(1:n_prim)
REAL(dbl)     :: orb_contribution(1:n_orbs)
REAL(dbl)     :: MO_r1(1:n_orbs),MO_r2(1:n_orbs),KSresponse(1:npts,1:2)


KSresponse(:,:) = 0.0_dbl

DO ipt=1,npts
   !Evaluate the values of the Gaussian primitives at these points
   call Gau_prim(prim_r1,R6d(1:3,ipt))
   call Gau_prim(prim_r2,R6d(4:6,ipt))
      
   !We will start by computing the alpha-spin density matrix and the alpha-spin
   !densities
   FORALL(j=1:n_orbs)
         MO_r1(j) = DOT_PRODUCT(MO_coeff(:,j,1),prim_r1(:))
         MO_r2(j) = DOT_PRODUCT(MO_coeff(:,j,1),prim_r2(:)) 
   ENDFORALL
   
   DO j=1,n_orbs
      DO k=1,n_orbs
         IF (ABS(MO_occ(j,1)-MO_occ(k,1)) > EPSILON(MO_occ(j,1))) THEN
            !There is a nonzero contribution from this orbital.
            KSresponse(ipt,1) = KSresponse(ipt,1)                               &
                                  +(MO_occ(k,1)-MO_occ(j,1))                    &
                                    /(MO_energy(k,1)-MO_energy(j,1))            &
                                       *MO_r1(j)*MO_r1(k)*MO_r2(j)*MO_r2(k)
         ENDIF
      ENDDO
   ENDDO
   
   
   IF (.not. spinorbs) THEN
      !The second component of the response function is the same as the first.
      KSresponse(ipt,2) = KSresponse(ipt,1)
   ELSE
      !Evaluate the beta-spin part of the KS response in the same way.
      FORALL(j=1:n_orbs)
            MO_r1(j) = DOT_PRODUCT(MO_coeff(:,j,2),prim_r1(:))
            MO_r2(j) = DOT_PRODUCT(MO_coeff(:,j,2),prim_r2(:)) 
      ENDFORALL
   
      DO j=1,n_orbs
         DO k=1,n_orbs
            IF (ABS(MO_occ(j,2)-MO_occ(k,2)) > EPSILON(MO_occ(j,2))) THEN
               !There is a nonzero contribution from this orbital.
               KSresponse(ipt,2) = KSresponse(ipt,2)                            &
                                  +(MO_occ(k,2)-MO_occ(j,2))                    &
                                    /(MO_energy(k,2)-MO_energy(j,2))            &
                                       *MO_r1(j)*MO_r1(k)*MO_r2(j)*MO_r2(k)
            ENDIF
         ENDDO
      ENDDO
   ENDIF
ENDDO

end subroutine KSresponse_wfn

!-------------------------------------------------------------------------------!
!                               orbs_wfn                                        !
!                                                                               !
! This subroutine gives the value of a molecular orbital at a set of grid points!
!                                                                               !
!      PSI(r) = SUM(j=1:n_prim) MO_coeff(j)*primitive(j,r)                      !
!-------------------------------------------------------------------------------!
! mo_coeff -- molecular orbital coefficients                                    !
! mo_occ -- molecular orbital occupation numbers                                !
! n_orbs -- the number of molecular orbitals read in.                           !
! n_prim -- the number of Gaussian primitives.                                  !
! spinorbs -- logical variable that tells us whether we have spin orbitals or   !
!             only spatial orbitals.                                            !
!-------------------------------------------------------------------------------!
! prim_r(1:n_prim) -- the value of the primitives at r.                         !
! orb_contribution(j) -- the contribution of the jth MO to the density matrix.  !
!-------------------------------------------------------------------------------!
! n_pts -- the number of points at which the density should be evaluated.       !
! XYZ(1:3,1:npts) -- the location of the point in Cartesian coordinates.        !
! PHIa(1:nlist,1:npts) -- alpha-spin molecular orbitals at the grid points.     !
! PHIb(1:nlist,1:npts) -- beta-spin MOs at the grid points.                     !
! nlist -- the number of orbitals in the list that one wishes to evaluate.      !
! iorbs(1:nlist) -- the list of MOs that one wish to evaluate at the points.    !
!-------------------------------------------------------------------------------!

subroutine orbs_wfn(PHIa,PHIb,nlist,iorbs,XYZ,npts)

USE utilities_wfn
USE grid_vars

IMPLICIT NONE

INTEGER(istd) :: ipt,npts,j,nlist,iorbs(1:nlist)
REAL(dbl)     :: PHIa(1:nlist,1:npts),PHIb(1:nlist,1:npts),XYZ(1:3,1:npts)
REAL(dbl)     :: prim_r(1:n_prim)

PHIa = 0.0_dbl; PHIb = 0.0_dbl

DO ipt=1,npts
   !Evaluate the values of the Gaussian primitives at these points
   call Gau_prim(prim_r,XYZ(1:3,ipt))
   
   !We will start by computing the alpha-spin orbitals
   FORALL(j=1:nlist)
      PHIa(j,ipt) = DOT_PRODUCT(MO_coeff(:,iorbs(j),1),prim_r(:)) 
   ENDFORALL

   IF (.not. spinorbs) THEN
      !The beta-spin orbitals are the same
      PHIb(:,ipt) = PHIa(:,ipt)
   ELSE
      !Evaluate the beta-spin orbitals
      FORALL(j=1:nlist)
         PHIb(j,ipt) = DOT_PRODUCT(MO_coeff(:,iorbs(j),2),prim_r(:)) 
      ENDFORALL
   ENDIF
ENDDO

end subroutine orbs_wfn

!-------------------------------------------------------------------------------!
!                           wavefunction_wfn                                    !
!                                                                               !
! This gives the value of a Slater determinant at a set of grid points. Only for!
! a Slater determinant.                                                         !
!                                                                               !
! A list of nlist_a alpha-spin orbitals and nlist_b beta-spin orbitals is input !
! A Slater determinant is then constructed, and the wavefunction is printed out !
! It is assumed that 3*(nlist_a + nlist_b) = dim, where dim is the dimension    !
! of the space, which is three times n_el_det, the number of electrons in the   !
! determinant.                                                                  !
!-------------------------------------------------------------------------------!
! dim_pt -- The dimension of the points.  It will normally be three times the   !
!           number of electrons but it is always 3 * (nlist_a + nlist_b)        !
!           This allows one to choose different N for testing purposes.         !
! nlist_a -- the number of alpha-spin orbitals in the determinant.              !
! nlist_b -- the number of beta-spin orbitals in the determinant.               !
! iorbs_a(1:nlist_a) -- the list of alpha-spin orbitals.                        !
! iorbs_b(1:nlist_b) -- the list of beta-spin orbitals.                         !
! XYZ_ND(1:3*nel_det,1:npts) -- the list of points to be evaluated.             !
! npts -- the number of points to be evaluated.                                 !
! Slaterdet(1:npts) -- the Slater determinant evaluated at the input points.    !
! k_elpos -- loops through the electronic coordinates.                          !
! elpos(1:3) -- an electronic position in the Slater determinant.               !
! Sl_matrix -- the matrix used to construct the Slater determinant.             !
! el_count -- a counter for the number of electrons.                            !
!-------------------------------------------------------------------------------!

subroutine wavefunction_wfn(SlaterDet,nlist_a,nlist_b,iorbs_a,iorbs_b,dim_pt, &
                            XYZ_Nd,npts)

USE utilities_wfn

IMPLICIT NONE

INTEGER(istd) :: ipt,npts,j,nlist_a,nlist_b,dim_pt,k_elpos,el_count
INTEGER(istd) :: iorbs_a(1:nlist_a),iorbs_b(1:nlist_b)
REAL(dbl)     :: SlaterDet(1:npts),XYZ_Nd(1:dim_pt,1:npts),elpos(1:3)
REAL(dbl)     :: prim_r(1:n_prim)
REAL(dbl), ALLOCATABLE :: Sl_matrix(:,:)

IF (dim_pt /= 3*(nlist_a+nlist_b)) THEN
   WRITE(*,*) 'Dimension mismatch in subroutine wavefunction_wfn for evaluating'
   WRITE(*,*) 'a Slater determinant.'
   STOP
ENDIF

ALLOCATE(Sl_matrix(1:dim_pt/3,1:dim_pt/3))

DO ipt=1,npts   
   !Loop through the different 3N-dimensional points.
   DO k_elpos=1,dim_pt/3 
      !Loop through the different electron positions.
      elpos(1:3) = XYZ_Nd(3*(k_elpos-1)+1,3*k_elpos)
      !Evaluate the values of the Gaussian primitives at these points
      call Gau_prim(prim_r,elpos(1:3))
      el_count = 0
      !We will start by computing the alpha-spin orbitals
      DO j=1,nlist_a
         el_count = el_count + 1
         Sl_matrix(k_elpos,el_count) = DOT_PRODUCT(MO_coeff(:,iorbs_a(j),1),prim_r(:)) 
      ENDDO
      DO j=1,nlist_b
         el_count = el_count + 1
         Sl_matrix(k_elpos,el_count) = DOT_PRODUCT(MO_coeff(:,iorbs_b(j),2),prim_r(:)) 
      ENDDO
   ENDDO
   !Now that the Slater matrix is constructed, we take its determinant
   call FindDet(SlaterDet(ipt),Sl_matrix,el_count)  
ENDDO

DEALLOCATE(Sl_matrix)

end subroutine wavefunction_wfn

!-------------------------------------------------------------------------------!
!                               dDM1_wfn                                        !
!                                                                               !
! Gives the gradient of the density matrix at a set of grid points.             !
!-------------------------------------------------------------------------------!
!      DM(r,r') = SUM(i=1:n_orbs) mo_occ(i) * PSI(r1)*PSI(r2)                   !
!      PSI(r) = SUM(j=1:n_prim) MO_coeff(j)*primitive(j,r)                      !
!-------------------------------------------------------------------------------!
! mo_coeff -- molecular orbital coefficients                                    !
! mo_occ -- molecular orbital occupation numbers                                !
! n_orbs -- the number of molecular orbitals read in.                           !
! n_prim -- the number of Gaussian primitives.                                  !
! spinorbs -- logical variable that tells us whether we have spin orbitals or   !
!             only spatial orbitals.                                            !
!-------------------------------------------------------------------------------!
! prim_r1(1:n_prim) -- the value of the primitives at r1.                       !
! prim_r2(1:n_prim) -- the value of the primitives at r2.                       !
! dORB_contribution(1:3,j) -- the contribution of the jth MO to the gradient of !
!                            the density matrix.                                !
! dDM(1:3,1:npts,1:2) -- the spin (1:2) components of the gradient of the first-!
!                        order density matrix at the given grid points (1:npts).!
!                        The components are:                                    !
!                         1 -- d/dx D(r',r)                                     !
!                         2 -- d/dy D(r',r)                                     !
!                         3 -- d/dz D(r',r)                                     !
!-------------------------------------------------------------------------------!

subroutine dDM1_wfn(dDM,R6d,npts)

USE utilities_wfn

IMPLICIT NONE

INTEGER(istd) :: ipt,npts,j,k_prim,l
REAL(dbl)     :: dDM(1:3,1:npts,1:2),R6d(1:6,1:npts)
REAL(dbl)     :: prim_r1(1:n_prim),prim_r2(1:n_prim),dprim_r2(1:3,1:n_prim)
REAL(dbl)     :: dORB_contribution(1:3,1:n_orbs)

DO ipt=1,npts
   !Evaluate the values of the Gaussian primitives at these points
   call Gau_prim(prim_r1,R6d(1:3,ipt))
   DO k_prim=1,n_prim
      call dGau_prim(prim_r2(k_prim),dprim_r2(1:3,k_prim),k_prim,R6d(4:6,ipt))
   ENDDO
         
   !We will start by computing the alpha-spin density matrix.
   dorb_contribution(:,:) = 0.0_dbl
   DO j=1,n_orbs
      IF (MO_occ(j,1) > EPSILON(MO_occ(j,1))) THEN
         FORALL(L=1:3)
               !For all components of the gradient:
               dorb_contribution(L,j) = DOT_PRODUCT(MO_coeff(:,j,1),dprim_r2(L,:)) 
         ENDFORALL
         dorb_contribution(:,j) = dorb_contribution(:,j)                               &
                                  *MO_occ(j,1)*DOT_PRODUCT(MO_coeff(:,j,1),prim_r1(:))   
         !This two-step method for computing the orbital contributions has the advantage
         !of avoiding repeated evaluation of the "second term", which does not have
         !any dependence on the component of the gradient being evaluated. In this
         !sense it is better than the more straightforward alternative:
         !         dorb_contribution(L,j) = MO_occ(j,1)                               &
         !                                 *DOT_PRODUCT(MO_coeff(:,j,2),prim_r1(:))   &
         !                                 *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r2(L,:))          
      ENDIF
   ENDDO
   FORALL(L=1:3)
         dDM(L,ipt,1) = SUM(dorb_contribution(L,:))
   ENDFORALL
   
   IF (.not. spinorbs) THEN
      !The second component of the density matrix is the same as the first.
      dDM(:,ipt,2) = dDM(:,ipt,1)
   ELSE
      !Evaluate the beta-spin part of the density matrix in the same way.
      DO j=1,n_orbs
         IF (MO_occ(j,2) > EPSILON(MO_occ(j,2))) THEN
            FORALL(L=1:3)
                  !For all components of the gradient:
                  dorb_contribution(L,j) = DOT_PRODUCT(MO_coeff(:,j,2),dprim_r2(L,:)) 
            ENDFORALL
            dorb_contribution(:,j) = dorb_contribution(:,j)                               &
                                     *MO_occ(j,2)*DOT_PRODUCT(MO_coeff(:,j,2),prim_r1(:))  
         ENDIF
      ENDDO
      FORALL(L=1:3)
            dDM(L,ipt,2) = SUM(dorb_contribution(L,:))
      ENDFORALL
   ENDIF
ENDDO

end subroutine dDM1_wfn


!-------------------------------------------------------------------------------!
!                               d2DM1_wfn                                       !
!                                                                               !
! Gives the second derivative components of the density matrix at a set of      !
! grid points.                                                                  !
!-------------------------------------------------------------------------------!
!      DM(r1,r2) = SUM(i=1:n_orbs) mo_occ(i) * PSI(r1)*PSI(r2)                  !
!      PSI(r) = SUM(j=1:n_prim) MO_coeff(j)*primitive(j,r)                      !
!-------------------------------------------------------------------------------!
! mo_coeff -- molecular orbital coefficients                                    !
! mo_occ -- molecular orbital occupation numbers                                !
! n_orbs -- the number of molecular orbitals read in.                           !
! n_prim -- the number of Gaussian primitives.                                  !
! spinorbs -- logical variable that tells us whether we have spin orbitals or   !
!             only spatial orbitals.                                            !
!-------------------------------------------------------------------------------!
! prim_r1(1:n_prim) -- the value of the primitives at r1.                       !
! prim_r2(1:n_prim) -- the value of the primitives at r2.                       !
! dprim_r1(1:n_prim) -- the gradient of the primitives at r1.                   !
! dprim_r2(1:n_prim) -- the gradient of the primitives at r2.                   !
! d2prim_r2(1:n_prim) -- the second derivative tensor of the primitives at r2.  !
! d2ORB_contribution(1:12,j) -- the contribution of the jth MO to the 2nd       !
!                              derivative tensor of the the density matrix.     !
! d2DM(1:12,1:npts,1:2) -- the spin (1:2) components of the second derivatives  !
!                         of the first-order density matrix at the given        !
!                         grid points (1:npts).  The components are:            !
!                         1 -- d^2/dx^2 D(r',r)                                 !
!                         2 -- d^2/dy^2 D(r',r)                                 !
!                         3 -- d^2/dz^2 D(r',r)                                 !
!                         4 -- d/dx d/dy D(r',r)                                !
!                         5 -- d/dx d/dz D(r',r)                                !
!                         6 -- d/dy d/dz D(r',r)                                !
!                         7 -- d/dx d/dx' D(r',r)                               !
!                         8 -- d/dy d/dy' D(r',r)                               !
!                         9 -- d/dz d/dz' D(r',r)                               !
!                        10 -- d/dx d/dy' D(r',r)                               !
!                        11 -- d/dx d/dz' D(r',r)                               !
!                        12 -- d/dy d/dz' D(r',r)                               !
!-------------------------------------------------------------------------------!
!       prim = G(x,y,z) = x^i * y^j * z^k * EXP(alpha(i)*(R-R0){dot}(R-R0))     !
!   dprim(1) = d/dx  G(x,y,z)                                                   !
!   dprim(2) = d/dy  G(x,y,z)                                                   !
!   dprim(3) = d/dz  G(x,y,z)                                                   !
!  d2prim(1) = d^2/dx^2  G(x,y,z)                                               !
!  d2prim(2) = d^2/dy^2  G(x,y,z)                                               !
!  d2prim(3) = d^2/dz^2  G(x,y,z)                                               !
!  d2prim(4) = d/dx d/dy G(x,y,z)                                               !
!  d2prim(5) = d/dx d/dz G(x,y,z)                                               !
!  d2prim(6) = d/dy d/dz G(x,y,z)                                               !
!-------------------------------------------------------------------------------!

subroutine d2DM1_wfn(d2DM,R6d,npts)

USE utilities_wfn

IMPLICIT NONE

INTEGER(istd) :: ipt,npts,j,k_prim,l
REAL(dbl)     :: d2DM(1:12,1:npts,1:2),R6d(1:6,1:npts)
REAL(dbl)     :: prim_r1(1:n_prim),prim_r2(1:n_prim)
REAL(dbl)     :: dprim_r1(1:3,1:n_prim),dprim_r2(1:3,1:n_prim)
REAL(dbl)     :: d2prim_r2(1:6,1:n_prim)
REAL(dbl)     :: d2ORB_contribution(1:12,1:n_orbs)

DO ipt=1,npts
   !Evaluate the values of the Gaussian primitives and gradients:
   DO k_prim=1,n_prim
      call dGau_prim(prim_r1(k_prim),dprim_r1(1:3,k_prim),k_prim,R6d(1:3,ipt))
      call d2Gau_prim(prim_r2(k_prim),dprim_r2(1:3,k_prim),d2prim_r2(1:6,k_prim), &
                      k_prim,R6d(4:6,ipt))
   ENDDO

   !We will start by computing the alpha-spin density matrix.
   d2orb_contribution(:,:) = 0.0_dbl
   DO j=1,n_orbs
      IF (MO_occ(j,1) > EPSILON(MO_occ(j,1))) THEN
         !We could do a two-step (well, now multi-step) method similar to what we
         !used in dDM1_wfn to evaluate the orbital contributions.  However, the 
         !complexity becomes daunting now, and it seems better to favor
         !readability over speed.  HOWEVER, THIS IS A PLACE THE PROGRAM COULD BE
         !SPED UP A LOT!!!
         ! 1 -- d^2/dx^2 D(r',r)   
         ! 2 -- d^2/dy^2 D(r',r)                                 !
         ! 3 -- d^2/dz^2 D(r',r)                                 !
         ! 4 -- d/dx d/dy D(r',r)                                !
         ! 5 -- d/dx d/dz D(r',r)                                !
         ! 6 -- d/dy d/dz D(r',r)                                !
         FORALL(L=1:6)
               !For all components of the gradient:
               d2orb_contribution(L,j) = DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(L,:)) 
         ENDFORALL
         d2orb_contribution(1:6,j) = d2orb_contribution(1:6,j)                            &
                                     *MO_occ(j,1)*DOT_PRODUCT(MO_coeff(:,j,1),prim_r1(:))   
                                     
         !The next section of code could be sped up by a factor of two but it 
         !would have a significant "readability price."  So I (PWA) didn't bother.
         !   7 -- d/dx d/dx' D(r',r)                               !
         !   8 -- d/dy d/dy' D(r',r)                               !
         !   9 -- d/dz d/dz' D(r',r)                               !
         !  10 -- d/dx d/dy' D(r',r)                               !
         !  11 -- d/dx d/dz' D(r',r)                               !
         !  12 -- d/dy d/dz' D(r',r)                               !
         d2orb_contribution(7,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(1,:))   &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r2(1,:)) 
         d2orb_contribution(8,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(2,:))   &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r2(2,:)) 
         d2orb_contribution(9,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(3,:))   &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r2(3,:))
         d2orb_contribution(10,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(2,:))    &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r2(1,:)) 
         d2orb_contribution(11,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(3,:))    &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r2(1,:)) 
         d2orb_contribution(12,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(3,:))    &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r2(2,:))
      ENDIF
   ENDDO
   FORALL(L=1:12)
         d2DM(L,ipt,1) = SUM(d2orb_contribution(L,:))
   ENDFORALL
   
   IF (.not. spinorbs) THEN
      !The second component of the density matrix is the same as the first.
      d2DM(:,ipt,2) = d2DM(:,ipt,1)
   ELSE
      !Evaluate the beta-spin part of the density matrix in the same way.
      DO j=1,n_orbs
         IF (MO_occ(j,2) > EPSILON(MO_occ(j,2))) THEN
            ! 1 -- d^2/dx^2 D(r',r)   
            ! 2 -- d^2/dy^2 D(r',r)                                 !
            ! 3 -- d^2/dz^2 D(r',r)                                 !
            ! 4 -- d/dx d/dy D(r',r)                                !
            ! 5 -- d/dx d/dz D(r',r)                                !
            ! 6 -- d/dy d/dz D(r',r)                                !
            FORALL(L=1:6)
                  !For all components of the gradient:
                  d2orb_contribution(L,j) = DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(L,:)) 
            ENDFORALL
            d2orb_contribution(1:6,j) = d2orb_contribution(1:6,j)                            &
                                        *MO_occ(j,2)*DOT_PRODUCT(MO_coeff(:,j,2),prim_r1(:))   
                                     
            !The next section of code could be sped up by a factor of two but it 
            !would have a significant "readability price."  So I (PWA) didn't bother.
            !   7 -- d/dx d/dx' D(r',r)                               !
            !   8 -- d/dy d/dy' D(r',r)                               !
            !   9 -- d/dz d/dz' D(r',r)                               !
            !  10 -- d/dx d/dy' D(r',r)                               !
            !  11 -- d/dx d/dz' D(r',r)                               !
            !  12 -- d/dy d/dz' D(r',r)                               !
            d2orb_contribution(7,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(1,:))   &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r2(1,:)) 
            d2orb_contribution(8,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(2,:))   &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r2(2,:)) 
            d2orb_contribution(9,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(3,:))   &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r2(3,:))
            d2orb_contribution(10,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(2,:))    &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r2(1,:)) 
            d2orb_contribution(11,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(3,:))    &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r2(1,:)) 
            d2orb_contribution(12,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(3,:))    &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r2(2,:))
         ENDIF
      ENDDO
      FORALL(L=1:12)
            d2DM(L,ipt,2) = SUM(d2orb_contribution(L,:))
      ENDFORALL
   ENDIF
ENDDO

end subroutine d2DM1_wfn

!-------------------------------------------------------------------------------!
!                               d3DM1_wfn                                       !
!                                                                               !
! This subroutine computes the spin-components of the third derivative tensor   !
! of the one-electron reduced density matrix.                                   !
!-------------------------------------------------------------------------------!
! mo_coeff -- molecular orbital coefficients                                    !
! mo_occ -- molecular orbital occupation numbers                                !
! n_orbs -- the number of molecular orbitals read in.                           !
! n_prim -- the number of Gaussian primitives.                                  !
! spinorbs -- logical variable that tells us whether we have spin orbitals or   !
!             only spatial orbitals.                                            !
!-------------------------------------------------------------------------------!
! prim_r1(1:n_prim) -- the value of the primitives at r1.                       !
! prim_r2(1:n_prim) -- the value of the primitives at r2.                       !
! dprim_r1(1:n_prim) -- the gradient of the primitives at r1.                   !
! dprim_r2(1:n_prim) -- the gradient of the primitives at r2.                   !
! d2prim_r1(1:n_prim) -- the second derivative tensor of the primitives at r1.  !
! d2prim_r2(1:n_prim) -- the second derivative tensor of the primitives at r2.  !
! d3prim_r2(1:n_prim) -- the third derivative tensor of the primitives at r3.   !
! d3ORB_contribution(1:28,j) -- the contribution of the jth MO to the 3rd       !
!                              derivative tensor of the the density matrix.     !
! d3DM(1:28,1:npts,1:2) -- the spin (1:2) components of the third derivatives   !
!                         of the first-order density matrix at the given        !
!                         grid points (1:npts).  The components are:            !
!                         1 -- d^3/dx^3 D(r',r)                                 !
!                         2 -- d^3/dy^3 D(r',r)                                 !
!                         3 -- d^3/dz^3 D(r',r)                                 !
!                         4 -- d^2/dx^2 d/dy D(r',r)                            !
!                         5 -- d^2/dx^2 d/dz D(r',r)                            !
!                         6 -- d^2/dy^2 d/dx D(r',r)                            !
!                         7 -- d^2/dy^2 d/dz D(r',r)                            !
!                         8 -- d^2/dz^2 d/dx D(r',r)                            !
!                         9 -- d^2/dz^2 d/dy D(r',r)                            !
!                        10 -- d/dx d/dy d/dz D(r',r)                           !
!                        11 -- d^2/dx^2 d/dx' D(r',r)                           !
!                        12 -- d^2/dy^2 d/dy' D(r',r)                           !
!                        13 -- d^2/dz^2 d/dz' D(r',r)                           !
!                        14 -- d^2/dx^2 d/dy' D(r',r)                           !
!                        15 -- d^2/dx^2 d/dz' D(r',r)                           !
!                        16 -- d^2/dy^2 d/dx' D(r',r)                           !
!                        17 -- d^2/dy^2 d/dz' D(r',r)                           !
!                        18 -- d^2/dz^2 d/dx' D(r',r)                           !
!                        19 -- d^2/dz^2 d/dy' D(r',r)                           !
!                        20 -- d/dx d/dy d/dx' D(r',r)                          !
!                        21 -- d/dx d/dz d/dx' D(r',r)                          !
!                        22 -- d/dy d/dx d/dy' D(r',r)                          !
!                        23 -- d/dy d/dx d/dy' D(r',r)                          !
!                        24 -- d/dz d/dx d/dz' D(r',r)                          !
!                        25 -- d/dz d/dy d/dz' D(r',r)                          !
!                        26 -- d/dx' d/dy d/dz D(r',r)                          !
!                        27 -- d/dx d/dy' d/dz D(r',r)                          !
!                        28 -- d/dx d/dy d/dz' D(r',r)                          !
!-------------------------------------------------------------------------------!
!       prim = G(x,y,z) = x^i * y^j * z^k * EXP(alpha(i)*(R-R0){dot}(R-R0))     !
!   dprim(1) = d/dx  G(x,y,z)                                                   !
!   dprim(2) = d/dy  G(x,y,z)                                                   !
!   dprim(3) = d/dz  G(x,y,z)                                                   !
!  d2prim(1) = d^2/dx^2  G(x,y,z)                                               !
!  d2prim(2) = d^2/dy^2  G(x,y,z)                                               !
!  d2prim(3) = d^2/dz^2  G(x,y,z)                                               !
!  d2prim(4) = d/dx d/dy G(x,y,z)                                               !
!  d2prim(5) = d/dx d/dz G(x,y,z)                                               !
!  d2prim(6) = d/dy d/dz G(x,y,z)                                               !
!  d3prim(1) = d^3/dx^3 G(x,y,z)                                                !
!  d3prim(2) = d^3/dy^3 G(x,y,z)                                                !
!  d3prim(3) = d^3/dz^3 G(x,y,z)                                                !
!  d3prim(4) = d^2/dx^2 d/dy G(x,y,z)                                           !
!  d3prim(5) = d^2/dx^2 d/dz G(x,y,z)                                           !
!  d3prim(6) = d^2/dy^2 d/dx G(x,y,z)                                           !
!  d3prim(7) = d^2/dy^2 d/dz G(x,y,z)                                           !
!  d3prim(8) = d^2/dz^2 d/dx G(x,y,z)                                           !
!  d3prim(9) = d^2/dz^2 d/dy G(x,y,z)                                           !
!  d3prim(10)= d/dx d/dy d/dz G(x,y,z)                                          !
!-------------------------------------------------------------------------------!

subroutine d3DM1_wfn(d3DM,R6d,npts)

USE utilities_wfn

IMPLICIT NONE

INTEGER(istd) :: ipt,npts,j,k_prim,l
REAL(dbl)     :: d3DM(1:28,1:npts,1:2),R6d(1:6,1:npts)
REAL(dbl)     :: prim_r1(1:n_prim),prim_r2(1:n_prim)
REAL(dbl)     :: dprim_r1(1:3,1:n_prim),dprim_r2(1:3,1:n_prim)
REAL(dbl)     :: d2prim_r2(1:6,1:n_prim)
REAL(dbl)     :: d3prim_r2(1:10,1:n_prim)
REAL(dbl)     :: d3ORB_contribution(1:28,1:n_orbs)

DO ipt=1,npts
   !Evaluate the values of the Gaussian primitives and gradients at this point
   DO k_prim=1,n_prim
      call dGau_prim(prim_r1(k_prim),dprim_r1(1:3,k_prim),k_prim,R6d(1:3,ipt))
      call d3Gau_prim(prim_r2(k_prim),dprim_r2(1:3,k_prim),d2prim_r2(1:6,k_prim), &
                      d3prim_r2(1:10,k_prim),k_prim,R6d(4:6,ipt))
   ENDDO

   !We will start by computing the alpha-spin density matrix.
   d3orb_contribution(:,:) = 0.0_dbl
   DO j=1,n_orbs
      IF (MO_occ(j,1) > EPSILON(MO_occ(j,1))) THEN
         !We could do a two-step (well, now multi-step) method similar to what we
         !used in dDM1_wfn to evaluate the orbital contributions.  However, the 
         !complexity becomes daunting now, and it seems better to favor
         !readability over speed.  HOWEVER, THIS IS A PLACE THE PROGRAM COULD BE
         !SPED UP A LOT!!!
         !                         1 -- d^3/dx^3 D(r',r)                                 !
         !                         2 -- d^3/dy^3 D(r',r)                                 !
         !                         3 -- d^3/dz^3 D(r',r)                                 !
         !                         4 -- d^2/dx^2 d/dy D(r',r)                            !
         !                         5 -- d^2/dx^2 d/dz D(r',r)                            !
         !                         6 -- d^2/dy^2 d/dx D(r',r)                            !
         !                         7 -- d^2/dy^2 d/dz D(r',r)                            !
         !                         8 -- d^2/dz^2 d/dx D(r',r)                            !
         !                         9 -- d^2/dz^2 d/dy D(r',r)                            !
         !                        10 -- d/dx d/dy d/dz D(r',r)                           !

         FORALL(L=1:10)
               !For all components of the gradient:
               d3orb_contribution(L,j) = DOT_PRODUCT(MO_coeff(:,j,1),d3prim_r2(L,:)) 
         ENDFORALL
         d3orb_contribution(1:10,j) = d3orb_contribution(1:10,j)                            &
                                     *MO_occ(j,1)*DOT_PRODUCT(MO_coeff(:,j,1),prim_r1(:))   
                                     
         !The next section of code could be sped up by a factor of two but it 
         !would have a significant "readability price."  So I (PWA) didn't bother.
         !                        11 -- d^2/dx^2 d/dx' D(r',r)                           !
         !                        12 -- d^2/dy^2 d/dy' D(r',r)                           !
         !                        13 -- d^2/dz^2 d/dz' D(r',r)                           !
         !                        14 -- d^2/dx^2 d/dy' D(r',r)                           !
         !                        15 -- d^2/dx^2 d/dz' D(r',r)                           !
         !                        16 -- d^2/dy^2 d/dx' D(r',r)                           !
         !                        17 -- d^2/dy^2 d/dz' D(r',r)                           !
         !                        18 -- d^2/dz^2 d/dx' D(r',r)                           !
         !                        19 -- d^2/dz^2 d/dy' D(r',r)                           !
         !                        20 -- d/dx d/dy d/dx' D(r',r)                          !
         !                        21 -- d/dx d/dz d/dx' D(r',r)                          !
         !                        22 -- d/dy d/dx d/dy' D(r',r)                          !
         !                        23 -- d/dy d/dz d/dy' D(r',r)                          !
         !                        24 -- d/dz d/dx d/dz' D(r',r)                          !
         !                        25 -- d/dz d/dy d/dz' D(r',r)                          !
         d3orb_contribution(11,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(1,:))   &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(1,:)) 
         d3orb_contribution(12,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(2,:))   &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(2,:)) 
         d3orb_contribution(13,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(3,:))   &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(3,:))
         d3orb_contribution(14,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(2,:))    &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(1,:)) 
         d3orb_contribution(15,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(3,:))    &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(1,:)) 
         d3orb_contribution(16,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(1,:))    &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(2,:))
         d3orb_contribution(17,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(3,:))   &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(2,:)) 
         d3orb_contribution(18,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(1,:))   &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(3,:)) 
         d3orb_contribution(19,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(2,:))   &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(3,:))
         d3orb_contribution(20,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(1,:))    &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(4,:)) 
         d3orb_contribution(21,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(1,:))    &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(5,:)) 
         d3orb_contribution(22,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(2,:))    &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(4,:))   
         d3orb_contribution(23,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(2,:))    &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(6,:)) 
         d3orb_contribution(24,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(3,:))    &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(5,:)) 
         d3orb_contribution(25,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(3,:))    &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(6,:)) 
         !                        26 -- d/dx' d/dy d/dz D(r',r)                          !
         !                        27 -- d/dx d/dy' d/dz D(r',r)                          !
         !                        28 -- d/dx d/dy d/dz' D(r',r)                          !
         d3orb_contribution(26,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(1,:))    &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(6,:)) 
         d3orb_contribution(27,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(2,:))    &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(5,:)) 
         d3orb_contribution(28,j) = MO_occ(j,1)                                           &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),dprim_r1(3,:))    &
                                          *DOT_PRODUCT(MO_coeff(:,j,1),d2prim_r2(4,:)) 
      ENDIF
   ENDDO
   FORALL(L=1:28)
         d3DM(L,ipt,1) = SUM(d3orb_contribution(L,:))
   ENDFORALL
   
   IF (.not. spinorbs) THEN
      !The second component of the density matrix is the same as the first.
      d3DM(:,ipt,2) = d3DM(:,ipt,1)
   ELSE
      !Evaluate the beta-spin part of the density matrix in the same way.
      DO j=1,n_orbs
         IF (MO_occ(j,2) > EPSILON(MO_occ(j,2))) THEN
            !                         1 -- d^3/dx^3 D(r',r)                                 !
            !                         2 -- d^3/dy^3 D(r',r)                                 !
            !                         3 -- d^3/dz^3 D(r',r)                                 !
            !                         4 -- d^2/dx^2 d/dy D(r',r)                            !
            !                         5 -- d^2/dx^2 d/dz D(r',r)                            !
            !                         6 -- d^2/dy^2 d/dx D(r',r)                            !
            !                         7 -- d^2/dy^2 d/dz D(r',r)                            !
            !                         8 -- d^2/dz^2 d/dx D(r',r)                            !
            !                         9 -- d^2/dz^2 d/dy D(r',r)                            !
            !                        10 -- d/dx d/dy d/dz D(r',r)                           !
            FORALL(L=1:10)
                  !For all components of the gradient:
                  d3orb_contribution(L,j) = DOT_PRODUCT(MO_coeff(:,j,2),d3prim_r2(L,:)) 
            ENDFORALL
            d3orb_contribution(1:10,j) = d3orb_contribution(1:10,j)                            &
                                        *MO_occ(j,2)*DOT_PRODUCT(MO_coeff(:,j,2),prim_r1(:))   
                                        
            !The next section of code could be sped up by a factor of two but it 
            !would have a significant "readability price."  So I (PWA) didn't bother.
            !                        11 -- d^2/dx^2 d/dx' D(r',r)                           !
            !                        12 -- d^2/dy^2 d/dy' D(r',r)                           !
            !                        13 -- d^2/dz^2 d/dz' D(r',r)                           !
            !                        14 -- d^2/dx^2 d/dy' D(r',r)                           !
            !                        15 -- d^2/dx^2 d/dz' D(r',r)                           !
            !                        16 -- d^2/dy^2 d/dx' D(r',r)                           !
            !                        17 -- d^2/dy^2 d/dz' D(r',r)                           !
            !                        18 -- d^2/dz^2 d/dx' D(r',r)                           !
            !                        19 -- d^2/dz^2 d/dy' D(r',r)                           !
            !                        20 -- d/dx d/dy d/dx' D(r',r)                          !
            !                        21 -- d/dx d/dz d/dx' D(r',r)                          !
            !                        22 -- d/dy d/dx d/dy' D(r',r)                          !
            !                        23 -- d/dy d/dz d/dy' D(r',r)                          !
            !                        24 -- d/dz d/dx d/dz' D(r',r)                          !
            !                        25 -- d/dz d/dy d/dz' D(r',r)                          !
            d3orb_contribution(11,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(1,:))   &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(1,:)) 
            d3orb_contribution(12,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(2,:))   &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(2,:)) 
            d3orb_contribution(13,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(3,:))   &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(3,:))
            d3orb_contribution(14,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(2,:))    &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(1,:)) 
            d3orb_contribution(15,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(3,:))    &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(1,:)) 
            d3orb_contribution(16,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(1,:))    &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(2,:))
            d3orb_contribution(17,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(3,:))   &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(2,:)) 
            d3orb_contribution(18,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(1,:))   &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(3,:)) 
            d3orb_contribution(19,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(2,:))   &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(3,:))
            d3orb_contribution(20,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(1,:))    &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(4,:)) 
            d3orb_contribution(21,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(1,:))    &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(5,:)) 
            d3orb_contribution(22,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(2,:))    &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(4,:))   
            d3orb_contribution(23,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(2,:))    &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(6,:)) 
            d3orb_contribution(24,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(3,:))    &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(5,:)) 
            d3orb_contribution(25,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(3,:))    &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(6,:)) 
            !                        26 -- d/dx' d/dy d/dz D(r',r)                          !
            !                        27 -- d/dx d/dy' d/dz D(r',r)                          !
            !                        28 -- d/dx d/dy d/dz' D(r',r)                          !
            d3orb_contribution(26,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(1,:))    &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(6,:)) 
            d3orb_contribution(27,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(2,:))    &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(5,:)) 
            d3orb_contribution(28,j) = MO_occ(j,2)                                           &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),dprim_r1(3,:))    &
                                             *DOT_PRODUCT(MO_coeff(:,j,2),d2prim_r2(4,:)) 
         ENDIF
      ENDDO
      FORALL(L=1:28)
            d3DM(L,ipt,2) = SUM(d3orb_contribution(L,:))
      ENDFORALL
   ENDIF
ENDDO

end subroutine d3DM1_wfn


END MODULE process_wfn




