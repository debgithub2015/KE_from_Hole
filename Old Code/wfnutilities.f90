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
! mo_energy(1:n_orb,1:2) -- molecular orbital energies for 1. alpha and 2. beta-!
!                           spin orbitals                                       !
! mo_coeff(1:n_prim,1:n_orb,1:2) -- molecular orbital coefficients of the spin  !
!                                   molecular orbitals (1 = alpha; 2 = beta)    !
! mo_occ(1:n_orb,1:2) -- molecular orbital occupation numbers for orbitals of   !
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

INTEGER(istd)   :: n_orb,n_prim,n_atoms,n_el,multiplicity
INTEGER(istd), ALLOCATABLE :: basis_center(:), basis_type(:), basis_exp(:)
INTEGER(istd), ALLOCATABLE :: basis_ijk(:,:)
INTEGER(istd), ALLOCATABLE :: atm_prim(:,:),Zatom(:)

REAL(dbl) :: Energy,virial,KE
REAL(dbl), ALLOCATABLE :: mo_occ(:,:),mo_coeff(:,:,:),mo_energy(:,:)
REAL(dbl), ALLOCATABLE :: atm_position(:,:),atm_charge(:)

LOGICAL :: spinorbs
 
CHARACTER(len=80) :: wfn_filename

END MODULE variables_wfn

!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!###############################################################################!
!*******************************************************************************!
!*******************************************************************************!
!*******************************************************************************!
!###############################################################################!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------!
!                               utilities_wfn                                   !
!                                                                               !
! This module contains "utility routines" needed to process at wfn file.        !
!                                                                               !
! INCLUDED SUBROUTINES:                                                         !
!       ajkfds;ajdf;                                                            !
!-------------------------------------------------------------------------------!

MODULE utilities_wfn

USE kinds

CONTAINS

!-------------------------------------------------------------------------------!
!                               Gau_prim                                        !
!                                                                               !
! This subroutine evaluates a Gaussian primitive as follows:                    !
!                   x^i * y^j * z^k * EXP(alpha(i)*(R-R0){dot}(R-R0))           !
!-------------------------------------------------------------------------------!
!                             DICTIONARY OF LOCAL VARIABLES                     !
!                                                                               !
! sqdist -- the square of the distance from the point where the Gaussian is to  !
!           be evaluated to the atom where it is centered.                      !
! Rvec -- the vector displacement of the point from the center of the primitive !
! b -- the number of the Gaussian primitive that is being evaluated.            !
! point(1:3) -- the point where the Gaussian is being evaluated.                !
! ijk(1:3) -- the exponents of x, y, and z in the Cartesian prefactor of the    !
!             Gaussian.                                                         !
! alpha -- the exponential prefactor of R**2 in the Gaussian.                   !
!-------------------------------------------------------------------------------!
!                           DICTIONARY OF VARIABLES_WFN                         !
!-------------------------------------------------------------------------------!
! atm_position(1:3,1:n_atoms) -- the position of each atom.                     !
! basis_center(1:n_prim) -- the center to which the gaussian primitive is       !
!                           assigned.                                           !
! basis_exp(1:n_prim) -- exponents of the gaussian primitives.                  !
! basis_ijk(1:3,1:nprim) -- assigns exponents to each Gaussian primitive        !
! n_prim -- number of Gaussian basis functions.                                 !
!-------------------------------------------------------------------------------!

subroutine Gau_prim(prim,point)

USE kinds
USE variables_wfn

IMPLICIT NONE

REAL(dbl), INTENT(OUT) :: prim(1:n_prim)
REAL(dbl), INTENT(IN)  :: point(1:3)
REAL(dbl)     :: alpha(1:n_prim)
INTEGER(istd) :: ijk(1:3,1:n_prim),b
REAL(dbl)     :: point(1:3),Rvec(1:3,1:n_prim),sqdist(1:n_prim)

FORALL(b=1:n_prim)
      Rvec(:,b) = point(:) - atm_position(:,basis_center(b))
      ijk(1:3,b) = basis_ijk(1:3,b)
      alpha(b) = -1*basis_exp(b)
ENDFORALL

!Compute square of distance from point to center of each basis function.      
FORALL(b=1:nprim)
      sqdist(b) = DOT_PRODUCT(Rvec(:,b),Rvec(:,b))
ENDFORALL
     
FORALL(b=1:nprim)     
      prim(b) = EXP(alpha(b)*sqdist(b))                                        &
                *Rvec(1,b)**ijk(1,b)*Rvec(2,b)**ijk(2,b)*Rvec(3,b)**ijk(3,b)
ENDFORALL

end subroutine Gau_prim

!-------------------------------------------------------------------------------!
!                               dGau_prim                                       !
!                                                                               !
! This subroutine evaluates the gradient of a Gaussian primitive as follows:    !
!      prim = G(x,y,z) = x^i * y^j * z^k * EXP(alpha(i)*(R-R0){dot}(R-R0))      !
!  dprim(1) = d/dx  G(x,y,z)                                                    !
!  dprim(2) = d/dy  G(x,y,z)                                                    !
!  dprim(3) = d/dz  G(x,y,z)                                                    !
!-------------------------------------------------------------------------------!
!                             DICTIONARY OF LOCAL VARIABLES                     !
!                                                                               !
! sqdist -- the square of the distance from the point where the Gaussian is to  !
!           be evaluated to the atom where it is centered.                      !
! Rvec -- the vector displacement of the point from the center of the primitive !
! i_prim -- the number of the Gaussian primitive that is being evaluated.       !
! point(1:3) -- the point where the Gaussian is being evaluated.                !
! exp_piece -- the exponential piece of the term.                               !
! dxyz(1:3) -- the polynomial piece of the term.                                !
! ider -- the term in the derivative (1,2, or 3) being evaluated.               !
! d -- a counter for the dimensions.                                            !
! ijk(1:3) -- the exponents of x, y, and z in the Cartesian prefactor of the    !
!             Gaussian.                                                         !
! alpha -- the exponential prefactor of R**2 in the Gaussian.                   !
! prim -- the value of the Gaussian primitive.                                  !
! dprim(1:3) -- the gradientof the Gaussian primitive.                          !
!-------------------------------------------------------------------------------!
!                           DICTIONARY OF VARIABLES_WFN                         !
!-------------------------------------------------------------------------------!
! atm_position(1:3,1:n_atoms) -- the position of each atom.                     !
! basis_center(1:n_prim) -- the center to which the gaussian primitive is       !
!                           assigned.                                           !
! basis_exp(1:n_prim) -- exponents of the gaussian primitives.                  !
! basis_ijk(1:3,1:nprim) -- assigns exponents to each Gaussian primitive        !
!-------------------------------------------------------------------------------!

subroutine dGau_prim(prim,dprim,iprim,point)

USE kinds; USE variables_wfn

IMPLICIT NONE

REAL(dbl) :: dprim(1:3),exp_piece,dxyz(1:3),alpha,prim
INTEGER(istd) :: i_prim,ider,d,ijk(1:3)
REAL(dbl)     :: point(1:3),Rvec(1:3)

Rvec(:) = point(:) - atm_position(:,basis_center(i_prim))
sqdist = DOT_PRODUCT(Rvec,Rvec)
ijk(1:3) = basis_ijk(1:3,i_prim)
alpha = -1*basis_exp(i_prim)
dxyz(:) = 1.0_dbl
d2xyz(:) = 1.0_dbl
exp_piece = EXP(alpha*sqdist)

prim = exp_piece * Rvec(1)**ijk(1) * Rvec(2)**ijk(2) * Rvec(3)**ijk(3)

!To avoid divide-by-zero errors, make sure Rvec is not too small.
IF (sqdist < EPSILON(sqdist)**4) THEN
   Rvec(1:3) = EPSILON(sqdist)**2
ENDIF

DO ider=1,3
   poly_piece = 1.0_dbl
   DO d=1,3
      IF (d /= ider) THEN
         dxyz(ider) = dxyz(ider)*(Rvec(d)**ijk(d))
      ELSE
         dxyz(ider) = dxyz(ider)                                               &
                      *(2*alpha*Rvec(ider)**(ijk(ider)+1)                      &
                         + ijk(ider)*Rvec(ider)**(ijk(ider)-1)
      ENDIF
   ENDDO
ENDDO

dprim(:) = exp_piece*dxyz(:)
         
end subroutine dprim

!-------------------------------------------------------------------------------!
!                               d2Gau_prim                                      !
!                                                                               !
! This function evaluates the second-derivative of a Gaussian primitive as follows:!
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
!                             DICTIONARY OF LOCAL VARIABLES                     !
!                                                                               !
! sqdist -- the square of the distance from the point where the Gaussian is to  !
!           be evaluated to the atom where it is centered.                      !
! Rvec -- the vector displacement of the point from the center of the primitive !
! i_prim -- the number of the Gaussian primitive that is being evaluated.       !
! point(1:3) -- the point where the Gaussian is being evaluated.                !
! exp_piece -- the exponential piece of the term.                               !
! poly_piece -- the polynomial piece of the term.                               !
! ider -- the term in the derivative (1,2, or 3) being evaluated.               !
! d -- a counter for the dimensions.                                            !
! ijk(1:3) -- the exponents of x, y, and z in the Cartesian prefactor of the    !
!             Gaussian.                                                         !
! alpha -- the exponential prefactor of R**2 in the Gaussian.                   !
! dxyz(1:3) -- the polynomial part of the first derivative.                     !
! d2xyz(1:3) -- the polynomial part of the second derivative.                   !
! prim -- the value of the Gaussian primitive.                                  !
! dprim(1:3) -- the first derivative of the Gaussian primitive.                 !
! d2prim(1:3) -- the second derivative of the Gaussian primitive.               !
!-------------------------------------------------------------------------------!
!                           DICTIONARY OF VARIABLES_WFN                         !
!-------------------------------------------------------------------------------!
! atm_position(1:3,1:n_atoms) -- the position of each atom.                     !
! basis_center(1:n_prim) -- the center to which the gaussian primitive is       !
!                           assigned.                                           !
! basis_exp(1:n_prim) -- exponents of the gaussian primitives.                  !
! basis_ijk(1:3,1:nprim) -- assigns exponents to each Gaussian primitive        !
!-------------------------------------------------------------------------------!

subroutine d2Gau_prim(prim,dprim,d2prim,iprim,point)

USE kinds; USE variables_wfn

IMPLICIT NONE

REAL(dbl) :: dprim(1:3),d2prim(1:6),exp_piece,dxyz(1:3),alpha,prim,d2xyz(1:3)
INTEGER(istd) :: i_prim,ider,d,ijk(1:3)
REAL(dbl)     :: point(1:3),Rvec(1:3)

Rvec(:) = point(:) - atm_position(:,basis_center(i_prim))
sqdist = DOT_PRODUCT(Rvec,Rvec)
ijk(1:3) = basis_ijk(1:3,i_prim)
dxyz(:) = 1.0_dbl
d2xyz(:) = 1.0_dbl
alpha = -1*basis_exp(i_prim)

exp_piece = EXP(alpha*sqdist)
prim = exp_piece * Rvec(1)**ijk(1) * Rvec(2)**ijk(2) * Rvec(3)**ijk(3)

!To avoid divide-by-zero errors, make sure Rvec is not too small.
IF (sqdist < EPSILON(sqdist)**4) THEN
   Rvec(1:3) = EPSILON(sqdist)**2
ENDIF

DO ider=1,3
   DO d=1,3
      IF (d /= ider) THEN
         dxyz(ider) = dxyz(ider)*(Rvec(d)**ijk(d))
      ELSE
         dxyz(ider) = dxyz(ider)                                               &
                      *(2*alpha*Rvec(ider)**(ijk(ider)+1)                      &
                         + ijk(ider)*Rvec(ider)**(ijk(ider)-1))
      ENDIF
   ENDDO
ENDDO

dprim(:) = exp_piece*dxyz(:)
         
DO ider=1,3
   DO d=1,3
      IF (d /= ider) THEN
         d2xyz(ider) = d2xyz(ider)*(Rvec(d)**ijk(d))
      ELSE
         d2xyz(ider) = d2xyz(ider)                                             &
                      *(4*alpha**2 * Rvec(ider)**(ijk(ider)+2)                 &
                         + 2*alpha(2*ijk(ider)+1)*Rvec(ider)**(ijk(ider))      &
                         + ijk(ider)*(ijk(ider)-1)*Rvec(ider)**(ijk(ider)-2))
      ENDIF
   ENDDO
ENDDO

d2prim(1:3) = d2xyz(1:3)*exp_piece
d2prim(4) = dxyz(1)*dxyz(2)*exp_piece
d2prim(5) = dxyz(1)*dxyz(3)*exp_piece
d2prim(6) = dxyz(2)*dxyz(3)*exp_piece

end subroutine d2Gau_prim

!-------------------------------------------------------------------------------!
!                               d3Gau_prim                                      !
!                                                                               !
! This function evaluates the 3rd derivative of a Gaussian primitive as follows:!
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
!-------------------------------------------------------------------------------!
!                             DICTIONARY OF LOCAL VARIABLES                     !
!                                                                               !
! sqdist -- the square of the distance from the point where the Gaussian is to  !
!           be evaluated to the atom where it is centered.                      !
! Rvec -- the vector displacement of the point from the center of the primitive !
! i_prim -- the number of the Gaussian primitive that is being evaluated.       !
! point(1:3) -- the point where the Gaussian is being evaluated.                !
! exp_piece -- the exponential piece of the term.                               !
! poly_piece -- the polynomial piece of the term.                               !
! ider -- the term in the derivative (1,2, or 3) being evaluated.               !
! d -- a counter for the dimensions.                                            !
! ijk(1:3) -- the exponents of x, y, and z in the Cartesian prefactor of the    !
!             Gaussian.                                                         !
! alpha -- the exponential prefactor of R**2 in the Gaussian.                   !
! dxyz(1:3) -- the polynomial part of the first derivative.                     !
! d2xyz(1:3) -- the polynomial part of the second derivative.                   !
! d3xyz(1:3) -- the polynomial part of the third derivative.                    !
! prim -- the value of the Gaussian primitive.                                  !
! dprim(1:3) -- the first derivative of the Gaussian primitive.                 !
! d2prim(1:3) -- the second derivative of the Gaussian primitive.               !
! d3prim(1:10) -- the third derivative of the Gaussian primitive.               !
!-------------------------------------------------------------------------------!
!                           DICTIONARY OF VARIABLES_WFN                         !
!-------------------------------------------------------------------------------!
! atm_position(1:3,1:n_atoms) -- the position of each atom.                     !
! basis_center(1:n_prim) -- the center to which the gaussian primitive is       !
!                           assigned.                                           !
! basis_exp(1:n_prim) -- exponents of the gaussian primitives.                  !
! basis_ijk(1:3,1:nprim) -- assigns exponents to each Gaussian primitive        !
!-------------------------------------------------------------------------------!

subroutine d3Gau_prim(prim,dprim,d2prim,d3prim,iprim,point)

USE kinds; USE variables_wfn

IMPLICIT NONE

REAL(dbl) :: prim,dprim(1:3),d2prim(1:6),d3prim(1:10),exp_piece
REAL(dbl) :: dxyz(1:3),d2xyz(1:3),d3xyz(1:3)
INTEGER(istd) :: i_prim,ider,d,ijk(1:3)
REAL(dbl)     :: point(1:3),Rvec(1:3)

Rvec(:) = point(:) - atm_position(:,basis_center(i_prim))
sqdist = DOT_PRODUCT(Rvec,Rvec)
ijk(1:3) = basis_ijk(1:3,i_prim)
dxyz(:) = 1.0_dbl
d2xyz(:) = 1.0_dbl
d3xyz(:) = 1.0_dbl
alpha = -1*basis_exp(i_prim)

exp_piece = EXP(alpha*sqdist)
prim = exp_piece * Rvec(1)**ijk(1) * Rvec(2)**ijk(2) * Rvec(3)**ijk(3)

!To avoid divide-by-zero errors, make sure Rvec is not too small.
IF (sqdist < EPSILON(sqdist)**4) THEN
   Rvec(1:3) = EPSILON(sqdist)**2
ENDIF

DO ider=1,3
   DO d=1,3
      IF (d /= ider) THEN
         dxyz(ider) = dxyz(ider)*(Rvec(d)**ijk(d))
      ELSE
         dxyz(ider) = dxyz(ider)                                               &
                      *(2*alpha*Rvec(ider)**(ijk(ider)+1)                      &
                         + ijk(ider)*Rvec(ider)**(ijk(ider)-1))
      ENDIF
   ENDDO
ENDDO

dprim(:) = exp_piece*dxyz(:)
         
DO ider=1,3
   DO d=1,3
      IF (d /= ider) THEN
         d2xyz(ider) = d2xyz(ider)*(Rvec(d)**ijk(d))
      ELSE
         d2xyz(ider) = d2xyz(ider)                                              &
                      *(4*alpha**2 * Rvec(ider)**(ijk(ider)+2)                  &
                         + 2*alpha(2*ijk(ider)+1) * Rvec(ider)**(ijk(ider))     &
                         + ijk(ider)*(ijk(ider)-1)*Rvec(ider)**(ijk(ider)-2))
      ENDIF
   ENDDO
ENDDO

d2prim(1:3) = d2xyz(1:3)*exp_piece
d2prim(4) = dxyz(1)*dxyz(2)*exp_piece
d2prim(5) = dxyz(1)*dxyz(3)*exp_piece
d2prim(6) = dxyz(2)*dxyz(3)*exp_piece

DO ider=1,3
   DO d=1,3
      IF (d /= ider) THEN
         d3xyz(ider) = d3xyz(ider)*(Rvec(d)**ijk(d))
      ELSE
         d3xyz(ider) = d3xyz(ider)                                              &
                      *(8*alpha**3 * Rvec(ider)**(ijk(ider)+3)                  &
                         + 4*alpha**2*(2*ijk(ider)+1)*(ijk(ider)+2)             &
                                        *Rvec(ider)**(ijk(ider)+1)              &
                         + 2*alpha*ijk(ider)*(2*ijk(ider)+1)*(ijk(ider)-1)      &
                                        *Rvec(ider)**(ijk(ider)-1)              &
                         + ijk(ider)*(ijk(ider)-1)*(ijk(ider)-2)                &
                                        *Rvec(ider)**(ijk(ider)-3))
      ENDIF
   ENDDO
ENDDO

d3prim(1:3) = d3xyz(1:3)*exp_piece
d3prim(4) = d2xyz(1)*dxyz(2)*exp_piece
d3prim(5) = d2xyz(1)*dxyz(3)*exp_piece
d3prim(6) = d2xyz(2)*dxyz(1)*exp_piece
d3prim(7) = d2xyz(2)*dxyz(3)*exp_piece
d3prim(8) = d2xyz(3)*dxyz(1)*exp_piece
d3prim(9) = d2xyz(3)*dxyz(2)*exp_piece
d3prim(10) = dxyz(1)*dxyz(2)*dxyz(3)*exp_piece

end subroutine d3Gau_prim

!-------------------------------------------------------------------------------!
!                           primitive_type_assignments                          !
! This subroutine assigns the exponents, i,j,k, for the Cartesian prefactors    !
! of a Gaussian primitive                                                       !
!                   x^i * y^j * z^k * EXP(alpha(i)*(R-R0){dot}(R-R0))           !
! Watch out that you wfn writer agrees with the order of a regular wfn file !!  !
! for each type N we store the exponents 1,2,3 for x,y,z                        !
!-------------------------------------------------------------------------------!
! type_exponents(1:3,1:20) -- the values of (i,j,k) for the Cartesian           !
!                         prefactor, x**i * y**j * z**k, of the Gausian         !
!                         function in the primitive basis functions.            !
!-------------------------------------------------------------------------------!


subroutine primitive_type_assignments(type_exponents)

USE kinds

implicit none

REAL(dbl) :: type_exponents(1:3,1:20)

type_exponents(1,1)=0.0
type_exponents(2,1)=0.0
type_exponents(3,1)=0.0

type_exponents(1,2)=1.0
type_exponents(2,2)=0.0
type_exponents(3,2)=0.0

type_exponents(1,3)=0.0
type_exponents(2,3)=1.0
type_exponents(3,3)=0.0

type_exponents(1,4)=0.0
type_exponents(2,4)=0.0
type_exponents(3,4)=1.0

type_exponents(1,5)=2.0
type_exponents(2,5)=0.0
type_exponents(3,5)=0.0

type_exponents(1,6)=0.0
type_exponents(2,6)=2.0
type_exponents(3,6)=0.0

type_exponents(1,7)=0.0
type_exponents(2,7)=0.0
type_exponents(3,7)=2.0

type_exponents(1,8)=1.0
type_exponents(2,8)=1.0
type_exponents(3,8)=0.0

type_exponents(1,9)=1.0
type_exponents(2,9)=0.0
type_exponents(3,9)=1.0

type_exponents(1,10)=0.0
type_exponents(2,10)=1.0
type_exponents(3,10)=1.0

type_exponents(1,11)=3.0
type_exponents(2,11)=0.0
type_exponents(3,11)=0.0

type_exponents(1,12)=0.0
type_exponents(2,12)=3.0
type_exponents(3,12)=0.0

type_exponents(1,13)=0.0
type_exponents(2,13)=0.0
type_exponents(3,13)=3.0

type_exponents(1,14)=2.0
type_exponents(2,14)=1.0
type_exponents(3,14)=0.0

type_exponents(1,15)=2.0
type_exponents(2,15)=0.0
type_exponents(3,15)=1.0

type_exponents(1,16)=0.0
type_exponents(2,16)=2.0
type_exponents(3,16)=1.0

type_exponents(1,17)=1.0
type_exponents(2,17)=2.0
type_exponents(3,17)=0.0

type_exponents(1,18)=1.0
type_exponents(2,18)=0.0
type_exponents(3,18)=2.0

type_exponents(1,19)=0.0
type_exponents(2,19)=1.0
type_exponents(3,19)=2.0

type_exponents(1,20)=1.0
type_exponents(2,20)=1.0
type_exponents(3,20)=1.0

end subroutine primitive_type_assigments

!Function to find the determinant of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Downloaded by Paul Ayers from:
!   http://www.dreamincode.net/code/snippet1273.htm
!and slightly modified (with kinds, turning a function into a subroutine, 
! using double precision, etc.)
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!
subroutine FindDet(determinant,matrix, n)
    USE kinds
	IMPLICIT NONE
	REAL(dbl), DIMENSION(n,n) :: matrix
	INTEGER(istd), INTENT(IN) :: n
	REAL(dbl) :: m, temp
	INTEGER(istd) :: i, j, k, l
	LOGICAL :: DetExists = .TRUE.
	l = 1
	!Convert to upper triangular form
	DO k = 1, n-1
		IF (matrix(k,k) == 0) THEN
			DetExists = .FALSE.
			DO i = k+1, n
				IF (matrix(i,k) /= 0) THEN
					DO j = 1, n
						temp = matrix(i,j)
						matrix(i,j)= matrix(k,j)
						matrix(k,j) = temp
					END DO
					DetExists = .TRUE.
					l=-l
					EXIT
				ENDIF
			END DO
			IF (DetExists .EQV. .FALSE.) THEN
				determinant = 0
				return
			END IF
		ENDIF
		DO j = k+1, n
			m = matrix(j,k)/matrix(k,k)
			DO i = k+1, n
				matrix(j,i) = matrix(j,i) - m*matrix(k,i)
			END DO
		END DO
	END DO
	
	!Calculate determinant by finding product of diagonal elements
	determinant = l
	DO i = 1, n
		determinant = determinant * matrix(i,i)
	END DO
	
END subroutine FindDet


END MODULE utilities_wfn
