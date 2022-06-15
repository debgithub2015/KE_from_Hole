!***************************************************************************************!
!									  MODULE Kinds                                      !
!                                                                                       !
!                                       Paul Ayers                                      !
!                                                                                       !
!                                       July 4, 2000                                    !
!                                                                                       !
!     Provenance; last modified July 6, 2000. (streamlined module)                      !
!                 renamed and revised, May 12, 2002.                                    !
!                 redefined integer types, August 14, 2010                              !
!                                                                                       !
! This module assigns various integer and real variable types to meaningful expressions.!
!***************************************************************************************!

MODULE kinds

IMPLICIT NONE

! First we define three levels of integer precision
  INTEGER, PARAMETER ::  ishrt  = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER ::  ihalf = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER ::  istd = SELECTED_INT_KIND(8) 
  INTEGER, PARAMETER ::  idbl = MAX(SELECTED_INT_KIND(11),istd)
  INTEGER, PARAMETER ::  ibig = MAX(SELECTED_INT_KIND(16),idbl)

! Now we define 3 levels of real precision.  
  INTEGER, PARAMETER ::  sngl = SELECTED_REAL_KIND(p=5)
  INTEGER, PARAMETER ::  dbl  = SELECTED_REAL_KIND(p=12)
! We assume that requiring 18 digits of precision is enough to force the program 
! into more the double precision; when that exists.  If it does not exist,
! then SELECTED_REAL_KIND < 0, so trpl <= dbl.
  INTEGER, PARAMETER ::  trpl = MAX(SELECTED_REAL_KIND(p=18),dbl)
! We try for quadruple precision:
  INTEGER, PARAMETER ::  qdpl = MAX(SELECTED_REAL_KIND(p=25),trpl)

END MODULE kinds

!-----------------------------------------------------------------------!
!                               CONSTANTS                               !
!                                                                       !
! This module contains the values of key constants like pi.             !
!                                                                       !
!***********************************************************************!
! Dependencies:                                                         !
!    kinds -- module containing variable types.                         !
!                                                                       !
!-----------------------------------------------------------------------!
!-----------------------------------------------------------------------!
!                                                                       !
! Contains:                                                             !
!     (none)                                                            !
!                                                                       !
!     OBITER DICTUM: Variable-storage modules should not contain any    !
!                    subroutines.                                       !
!                                                                       !
!-----------------------------------------------------------------------!
! Authors:                                                              !
!     Juan I. Rodriguez (rodrigji@mcmaster.ca)                          !
!     Paul W. Ayers (ayers@mcmaster.ca)                                 !
!                                                                       !
! Provenance:                                                           !
!     May 1, 2007 (PWA) documentation of existing module.               !
!                                                                       !
!-----------------------------------------------------------------------!

MODULE constants

USE kinds

IMPLICIT NONE

SAVE

REAL(dbl) :: pi = 3.14159265358979323846_dbl
REAL(dbl) :: Ang2Bohr = .529177249_dbl
REAL(dbl) :: Bohr2Ang = 1.88972599_dbl

END MODULE constants



