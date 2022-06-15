! This MODULE contains the following SUBROUTINES:
!
! SUBROUTINE EigenMatrix(m,A,EigenValues,EigenVectors)
!
! This subroutine computes the eigenvalues and eigenvectors of a 
! given square matrix. The matrix has to be symmetric.
! This subrotuine uses the following LAPACK routine 
! SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
! For more detailed information about this subroutine please follow
! the link http://netlib.org/lapack/double/dsyev.f
! 

MODULE CalcMatrix

CONTAINS


!=============================================================================
!=============================================================================

!
! This subroutine computes the eigenvalues and eigenvectors of a 
! given square matrix. The matrix has to be symmetric
! 

SUBROUTINE EigenMatrix(m,A,EigenValues,EigenVectors)


USE nrtype			! It contains the types of used variables


IMPLICIT NONE


INTEGER, INTENT(IN) :: m						!(input) Order of the matrix
REAL(dp), INTENT(IN) :: A(1:m,1:m)				!(input) The given square matrix
REAL(dp), INTENT(OUT) :: EigenValues(1:m)		!(output) Array containing the eigenvalues of the matrix
REAL(dp), INTENT(OUT) :: EigenVectors(1:m,1:m)	!(output) Matrix containing the eigenvectors of the matrix


REAL(dp)   :: a_copy(m,m)	! Auxiliary square matrix. Just a copy of the original given matrix
INTEGER(idbl) :: i			! Just a counter
INTEGER(idbl) :: info		! Integer requiered by the LAPACK subroutine
INTEGER(idbl) :: lda		! Integer requiered by the LAPACK subroutine
character :: jobz			! Character requiered by the LAPACK subroutine
character :: uplo			! Character requiered by the LAPACK subroutine
INTEGER(idbl) :: lwork		! Integer requiered by the LAPACK subroutine

real(dp)  , allocatable, dimension ( : ) :: work  ! Allocatable array required by the LAPACK subroutine

lwork = 5*m ! The value of lwork is set as sugested by the 
allocate ( work(1:lwork) )	! The array work is allocated to be of dimension lwork
!
!
!

jobz = 'V' !Character is set according to the LAPACK routine specification
uplo = 'U' !Character is set according to the LAPACK routine specification

! Value of lda is set to be same as the matrix order
lda = m 

! The copy of the matrix is obtained
a_copy(1:m,1:m) = a(1:m,1:m)


!!===============================================================================================
!!===============================================================================================
!!
!!
!!  LAPACK IS CALLED OUT 
!!

   call dsyev( jobz, uplo, m, a_copy, lda, EigenValues, work, lwork, info)

!!
!!
!!===============================================================================================
!!===============================================================================================

 !info = 0

! This IF  condition is coded so that some info is displayed in case
! the LAPACK routine fails
  if ( info /= 0 ) then
    write ( *,* ) ' '
    write ( *,* ) '  Eigen_lapack - Failure!'
    write ( *,* ) '  The Eigenvalues could not be calculated.'
    write ( *,* ) '  LAPACK routine DSYEV returned a nonzero'
    write ( *,* ) '  value of the error flag, INFO = ', info
    return
  end if

! Since on exit the eigenvectors are storaged in the original input matrix A, its value
! is transfered to the matrix of interest
EigenVectors(1:m,1:m) =  a_copy(1:m,1:m)

 
! The array work is deallocated
deallocate ( work )


! The file "Check_EigenVectors.out" is open.
! In this file the eigenvalue label, the eigenvalue and the error asociated to that
! eigenvalue and its eigenvector are printed out.
! The error is just the Sqrt(Norm(A*EigenVector-EigenValue*EigenVector))
! This file is recognized as UNIT 500
OPEN(UNIT=500,STATUS='UNKNOWN',FILE='Check_Eigenvectors.out', &
ACTION='write')
write(500,*) ' No. Eigenvalue --- EigenValue --- Error --- '
write(500,*) '														 '
DO i = 1, m
	write(500,505) i , EigenValues(i) , &
	Sqrt(Dot_Product(MatMul(A,EigenVectors(:,i))-(EigenValues(i)*EigenVectors(:,i)), &
	MatMul(A,EigenVectors(:,i))-(EigenValues(i)*EigenVectors(:,i))))
	505 FORMAT (I7,2E27.17,2E27.17)
END DO
CLOSE(500)
! The file is closed

  return





END SUBROUTINE EigenMatrix

!=============================================================================
!=============================================================================



END MODULE CalcMatrix