module diagonalize
  use iso_fortran_env, only: real64
  implicit none
  private
  public :: diagonalize_hamiltonian ! zheev_wrapper

  integer, parameter :: dp = real64
  integer, parameter :: cp = dp

  ! Explicit interface for the external LAPACK routine ZHEEV
  interface
    subroutine zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
      import :: dp, cp
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      complex(kind=cp), intent(inout) :: a(lda,*)
      real(kind=dp), intent(out) :: w(*)
      complex(kind=cp), intent(out) :: work(*)
      real(kind=dp), intent(out) :: rwork(*)
      integer, intent(out) :: info
    end subroutine zheev
  end interface

contains

  ! Wrapper for LAPACK's ZHEEV routine
  subroutine zheev_wrapper(A_in, n, eigenvalues, eigenvectors, info)
    ! Use the defined kinds
    complex(kind=cp), intent(in) :: A_in(:,:)           ! Matrix to diagonalize
    integer, intent(in) :: n                            ! Matrix dimension
    real(kind=dp), intent(out) :: eigenvalues(:)        ! Array to store eigenvalues
    complex(kind=cp), intent(out) :: eigenvectors(:,:)  ! Matrix to store eigenvectors
    integer, intent(out) :: info                        ! Success indicator

    integer :: lwork
    complex(kind=cp), allocatable :: work(:)
    real(kind=dp), allocatable :: rwork(:)
    character(1) :: jobz = 'V', uplo = 'U'
    complex(kind=cp), allocatable :: matrix_for_zheev(:,:) ! Local copy for zheev

    ! Allocate matrix that zheev will use and overwrite
    allocate(matrix_for_zheev(n,n))
    matrix_for_zheev = A_in ! Copy original data here

    ! Determine optimal workspace size using the copy
    lwork = -1
    allocate(work(1), rwork(max(1, 3*n-2)))
    ! Pass matrix_for_zheev to query call
    ! NOTE: Pass n as LDA since matrix_for_zheev has dimensions (n,n)
    call zheev(jobz, uplo, n, matrix_for_zheev, n, eigenvalues, work, lwork, rwork, info)

    ! Check if query itself failed
    if (info /= 0) then
        write(*,*) "ZHEEV workspace query failed, info = ", info
        deallocate(work, rwork, matrix_for_zheev)
        return ! Exit subroutine on query failure
    endif

    ! Allocate workspace with optimal size
    lwork = int(real(work(1), kind=dp))
    deallocate(work)
    allocate(work(lwork))

    ! Re-copy original matrix A_in into matrix_for_zheev
    matrix_for_zheev = A_in

    ! Diagonalize the matrix using the local copy
    ! zheev will overwrite matrix_for_zheev with eigenvectors
    ! NOTE: Pass n as LDA
    call zheev(jobz, uplo, n, matrix_for_zheev, n, eigenvalues, work, lwork, rwork, info)

    ! If successful, copy the eigenvectors to the output argument
    if (info == 0) then
       eigenvectors = matrix_for_zheev
    endif

    ! Clean up
    deallocate(work, rwork, matrix_for_zheev)

  end subroutine zheev_wrapper

  ! Diagonalize Hamiltonian matrix and return eigenvalues and eigenvectors
  subroutine diagonalize_hamiltonian(H, eigenvalues, eigenvectors, info)
    ! Use the defined kinds
    complex(kind=cp), intent(in) :: H(:,:)
    real(kind=dp), allocatable, intent(out) :: eigenvalues(:)
    complex(kind=cp), allocatable, intent(out) :: eigenvectors(:,:)
    integer, intent(out) :: info

    integer :: n

    n = size(H, 1)
    if (size(H, 2) /= n) then
        write(*,*) "Error: Input Hamiltonian matrix H must be square."
        info = -99 ! Assign a custom error code
        if (n > 0) then
           allocate(eigenvalues(n), eigenvectors(n,n))
        else
           allocate(eigenvalues(0), eigenvectors(0,0))
        endif
        eigenvalues = -huge(1.0_dp)
        eigenvectors = cmplx(0.0_dp, 0.0_dp, kind=cp)
        return
    endif


    ! Allocate output arrays
    allocate(eigenvalues(n), eigenvectors(n,n))

    ! Call LAPACK wrapper (passing H directly)
    call zheev_wrapper(H, n, eigenvalues, eigenvectors, info)

    ! Handle errors from zheev_wrapper
    if (info /= 0) then
      write(*,*) "Error during diagonalization: ZHEEV info = ", info
    end if
  end subroutine diagonalize_hamiltonian

end module diagonalize