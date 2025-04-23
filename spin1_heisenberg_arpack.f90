PROGRAM spin1_heisenberg_arpack
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0D0)

  ! --- ALL Declarations MUST come first ---
  INTEGER :: N, nnz, i, info
  REAL(dp), ALLOCATABLE :: H_data(:)
  INTEGER, ALLOCATABLE :: H_indices(:), H_indptr(:)

  ! ARPACK setup parameters (declarations with initialization are fine here)
  INTEGER :: nev = 1, ncv = 20, max_iter = 1000
  REAL(dp) :: tol = 1.0D-10, sigma = 0.0D0  ! sigma not strictly needed for mode 1 'SM', but harmless
  REAL(dp), ALLOCATABLE :: eigenvalues(:), eigenvectors(:,:)
  ! --- End of Declarations ---


  ! --- Executable statements start here ---
  PRINT *, "Reading CSR matrix from files..."
  CALL read_csr_matrix('spin1_meta.bin', 'spin1_data.bin', 'spin1_indices.bin', 'spin1_indptr.bin', &
                       N, nnz, H_data, H_indices, H_indptr)
  PRINT *, "Matrix read complete. N =", N, " NNZ =", nnz

  PRINT *, "Calling ARPACK..."
  CALL arpack_diagonalize(N, nnz, H_data, H_indices, H_indptr, &
                          'SA', nev, ncv, tol, max_iter, eigenvalues, eigenvectors, info)

  IF (info == 0) THEN
    PRINT *, "ARPACK converged."
    PRINT *, "Ground state energy:", eigenvalues(1)
    ! If you wanted to print the eigenvector (can be long):
    ! PRINT *, "Ground state vector (first few elements):"
    ! PRINT *, eigenvectors(1:MIN(10,N), 1)
  ELSE
    PRINT *, "ARPACK error. Info code:", info
    ! Add more detailed error messages based on ARPACK documentation for info values if needed
  END IF

CONTAINS

  SUBROUTINE read_csr_matrix(meta_file, data_file, indices_file, indptr_file, N, nnz, data, indices, indptr)
    CHARACTER(*), INTENT(IN) :: meta_file, data_file, indices_file, indptr_file
    INTEGER, INTENT(OUT) :: N, nnz
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: data(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: indices(:), indptr(:)
    INTEGER :: unit, io_stat

    ! Read metadata
    OPEN(NEWUNIT=unit, FILE=meta_file, FORM='UNFORMATTED', ACCESS='STREAM', STATUS='OLD', IOSTAT=io_stat)
    IF (io_stat /= 0) THEN
        PRINT *, "Error opening meta file: ", meta_file, " IOSTAT=", io_stat
        STOP 1
    END IF
    READ(unit, IOSTAT=io_stat) N, nnz
    IF (io_stat /= 0) THEN
        PRINT *, "Error reading from meta file: ", meta_file, " IOSTAT=", io_stat
        CLOSE(unit)
        STOP 1
    END IF
    CLOSE(unit)

    ! Allocate arrays based on metadata
    ALLOCATE(data(nnz), indices(nnz), indptr(N+1), STAT=io_stat)
     IF (io_stat /= 0) THEN
        PRINT *, "Error allocating arrays. N=", N, " NNZ=", nnz, " STAT=", io_stat
        STOP 1
    END IF

    ! Read data array
    OPEN(NEWUNIT=unit, FILE=data_file, FORM='UNFORMATTED', ACCESS='STREAM', STATUS='OLD', IOSTAT=io_stat)
     IF (io_stat /= 0) THEN
        PRINT *, "Error opening data file: ", data_file, " IOSTAT=", io_stat
        STOP 1
    END IF
    READ(unit, IOSTAT=io_stat) data
    IF (io_stat /= 0) THEN
        PRINT *, "Error reading from data file: ", data_file, " IOSTAT=", io_stat
        CLOSE(unit)
        STOP 1
    END IF
    CLOSE(unit)

    ! Read indices array
    OPEN(NEWUNIT=unit, FILE=indices_file, FORM='UNFORMATTED', ACCESS='STREAM', STATUS='OLD', IOSTAT=io_stat)
     IF (io_stat /= 0) THEN
        PRINT *, "Error opening indices file: ", indices_file, " IOSTAT=", io_stat
        STOP 1
    END IF
    READ(unit, IOSTAT=io_stat) indices
    IF (io_stat /= 0) THEN
        PRINT *, "Error reading from indices file: ", indices_file, " IOSTAT=", io_stat
        CLOSE(unit)
        STOP 1
    END IF
    CLOSE(unit)
    indices = indices + 1  ! Convert to 1-based

    ! Read indptr array
    OPEN(NEWUNIT=unit, FILE=indptr_file, FORM='UNFORMATTED', ACCESS='STREAM', STATUS='OLD', IOSTAT=io_stat)
     IF (io_stat /= 0) THEN
        PRINT *, "Error opening indptr file: ", indptr_file, " IOSTAT=", io_stat
        STOP 1
    END IF
    READ(unit, IOSTAT=io_stat) indptr
     IF (io_stat /= 0) THEN
        PRINT *, "Error reading from indptr file: ", indptr_file, " IOSTAT=", io_stat
        CLOSE(unit)
        STOP 1
    END IF
    CLOSE(unit)
    indptr = indptr + 1 ! Convert to 1-based

  END SUBROUTINE read_csr_matrix

  ! --- arpack_diagonalize subroutine remains the same ---
      SUBROUTINE arpack_diagonalize(N, nnz, data, indices, indptr, which, nev, ncv_in, tol, maxit, evals, evecs, info)
    ! --- Arguments ---
    INTEGER, INTENT(IN) :: N, nnz, nev, ncv_in, maxit  ! Renamed ncv to ncv_in
    REAL(dp), INTENT(IN) :: data(nnz), tol
    INTEGER, INTENT(IN) :: indices(nnz), indptr(N+1)
    CHARACTER(*), INTENT(IN) :: which
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: evals(:), evecs(:,:)
    INTEGER, INTENT(OUT) :: info

    ! --- Local Variables ---
    INTEGER :: ido, iparam(11), ipntr(11), lworkl
    REAL(dp), ALLOCATABLE :: resid(:), v(:,:), workd(:), workl(:)
    LOGICAL, ALLOCATABLE :: select(:)
    CHARACTER(1) :: bmat = 'I'
    CHARACTER(2) :: whch
    REAL(dp) :: sigma_dummy = 0.0_dp
    REAL(dp) :: tol_internal
    INTEGER :: ncv ! Adjusted ncv to be used internally

    ! --- Executable Statements ---

    ! 1. Check basic nev requirement
    IF (nev >= N) THEN
        PRINT *, "Error: nev (", nev, ") must be >= 1 and < N (", N, ")"
        info = -999 ! Assign a custom error code
        RETURN
    END IF
    IF (nev < 1) THEN
       PRINT *, "Error: nev must be at least 1."
       info = -997
       RETURN
    END IF

    ! 2. Determine appropriate internal ncv
    !    Rule: nev < ncv <= N. Try to use ncv_in, but adjust if needed.
    !    A common heuristic is ncv >= 2*nev + 1 for stability.
    ncv = MAX(nev + 1, 2 * nev + 1) ! Minimum sensible ncv
    ncv = MAX(ncv, ncv_in)          ! Try to respect user's request if larger
    ncv = MIN(ncv, N)               ! Ensure ncv <= N

    ! Re-check if adjustment made it invalid (can happen if N is very small, e.g., N=2, nev=1)
    IF (ncv <= nev) THEN
        PRINT *, "Error: Cannot satisfy nev < ncv <= N. N=", N, " nev=", nev, " Adjusted ncv=", ncv
        info = -998 ! Another custom error code
        RETURN
    END IF

    ! Print a warning if the requested ncv was adjusted
    IF (ncv /= ncv_in) THEN
         PRINT *, "Warning: Input ncv (", ncv_in, ") adjusted to ", ncv, " to satisfy nev < ncv <= N."
    END IF

    ! 3. Proceed with adjusted ncv
    whch = which(1:2)
    IF (tol <= 0.0_dp) THEN
      tol_internal = 0.0_dp
    ELSE
      tol_internal = tol
    END IF

    lworkl = ncv * (ncv + 8)
    ALLOCATE(resid(N), v(N, ncv), workd(3*N), workl(lworkl), select(ncv), STAT=info) ! Use adjusted ncv
    IF (info /= 0) THEN
        PRINT *, "Error allocating ARPACK workspace arrays. STAT=", info
        RETURN
    END IF

    ido = 0
    iparam = 0
    iparam(1) = 1
    iparam(3) = maxit
    iparam(7) = 1

    DO WHILE (ido /= 99)
      CALL dsaupd(ido, bmat, N, whch, nev, tol_internal, resid, ncv, v, N, &  ! Use adjusted ncv
                  iparam, ipntr, workd, workl, lworkl, info)

      IF (ido == -1 .OR. ido == 1) THEN
        CALL spmv_csr(N, data, indices, indptr, workd(ipntr(1)), workd(ipntr(2)))
      ELSE IF (ido /= 99) THEN
         PRINT *, "Unexpected IDO value from dsaupd:", ido
         GOTO 100
      END IF
    END DO

100 CONTINUE

    IF (info == 0 .OR. info == 1) THEN
      ALLOCATE(evals(nev), evecs(N, nev), STAT=iparam(1))
      IF (iparam(1) /= 0) THEN
          PRINT *, "Error allocating eigenvalue/vector arrays. STAT=", iparam(1)
          info = -99 ! Assign a specific error code for allocation failure
          GOTO 200
      END IF

      CALL dseupd(.TRUE., 'A', select, evals, evecs, N, sigma_dummy, bmat, N, whch, &
                  nev, tol_internal, resid, ncv, v, N, iparam, ipntr, workd, workl, lworkl, info) ! Use adjusted ncv
      IF (info /= 0) THEN
         PRINT *, "Error during dseupd (post-processing). Info code:", info
      END IF
    ELSE
       PRINT *, "dsaupd did not complete successfully. Info code:", info
       IF (ALLOCATED(evals)) DEALLOCATE(evals)
       IF (ALLOCATED(evecs)) DEALLOCATE(evecs)
       ALLOCATE(evals(0), evecs(N,0))
    END IF

200 CONTINUE

    DEALLOCATE(resid, v, workd, workl, select)

  END SUBROUTINE arpack_diagonalize

  ! --- spmv_csr subroutine remains the same ---
  SUBROUTINE spmv_csr(N, data, indices, indptr, x, y)
    INTEGER, INTENT(IN) :: N
    REAL(dp), INTENT(IN) :: data(:), x(N)
    INTEGER, INTENT(IN) :: indices(:), indptr(N+1) ! indptr is size N+1
    REAL(dp), INTENT(OUT) :: y(N)
    INTEGER :: i, k, start_ptr, end_ptr, col
    REAL(dp) :: yi
    ! Optimized slightly to reduce memory access inside inner loop
    DO i = 1, N
      yi = 0.0D0
      start_ptr = indptr(i)
      end_ptr = indptr(i+1) - 1
      DO k = start_ptr, end_ptr
        col = indices(k)
        yi = yi + data(k) * x(col)
      END DO
      y(i) = yi
    END DO
  END SUBROUTINE spmv_csr

END PROGRAM spin1_heisenberg_arpack