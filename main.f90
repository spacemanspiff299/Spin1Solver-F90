program main
  ! Use modules, specifying ONLY the needed components
  use iso_fortran_env, only: real64
  use basis_generator, only: generate_basis, get_basis_size
  use hamiltonian_builder, only: build_heisenberg_hamiltonian
  use diagonalize, only: diagonalize_hamiltonian
  implicit none

  integer, parameter :: dp = real64

  ! System parameters
  integer :: N                  ! Number of sites
  real(dp) :: J                 ! Coupling constant
  logical :: periodic           ! Periodic boundary conditions
  integer :: basis_size         ! Size of the Hilbert space (integer is fine)

  ! Variables for calculations
  integer, allocatable :: basis(:,:)            ! Integer basis states
  complex(dp), allocatable :: H(:,:)            ! Hamiltonian
  real(dp), allocatable :: eigenvalues(:)       ! Eigenvalues
  complex(dp), allocatable :: eigenvectors(:,:) ! Eigenvectors
  integer :: info                               ! LAPACK info flag
  integer :: i, j_idx

  ! File handling
  integer :: outfile
  character(len=256) :: filename

  ! --- Set system parameters ---
  N = 4                 ! System size
  J = 1.0_dp            ! Antiferromagnetic coupling
  periodic = .true.     ! Use periodic boundary conditions
  if (N .eq. 2) then
    periodic = .false.
  end if

  ! Calculate basis size
  basis_size = get_basis_size(N)

  write(*,*) "Spin-1 Heisenberg chain calculator (Double Precision)"
  write(*,*) "Number of sites (N): ", N
  write(*,*) "Coupling constant (J): ", J
  write(*,*) "Periodic BC: ", periodic
  write(*,*) "Hilbert space dimension: ", basis_size

  ! --- Allocate arrays ---
  allocate(basis(basis_size, N))
  allocate(H(basis_size, basis_size))
  allocate(eigenvalues(basis_size))
  allocate(eigenvectors(basis_size, basis_size))

  ! --- Calculations ---
  write(*,*) "Generating basis..."
  call generate_basis(N, basis)

  write(*,*) "Building Hamiltonian..."
  call build_heisenberg_hamiltonian(N, J, periodic, H)

  write(*,*) "Diagonalizing Hamiltonian..."
  call diagonalize_hamiltonian(H, eigenvalues, eigenvectors, info)

  ! Check if diagonalization was successful
  if (info /= 0) then
    write(*,'(A,I0)') "Diagonalization failed with LAPACK info = ", info
    stop "Diagonalization Error"
  end if
  write(*,*) "Diagonalization successful."
  write(*,*) "Ground state energy: ", eigenvalues(1)

  ! Write eigenvalues and eigenvectors to file
  filename = "eigenspectrum_N" // trim(adjustl(write_int(N))) // "_J" // &
             trim(adjustl(write_real(J))) // "_PBC" // trim(adjustl(write_logical(periodic))) // ".dat"
  
  write(*,*) "Writing eigenvalues and eigenvectors to file: ", trim(filename)
  
  open(newunit=outfile, file=trim(filename), status='replace', action='write')
  
  ! Write header with system parameters
  write(outfile, '(A)') "# Heisenberg Model Eigenspectrum"
  write(outfile, '(A,I0)') "# Number of sites (N): ", N
  write(outfile, '(A,F10.4)') "# Coupling constant (J): ", J
  write(outfile, '(A,L1)') "# Periodic boundary conditions: ", periodic
  write(outfile, '(A,I0)') "# Hilbert space dimension: ", basis_size
  write(outfile, '(A)') "#"
  
  ! Write eigenvalues
  write(outfile, '(A)') "# Eigenvalues:"
  do i = 1, basis_size
    write(outfile, '(I6,A,ES20.12)') i, ": ", eigenvalues(i)
  end do
  
  ! Write eigenvectors
  write(outfile, '(A)') "#"
  write(outfile, '(A)') "# Eigenvectors (columns):"
  
  ! Write column headers for eigenvectors (one header per column)
  write(outfile, '(A7)', advance='no') "# State"
  do j_idx = 1, basis_size
    write(outfile, '(A12)', advance='no') "Re[Psi_" // trim(adjustl(write_int(j_idx))) // "]"
  end do
  write(outfile, '(A)') ""
  
  ! Write eigenvector values
  do i = 1, basis_size
    write(outfile, '(I7)', advance='no') i
    do j_idx = 1, basis_size
      write(outfile, '(ES12.4)', advance='no') real(eigenvectors(i,j_idx))
    end do
    write(outfile, '(A)') ""
  end do
  
  close(outfile)
  
  write(*,*) "Results written to file successfully."
  
  ! --- Clean up ---
  deallocate(basis, H, eigenvalues, eigenvectors)

contains

  ! Helper function to write integer to string
  character(len=20) function write_int(ival)
    implicit none
    integer, intent(in) :: ival
    write(write_int, '(I0)') ival
  end function write_int

  ! Helper function to write real to string
  character(len=20) function write_real(rval)
    implicit none
    real(dp), intent(in) :: rval
    write(write_real, '(F10.4)') rval 
  end function write_real

  ! Helper function to write logical to string
  character(len=5) function write_logical(lval)
    implicit none
    logical, intent(in) :: lval
    if (lval) then
      write_logical = "true"
    else
      write_logical = "false"
    end if
  end function write_logical

end program main