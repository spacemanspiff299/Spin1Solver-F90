module hamiltonian_builder

  use iso_fortran_env, only: real64
  use spin_matrices, only: Sx, Sy, Sz
  use basis_generator, only: get_basis_size, generate_basis
  implicit none

  integer, parameter :: dp = real64

  private

  ! Public interfaces
  public :: build_heisenberg_hamiltonian
  public :: apply_operator_multi 
  public :: calc_bond_element

  ! --- Constants ---
  real(dp), parameter :: TOLERANCE = 1.0e-10_dp

contains

  !--------------------------------------------------------------------
  ! Build the full Heisenberg Hamiltonian matrix for N spin-1 sites
  ! H = J ∑_i (S_i · S_{i+1})
  ! Exploits Hermiticity and conservation of total Sz
  !--------------------------------------------------------------------
  subroutine build_heisenberg_hamiltonian(N, J, periodic, H)
    integer, intent(in) :: N
    real(dp), intent(in) :: J
    logical, intent(in) :: periodic
    complex(dp), allocatable, intent(out) :: H(:,:)

    integer :: basis_size, i, bra_idx, ket_idx
    integer, allocatable :: basis(:,:)
    complex(dp) :: matrix_element
    integer, allocatable :: total_sz(:)

    basis_size = get_basis_size(N)
    allocate(H(basis_size, basis_size))
    H = (0.0_dp, 0.0_dp)

    ! Generate the basis
    allocate(basis(basis_size, N))
    call generate_basis(N, basis)

    ! Calculate total Sz for each basis state for block diagonal optimization
    allocate(total_sz(basis_size))
    do i = 1, basis_size
      total_sz(i) = sum(basis(i,:))
    end do

    ! Only compute upper triangular part of the Hamiltonian (bra_idx <= ket_idx)
    ! This exploits Hermiticity
    do ket_idx = 1, basis_size
      do bra_idx = 1, ket_idx  ! Compute upper triangle including diagonal

        ! Skip calculation if states are in different Sz sectors
        if (total_sz(bra_idx) /= total_sz(ket_idx)) then
          cycle ! Element must be zero
        end if

        ! Initialize the matrix element <bra| H |ket>
        matrix_element = (0.0_dp, 0.0_dp)

        ! Loop over bonds (i, i+1)
        do i = 1, N - 1
          matrix_element = matrix_element + calc_bond_element(basis(bra_idx,:), basis(ket_idx,:), i, i+1)
        end do

        ! Add periodic boundary term if needed (N, 1)
        if (periodic .and. N > 1) then
          matrix_element = matrix_element + calc_bond_element(basis(bra_idx,:), basis(ket_idx,:), N, 1)
        end if

        ! Scale by coupling constant
        matrix_element = J * matrix_element

        ! Set the matrix element in the upper triangle
        H(bra_idx, ket_idx) = matrix_element

        ! Set the conjugate element in the lower triangle (Hermiticity)
        if (bra_idx /= ket_idx) then
          H(ket_idx, bra_idx) = conjg(matrix_element)
        end if
      end do
    end do

    deallocate(total_sz)
    deallocate(basis)

  end subroutine build_heisenberg_hamiltonian

  !--------------------------------------------------------------------
  ! Calculate matrix element for S_i · S_j between two basis states
  ! element = <bra| S_i · S_j |ket>
  !         = <bra| S_i^x S_j^x |ket> + <bra| S_i^y S_j^y |ket> + <bra| S_i^z S_j^z |ket>
  !--------------------------------------------------------------------
  function calc_bond_element(bra_state, ket_state, i, j) result(element)
    integer, intent(in) :: bra_state(:) ! Vector of size N with spins (-1, 0, 1)
    integer, intent(in) :: ket_state(:) ! Vector of size N
    integer, intent(in) :: i, j         ! Site indices (1-based)
    complex(dp) :: element

    element = (0.0_dp, 0.0_dp)

    ! S_i^x S_j^x term
    element = element + calculate_two_op_term_element(Sx, Sx, i, j, bra_state, ket_state)

    ! S_i^y S_j^y term
    element = element + calculate_two_op_term_element(Sy, Sy, i, j, bra_state, ket_state)

    ! S_i^z S_j^z term
    element = element + calculate_two_op_term_element(Sz, Sz, i, j, bra_state, ket_state)

  end function calc_bond_element

  !--------------------------------------------------------------------
  ! Helper Function: Calculate <bra| op_i op_j |ket>
  ! Applies op_j to ket_state, then applies op_i to the results,
  ! and sums amplitudes where the final state matches bra_state.
  !--------------------------------------------------------------------
  function calculate_two_op_term_element(op_i, op_j, i, j, bra_state, ket_state) result(term_element)
    complex(dp), dimension(3,3), intent(in) :: op_i, op_j
    integer, intent(in) :: i, j
    integer, intent(in) :: bra_state(:), ket_state(:)
    complex(dp) :: term_element

    integer :: N, k, l
    integer :: num_intermediate_states, num_final_states
    integer, allocatable :: intermediate_states(:,:), final_states(:,:)
    complex(dp), allocatable :: intermediate_amplitudes(:), final_amplitudes(:)

    N = size(ket_state)
    term_element = (0.0_dp, 0.0_dp)

    ! Maximum 3 intermediate states from first operator
    allocate(intermediate_states(3, N), intermediate_amplitudes(3))
    ! Maximum 3*3 = 9 final states from second operator (before checking overlap)
    allocate(final_states(3, N), final_amplitudes(3)) ! Reuse buffer inside loop

    ! Step 1: Apply op_j to |ket_state> at site j
    call apply_operator_multi(op_j, j, ket_state, intermediate_states, &
                              intermediate_amplitudes, num_intermediate_states)

    ! Step 2: For each intermediate state, apply op_i at site i and check overlap with <bra_state|
    do k = 1, num_intermediate_states
      if (abs(intermediate_amplitudes(k)) > TOLERANCE) then
        ! Apply op_i to intermediate_states(k,:) at site i
        call apply_operator_multi(op_i, i, intermediate_states(k,:), final_states, &
                                  final_amplitudes, num_final_states)

        ! Step 3: Check if any final state matches bra_state
        do l = 1, num_final_states
          if (abs(final_amplitudes(l)) > TOLERANCE) then
            if (all(final_states(l,:) == bra_state)) then
              ! Found a contribution: amplitude_k * amplitude_l
              term_element = term_element + intermediate_amplitudes(k) * final_amplitudes(l)
            end if
          end if
        end do ! loop l over final states
      end if
    end do ! loop k over intermediate states

    deallocate(intermediate_states, intermediate_amplitudes)
    deallocate(final_states, final_amplitudes)

  end function calculate_two_op_term_element

  !--------------------------------------------------------------------
  ! Apply a spin operator op to site 'site' in 'in_state'.
  ! Returns ALL possible resulting states in 'out_states' (up to 3 for S=1),
  ! their corresponding 'amplitudes', and the 'num_states' found.
  !--------------------------------------------------------------------
  subroutine apply_operator_multi(op, site, in_state, out_states, amplitudes, num_states)
    complex(dp), dimension(3,3), intent(in) :: op
    integer, intent(in) :: site              ! Site index (1 to N)
    integer, intent(in) :: in_state(:)       ! Input state vector
    integer, intent(out) :: out_states(:,:)  ! Array(3, N) to store output states
    complex(dp), intent(out) :: amplitudes(:)! Array(3) for amplitudes
    integer, intent(out) :: num_states       ! Number of non-zero transitions

    integer :: in_spin_idx, out_spin_idx ! Matrix indices (1, 2, 3)
    integer :: out_spin_val              ! Spin value (-1, 0, 1)
    complex(dp) :: amp

    num_states = 0
    ! Convert input spin value (-1, 0, 1) at target site to matrix index (1, 2, 3)
    in_spin_idx = in_state(site) + 2

    ! Iterate through all possible output spin *indices* (1, 2, 3)
    do out_spin_idx = 1, 3
      amp = op(out_spin_idx, in_spin_idx)

      ! If the transition amplitude is non-zero
      if (abs(amp) > TOLERANCE) then
        num_states = num_states + 1

        ! Copy input state to the next available output state slot
        out_states(num_states, :) = in_state

        ! Modify the spin at the target site in the output state
        out_spin_val = out_spin_idx - 2 ! Convert index back to spin value
        out_states(num_states, site) = out_spin_val

        ! Store the corresponding amplitude
        amplitudes(num_states) = amp
      end if
    end do

    ! Zero out remaining amplitude slots if fewer than 3 states were generated
    do out_spin_idx = num_states + 1, 3
        amplitudes(out_spin_idx) = (0.0_dp, 0.0_dp)
    end do


  end subroutine apply_operator_multi

end module hamiltonian_builder