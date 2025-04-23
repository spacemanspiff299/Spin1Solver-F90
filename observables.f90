module observables
  use iso_fortran_env, only: real64
  use spin_matrices, only: Sx, Sy, Sz  ! Assuming this provides S^2, Sx, Sy, Sz operators if needed elsewhere
  use basis_generator ! Assuming this generates the basis states
  use hamiltonian_builder ! Assuming this provides related functionality (e.g., calc_bond_element)
  implicit none
  private
  public :: calculate_magnetization, calculate_correlation_function, &
           calculate_thermal_average, calculate_specific_heat, &
           calculate_susceptibility

  ! Define precision kind parameter
  integer, parameter :: dp = real64
  ! Define Boltzmann constant (often set to 1 in theoretical calculations)
  real(dp), parameter :: kB = 1.0_dp

contains

  !> Calculate expectation value of total magnetization per site for a given state (eigenvector).
  !> <S_z> = (1/N) * <psi| sum_j S_z^j |psi>
  !> Assumes basis(i,j) stores the eigenvalue of S_z^j (+/- 0.5) for the i-th basis state.
  function calculate_magnetization(eigenvector, basis, N) result(magnetization)
    complex(dp), intent(in) :: eigenvector(:) ! State vector (|psi>) in the computational basis
    integer, intent(in) :: basis(:,:)       ! Computational basis states [basis_size, N], storing Sz values per site
    integer, intent(in) :: N                 ! Number of sites
    real(dp) :: magnetization                ! Result: <S_z> per site

    integer :: i, j, basis_size
    real(dp) :: total_sz_expectation

    basis_size = size(eigenvector) ! Should match size(basis, 1)
    if (basis_size /= size(basis, 1)) then
        print *, "Error (calculate_magnetization): eigenvector size does not match basis size."
        stop 1
    end if
    if (N /= size(basis, 2)) then
        print *, "Error (calculate_magnetization): N does not match basis dimension."
        stop 2
    end if

    total_sz_expectation = 0.0_dp

    ! Sum over basis states |phi_i>: <psi| M |psi> = sum_i <psi|phi_i> <phi_i| M |phi_i> <phi_i|psi>
    ! Since M = sum_j S_z^j is diagonal in this basis: <phi_i| M |phi_i> = sum_j <phi_i| S_z^j |phi_i> = sum_j basis(i,j)
    do i = 1, basis_size
      ! Calculate total Sz for the i-th basis state
      ! NOTE: Assumes basis(i, j) directly contains the Sz eigenvalue (e.g., +0.5 or -0.5)
      total_sz_expectation = total_sz_expectation + &
           real(conjg(eigenvector(i)) * eigenvector(i), kind=dp) * sum(real(basis(i,:), kind=dp))
           ! abs(eigenvector(i))**2 would also work if complex arithmetic efficiency is a concern
    end do

    ! Divide by number of sites to get per-site magnetization
    if (N > 0) then
      magnetization = total_sz_expectation / real(N, kind=dp)
    else
      magnetization = 0.0_dp
    end if
  end function calculate_magnetization

  !> Calculate spin-spin correlation function <psi| S_i · S_j |psi> for a given state.
  !> Note: Correctness depends on the external function 'calc_bond_element'.
  function calculate_correlation_function(eigenvector, basis, i, j) result(correlation)
    complex(dp), intent(in) :: eigenvector(:) ! State vector (|psi>) in the computational basis
    integer, intent(in) :: basis(:,:)       ! Computational basis states [basis_size, N]
    integer, intent(in) :: i, j             ! Site indices (1-based)
    real(dp) :: correlation                  ! Result: <S_i . S_j>

    integer :: bra_idx, ket_idx, basis_size, N
    complex(dp) :: matrix_element
    integer :: bra_state(size(basis,2)), ket_state(size(basis,2))
    complex(dp) :: total_correlation_complex

    N = size(basis, 2)
    basis_size = size(eigenvector)
    if (basis_size /= size(basis, 1)) then
        print *, "Error (calculate_correlation_function): eigenvector size does not match basis size."
        stop 1
    end if
    if (i < 1 .or. i > N .or. j < 1 .or. j > N) then
        print *, "Error (calculate_correlation_function): site indices out of bounds."
        stop 2
    end if

    total_correlation_complex = (0.0_dp, 0.0_dp)

    ! Loop over all basis states for |ket> and <bra|
    ! <psi|O|psi> = sum_{bra,ket} <psi|bra><bra|O|ket><ket|psi>
    do ket_idx = 1, basis_size
      ket_state = basis(ket_idx, :)
      do bra_idx = 1, basis_size
        bra_state = basis(bra_idx, :)

        ! Calculate <bra|S_i · S_j|ket> using an external function
        ! This function needs to be provided and correctly implement the S.S operator.
        ! It should return (0,0) if the operator doesn't connect bra_state and ket_state.
        matrix_element = calc_bond_element(bra_state, ket_state, i, j)

        ! Add contribution only if matrix element is non-zero
        if (abs(matrix_element) > epsilon(0.0_dp)) then ! Use epsilon for floating point comparison
            total_correlation_complex = total_correlation_complex + &
                 conjg(eigenvector(bra_idx)) * matrix_element * eigenvector(ket_idx)
        end if
      end do
    end do

    ! Expectation value of a Hermitian operator is real
    correlation = real(total_correlation_complex, kind=dp)
    ! Small imaginary parts might arise due to numerical precision, could add a check here if needed.
    ! if (abs(aimag(total_correlation_complex)) > some_threshold) print *, "Warning: Non-zero imaginary part in correlation function."

  end function calculate_correlation_function

  !> Calculate thermal average of an observable O at temperature T.
  !> <O> = Tr(O*exp(-H/kT))/Tr(exp(-H/kT))
  !> Assumes the input 'observable' matrix is already in the energy eigenbasis.
  !> i.e., observable(n,n) = <E_n| O |E_n>.
  function calculate_thermal_average(observable_diag, eigenvalues, temperature) result(avg)
    real(dp), intent(in) :: observable_diag(:) ! Diagonal elements of Observable matrix in the energy basis <E_n|O|E_n>
    real(dp), intent(in) :: eigenvalues(:)     ! Eigenvalues of the Hamiltonian (E_n)
    real(dp), intent(in) :: temperature        ! Temperature (in energy units, kB=1)
    real(dp) :: avg                            ! Result: <O>_T

    integer :: i, n
    real(dp) :: partition_function, sum_o_boltzmann, boltzmann_factor
    real(dp) :: min_eigenvalue ! For numerical stability

    n = size(eigenvalues)
    if (n /= size(observable_diag)) then
        print *, "Error (calculate_thermal_average): observable_diag size does not match eigenvalues size."
        stop 1
    end if
    if (temperature <= 0.0_dp) then
        print *, "Error (calculate_thermal_average): Temperature must be positive."
        ! Handle T=0 case? Return ground state expectation? For now, stop.
        stop 2
    end if

    partition_function = 0.0_dp
    sum_o_boltzmann = 0.0_dp

    ! Find minimum eigenvalue for numerical stability (optional but recommended)
    if (n > 0) then
      min_eigenvalue = minval(eigenvalues)
    else
      avg = 0.0_dp ! Or NaN? Or error?
      return
    end if

    ! Calculate partition function Z = sum_n exp(-(E_n - E_min)/kT)
    ! And sum_n <E_n|O|E_n> * exp(-(E_n - E_min)/kT)
    do i = 1, n
      ! Avoid overflow/underflow by calculating relative Boltzmann factors
      boltzmann_factor = exp(-(eigenvalues(i) - min_eigenvalue) / (kB * temperature))
      partition_function = partition_function + boltzmann_factor
      sum_o_boltzmann = sum_o_boltzmann + observable_diag(i) * boltzmann_factor
    end do

    ! Normalize by the partition function
    if (partition_function > epsilon(0.0_dp)) then
      avg = sum_o_boltzmann / partition_function
    else
      ! Handle case where partition function is zero (e.g., T=0 and degenerate ground state?)
      ! Or if all Boltzmann factors underflowed.
      print *, "Warning (calculate_thermal_average): Partition function is close to zero."
      ! Depending on context, could return ground state average if T is very low
      ! Or NaN or 0. For now, return 0.
      avg = 0.0_dp
      ! If T->0, avg should approach the ground state expectation value(s).
      ! A more robust implementation might check temperature and handle T=0 separately.
    end if

  end function calculate_thermal_average

  !> Calculate specific heat from energy eigenvalues.
  !> C = (1/N_sites * k_B T^2) * (<E^2> - <E>^2)
  !> The result is specific heat *per site*. Add N_sites if total C is needed.
  function calculate_specific_heat(eigenvalues, N, temperature) result(specific_heat_per_site)
    real(dp), intent(in) :: eigenvalues(:) ! Eigenvalues of the Hamiltonian (E_n)
    integer, intent(in) :: N              ! Number of sites (for per-site normalization)
    real(dp), intent(in) :: temperature    ! Temperature (in energy units, kB=1)
    real(dp) :: specific_heat_per_site     ! Result: C / N

    integer :: i, n_states
    real(dp) :: partition_function, energy_avg_num, energy_sqr_avg_num
    real(dp) :: boltzmann_factor, energy_avg, energy_sqr_avg
    real(dp) :: min_eigenvalue ! For numerical stability

    n_states = size(eigenvalues)
    if (temperature <= 0.0_dp) then
        print *, "Error (calculate_specific_heat): Temperature must be positive."
        stop 1
    end if
    if (N <= 0) then
        print *, "Error (calculate_specific_heat): Number of sites N must be positive."
        stop 2
    end if

    partition_function = 0.0_dp
    energy_avg_num = 0.0_dp
    energy_sqr_avg_num = 0.0_dp

    ! Find minimum eigenvalue for numerical stability
    if (n_states > 0) then
      min_eigenvalue = minval(eigenvalues)
    else
      specific_heat_per_site = 0.0_dp
      return
    end if

    ! Calculate partition function Z = sum_n exp(-(E_n - E_min)/kT)
    ! Calculate numerator for <E> = sum_n E_n * exp(-(E_n - E_min)/kT)
    ! Calculate numerator for <E^2> = sum_n E_n^2 * exp(-(E_n - E_min)/kT)
    do i = 1, n_states
      boltzmann_factor = exp(-(eigenvalues(i) - min_eigenvalue) / (kB * temperature))
      partition_function = partition_function + boltzmann_factor
      energy_avg_num = energy_avg_num + eigenvalues(i) * boltzmann_factor
      energy_sqr_avg_num = energy_sqr_avg_num + eigenvalues(i)**2 * boltzmann_factor
    end do

    ! Calculate thermal averages <E> and <E^2>
    if (partition_function > epsilon(0.0_dp)) then
      energy_avg = energy_avg_num / partition_function
      energy_sqr_avg = energy_sqr_avg_num / partition_function
    else
      print *, "Warning (calculate_specific_heat): Partition function is close to zero."
      specific_heat_per_site = 0.0_dp
      ! Again, a T->0 limit might need special handling. C->0 as T->0.
      return
    end if

    ! Calculate specific heat per site
    specific_heat_per_site = (energy_sqr_avg - energy_avg**2) / (kB * temperature**2 * real(N, kind=dp))

    ! Handle potential NaN if T=0 was allowed (T**2 in denominator)
    if (specific_heat_per_site /= specific_heat_per_site) specific_heat_per_site = 0.0_dp ! Check for NaN

  end function calculate_specific_heat

  !> Calculate magnetic susceptibility per site.
  !> χ = (1 / N_sites * k_B T) * (<M^2> - <M>^2), where M = sum_j S_z^j.
  !> Assumes basis(k,j) stores the eigenvalue of S_z^j (+/- 0.5) for the k-th basis state.
  function calculate_susceptibility(eigenvalues, eigenvectors, basis, N, temperature) result(susceptibility_per_site)
    real(dp), intent(in) :: eigenvalues(:)       ! Hamiltonian eigenvalues E_n
    complex(dp), intent(in) :: eigenvectors(:,:)  ! Eigenvectors V [basis_size, n_states]. V(:,n) is the n-th eigenvector
    integer, intent(in) :: basis(:,:)            ! Computational basis states [basis_size, N]
    integer, intent(in) :: N                      ! Number of sites
    real(dp), intent(in) :: temperature          ! Temperature (in energy units, kB=1)
    real(dp) :: susceptibility_per_site          ! Result: Chi / N

    integer :: i, k, site, n_states, basis_size
    real(dp) :: partition_function, mag_avg_num, mag_sqr_avg_num
    real(dp) :: boltzmann_factor
    real(dp) :: mag_n              ! <E_n| M |E_n>
    real(dp) :: mag_sqr_n          ! <E_n| M^2 |E_n>
    real(dp) :: mag_avg, mag_sqr_avg
    real(dp) :: Mk                 ! Total Sz of k-th basis state
    real(dp) :: min_eigenvalue     ! For numerical stability

    n_states = size(eigenvalues)
    basis_size = size(basis, 1)

    ! Input validation
    if (n_states /= size(eigenvectors, 2)) then
        print *, "Error (calculate_susceptibility): Eigenvalue count mismatch with eigenvector columns."
        stop 1
    end if
    if (basis_size /= size(eigenvectors, 1)) then
        print *, "Error (calculate_susceptibility): Basis size mismatch with eigenvector rows."
        stop 2
    end if
     if (N /= size(basis, 2)) then
        print *, "Error (calculate_susceptibility): N does not match basis dimension."
        stop 3
    end if
    if (temperature <= 0.0_dp) then
        print *, "Error (calculate_susceptibility): Temperature must be positive."
        stop 4
    end if
    if (N <= 0) then
        print *, "Error (calculate_susceptibility): Number of sites N must be positive."
        stop 5
    end if

    partition_function = 0.0_dp
    mag_avg_num = 0.0_dp
    mag_sqr_avg_num = 0.0_dp

    ! Find minimum eigenvalue for numerical stability
    if (n_states > 0) then
      min_eigenvalue = minval(eigenvalues)
    else
      susceptibility_per_site = 0.0_dp
      return
    end if

    ! Calculate numerators for <M> and <M^2> thermal averages
    do i = 1, n_states  ! Loop over energy eigenstates |E_i>
      mag_n = 0.0_dp      ! Calculate <E_i| M |E_i>
      mag_sqr_n = 0.0_dp  ! Calculate <E_i| M^2 |E_i>

      ! Calculate expectation values in the energy eigenstate |E_i>
      ! using the computational basis expansion: |E_i> = sum_k V(k,i) |phi_k>
      ! M = sum_j S_z^j is diagonal in |phi_k> basis: M |phi_k> = Mk |phi_k>
      ! Mk = sum_j basis(k,j) (Total Sz of basis state k)
      ! <E_i| M |E_i> = sum_k |V(k,i)|^2 Mk
      ! <E_i| M^2 |E_i> = sum_k |V(k,i)|^2 Mk^2
      do k = 1, basis_size ! Loop over computational basis states |phi_k>
        ! Calculate Mk = total Sz for basis state k
        ! NOTE: Assumes basis(k, site) directly contains the Sz eigenvalue (e.g., +0.5 or -0.5)
        Mk = sum(real(basis(k,:), kind=dp))

        mag_n = mag_n + real(conjg(eigenvectors(k,i)) * eigenvectors(k,i), kind=dp) * Mk
        mag_sqr_n = mag_sqr_n + real(conjg(eigenvectors(k,i)) * eigenvectors(k,i), kind=dp) * Mk**2
           ! abs(eigenvectors(k,i))**2 is equivalent
      end do

      ! Calculate Boltzmann factor relative to E_min
      boltzmann_factor = exp(-(eigenvalues(i) - min_eigenvalue) / (kB * temperature))

      ! Add contributions to partition function and thermal average numerators
      partition_function = partition_function + boltzmann_factor
      mag_avg_num = mag_avg_num + mag_n * boltzmann_factor
      mag_sqr_avg_num = mag_sqr_avg_num + mag_sqr_n * boltzmann_factor
    end do

    ! Calculate thermal averages <M> and <M^2>
    if (partition_function > epsilon(0.0_dp)) then
      mag_avg = mag_avg_num / partition_function
      mag_sqr_avg = mag_sqr_avg_num / partition_function
    else
      print *, "Warning (calculate_susceptibility): Partition function is close to zero."
      susceptibility_per_site = 0.0_dp
      ! Handle T->0 limit? Chi should go to 0 or a constant depending on degeneracy.
      return
    end if

    ! Calculate susceptibility per site
    susceptibility_per_site = (mag_sqr_avg - mag_avg**2) / (kB * temperature * real(N, kind=dp))

    ! Handle potential NaN if T=0 was allowed
    if (susceptibility_per_site /= susceptibility_per_site) susceptibility_per_site = 0.0_dp

  end function calculate_susceptibility

end module observables