module basis_generator
  implicit none
  private
  public :: generate_basis, get_basis_size

contains
  
  ! Calculate the size of the full Hilbert space for N spin-1 sites
  ! Each site has 3 possible states (-1, 0, 1), so total is 3^N
  function get_basis_size(N) result(basis_size)
    integer, intent(in) :: N
    integer :: basis_size
    
    basis_size = 3**N
  end function get_basis_size
  
  ! Generate the full basis for a system of N spin-1 sites
  ! Each row of the basis array represents one basis state
  ! Each element of a basis state is -1, 0, or 1 (corresponding to Sz eigenvalue)
  subroutine generate_basis(N, basis)
    integer, intent(in) :: N
    integer, allocatable, intent(out) :: basis(:,:)
    integer :: i, j, state, remainder
    integer :: basis_size
    
    basis_size = get_basis_size(N)
    allocate(basis(basis_size, N))
    
    do i = 1, basis_size
      state = i - 1  ! Start from 0 for easier conversion
      
      do j = N, 1, -1
        remainder = mod(state, 3)
        basis(i, j) = remainder - 1  ! Convert {0,1,2} to {-1,0,1}
        state = state / 3
      end do
    end do
  end subroutine generate_basis
  
end module basis_generator 