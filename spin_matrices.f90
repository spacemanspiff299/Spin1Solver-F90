module spin_matrices
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: dp = real64

  complex(dp), parameter :: iunit = (0.0_dp, 1.0_dp)

  ! Define the spin value (although not directly used in matrices below, keep for consistency)
  real(dp), parameter :: s = 1.0_dp

  ! Sx = 1/sqrt(2) * [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
  complex(dp), dimension(3,3), parameter :: &
  Sx = (1.0_dp / sqrt(2.0_dp)) * reshape([ 0.0_dp, 1.0_dp, 0.0_dp, & ! Col 1
                                           1.0_dp, 0.0_dp, 1.0_dp, & ! Col 2
                                           0.0_dp, 1.0_dp, 0.0_dp ],& ! Col 3
                                          [3,3], order=[2,1]), & ! Ensure Fortran column-major order

  ! Sy = (1/sqrt(2)) * [[0, -i, 0], [i, 0, -i], [0, i, 0]]
  Sy = (1.0_dp / sqrt(2.0_dp)) * reshape([ (0.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), (0.0_dp, 0.0_dp), & ! Col 1
                                           (0.0_dp,-1.0_dp), (0.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), & ! Col 2
                                           (0.0_dp, 0.0_dp), (0.0_dp,-1.0_dp), (0.0_dp, 0.0_dp) ],& ! Col 3
                                          [3,3], order=[2,1]), & ! Ensure Fortran column-major order

  ! Sz = [[1, 0, 0], [0, 0, 0], [0, 0, -1]]
  Sz = reshape([ 1.0_dp,  0.0_dp,  0.0_dp, & ! Col 1
                 0.0_dp,  0.0_dp,  0.0_dp, & ! Col 2
                 0.0_dp,  0.0_dp, -1.0_dp ],& ! Col 3
                [3,3], order=[2,1]), & ! Ensure Fortran column-major order

  ! S+ = [[0, sqrt(2), 0], [0, 0, sqrt(2)], [0, 0, 0]]
  Sp = reshape([ 0.0_dp, 0.0_dp       , 0.0_dp, & ! Col 1
                 sqrt(2.0_dp), 0.0_dp, 0.0_dp, & ! Col 2
                 0.0_dp, sqrt(2.0_dp), 0.0_dp ],& ! Col 3
                [3,3], order=[2,1]), & ! Ensure Fortran column-major order

  ! S- = (S+)^dagger = [[0, 0, 0], [sqrt(2), 0, 0], [0, sqrt(2), 0]]
  Sm = reshape([ 0.0_dp, sqrt(2.0_dp), 0.0_dp, & ! Col 1: S_- |1>, S_- |0>, S_- |-1>
                 0.0_dp, 0.0_dp      , sqrt(2.0_dp), & ! Col 2
                 0.0_dp, 0.0_dp      , 0.0_dp ],& ! Col 3
                [3,3], order=[2,1])
                
end module spin_matrices