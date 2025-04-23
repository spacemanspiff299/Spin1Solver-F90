# Spin-1 Heisenberg Chain Simulation

This package provides tools for simulating spin-1 Heisenberg chains using exact diagonalization and the Lanczos algorithm via ARPACK.

## Project Structure

- `main.f90` - Main program for exact diagonalization of the spin-1 Heisenberg model
- `spin_matrices.f90` - Module containing spin matrices definitions
- `basis_generator.f90` - Module to generate the basis states
- `hamiltonian_builder.f90` - Module to build the Hamiltonian matrix
- `diagonalize.f90` - Module for diagonalization routines
- `observables.f90` - Module for calculating physical observables
- `sparse_matrix.py` - Python script to generate sparse matrix representation for Lanczos method
- `spin1_heisenberg_arpack.f90` - Fortran program using ARPACK for Lanczos algorithm

## Getting Started

### Prerequisites

- gfortran compiler
- OpenBLAS library
- Python 3 with NumPy and SciPy (for sparse matrix generation)

## Compilation and Usage

### Method 1: Exact Diagonalization

1. Compile the program:
   ```
   make
   ```

2. Run the program:
   ```
   ./spin_program
   ```

3. To clean up object files and executables:
   ```
   make clean
   ```

### Using Multithreading with OpenBLAS

You can control the number of threads used by OpenBLAS by setting the `OPENBLAS_NUM_THREADS` environment variable:

```bash
export OPENBLAS_NUM_THREADS=4  # Use 4 threads for OpenBLAS
./spin_program
```

Or directly on the command line:

```bash
OPENBLAS_NUM_THREADS=4 ./spin_program
```

### Method 2: Lanczos Algorithm (ARPACK)

For larger systems, the Lanczos algorithm is more efficient as it finds only a few eigenvalues/eigenvectors:

1. First, generate the sparse matrix representation using Python:
   ```bash
   python sparse_matrix.py N
   ```
   Where `N` is the number of sites in your spin chain. This will create the binary files needed for the ARPACK solver.

2. Compile the ARPACK solver:
   ```bash
   gfortran -O2 spin1_heisenberg_arpack.f90 -o spin1_arpack -lopenblas
   ```

3. Run the ARPACK solver:
   ```bash
   ./spin1_arpack
   ```

4. For multithreading with ARPACK:
   ```bash
   OPENBLAS_NUM_THREADS=4 ./spin1_arpack
   ```

## Output Files

- Exact diagonalization: `eigenspectrum_N{N}_J{J}_PBC{periodic}.dat`
- ARPACK method: Results are printed to standard output

## System Parameters

You can modify the system parameters in `main.f90` or by editing the Python script for the ARPACK method:
- `N`: Number of sites
- `J`: Coupling constant (J > 0 is antiferromagnetic)
- `periodic`: Whether to use periodic boundary conditions 