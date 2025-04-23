import numpy as np
from scipy.sparse import csr_matrix
import sys

def build_spin1_heisenberg(N, J=1.0, periodic=True):
    """
    Build the Hamiltonian for a spin-1 Heisenberg chain with N sites.
    
    Parameters:
    -----------
    N : int
        Number of sites
    J : float
        Coupling constant (J > 0 is antiferromagnetic)
    periodic : bool
        Whether to use periodic boundary conditions
        
    Returns:
    --------
    H : scipy.sparse.csr_matrix
        The Hamiltonian in CSR format
    """
    dim = 3**N
    
    # Construct the spin matrices for spin-1
    Sp = np.array([[0, np.sqrt(2), 0], [0, 0, np.sqrt(2)], [0, 0, 0]])  # S+
    Sm = np.array([[0, 0, 0], [np.sqrt(2), 0, 0], [0, np.sqrt(2), 0]])  # S-
    Sz = np.array([[1, 0, 0], [0, 0, 0], [0, 0, -1]])                   # Sz
    
    # Initialize data structure for sparse matrix
    row_indices = []
    col_indices = []
    values = []

    if N == 2:
        periodic = False
    
    # Function to get the state representation
    def state_to_spins(state_idx):
        spins = []
        temp = state_idx
        for _ in range(N):
            spins.append(temp % 3)  # 0 = +1, 1 = 0, 2 = -1
            temp //= 3
        return spins
    
    # Function to get state index from spins
    def spins_to_state(spins):
        state_idx = 0
        for i, s in enumerate(spins):
            state_idx += s * (3**i)
        return state_idx
    
    # Process all basis states
    for state_idx in range(dim):
        spins = state_to_spins(state_idx)
        
        # Loop over all bonds
        for i in range(N):
            j = (i + 1) % N if periodic else i + 1
            if j == N:  # Skip if non-periodic and at boundary
                continue
                
            si, sj = spins[i], spins[j]
            
            # Diagonal term: Sz_i * Sz_j
            # Map 0->+1, 1->0, 2->-1
            sz_i = 1 if si == 0 else (-1 if si == 2 else 0)
            sz_j = 1 if sj == 0 else (-1 if sj == 2 else 0)
            
            # Add diagonal contribution
            row_indices.append(state_idx)
            col_indices.append(state_idx)
            values.append(J * sz_i * sz_j)
            
            # Off-diagonal terms: 0.5 * (S+_i * S-_j + S-_i * S+_j)
            
            # S+_i * S-_j term
            if si < 2 and sj > 0:  # Can apply S+_i and S-_j
                new_spins = spins.copy()
                new_spins[i] = si + 1  # Apply S+
                new_spins[j] = sj - 1  # Apply S-
                new_state_idx = spins_to_state(new_spins)
                
                coeff = 0.5 * J * Sp[si, si+1] * Sm[sj, sj-1]
                
                row_indices.append(state_idx)
                col_indices.append(new_state_idx)
                values.append(coeff)
            
            # S-_i * S+_j term
            if si > 0 and sj < 2:  # Can apply S-_i and S+_j
                new_spins = spins.copy()
                new_spins[i] = si - 1  # Apply S-
                new_spins[j] = sj + 1  # Apply S+
                new_state_idx = spins_to_state(new_spins)
                
                coeff = 0.5 * J * Sm[si, si-1] * Sp[sj, sj+1]
                
                row_indices.append(state_idx)
                col_indices.append(new_state_idx)
                values.append(coeff)
    
    # Create CSR matrix
    H = csr_matrix((values, (row_indices, col_indices)), shape=(dim, dim), dtype=np.float64)
    
    # Verify Hermiticity
    if not np.allclose((H - H.T).data, 0):
        print("WARNING: Hamiltonian is not Hermitian!")
        
    return H

# Example usage
if __name__ == "__main__":
    if len(sys.argv) > 1:
        N = int(sys.argv[1])
    else:
        N = 4  # Default
    
    print(f"Building spin-1 Heisenberg chain with N = {N}")
    H = build_spin1_heisenberg(N, J=1.0, periodic=True)
    
    # Print matrix info
    print(f"Hamiltonian dimension: {H.shape[0]}x{H.shape[1]}")
    print(f"Number of non-zeros: {H.nnz}")
    print(f"Sparsity: {H.nnz/(H.shape[0]*H.shape[1])*100:.4f}%")
    
    # Check Hermiticity
    hermitian_diff = (H - H.T).data
    print(f"Hermiticity check: max|H-H†| = {np.max(np.abs(hermitian_diff)) if len(hermitian_diff) > 0 else 0}")
    
    # For small systems, we can use dense eigenvalue solver
    if N <= 5:
        from scipy.sparse.linalg import eigsh
        print("Computing ground state energy...")
        eigenvalues, _ = eigsh(H, k=1, which='SA')
        print(f"Ground state energy = {eigenvalues[0]}")
    
    # Save matrix to binary files
    print("Saving matrix to binary files...")
    H.data.astype(np.float64).tofile('spin1_data.bin')
    H.indices.astype(np.int32).tofile('spin1_indices.bin')
    H.indptr.astype(np.int32).tofile('spin1_indptr.bin')
    with open('spin1_meta.bin', 'wb') as f:
        np.array([H.shape[0], H.nnz], dtype=np.int32).tofile(f)
    
    print("Done.")