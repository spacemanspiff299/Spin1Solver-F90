FC = gfortran
FFLAGS = -O2
LDFLAGS = -lopenblas

# Object files
OBJS = spin_matrices.o basis_generator.o hamiltonian_builder.o diagonalize.o main.o

# Default target
all: spin_program

# Link the executable
spin_program: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDFLAGS)

# Compile the source files
spin_matrices.o: spin_matrices.f90
	$(FC) $(FFLAGS) -c $<

basis_generator.o: basis_generator.f90 spin_matrices.o
	$(FC) $(FFLAGS) -c $<

hamiltonian_builder.o: hamiltonian_builder.f90 spin_matrices.o basis_generator.o
	$(FC) $(FFLAGS) -c $<

diagonalize.o: diagonalize.f90
	$(FC) $(FFLAGS) -c $<

main.o: main.f90 spin_matrices.o basis_generator.o hamiltonian_builder.o diagonalize.o
	$(FC) $(FFLAGS) -c $<

# Clean up
clean:
	rm -f *.o *.mod *.dat spin_program

# Run the program
run: spin_program
	./spin_program 