## Ising-1D

/2014-10-01/

C code;

Solves 1D Ising model with N spins,
using Metropolis Monte Carlo method.

main.cpp -- Simple (slower) solver, spins stored as integers.

main-bitwise.cpp -- Faster solver, spins stored as bits.

Number of iterations can be set through nsteps.

Outputs 2 text files with magnetization and energy values,
at several temperatures and magnetic fields.
