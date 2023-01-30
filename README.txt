<About>

This is the program to Bott index, a topological number of insulators or superconductors in real space.

The Hamiltonian is of 2-dimentional superconductor when adding magnetic field h_z, and Rashba spin-orbit interaction.[1]



<Hop file>

You have to prepare the hopping indices and the unit vector in the direction.

The included hop.txt is the example. When the electron hopps from site j to site i, and the direction are x_ij, and y_ij, you write the file as follows:

i    j    x_ij    y_ij


The example file is for one of Penrose approximants.



<Pair potential and particle number file>

You have to prepare the pair potentials and the particle numbers.

The included pair_potential.txt and particle_number.txt are the example. When the site i has a pair potential Delta_i, you write the file as follows:

Delta_i


When the partilce number at the site i and spin up is n_iup, and at the site i and spin down is n_idown, you write the file as follows:

n_iup    n_idown



The example file is for one of Penrose approximants.

<Referrence>

[1] R. Ghadimi, T. Sugimoto, K. Tanaka, and T. Tohyama. (2021). Phys. Rev. B. 104, 144511.