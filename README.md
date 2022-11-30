# CP-AFQMC
The method implements the imaginary-time evolution of a trial state towards the ground state. Here we simulate the Triangular Fermi-Hubbard model projecting the ground state from a Generalized-Hartree Fock wave function.

The CP-AFQMC method implements the imaginary-time evolution of a trial state towards the ground state [1]. Here we consider the Triangular Fermi-Hubbard model. We implement the spin-discrete Hubbard-Stratonovich decomposition to construct the imaginary-time evolution operator [2]. A constraint in the auxiliary-field paths is imposed to deal with the sign problem at the cost of a systematic error. The constraint can be released if the model has particle-hole symmetry. We project the ground state from a Generalized-Hartree Fock (GHF) wave function constructed by the self-consistent diagonalization of the GHF mean-field Hamiltonian [3].

## INSTRUCTIONS:
1. To compile change the makefile according to your machine and type make on the terminal;
2. To run you need two input files:
 * a) cpmc.in contains the relevant physical parameters and simulation instructions.
 * b) psi.in contains the initial wave function to be evolved in imaginary time.
3. The simulation produces the following output files:
 * a) cpmc.out contains the average total, kinetic and potential energies, the doublon density
and the chiral order parameter [4].
 * b) nn(ss).out contains charge (spin) correlations between all pairs of sites as a function of
the relative position.
 * c) t(s)pair.out contains the triplet (singlet) Cooper-pair correlations as a function of the
distance between lattice sites.
 * d) nk.out contains $\langle c^\dagger_{i \sigma} c_{j \sigma} \rangle$ displayed as in b). This is relevant to the computation of the
momentum distribution [5].
4. To create psi.in run GHF_wf.m script on MATLAB or OCTAVE after updating the relevant
Physical parameters (Lx, Ly, N_up, N_dwn and U).

## TO DO:
1. Implement twist boundary conditions.

## Acknowledgements
We acknowledge the didactic CPMC-Lab package which was used as a starting point for our
code [6].

## Giving Credit
If you use this code in your work, please cite the associated papers.
The Arxiv of "Chiral superconductivity in the doped triangular-lattice Fermi-Hubbard model in two dimensions" can be found on https://arxiv.org/abs/2210.13551.

This repository can be cited using:
```
@software{cp-afqmc,
  author = {Vinicius Zampronio},
  title = {{CP-AFQMC}},
  url = {https://github.com/quantechsimulations/CP-AFQMC},
  year = {2022},
}
```

## Bibliography
[1] Shiwei Zhang et al., Phys. Rev. B 55, 7464 (1997).

[2] J. E. Hirsch, Phys. Rev. B 31, 4403 (1985).

[3] Mingpu Qin et al., Phys. Rev. B 94, 085103 (2016).

[4] Aaron Szasz et al., Phys. Rev. X 10, 021042 (2020).

[5] Zheng Zhu et al., Phys. Rev. B 105, 205110 (2022).

[6] Huy Nguyen et al., Computer Physics Communications 185, 3344-3357 (2014).
