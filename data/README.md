# Runtime data

Each thermal-feedback cycle creates an `eval_NNNN` directory here. It contains
the initial and compressed LAMMPS data files, compression and NEMD inputs and
logs, temperature profiles, heat-flux samples, and the final NEMD structure.

The generated evaluation directories are ignored by Git because trajectory and
data files can be large. This README keeps the directory and its purpose visible
in a fresh clone.
