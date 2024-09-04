Recompiling Namd
________________________________

If careful comparisons of the energies in NAMD and Amber are desired, the NAMD source code must be modified to use the same conversion factors, particularly to convert the electrostatic interactions to kcal/mol. These constants depend on the program:

AMBER = 332.0522173
CHARMM = 332.054
NAMD = 332.0636

To make the change, edit the line of the src/common.h file in the NAMD source code that contain

#define COULOMB 332.0636

to

#define COULOMB 332.0522173

Once this modification is made, recompile it and the program will be ready to run.