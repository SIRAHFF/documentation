This tutorial shows how to apply the hybrid solvation approach of SIRAH force field to speed up the simulation of an atomistic solute. The example system contains a DNA fiber surrounded by a shell of atomistic waters, which are embedded in coarse-grained (CG) molecules, called WT4, to represent
bulk water. The general procedure is extensible to any other solute. The hybrid solvation methodology is well tested for SPC, SPC/E, TIP3P atomistic water models. The main references for this tutorial are: `WAT4 <https://pubs.acs.org/doi/abs/10.1021/ct100379f>`_, `All-atoms/CG solvation <https://pubs.acs.org/doi/abs/10.1021/ct3001816>`_, `Transferable All-atoms/CG solvation <https://pubs.acs.org/doi/abs/10.1021/jp4079579>`_, and `SIRAH Tools <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_. We strongly advise you to read these articles before starting the tutorial.

.. important::

    Check :ref:`download <download gromacs>` section for download and set up details before to start this tutorial.
    Since this is **tutorial 2**, remember to replace ``X.X``, the files corresponding to this tutorial can be found in: ``sirah.ff/tutorial/2/``


2.1. PDB to GROMACS format
__________________________

Use pdb2gmx to convert the PDB file of the solute ``dna.pdb`` into GROMACS format: 

.. code-block:: bash

	gmx pdb2gmx -f sirah.ff/tutorial/2/dna.pdb -o dna.gro

.. caution::
	
	In GROMACS versions prior to 5.x, the "gmx" command should not be used.
	
When prompted, choose *AMBER99SB* and then *TIP3P* as water model.

For GROMACS to recognize SIRAH, edit your topology file ``topol.top`` adding the following lines **after** the force field definition:  

.. list-table::
   :align: center
   :widths: 50 50
   :header-rows: 1

   * - Topology before editing
     - Topology after editing
   * - | ; Include forcefield parameters 
       | #include “amber99sb.ff/forcefield.itp” 
       |     
       |     
              
     - | ; Include forcefield parameters 
       | #include “amber99sb.ff/forcefield.itp” 
       | #include “./sirah.ff/hybsol_comb2.itp” 
       | #include “./sirah.ff/solv.itp” 


.. important::

	``hybsol_comb2.itp`` is a parameter file, while ``solv.itp`` links to the topologies of WT4 and SIRAH ions. The choice of the right parameter file depends on the chosen atomistic force field (see :ref:`Table 1 <tableg>`). The plug-in works smoothly with the implemented force fields of the GROMACS distribution. When using a customized force field (e.g. for lipids) choose the parameter file according to the combination rule and check that the atom type you are using for the atomistic water match the DEFAULT in the SIRAH parameter file, which is OW. If they don't match then rename OW accordingly to your definition.

.. _tableg:

Table 1. SIRAH parameter file to include in the topology according to the chosen atomistic force field.

+------------------+------------------+-----------------------+--------+-------+--------+------+
|                                     | .. centered::   Atomistic force field                  |
+==================+==================+=======================+========+=======+========+======+
|  Parameter file  | Combination rule |          GMX          | GROMOS | AMBER | CHARMM | OPLS |
+------------------+------------------+-----------------------+--------+-------+--------+------+
| hybsol_comb1.itp | .. centered::  1 |           X           |   X    |       |        |      |
+------------------+------------------+-----------------------+--------+-------+--------+------+
| hybsol_comb2.itp | .. centered::  2 |                       |        |   X   |   X    |      |
+------------------+------------------+-----------------------+--------+-------+--------+------+
| hybsol_comb3.itp | .. centered::  3 |                       |        |       |        |  X   |
+------------------+------------------+-----------------------+--------+-------+--------+------+



2.2. Solvate the system
_______________________


Define the simulation box of the system:

.. code-block:: bash 
	
	gmx editconf -f dna.gro -o dna_box.gro -bt octahedron -d 2 -c

Then add an atomistic water shell:

.. code-block:: bash 

	gmx solvate -cp dna_box.gro -cs spc216.gro -o dna_shell.gro -shell 1

.. note:: 

	Before GROMACS version 5.x, the command *gmx solvate* was called *genbox*.

Finally add the CG solvent:

.. code-block:: bash 

	gmx solvate -cp dna_shell.gro -cs ./sirah.ff/wt4tip3p.gro -o dna_sol.gro


Edit the [ molecules ] section in ``topol.top`` to include the number of SOL and WT4 molecules:
	
.. hint::
	
	If you forget to look at the result of *solvate* to see how many SOL and WT4 molecules were added, you can use the following command line to get these numbers: 

	.. code-block:: console

		grep -c OW dna_sol.gro; grep -c WP1 dna_sol.gro

.. list-table::
   :align: center
   :widths: 50 50
   :header-rows: 1

   * - Topology before editing
     - Topology after editing
   * - | [ molecules ] 
       | ; Compound #mols 
       | DNA_chain_A 1    
       | DNA_chain_B 1    
              
     - | [ molecules ] 
       | ; Compound #mols 
       | DNA_chain_A 1 
       | DNA_chain_B 1  
       | SOL 3580 
       | WT4 2659  



.. caution::
	
	The number of added SOL and WT4 molecules (3580, 2659) may change according to the software version.

Add CG counterions:

.. code-block:: bash

	gmx grompp -f sirah.ff/tutorial/2/GPU/em_HYBSOL.mdp -p topol.top -c dna_sol.gro -o dna_sol.tpr

.. code-block:: bash

	gmx genion -s dna_sol.tpr -o dna_ion.gro -np 38 -pname NaW


When prompted, choose to substitute WT4 molecules by NaW ions.

.. note:: 

	The available ionic species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW). One
	ion pair (e.g. NaW-ClW) each 34 WT4 molecules renders a salt concentration of ~0.15M (see :ref:`Appendix <Appendix>` for details). 
	Be aware that SIRAH ions remain within the CG phase. So, if the presence of atomistic electrolytes in close contact with the solute is important to describe the physics of the system you will have to add them.

Edit the [ molecules ] section in ``topol.top`` to include the 38 NaW ions and the correct number of WT4 molecules.

Before running the simulation it may be a good idea to visualize your molecular system. CG molecules
are not recognized by molecular visualizers and will not display correctly. To fix this problem you may
generate a PSF file of the system using the script *g_top2psf.pl*:

.. code-block:: bash

	./sirah.ff/tools/g_top2psf.pl -i topol.top -o dna_ion.psf

.. note::

	This is the basic usage of the script ``g_top2psf.pl``, you can learn other capabilities from its help:
	``./sirah.ff/tools/g_top2psf.pl -h``


Use VMD to check how the hybrid system looks like:

.. code-block::

	vmd dna_ion.psf dna_ion.gro -e sirah.ff/tools/sirah_vmdtk.tcl

.. note::

	VMD assigns default radius to unknown atom types, the script sirah_vmdtk.tcl sets the right
	ones. It also provides a kit of useful selection macros, coloring methods and a backmapping utility.
	Use the command sirah_help in the Tcl/Tk console of VMD to access the manual pages.

2.3. Run the simulation
________________________

.. important:: 

	By default in this tutorial we will use input files for GROMACS on GPU (``sirah.ff/tutorial/2/GPU``). Example input files for using GROMACS on CPU can be found at: ``sirah.ff/tutorial/2/CPU``.

The folder ``sirah.ff/tutorial/2/GPU/`` contains typical input files for energy minimization
(``em_HYBSOL.mdp``), equilibration (``eq_HYBSOL.mdp``) and production (``md_HYBSOL.mdp``) runs. Please
check carefully the input flags therein.

Create an index files adding a group for WT4 and NaW:

.. code-block:: bash

	echo -e "r WT4 | r NaW\nq\n" | make_ndx -f dna_ion.gro -o dna_ion.ndx

.. note::

	WT4 and CG ions (NaW and ClW) are automatically set to the group “SIRAH-Solvent”.

Make a new folder for the run:

.. code-block:: bash

	mkdir run; cd run

Energy Minimization:

.. code-block:: bash

	gmx grompp -f ../sirah.ff/tutorial/2/GPU/em_HYBSOL.mdp -p ../topol.top -po em.mdp -n ../dna_ion.ndx -c ../dna_ion.gro -o dna_em.tpr mdrun -deffnm dna_em &> EM.log &

Equilibration:

.. code-block:: bash 

	grompp -f ../sirah.ff/tutorial/2/GPU/eq_HYBSOL.mdp -p ../topol.top -po eq.mdp -n ../dna_ion.ndx -c dna_em.gro -o dna_eq.tpr mdrun -deffnm dna_eq &> EQ.log &

Production (100ns):

.. code-block:: bash

	gmx grompp -f ../sirah.ff/tutorial/2/GPU/md_HYBSOL.mdp -p ../topol.top -po md.mdp -n ../dna_ion.ndx -c dna_eq.gro -o dna_md.tpr mdrun -deffnm dna_md &> MD.log &

.. note::

	GPU flags were set for GROMACS 4.6.7, different versions may complain about some specifications.


2.4. Visualizing the simulation
________________________________


That’s it! Now you can analyze the trajectory.

Process the output trajectory at folder ``run/`` to account for the Periodic Boundary Conditions (PBC):

.. code-block:: bash

	gmx trjconv -s dna_em.tpr -f dna_md.xtc -o dna_md_pbc.xtc -n ../dna_ion.ndx -ur compact -center -pbc mol

When prompted, choose *DNA* for centering and *System* for output.

Now you can check the simulation using VMD:

.. code-block:: bash

	vmd ../dna_ion.psf ../dna_ion.gro dna_md_pbc.xtc -e ../sirah.ff/tools/sirah_vmdtk.tcl


