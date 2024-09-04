.. note::

   Please report bugs, errors or enhancement requests through `Issue Tracker <https://github.com/SIRAHFF/documentation/issues>`_ or if you have a question about SIRAH open a `New Discussion <https://github.com/SIRAHFF/documentation/discussions>`_.
   
This tutorial shows how to use the SIRAH force field to perform a coarse grained (CG) simulation of a
double stranded DNA in explicit solvent (called WatFour, WT4). The main references for
this tutorial are: `Dans et al. <https://pubs.acs.org/doi/abs/10.1021/ct900653p>`_ (latest parameters are those reported `here <https://pubs.acs.org/doi/abs/10.1021/ct100379f>`_), `Darré et al. <https://pubs.acs.org/doi/abs/10.1021/ct100379f>`_, and `Machado & Pantano <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_.
We strongly advise you to read these articles before starting the tutorial.

.. important::

    Check the :ref:`Setting up SIRAH <download gromacs>` section for download and set up details before starting this tutorial.
    Since this is **Tutorial 1**, remember to replace ``X.X``, the files corresponding to this tutorial can be found in: ``sirah.ff/tutorial/1/``


1.1. Build CG representations
______________________________

Map the atomistic structure of a 20-mer DNA to its CG representation:

.. code-block:: bash

	./sirah.ff/tools/CGCONV/cgconv.pl -i sirah.ff/tutorial/1/dna.pdb -o dna_cg.pdb

The input file ``-i`` dna.pdb contains all the heavy atoms composing the DNA molecule, while the  output ``-o`` dna_cg.pdb preserves a few of them.

Please check both PDB structures using VMD:

.. code-block:: bash

    vmd -m ./sirah.ff/tutorial/1/dna.pdb dna_cg.pdb

.. tip::

    This is the basic usage of the script **cgconv.pl**, you can learn other capabilities from its help by typing:

    .. code-block:: bash

    	./sirah.ff/tools/CGCONV/cgconv.pl -h	

From now on it is just normal GROMACS stuff!

.. caution::
	
	In GROMACS versions prior to 5.x, the "gmx" command should not be used.

1.2. PDB to GROMACS format
__________________________

Use ``pdb2gmx`` to convert your PDB file into GROMACS format: 

.. code-block:: bash

	gmx pdb2gmx -f dna_cg.pdb -o dna_cg.gro -merge all

When prompted, choose *SIRAH force field* and then *SIRAH solvent models*.

.. note::

	Warning messages about long, triangular or square bonds are fine and expected due to the CG topology of some residues.

.. caution::

	However, missing atom messages are errors which probably trace back to the
	mapping step. In that case, check your atomistic and mapped structures and do not carry on the
	simulation until the problem is solved.

.. note::

	Merging both DNA chains is convenient when planning to apply restraints between them.


During long simulations of DNA, capping residues may eventually separate. If you want to avoid this
effect, which is called helix fraying, add Watson-Crick (WC) restraints at terminal base pairs. To add
these restraints edit ``topol.top`` to include the file ``WC_RST.itp`` at the end of the [ moleculetype ] section:

.. list-table::
   :align: center
   :widths: 50 50
   :header-rows: 1

   * - Topology without WC restraints
     - Topology with WC restraints
   * - | ; Include Position restraint file
       | #ifdef POSRES
       | #include \"posre.itp\"
       | #endif
       | 
       | 
       |
       
     - | ; Include Position restraint file
       | #ifdef POSRES
       | #include \"posre.itp\"
       | #endif
       |
       | ; Watson-Crick restraints
       | #include \"./sirah.ff/tutorial/1/WC_RST.itp\"


1.3. Solvate the system
_______________________


Define the simulation box of the system

.. code-block:: bash 
	
	gmx editconf -f dna_cg.gro -o dna_cg_box.gro -bt octahedron -d 2 -c

Add WT4 molecules:

.. code-block:: bash 

	gmx solvate -cp dna_cg_box.gro -cs sirah.ff/wt416.gro -o dna_cg_sol.gro

.. note:: 

	In GROMACS versions prior to 5.x, the command *gmx solvate* was called *genbox*.

Edit the [ molecules ] section in ``topol.top`` to include the number of added WT4 molecules:

.. list-table::
   :align: center
   :widths: 50 50
   :header-rows: 1

   * - Topology before editing
     - Topology after editing
   * - | [ molecules ]
       | ; Compound        #mols
       | DNA_chain_A         1
       | 
              
     - | [ molecules ]
       | ; Compound        #mols
       | DNA_chain_A         1
       | WT4              3179

.. hint::
	
	If you forget to read the number of added WT4 molecules from the output of *solvate*, then use the following command line to get it 

	.. code-block:: console

		grep -c WP1 dna_cg_sol.gro

.. caution::
	
	The number of added WT4 molecules, **3179**, may change according to the software version.

Add CG counterions and 0.15M NaCl:

.. code-block:: bash

	gmx grompp -f sirah.ff/tutorial/1/GPU/em_CGDNA.mdp -p topol.top -c dna_cg_sol.gro -o dna_cg_sol.tpr

.. code-block:: bash

	gmx genion -s dna_cg_sol.tpr -o dna_cg_ion.gro -np 113 -pname NaW -nn 75 -nname ClW


When prompted, choose to substitute *WT4* molecules by *ions*.

.. note:: 

	The available electrolyte species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW) which represent solvated ions in solution. One ion pair (e.g., NaW-ClW) each 34 WT4 molecules results in a salt concentration of ~0.15M (see :ref:`Appendix <Appendix>` for details). Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`_.

Edit the [ molecules ] section in ``topol.top`` to include the CG ions and the correct number of WT4.

Before running the simulation it may be a good idea to visualize your molecular system. CG molecules
are not recognized by molecular visualizers and will not display correctly. To fix this problem you may
generate a PSF file of the system using the script ``g_top2psf.pl``:

.. code-block:: bash

	./sirah.ff/tools/g_top2psf.pl -i topol.top -o dna_cg_ion.psf

.. note::

	This is the basic usage of the script ``g_top2psf.pl``, you can learn other capabilities from its help:

	.. code-block:: bash

		./sirah.ff/tools/g_top2psf.pl -h


Use VMD to check how the CG system looks like:

.. code-block::

	vmd dna_cg_ion.psf dna_cg_ion.gro -e sirah.ff/tools/sirah_vmdtk.tcl

.. tip::

	VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right ones, according to the CG representation. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
	Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.

1.4. Run the simulation
________________________

.. important:: 

	By default, in this tutorial we will use input files for GROMACS on GPU (``sirah.ff/tutorial/1/GPU``). Example input files for using GROMACS on CPU can be found at: ``sirah.ff/tutorial/1/CPU``.

The folder ``sirah.ff/tutorial/1/GPU/`` contains typical input files for energy minimization
(``em_CGDNA.mdp``), equilibration (``eq_CGDNA.mdp``) and production (``md_CGDNA.mdp``) runs. Please
check carefully the input flags therein.

Create an index file:

.. code-block:: bash

	echo "q" | gmx make_ndx -f dna_cg_ion.gro -o dna_cg_ion.ndx

.. note::

	WT4 and CG ions (NaW and ClW) are automatically set to the group *SIRAH-Solvent*.

Make a new folder for the run:

.. code-block:: bash

	mkdir run; cd run

**Energy Minimization**:

.. code-block:: bash

	gmx grompp -f ../sirah.ff/tutorial/1/GPU/em_CGDNA.mdp -p ../topol.top -po em.mdp -n ../dna_cg_ion.ndx -c ../dna_cg_ion.gro -o dna_cg_em.tpr 

.. code-block:: bash

	mdrun -deffnm dna_cg_em &> EM.log &

**Equilibration**:

.. code-block:: bash 

	gmx grompp -f ../sirah.ff/tutorial/1/GPU/eq_CGDNA.mdp -p ../topol.top -po eq.mdp -n ../dna_cg_ion.ndx -c dna_cg_em.gro -r dna_cg_em.gro -o dna_cg_eq.tpr 

.. code-block:: bash 

	mdrun -deffnm dna_cg_eq &> EQ.log &

**Production (100ns)**:

.. code-block:: bash

	gmx grompp -f ../sirah.ff/tutorial/1/GPU/md_CGDNA.mdp -p ../topol.top -po md.mdp -n ../dna_cg_ion.ndx -c dna_cg_eq.gro -o dna_cg_md.tpr 

.. code-block:: bash

	mdrun -deffnm dna_cg_md &> MD.log &

.. note::

	GPU flags have been set for GROMACS 4.6.7; however, different versions may object to certain specifications.


1.5. Visualizing the simulation
________________________________

That’s it! Now you can analyze the trajectory.

Process the output trajectory at folder ``run/`` to account for the Periodic Boundary Conditions (PBC):

.. code-block:: bash

	gmx trjconv -s dna_cg_em.tpr -f dna_cg_md.xtc -o dna_cg_md_pbc.xtc -n ../dna_cg_ion.ndx -ur compact -center -pbc mol

When prompted, choose *DNA* for centering and *System* for output.

Now you can check the simulation using VMD:

.. code-block:: bash

	vmd ../dna_cg_ion.psf ../dna_cg_ion.gro dna_cg_md_pbc.xtc -e ../sirah.ff/tools/sirah_vmdtk.tcl


.. note::

    The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. Use the command ``sirah-help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.
