.. note::

   Please report bugs, errors or enhancement requests through `Issue Tracker <https://github.com/SIRAHFF/documentation/issues>`_ or if you have a question about SIRAH open a `New Discussion <https://github.com/SIRAHFF/documentation/discussions>`_.
   
This tutorial shows how to use the SIRAH force field to perform a coarse grained (CG) simulation of a
protein in explicit solvent (called WatFour, WT4). The main references for this tutorial are:  `Dans et al. <https://pubs.acs.org/doi/abs/10.1021/ct900653p>`__, `Darré et al. <https://pubs.acs.org/doi/abs/10.1021/ct100379f>`_, `Machado et al. <https://doi.org/10.1021/acs.jctc.9b00006>`__,  and `Machado & Pantano <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_. We strongly advise you to read these articles before starting the tutorial.

.. important::

    Check the :ref:`Setting up SIRAH <download gromacs>` section for download and set up details before starting this tutorial.
    Since this is **Tutorial 3**, remember to replace ``X.X``, the files corresponding to this tutorial can be found in: ``sirah.ff/tutorial/3/``

3.1. Build CG representations
______________________________

.. caution::

	The mapping to CG requires the correct protonation state of each residue at a given pH. We recommend using the `CHARMM-GUI server <https://www.charmm-gui.org/>`_ and use the **PDB Reader & Manipulator** to prepare your system. An account is required to access any of the CHARMM-GUI Input Generator modules, and it can take up to 24 hours to obtain one. 
	
	Other option is the `PDB2PQR server <https://server.poissonboltzmann.org/pdb2pqr>`_ and choosing the output naming scheme of AMBER for best compatibility. This server was utilized to generate the *PQR* file featured in this tutorial. Be aware that modified residues lacking parameters such as: MSE (seleno MET), TPO (phosphorylated THY), SEP (phosphorylated SER) or others are deleted from the PQR file by the server. In that case, mutate the residues to their unmodified form before submitting the structure to the server.

Map the protonated atomistic structure of protein `1CRN <https://www.rcsb.org/structure/1CRN>`_ to its CG representation:

.. code-block:: bash 
	
	./sirah.ff/tools/CGCONV/cgconv.pl -i sirah.ff/tutorial/3/1CRN.pqr -o 1CRN_cg.pdb

The input file ``-i`` 1CRN.pqr contains the atomistic representation of `1CRN <https://www.rcsb.org/structure/1CRN>`_ structure at pH **7.0**, while the output ``-o`` 1CRN_cg.pdb is its SIRAH CG representation.

.. tip::

	This is the basic usage of the script **cgconv.pl**, you can learn other capabilities from its help by typing:

	.. code-block:: bash

		./sirah.ff/tools/CGCONV/cgconv.pl -h	

.. note:: 

	**Pay attention to residue names when mapping structures from other atomistic force fields or experimental structures.** Although we provide compatibility for naming schemes in PDB, GMX, GROMOS, CHARMM and OPLS, there might always be some ambiguity in the residue naming, specially regarding protonation states, that may lead to a wrong mapping. For example, SIRAH Tools always maps the residue name “HIS” to a Histidine protonated at the epsilon nitrogen (:math:`N_{\epsilon}`) regardless the actual proton placement. Similarly, protonated Glutamic and Aspartic acid residues must be named “GLH” and “ASH”, otherwise they will be treated as negative charged residues. In addition, protonated and disulfide bonded Cysteines must be named “CYS” and “CYX” respectively. These kind of situations need to be carefully checked by the users. In all cases the residues preserve their identity when mapping and back-mapping the structures. Hence, the total charge of the protein should be the same at atomistic and SIRAH levels. You can check the following mapping file to be sure of the compatibility: ``sirah.ff/tools/CGCONV/maps/sirah_prot.map``.    

Please check both PDB and PQR structures using VMD:

.. code-block:: bash 
	
	vmd -m sirah.ff/tutorial/3/1CRN.pqr 1CRN_cg.pdb

From now on it is just normal GROMACS stuff!

.. caution::
	
	In GROMACS versions prior to 5.x, the "gmx" command should not be used.


3.2. PDB to GROMACS format
__________________________

Use ``pdb2gmx`` to convert your PDB file into GROMACS format: 

.. code-block:: bash

	gmx pdb2gmx -f 1CRN_cg.pdb -o 1CRN_cg.gro

When prompted, choose *SIRAH force field* and then *SIRAH solvent models*.

.. note:: 

	By default, charged terminal are used, but it is possible to set them neutral with option ``-ter``.

.. note::

	Warning messages about long, triangular or square bonds are fine and expected due to the CG topology of some residues.

.. caution::

	However, missing atom messages are errors which probably trace back to the
	mapping step. In that case, check your atomistic and mapped structures and do not carry on the
	simulation until the problem is solved.

3.3. Solvate the system
________________________

Define the simulation box of the system

.. code-block:: bash 

	gmx editconf -f 1CRN_cg.gro -o 1CRN_cg_box.gro -bt octahedron -d 2.0 -c

Add WT4 molecules:

.. code-block:: bash 
	
	gmx solvate -cp 1CRN_cg_box.gro -cs sirah.ff/wt416.gro -o 1CRN_cg_sol1.gro

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
       | Protein_chain_A     1
       |  
              
     - | [ molecules ]
       | ; Compound        #mols
       | Protein_chain_A     1
       | WT4               756

.. hint::
	
	If you forget to read the number of added WT4 molecules from the output of *solvate*, then use the following command line to get it

	.. code-block:: console

		grep -c WP1 1CRN_cg_sol1.gro

.. caution::
	
	The number of added WT4 molecules, **756**, may change according to the software version.

Remove WT4 molecules within 0.3 nm of protein:

.. code-block:: bash

	echo q | gmx make_ndx -f 1CRN_cg_sol1.gro -o 1CRN_cg_sol1.ndx

.. code-block:: bash 

	gmx grompp -f sirah.ff/tutorial/3/GPU/em1_CGPROT.mdp -p topol.top -po delete1.mdp -c 1CRN_cg_sol1.gro -o 1CRN_cg_sol1.tpr

.. code-block:: bash 

	gmx select -f 1CRN_cg_sol1.gro -s 1CRN_cg_sol1.tpr -n 1CRN_cg_sol1.ndx -on rm_close_wt4.ndx -select 'not (same residue as (resname WT4 and within 0.3 of group Protein))'

.. code-block:: bash 

	gmx editconf -f 1CRN_cg_sol1.gro -o 1CRN_cg_sol2.gro -n rm_close_wt4.ndx

.. note:: 
	
	In GROMACS versions prior to 5.x, the command *gmx select* was called *g_select*.

Edit the [ molecules ] section in ``topol.top`` to include the correct number of WT4 molecules:

.. code-block:: bash

	grep -c WP1 1CRN_cg_sol2.gro

Add CG counterions and 0.15M NaCl:

.. code-block:: bash

	gmx grompp -f sirah.ff/tutorial/3/GPU/em1_CGPROT.mdp -p topol.top -po delete2.mdp -c 1CRN_cg_sol2.gro -o 1CRN_cg_sol2.tpr

.. code-block:: bash

	gmx genion -s 1CRN_cg_sol2.tpr -o 1CRN_cg_ion.gro -np 22 -pname NaW -nn 22 -nname ClW

When prompted, choose to substitute *WT4* molecules by *ions*.

.. note:: 

	The available electrolyte species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW) which represent solvated ions in solution. One ion pair (e.g., NaW-ClW) each 34 WT4 molecules results in a salt concentration of ~0.15M (see :ref:`Appendix <Appendix>` for details). Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`_.

Edit the [ molecules ] section in ``topol.top`` to include the CG ions and the correct number of WT4.

Before running the simulation it may be a good idea to visualize your molecular system. CG molecules
are not recognized by molecular visualizers and will not display correctly. To fix this problem you may
generate a PSF file of the system using the script ``g_top2psf.pl``:

.. code-block:: bash

	./sirah.ff/tools/g_top2psf.pl -i topol.top -o 1CRN_cg_ion.psf

.. note::
	
	This is the basic usage of the script ``g_top2psf.pl``, you can learn other capabilities from its help:

	.. code-block:: bash

		./sirah.ff/tools/g_top2psf.pl -h


Use VMD to check how the CG system looks like:

.. code-block::

	vmd 1CRN_cg_ion.psf 1CRN_cg_ion.gro -e sirah.ff/tools/sirah_vmdtk.tcl

.. note::

	VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right ones, according to the CG representation. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
	Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.

Create an index file including a group for the backbone GN and GO beads:

.. code-block:: bash 

	echo -e "a GN GO\n\nq" | gmx make_ndx -f 1CRN_cg_ion.gro -o 1CRN_cg_ion.ndx

.. note::

	WT4 and CG ions (NaW, ClW) are automatically set to the group *SIRAH-Solvent*.

Generate restraint files for the backbone *GN* and *GO* beads:

.. code-block:: bash

	gmx genrestr -f 1CRN_cg.gro -n 1CRN_cg_ion.ndx -o bkbres.itp

.. code-block:: bash

	gmx genrestr -f 1CRN_cg.gro -n 1CRN_cg_ion.ndx -o bkbres_soft.itp -fc 100 100 100

When prompted, choose the group *GN_GO*

Add the restraints to ``topol.top``:

.. list-table:: 
   :align: center
   :widths: 50 50
   :header-rows: 1

   * - Topology before editing
     - Topology after editing
   * - | ; Include Position restraint file
       | #ifdef POSRES
       | #include \"posre.itp\"
       | #endif
       
       | 
       | 
       | 
       | 
       
       | 
       | 
       | 
       |
       
                  
     - | ; Include Position restraint file
       | #ifdef POSRES
       | #include \"posre.itp\"
       | #endif
        
       | ; Backbone restraints
       | #ifdef GN_GO
       | #include \"bkbres.itp\"
       | #endif
	   
       | ; Backbone soft restrains
       | #ifdef GN_GO_SOFT
       | #include \"bkbres_soft.itp\"
       | #endif

3.4. Run the simulation
________________________

.. important:: 

	By default, in this tutorial we will use input files for GROMACS on GPU (``sirah.ff/tutorial/3/GPU``). Example input files for using GROMACS on CPU can be found at: ``sirah.ff/tutorial/3/CPU``.

The folder ``sirah.ff/tutorial/3/GPU/`` contains typical input files for energy minimization
(``em1_CGPROT.mdp``, ``em2_CGPROT.mdp``), equilibration (``eq1_CGPROT.mdp``, ``eq2_CGPROT.mdp``)
and production (``md_CGPROT.mdp``) runs. Please check carefully the input flags therein.

Make a new folder for the run:

.. code-block:: bash 

	mkdir -p run; cd run

**Energy Minimization of side chains by restraining the backbone**:

.. code-block:: bash 

	gmx grompp -f ../sirah.ff/tutorial/3/GPU/em1_CGPROT.mdp -p ../topol.top -po em1.mdp -n ../1CRN_cg_ion.ndx -c ../1CRN_cg_ion.gro -r ../1CRN_cg_ion.gro -o 1CRN_cg_em1.tpr 

.. code-block:: bash 
	
	gmx mdrun -deffnm 1CRN_cg_em1 &> EM1.log &


**Energy Minimization of whole system**:

.. code-block:: bash 

	gmx grompp -f ../sirah.ff/tutorial/3/GPU/em2_CGPROT.mdp -p ../topol.top -po em2.mdp -n ../1CRN_cg_ion.ndx -c 1CRN_cg_em1.gro -o 1CRN_cg_em2.tpr

.. code-block:: bash 

	gmx mdrun -deffnm 1CRN_cg_em2 &> EM2.log &

**Solvent equilibration**:

.. code-block:: bash 

	gmx grompp -f ../sirah.ff/tutorial/3/GPU/eq1_CGPROT.mdp -p ../topol.top -po eq1.mdp -n ../1CRN_cg_ion.ndx -c 1CRN_cg_em2.gro -r 1CRN_cg_em2.gro -o 1CRN_cg_eq1.tpr

.. code-block:: bash 

	gmx mdrun -deffnm 1CRN_cg_eq1 &> EQ1.log &

**Soft equilibration to improve side chain solvation**:

.. code-block:: bash

	gmx grompp -f ../sirah.ff/tutorial/3/GPU/eq2_CGPROT.mdp -p ../topol.top -po eq2.mdp -n ../1CRN_cg_ion.ndx -c 1CRN_cg_eq1.gro -r 1CRN_cg_eq1.gro -o 1CRN_cg_eq2.tpr

.. code-block:: bash

	gmx mdrun -deffnm 1CRN_cg_eq2 &> EQ2.log &

**Production (1000 ns)**:

.. code-block:: bash

	gmx grompp -f ../sirah.ff/tutorial/3/GPU/md_CGPROT.mdp -p ../topol.top -po md.mdp -n ../1CRN_cg_ion.ndx -c 1CRN_cg_eq2.gro -o 1CRN_cg_md.tpr

.. code-block:: bash

	gmx mdrun -deffnm 1CRN_cg_md &> MD.log &

.. note::

	GPU flags have been set for GROMACS 4.6.7; however, different versions may object to certain specifications.

3.5. Visualizing the simulation
________________________________

That’s it! Now you can analyze the trajectory.

Process the output trajectory at folder ``run/`` to account for the Periodic Boundary Conditions (PBC):

.. code-block:: bash

	gmx trjconv -s 1CRN_cg_em1.tpr -f 1CRN_cg_md.xtc -o 1CRN_cg_md_pbc.xtc -n ../1CRN_cg_ion.ndx -ur compact -center -pbc mol

When prompted, choose *Protein* for centering and *System* for output.

Now you can check the simulation using VMD:

.. code-block:: bash 

	vmd ../1CRN_cg_ion.psf ../1CRN_cg_ion.gro 1CRN_cg_md_pbc.xtc -e ../sirah.ff/tools/sirah_vmdtk.tcl

.. note::

    The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. Use the command ``sirah-help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.


