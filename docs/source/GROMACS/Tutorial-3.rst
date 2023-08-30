This tutorial shows how to use the SIRAH force field to perform a coarse grained (CG) simulation of a
protein in explicit solvent (called WatFour, WT4) in four simple steps: 1 download; 2 map; 3 solvate
and; 4 run. The main references for this tutorial are:  `WAT4 <https://pubs.acs.org/doi/abs/10.1021/ct100379f>`_, `SIRAH 2.0 <https://doi.org/10.1021/acs.jctc.9b00006>`_ and `SIRAH Tools <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_. We strongly advise you to read these articles before starting the tutorial.

.. important::

    Check :ref:`download <download gromacs>` section for download and set up details before to start this tutorial.
    Since this is **tutorial 3**, remember to replace ``X.X``, the files corresponding to this tutorial can be found in: ``sirah.ff/tutorial/3/``

3.1. Build CG representations
______________________________

Map the protonated atomistic structure of protein 1CRN to its CG representation:

.. code-block:: bash 
	
	./sirah.ff/tools/CGCONV/cgconv.pl -i sirah.ff/tutorial/3/1CRN.pqr -o 1CRN_cg.pdb

.. note:: 

	The mapping to CG requires the correct protonation state of each residue at a given pH. We
	recommend using the `PDB2PQR server <https://server.poissonboltzmann.org/pdb2pqr>`_ and choosing the output
	naming scheme of AMBER for best compatibility. Be aware that modified residues lacking parameters
	such as: MSE (seleno MET), TPO (phosphorylated THY), SEP (phosphorylated SER) or others are
	deleted from the PQR file by the server. In that case, mutate the residues to their unmodified form
	before submitting the structure to the server.

.. warning:: 

	Pay attention to residue names when mapping structures from other atomistic force fields or
	experimental structures. Although we provide compatibility for naming schemes in PDB, GMX,
	GROMOS, CHARMM and OPLS, there always may be some ambiguity in the residue naming,
	specially regarding protonation states, that may lead to a wrong mapping. For example, SIRAH Tools
	always maps the residue name “HIS” to a Histidine protonated at Ne regardless the actual proton
	placement. Similarly, protonated Glutamic and Aspartic acid residues must be named “GLH” and
	“ASH”, otherwise they will be treated as negative charged residues. In addition, protonated and
	disulfide bonded Cysteines must be named “CYS” and “CYX” respectively. These kind of situations
	need to be carefully checked by the users. In all cases the residues preserve their identity when
	mapping and back-mapping the structures. Hence, the total charge of the protein should be the same
	at atomistic and SIRAH level. You can check the following mapping file to be sure of the compatibility:``sirah.ff/tools/CGCONV/maps/sirah_prot.map``.

The input file 1CRN.pqr contains all the heavy atoms composing the protein, while the output
1CRN_cg.pdb preserves a few of them. Please check both PDB and PQR structures using VMD:

.. code-block:: bash 
	
	vmd -m sirah.ff/tutorial/3/1CRN.pqr 1CRN_cg.pdb

.. tip::

    This is the basic usage of the script **cgconv.pl**, you can learn other capabilities from its help:
    ``./sirah.amber/tools/CGCONV/cgconv.pl -h``

From now on it is just normal GROMACS stuff!

.. caution::
	
	In GROMACS versions prior to 5.x, the "gmx" command should not be used.


3.2. PDB to GROMACS format
__________________________

Use pdb2gmx to convert your PDB file into GROMACS format: 

.. code-block:: bash

	gmx pdb2gmx -f 1CRN_cg.pdb -o 1CRN_cg.gro

When prompted, choose **SIRAH force field** and then **SIRAH solvent models**.

.. note:: 

	By default charged terminal are used but it is possible to set them neutral with option ``-ter``

.. caution::

	Getting warning messages of long bonds is fine and expected due to the CG nature of the
	residue topologies. However missing atom messages are errors which probably trace back to the
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

	In GROMACS in versions earlier than 5.x the command  *gmx solvate* was named to *genbox*.

Edit the [ molecules ] section in topol.top to include the number of added WT4 molecules:

.. list-table::
   :align: center
   :widths: 50 50
   :header-rows: 1

   * - Topology before editing
     - Topology after editing
   * - | [ molecules ]
       | ; Compound #mols
       | Protein_chain_A 1
       | 
              
     - | [ molecules ]
       | ; Compound #mols
       | Protein_chain_A 1
       | WT4 756

.. hint::
	
	If you forget to read the number of added WT4 molecules from the output of *solvate*, then use the following command line to get it

	.. code-block:: console

		grep -c WP1 1CRN_cg_sol1.gro

.. caution::
	
	The number of added WT4 molecules **756** may change according to the software version.

Remove WT4 molecules within 0.3 nm of protein:

.. code-block:: bash

	echo q | make_ndx -f 1CRN_cg_sol1.gro -o 1CRN_cg_sol1.ndx

.. code-block:: bash 

	gmx grompp -f sirah.ff/tutorial/3/CPU/em1_CGPROT.mdp -p topol.top -po delete1.mdp -c 1CRN_cg_sol1.gro -o 1CRN_cg_sol1.tpr

.. code-block:: bash 

	gmx select -f 1CRN_cg_sol1.gro -s 1CRN_cg_sol1.tpr -n 1CRN_cg_sol1.ndx -on rm_close_wt4.ndx -select 'not (same residue as (resname WT4 and within 0.3 of group Protein))'

.. code-block:: bash 

	gmx editconf -f 1CRN_cg_sol1.gro -o 1CRN_cg_sol2.gro -n rm_close_wt4.ndx

.. note:: 
	
	In GROMACS in versions earlier than 5.x the command  *gmx select* was named to *g_select*.

Edit the [ molecules ] section in topol.top to include the correct number of WT4 molecules:

.. code-block:: bash

	grep -c WP1 1CRN_cg_sol2.gro

Add CG counterions and 0.15M NaCl:

.. code-block:: bash

	gmx grompp -f sirah.ff/tutorial/3/CPU/em1_CGPROT.mdp -p topol.top -po delete2.mdp -c 1CRN_cg_sol2.gro -o 1CRN_cg_sol2.tpr

.. code-block:: bash

	genion -s 1CRN_cg_sol2.tpr -o 1CRN_cg_ion.gro -np 22 -pname NaW -nn 22 -nname ClW

When prompted, choose to substitute **WT4** molecules by **ions**.

.. note:: 

	The available ionic species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW). One
	ion pair (e.g. NaW-ClW) each 34 WT4 molecules renders a salt concentration of ~0.15M (see :ref:`Appendix <Appendix>` for details). 
	Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`_.

Edit the [ molecules ] section in ``topol.top`` to include the CG ions and the correct number of WT4.

Before running the simulation it may be a good idea to visualize your molecular system. CG molecules
are not recognized by molecular visualizers and will not display correctly. To fix this problem you may
generate a PSF file of the system using the script *g_top2psf.pl*:

.. code-block:: bash

	./sirah.ff/tools/g_top2psf.pl -i topol.top -o 1CRN_cg_ion.psf

.. note::

	This is the basic usage of the script ``g_top2psf.pl``, you can learn other capabilities from its help:
	``./sirah.ff/tools/g_top2psf.pl -h``


Use VMD to check how the CG system looks like:

.. code-block::

	vmd 1CRN_cg_ion.psf 1CRN_cg_ion.gro -e sirah.ff/tools/sirah_vmdtk.tcl

.. note::

	VMD assigns default radius to unknown atom types, the script sirah_vmdtk.tcl sets the right
	ones. It also provides a kit of useful selection macros, coloring methods and a backmapping utility.
	Use the command sirah_help in the Tcl/Tk console of VMD to access the manual pages.

Create an index file including a group for the backbone GN and GO beads:

.. code-block:: bash 

	echo -e "a GN GO\n\nq" | gmx make_ndx -f 1CRN_cg_ion.gro -o 1CRN_cg_ion.ndx

.. note::

	WT4 and CG ions (NaW, ClW) are automatically set to the group **SIRAH-Solvent**.

Generate restraint files for the backbone GN and GO beads:

.. code-block:: bash

	genrestr -f 1CRN_cg.gro -n 1CRN_cg_ion.ndx -o bkbres.itp

.. code-block:: bash

	genrestr -f 1CRN_cg.gro -n 1CRN_cg_ion.ndx -o bkbres_soft.itp -fc 100 100 100

When prompted, choose the group **GN_GO**

Add restraints to topol.top

.. list-table:: restraints to topol.top
   :align: center
   :widths: 50 50
   :header-rows: 1

   * - Topology before editing
     - Topology after editing
   * - | ; Include Position restraint file
       | #ifdef POSRES
       | #include "posre.itp"
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
       | #include "posre.itp"
       | #endif
        
       | ; Backbone restraints
	   | #ifdef GN_GO
       | #include "bkbres.itp"
       | #endif
	   
       | ; Backbone soft restrains
	   | #ifdef GN_GO_SOFT
       | #include "bkbres_soft.itp"
       | #endif

3.4. Run the simulation
________________________

.. important:: 

	By default in this tutorial we will use input files for GROMACS on GPU (``sirah.ff/tutorial/3/GPU``). Example input files for using GROMACS on CPU can be found at: ``sirah.ff/tutorial/3/CPU``.

The folder ``sirah.ff/tutorial/3/GPU/`` contains typical input files for energy minimization
(``em1_CGPROT.mdp``, ``em2_CGPROT.mdp``), equilibration (``eq1_CGPROT.mdp``, ``eq2_CGPROT.mdp``)
and production (``md_CGPROT.mdp``) runs. Please check carefully the input flags therein.

Make a new folder for the run:

.. code-block:: bash 

	mkdir -p run; cd run

**Energy Minimization of side chains by restraining the backbone:**

.. code-block:: bash 

	gmx grompp -f ../sirah.ff/tutorial/3/CPU/em1_CGPROT.mdp -p ../topol.top -po em1.mdp -n ../1CRN_cg_ion.ndx -c ../1CRN_cg_ion.gro -r ../1CRN_cg_ion.gro -o 1CRN_cg_em1.tpr 

.. code-block:: bash 
	
	gmx mdrun -deffnm 1CRN_cg_em1 &> EM1.log &



**Energy Minimization of whole system:**

.. code-block:: bash 

	gmx grompp -f ../sirah.ff/tutorial/3/CPU/em2_CGPROT.mdp -p ../topol.top -po em2.mdp -n ../1CRN_cg_ion.ndx -c 1CRN_cg_em1.gro -o 1CRN_cg_em2.tpr

.. code-block:: bash 

	gmx mdrun -deffnm 1CRN_cg_em2 &> EM2.log &

**Solvent equilibration:**

.. code-block:: bash 

	gmx grompp -f ../sirah.ff/tutorial/3/CPU/eq1_CGPROT.mdp -p ../topol.top -po eq1.mdp -n ./1CRN_cg_ion.ndx -c 1CRN_cg_em2.gro -r 1CRN_cg_em2.gro -o 1CRN_cg_eq1.tpr

.. code-block:: bash 

	gmx mdrun -deffnm 1CRN_cg_eq1 &> EQ1.log &

**Soft equilibration to improve side chain solvation:**

.. code-block:: bash

	gmx grompp -f ../sirah.ff/tutorial/3/CPU/eq2_CGPROT.mdp -p ../topol.top -po eq2.mdp -n ../1CRN_cg_ion.ndx -c 1CRN_cg_eq1.gro -r 1CRN_cg_eq1.gro -o 1CRN_cg_eq2.tpr

.. code-block:: bash

	gmx mdrun -deffnm 1CRN_cg_eq2 &> EQ2.log &

**Production** (1 :math:`\mu s`):

.. code-block:: bash

	gmx grompp -f ../sirah.ff/tutorial/3/CPU/md_CGPROT.mdp -p ../topol.top -po md.mdp -n ../1CRN_cg_ion.ndx -c 1CRN_cg_eq2.gro -o 1CRN_cg_md.tpr

.. code-block:: bash

	gmx mdrun -deffnm 1CRN_cg_md &> MD.log &

.. note::

	GPU flags were set for GROMACS 4.6.7, different versions may complain about some specifications.

3.5. Visualizing the simulation
________________________________

That’s it! Now you can analyze the trajectory.

Process the output trajectory at folder run/ to account for the Periodic Boundary Conditions (PBC):

.. code-block:: bash

	gmx trjconv -s 1CRN_cg_em1.tpr -f 1CRN_cg_md.xtc -o 1CRN_cg_md_pbc.xtc -n ../1CRN_cg_ion.ndx -ur compact -center -pbc mol

When prompted, choose **Protein** for centering and **System** for output.

Now you can check the simulation using VMD:

.. code-block:: bash 

	vmd ../1CRN_cg_ion.psf ../1CRN_cg_ion.gro 1CRN_cg_md_pbc.xtc -e ../sirah.ff/tools/sirah_vmdtk.tcl


3.6. Calculate the solvent accessible surface (SAS)
____________________________________________________

Create the following symbolic link in the folder run/:

.. code-block:: bash 

	ln -s ../sirah.ff/vdwradii.dat

Calculate the SAS of the protein along the trajectory:

.. code-block:: bash 

	g_sas -s 1CRN_cg_md.tpr -f 1CRN_cg_md_pbc.xtc -n ../1CRN_cg_ion.ndx -qmax 0 -probe 0.21 -o area.xvg

When prompted, choose **Protein** as both the group for calculation and the output.

.. note:: 

	The solvent probe radius corresponds to a WT4 bead while a charge of 0e refers to any
	hydrophobic bead. The file vdwradii.dat must be placed at the same folder where *gmx sasa* is executed to 	assure that the correct van der Waals radii of SIRAH beads are used in the calculation.

.. important:: 

	g_sas is deprecated, the tool no longer automatically divides the surface into hydrophobic and hydrophilic areas, and there is no -f_index option. The same effects can be obtained by defining suitable selections for -output. If you want output that contains the same numbers as with the old tool for a calculation group A and output group B, you can use `[1] <https://manual.gromacs.org/current/user-guide/cmdline.html>`_. ::

	 gmx sasa -surface 'group "A"' -output '"Hydrophobic" group "A" and charge {-0.2 to 0.2}; "Hydrophilic" group "B" and not charge {-0.2 to 0.2}; "Total" group "B"'


You can use Grace to plot the results::
	
	xmgrace -nxy area.xvg

3.7. Visualize the secondary structure
________________________________________


Load the processed trajectory in VMD::

	vmd ../1CRN_cg_ion.psf ../1CRN_cg_ion.gro 1CRN_cg_md_pbc.xtc -e ../sirah.ff/tools/sirah_vmdtk.tcl

.. note::

    The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. Use the command ``sirah-help`` in the Tcl/Tk console of VMD to access the manual pages.

At the **Tk/Tcl console** run the command ``sirah_ss`` to get the secondary structure of the CG protein.

.. note:: 
	
	After assigning the secondary structure it is possible to represent a-helices with Bendix in VMD
	1.9.2 or upper by setting the backbone particle name to GC (do not check the CG box).

To analyze the output files from sirah_ss, go back at the shell command line and execute::

	xmgrace -nxy ss_by_frame.xvg

.. code-block:: bash 

	xmgrace -nxy ss_by_res.xvg

The file ss.mtx can be processed to visualize the time evolution of the secondary structure by residue::

	../sirah.ff/tools/ssmtx2png.R --mtx=ss.mtx

.. code-block:: bash

	display ssmtx.png

.. hint::

	The usage of ssmtx2png.R can be accessed through::

	../sirah.ff/tools/ssmtx2png.R --help
