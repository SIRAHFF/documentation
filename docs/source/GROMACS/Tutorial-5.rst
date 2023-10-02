This tutorial shows how to use the SIRAH force field to perform a coarse grained (CG) simulation of a
DMPC bilayer in explicit solvent (called WatFour, WT4). The main references for
this tutorial are: `Barrera et al. <https://doi.org/10.1021/acs.jctc.9b00435>`_ and `Machado & Pantano <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_.
We strongly advise you to read these articles before starting the tutorial. You may also find interesting `this book chapter <https://pubs.aip.org/books/monograph/137/chapter-abstract/58880922/Simulating-Transmembrane-Proteins-with-the-Coarse?redirectedFrom=fulltext>`_.

.. important::

    Check :ref:`download <download gromacs>` section for download and set up details before to start this tutorial.
    Since this is **tutorial 5**, remember to replace ``X.X``, the files corresponding to this tutorial can be found in: ``sirah.ff/tutorial/5/``


5.1. Build CG representations
______________________________

Map the atomistic structure of the preassembled DMPC bilayer to its CG representation:

.. code-block:: bash

	./sirah.ff/tools/CGCONV/cgconv.pl -i sirah.ff/tutorial/5/DMPC64.pdb -o DMPC64_cg.pdb -a sirah.ff/tools/CGCONV/maps/tieleman_lipid.map

The input file ``-i`` DMPC64.pdb contains the atomistic representation of the DMPC bilayer, while the output ``-o`` DMPC64_cg.pdb is its SIRAH CG representation. The flag ``-a`` is a mapping scheme for lipids.

.. important::
	
	5' end residues mast be renamed to AW5, TW5, GW5 or CW5 to represent the corresponding Adenine, Thymine, Guanine or Cytosine extremes in a closed circular DNA.

	By default, no mapping is applied to lipids, as there is no standard naming convention for them. So users are requested to append a MAP file from the list in :ref:`Table 1 <table>`, by setting the flag ``-a`` in ``cgconv.pl``. We recommend using `PACKMOL <https://m3g.github.io/packmol/>`_ for building the system. Reference building-block structures are provided at folder ``sirah.ff/PDB/``, which agree with the mapping scheme in ``sirah.ff/tools/CGCONV/maps/tieleman_lipid.map``. The provided DMPC bilayer contains 64 lipid molecules per leaflet distributed in a 6.4 \* 6.4 nm surface, taking into account an approximate area per lipid of 0.64 nm\ :sup:`2` \ at 333 K. The starting configuration was created with the input file ``sirah.ff/tutorial/5/DMPC_bilayer.pkm``. See :doc:`FAQs <../FAQ>` for cautions on mapping lipids to SIRAH and tips on using fragment-based topologies. If you do not find your issue please start a discussion in our `github discussion page F&Q <https://github.com/SIRAHFF/documentation/discussions>`_.

.. tip::

  This an advanced usage of the script **cgconv.pl**, you can learn other capabilities from its help by typing:

  .. code-block:: bash

    ./sirah.ff/tools/CGCONV/cgconv.pl -h

The input file ``DMPC64.pdb`` contains all the heavy atoms composing the lipids, while the output ``DMPC64_cg.pdb`` preserves a few of them. Please check both PDB structures using VMD:	

.. code-block:: bash

  vmd -m sirah.ff/tutorial/5/DMPC64.pdb DMPC64_cg.pdb


From now on it is just normal GROMACS stuff!

.. caution::
	
	In GROMACS versions prior to 5.x, the "gmx" command should not be used.

5.2. PDB to GROMACS format
__________________________

Use pdb2gmx to convert your PDB file into GROMACS format: 

.. code-block:: bash

	gmx pdb2gmx -f DMPC64_cg.pdb -o DMPC64_cg.gro

When prompted, choose *SIRAH force field* and then *SIRAH solvent models*.

.. caution::

	Getting warning messages of long bonds is fine and expected due to the CG nature of the
	residue topologies. However missing atom messages are errors which probably trace back to the
	mapping step. In that case, check your atomistic and mapped structures and do not carry on the
	simulation until the problem is solved.


5.3. Solvate the system
_______________________


Define the simulation box of the system

.. code-block:: bash 
	
	gmx editconf -f DMPC64_cg.gro -o DMPC64_cg_box.gro -box 6.6 6.6 10 -c

.. note:: 

	As PACKMOL does not consider periodicity while building up the system, increasing the box size a few Angstroms may be required to avoid bad contacts between images.

Add WT4 molecules:

.. code-block:: bash 

	gmx solvate -cp DMPC64_cg_box.gro -cs sirah.ff/wt416.gro -o DMPC64_cg_sol1.gro

.. note:: 

	Before GROMACS version 5.x, the command *gmx solvate* was called *genbox*.

Edit the [ molecules ] section in ``topol.top`` to include the number of added WT4 molecules:

.. list-table::
   :align: center
   :widths: 50 50
   :header-rows: 1

   * - Topology before editing
     - Topology after editing
   * - | [ molecules ] 
       | ; Compound		#mols 
       | Lipid_chain_A		1    
       | Lipid_chain_B		1    
              
     - | [ molecules ] 
       | ; Compound		#mols 
       | Lipid_chain_A		1
       | Lipid_chain_B		1
       | WT4 			  850

.. hint::
	
	If you forget to read the number of added WT4 molecules from the output of *solvate*, then use the following command line to get it 

	.. code-block:: console

		grep -c WP1 DMPC64_cg_sol1.gro

.. caution::
	
	The number of added WT4 molecules **850** may change according to the software version.

Remove misplaced WT4 molecules inside the bilayer (read sirah.ff/0ISSUES and `FAQs <../FAQ>` for
comments on known solvation problems):

.. code-block:: bash
	
	gmx grompp -f sirah.ff/tutorial/5/CPU/em_CGLIP.mdp -p topol.top -c DMPC64_cg_sol1.gro -o DMPC64_cg_sol1.tpr

.. code-block:: bash
	
	echo -e "a BE* BC1* BC2* BCT*\n name 7 tail\n\nq" | gmx make_ndx -f DMPC64_cg_sol1.gro -o DMPC64_cg_sol1.ndx

.. code-block:: bash
	
	gmx select -f DMPC64_cg_sol1.gro -s DMPC64_cg_sol1.tpr -n DMPC64_cg_sol1.ndx -on rm_close_wt4.ndx -select "not (same residue as (resname WT4 and within 0.4 of group tail))"

.. code-block:: bash
	
	gmx editconf -f DMPC64_cg_sol1.gro -o DMPC64_cg_sol2.gro -n rm_close_wt4.ndx

Edit the [ molecules ] section in ``topol.top`` to correct the number of WT4 molecules:

.. hint::
	
	If you forget to read the number of added WT4 molecules from the output of *solvate*, then use the following command line to get it 

	.. code-block:: console

		grep -c WP1 DMPC64_cg_sol2.gro

Add CG counterions and 0.15M NaCl:

.. code-block:: bash

	gmx grompp -f sirah.ff/tutorial/5/GPU/em_CGLIP.mdp -p topol.top -c DMPC64_cg_sol2.gro -o DMPC64_cg_sol2.tpr

.. code-block:: bash

	gmx genion -s DMPC64_cg_sol2.tpr -o DMPC64_cg_ion.gro -np 21 -pname NaW -nn 21 -nname ClW


When prompted, choose to substitute **WT4** molecules by **ions**.

.. note:: 

	The available ionic species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW). One
	ion pair (e.g. NaW-ClW) each 34 WT4 molecules renders a salt concentration of ~0.15M (see :ref:`Appendix <Appendix>` for details). 
	Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`_.

Edit the [ molecules ] section in ``topol.top`` to include the CG ions and the correct number of WT4.

Before running the simulation it may be a good idea to visualize your molecular system. CG molecules are not recognized by molecular visualizers and will not display correctly. To fix this problem you may
generate a PSF file of the system using the script *g_top2psf.pl*:

.. code-block:: bash

	./sirah.ff/tools/g_top2psf.pl -i topol.top -o DMPC64_cg_ion.psf

.. note::

	This is the basic usage of the script ``g_top2psf.pl``, you can learn other capabilities from its help:
	``./sirah.ff/tools/g_top2psf.pl -h``


Use VMD to check how the CG system looks like:

.. code-block::

	vmd DMPC64_cg_ion.psf DMPC64_cg_ion.gro -e sirah.ff/tools/sirah_vmdtk.tcl

.. tip::

    VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right
    ones, according to the CG representation. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
    Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.

5.4. Run the simulation
________________________

.. important:: 

	By default in this tutorial we will use input files for GROMACS on GPU (``sirah.ff/tutorial/5/GPU``). Example input files for using GROMACS on CPU can be found at: ``sirah.ff/tutorial/5/CPU``.

The folder ``sirah.ff/tutorial/5/GPU/`` contains typical input files for energy minimization
(``em_CGLIP.mdp``), equilibration (``eq_CGLIP.mdp``) and production (``md_CGLIP.mdp``) runs. Please
check carefully the input flags therein.

Create an index file:

.. code-block:: bash

	echo "q" | gmx make_ndx -f DMPC64_cg_ion.gro -o DMPC64_cg_ion.ndx

.. note::

	WT4 and CG ions (NaW, ClW) are automatically set to the group “SIRAH-Solvent” while DMPC (named CMM at CG level) is assigned to group “Lipid”.

Make a new folder for the run:

.. code-block:: bash

	mkdir run; cd run

Energy Minimization:

.. code-block:: bash

	gmx grompp -f ../sirah.ff/tutorial/5/GPU/em_CGLIP.mdp -p ../topol.top -po em.mdp -n ../DMPC64_cg_ion.ndx -c ../DMPC64_cg_ion.gro -o DMPC64_cg_em.tpr 

.. code-block:: bash

	gmx mdrun -deffnm DMPC64_cg_em &> EM.log &

Equilibration:

.. code-block:: bash 

	gmx grompp -f ../sirah.ff/tutorial/5/GPU/eq_CGLIP.mdp -p ../topol.top -po eq.mdp -n ../DMPC64_cg_ion.ndx -c DMPC64_cg_em.gro -o DMPC64_cg_eq.tpr

.. code-block:: bash

	gmx mdrun -deffnm DMPC_cg_eq &> EQ.log &

Production (100ns):

.. code-block:: bash

	gmx grompp -f ../sirah.ff/tutorial/5/GPU/md_CGLIP.mdp -p ../topol.top -po md.mdp -n ../DMPC64_cg_ion.ndx -c DMPC64_cg_eq.gro -o DMPC64_cg_md.tpr

.. code-block:: bash

	gmx mdrun -deffnm DMPC64_cg_md &> MD.log &

.. note::

	GPU flags were set for GROMACS 4.6.7, different versions may complain about some specifications.

5.5. Visualizing the simulation
________________________________

That’s it! Now you can analyze the trajectory.

Process the output trajectory at folder ``run/`` to account for the Periodic Boundary Conditions (PBC):

.. code-block:: bash

	gmx trjconv -s DMPC64_cg_em.tpr -f DMPC64_cg_md.xtc -o DMPC64_cg_md_pbc.xtc -n ../DMPC64_cg_ion.ndx -pbc mol

When prompted, choose *System* for output.

Now you can check the simulation using VMD:

.. code-block:: bash

	vmd ../DMPC64_cg_ion.psf ../DMPC64_cg_ion.gro DMPC64_cg_md_pbc.xtc -e ../sirah.ff/tools/sirah_vmdtk.tcl

.. note::

    The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. Use the command ``sirah-help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.

Calculate the area per lipid

.. code-block:: bash

	gmx energy -f DMPC64_cg_md.edr -o box_XY.xvg

When prompted, choose "Box-X" and "Box-Y" for output. End your selection with zero.

.. note::

    To calculate the area per lipid, divide the membrane's area by the DMPC molecules per leaflet:   
	
	.. math::
		\frac{Area}{Lipid} = \frac{Box(x) * Box(y)}{64} 

Density profiles and bilayer thickness
Include a new group for phosphate beads in the index file: 

.. code-block:: bash

	echo -e "a BFO\nq\n" | gmx make_ndx -f DMPC64_cg_em.gro -n ../DMPC64_cg_ion.ndx -o DMPC64_cg_ion.ndx

.. code-block:: bash

	gmx density -sl 1000 -ng 5 -f DMPC64_cg_md_pbc.xtc -s DMPC64_cg_em.tpr -n DMPC64_cg_ion.ndx -o density_profile.xvg

When prompted, choose "Lipid", "WT4", "NaW", "ClW", and "BFO"

Use Grace to plot the results:

.. code-block:: bash

	xmgrace -nxy density_profile.xvg

.. note::

    The thickness of the bilayer is the distance between the two peaks corresponding to the position of phosphate beads (BFO) along the z-axis.