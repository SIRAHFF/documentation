.. note::

   Please report bugs, errors or enhancement requests through `Issue Tracker <https://github.com/SIRAHFF/documentation/issues>`_ or if you have a question about SIRAH open a `New Discussion <https://github.com/SIRAHFF/documentation/discussions>`_.
   
This tutorial shows how to use the SIRAH force field to perform a coarse grained (CG) simulation of a protein embedded in a lipid bilayer using explicit solvent (called WatFour, WT4). The main references for
this tutorial are: `Machado et al. <https://doi.org/10.1021/acs.jctc.9b00006>`__, `Barrera et al. <https://doi.org/10.1021/acs.jctc.9b00435>`_, and `Machado & Pantano <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_.
We strongly advise you to read these articles before starting the tutorial.

.. note::

	We strongly advise you to read and complete :ref:`Tutorial 3 <Tutorial 3 G>` and :ref:`Tutorial 5 <Tutorial 5 G>` before starting.
	

.. important::

    Check the :ref:`Setting up SIRAH <download gromacs>` section for download and set up details before starting this tutorial.
    Since this is **Tutorial 6**, remember to replace ``X.X``, the files corresponding to this tutorial can be found in: ``sirah.ff/tutorial/6/``
	
6.1. Build CG representations
______________________________

.. caution::

  The mapping to CG requires the correct protonation state of each residue at a given pH. We recommend using the `CHARMM-GUI server <https://www.charmm-gui.org/>`_ and use the **PDB Reader & Manipulator** to prepare your system. An account is required to access any of the CHARMM-GUI Input Generator modules, and it can take up to 24 hours to obtain one. 
  
  Other option is the `PDB2PQR server <https://server.poissonboltzmann.org/pdb2pqr>`_ and choosing the output naming scheme of AMBER for best compatibility. This server was utilized to generate the *PQR* file featured in this tutorial. Be aware that modified residues lacking parameters such as: MSE (seleno MET), TPO (phosphorylated THY), SEP (phosphorylated SER) or others are deleted from the PQR file by the server. In that case, mutate the residues to their unmodified form before submitting the structure to the server.

  See :ref:`Tutorial 3 <Tutorial 3 G>` for cautions while preparing and mapping atomistic proteins to SIRAH.
	
Map the atomistic structure of protein `2KYV <https://www.rcsb.org/structure/2KYV>`__ to its CG representation:  

.. code-block:: bash

  ./sirah.ff/tools/CGCONV/cgconv.pl -i sirah.ff/tutorial/6/2kyv.pqr -o 2kyv_cg.pdb 
  
The input file ``-i`` 2kyv.pqr contains the atomistic representation of `2KYV <https://www.rcsb.org/structure/2KYV>`__ structure at pH **7.0**, while the output ``-o`` 2kyv_cg.pdb is its SIRAH CG representation. 

.. important::

	If you already have an atomistic protein within a membrane, then you can simply map the entire system to SIRAH (this is highly recommended) and skip the step of embedding the protein into a lipid bilayer, however clipping the membrane patch may be required to set a correct solvation box (see bellow). By default, no mapping is applied to lipids, as there is no standard naming convention for them. So users are requested to append a MAP file from the list in :ref:`Table 1 <table>`, by setting the flag ``-a`` in ``cgconv.pl``. We recommend using `PACKMOL <https://m3g.github.io/packmol/>`__ for building the system. Reference building-block structures are provided at folder ``sirah.ff/PDB/``, which agree with the mapping scheme ``sirah.ff/tools/CGCONV/maps/tieleman_lipid.map``. See :doc:`FAQs <../FAQ>` for cautions on mapping lipids to SIRAH and tips on using fragment-based topologies.  

.. tip::

  This an advanced usage of the script **cgconv.pl**, you can learn other capabilities from its help by typing:

  .. code-block:: bash

    ./sirah.ff/tools/CGCONV/cgconv.pl -h
	

6.1.1. Embed the protein in a lipid bilayer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We will show one possible way to do it by starting from a pre-equilibrated CG membrane patch.

In ``sirah_[version].ff/tutorial/6/`` you will find a pre-stabilized CG DMPC bilayer patch, concatenate it with the previously generated CG representation of the phospholamban (PLN) pentamer (PDB code: `2KYV <https://www.rcsb.org/structure/2KYV>`__):

.. code-block:: bash

	head -qn -1 2kyv_cg.pdb ./sirah.ff/tutorial/6/DMPC_cg.pdb > 2kyv_DMPC_cg_init.pdb

Luckily, we already oriented the protein inside the membrane. For setting up your own system you can go to `Orientations of Proteins in Membranes (OPM) database <https://opm.phar.umich.edu/>`__ and (if your structure is available) use the dummy atoms provided there to make them match with your membrane model (see **Figure 1**).

.. figure:: /../images/Tuto7.png
   :align: center
   :width: 100%

   **Figure 1.** Protein oriented inside the membrane from the OPM database with dummy atoms represented as orange spheres.
   

6.1.2. Delete close contact lipid molecules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We need to use VMD to to delete lipid molecules in close contact with the protein. For a proper treatment and visualization of the system in VMD you must first generate the molecular topology and initial coordinate files.

Use ``pdb2gmx`` to convert your PDB file into GROMACS format: 

.. code-block:: bash

  gmx pdb2gmx -f 2kyv_DMPC_cg_init.pdb -o 2kyv_DMPC_cg_init.gro -p init_topol -i init_posre

When prompted, choose *SIRAH force field* and then *SIRAH solvent models*.

.. caution::
  
  In GROMACS versions prior to 5.x, the "gmx" command should not be used.

CG molecules are not recognized by molecular visualizers and will not display correctly. To fix this problem you may
generate a PSF file of the system using the script ``g_top2psf.pl``:

.. code-block:: bash

  ./sirah.ff/tools/g_top2psf.pl -i init_topol.top -o 2kyv_DMPC_cg_init.psf

.. note::

  This is the basic usage of the script ``g_top2psf.pl``, you can learn other capabilities from its help:
  
  .. code-block:: bash

    ./sirah.ff/tools/g_top2psf.pl -h


Use VMD to check how the CG system looks like:

.. code-block::

  vmd 2kyv_DMPC_cg_init.psf 2kyv_DMPC_cg_init.gro -e sirah.ff/tools/sirah_vmdtk.tcl

.. tip::

    VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right
    ones, according to the CG representation. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
    Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.

In the VMD main window, select *Graphics* > *Representations*. In the *Selected Atoms* box, type:

.. code-block:: text

	not (same residue as (sirah_membrane within 3.5 of sirah_protein))


To save the refined protein-membrane system, in the VMD main window click on ``2kyv_DMPC_cg_init.psf``, then select *File* > *Save Coordinates*. In the *Selected atoms* option choose the selection you have just created and Save as ``2kyv_DMPC_cg.pdb``.
	
From now on it is just normal GROMACS stuff!

6.2. PDB to GROMACS format
__________________________

Use ``pdb2gmx`` to convert your PDB file into GROMACS format: 

.. code-block:: bash

  gmx pdb2gmx -f 2kyv_DMPC_cg.pdb -o 2kyv_DMPC_cg.gro

When prompted, choose *SIRAH force field* and then *SIRAH solvent models*.

.. note:: 

  By default, charged terminal are used but it is possible to set them neutral with option ``-ter``

.. note::

  Warning messages about long, triangular or square bonds are fine and expected due to the CG topology of some residues.

.. caution::

  However, missing atom messages are errors which probably trace back to the
  mapping step. In that case, check your atomistic and mapped structures and do not carry on the
  simulation until the problem is solved.

6.3. Solvate the system
_______________________


Define the simulation box of the system

.. code-block:: bash 
  
  gmx editconf -f 2kyv_DMPC_cg.gro -o 2kyv_DMPC_cg_box.gro -box 12 12 12 -c


Add WT4 molecules:

.. code-block:: bash 

  gmx solvate -cp 2kyv_DMPC_cg_box.gro -cs sirah.ff/wt416.gro -o 2kyv_DMPC_cg_solv1.gro

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
       | ; Compound #mols 
       | Protein_chain_A    1    
       | Protein_chain_B    1 
       | Protein_chain_C    1
       | Protein_chain_D    1
       | Protein_chain_E    1
       | Lipid_chain_F      1
       | 
    
     - | [ molecules ] 
       | ; Compound #mols
       | Protein_chain_A    1    
       | Protein_chain_B    1 
       | Protein_chain_C    1
       | Protein_chain_D    1
       | Protein_chain_E    1
       | Lipid_chain_F      1 
       | WT4             3875

.. hint::
  
  If you forget to read the number of added WT4 molecules from the output of *solvate*, then use the following command line to get it 

  .. code-block:: console

    grep -c WP1 2kyv_DMPC_cg_solv1.gro

.. caution::
  
  The number of added WT4 molecules, **3875**, may change according to the software version.

Remove misplaced WT4 molecules inside the bilayer:

.. code-block:: bash
  
  gmx grompp -f sirah.ff/tutorial/6/GPU/em1_CGLIPROT.mdp -p topol.top -c 2kyv_DMPC_cg_solv1.gro -o 2kyv_DMPC_cg_solv1.tpr

.. code-block:: bash
  
  echo "q" | gmx make_ndx -f 2kyv_DMPC_cg_solv1.gro -o 2kyv_DMPC_cg_solv1.ndx

.. code-block:: bash
  
  gmx select -f 2kyv_DMPC_cg_solv1.gro -s 2kyv_DMPC_cg_solv1.tpr -n 2kyv_DMPC_cg_solv1.ndx -on rm_close_wt4.ndx -sf sirah.ff/tutorial/6/rm_close_wt4.dat

.. code-block:: bash
  
  gmx editconf -f 2kyv_DMPC_cg_solv1.gro -o 2kyv_DMPC_cg_solv2.gro -n rm_close_wt4.ndx

.. caution::
  
  New GROMACS versions may complain about the non-neutral charge of the system, aborting the generation of the TPR file by command grompp. We will neutralize the system later, so to overcame this issue, just allow a warning message by adding the following keyword to the grompp command line: ``-maxwarn 1``

.. note::
  
  Consult ``sirah.ff/0ISSUES`` and :doc:`FAQs <../FAQ>` for information on known solvation issues.

Edit the [ molecules ] section in ``topol.top`` to correct the number of WT4 molecules:

.. hint::
  
  If you forget to read the number of added WT4 molecules from the output of *solvate*, then use the following command line to get it 

  .. code-block:: console

    grep -c WP1 2kyv_DMPC_cg_solv2.gro

Add CG counterions and 0.15M NaCl:

.. code-block:: bash

  gmx grompp -f sirah.ff/tutorial/6/GPU/em1_CGLIPROT.mdp -p topol.top -c 2kyv_DMPC_cg_solv2.gro -o 2kyv_DMPC_cg_solv2.tpr

.. code-block:: bash

  gmx genion -s 2kyv_DMPC_cg_solv2.tpr -o 2kyv_DMPC_cg_ion.gro -np 95 -pname NaW -nn 110 -nname ClW


When prompted, choose to substitute *WT4* molecules by *ions*.

.. note:: 

  The available electrolyte species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW) which represent solvated ions in solution. One ion pair (e.g., NaW-ClW) each 34 WT4 molecules results in a salt concentration of ~0.15M (see :ref:`Appendix <Appendix>` for details). Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`_.

Edit the [ molecules ] section in ``topol.top`` to include the CG ions and the correct number of WT4.

Before running the simulation it may be a good idea to visualize your molecular system. CG molecules are not recognized by molecular visualizers and will not display correctly. To fix this problem you may
generate a PSF file of the system using the script ``g_top2psf.pl``:

.. code-block:: bash

  ./sirah.ff/tools/g_top2psf.pl -i topol.top -o 2kyv_DMPC_cg_ion.psf

.. note::

  This is the basic usage of the script ``g_top2psf.pl``, you can learn other capabilities from its help:
  
  .. code-block:: bash

    ./sirah.ff/tools/g_top2psf.pl -h


Use VMD to check how the CG system looks like:

.. code-block::

  vmd 2kyv_DMPC_cg_ion.psf 2kyv_DMPC_cg_ion.gro -e sirah.ff/tools/sirah_vmdtk.tcl

.. tip::

  VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right ones, according to the CG representation. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
  Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.

6.4. Generate position restraint files
_______________________________________

To achive a proper interaction between the protein and bilayer, we will perform a equilibration step applying restraints over the protein backbone and lipids' phosphate groups.

.. important:: 

  GROMACS enumerates atoms starting from 1 in each topology file (e.i. ``topol_*.itp``). In order to avoid atom index discordance between the .gro file of the system and each topology, we will need independent .gro files for the PLN monomer and the DMPC bilayer.

Crate an index file of the system with a group for PLN monomer, then generate the .gro file of it.

.. note::

  WT4 and CG ions (NaW, ClW) are automatically set to the group *SIRAH-Solvent* while DMPC (named CMM at CG level) is assigned to group *Lipid*.

.. code-block:: bash
  
  echo -e "ri 1-52\nq" | gmx make_ndx -f 2kyv_DMPC_cg_ion.gro -o 2kyv_DMPC_cg_ion.ndx

.. code-block:: bash
  
  gmx editconf -f 2kyv_DMPC_cg_ion.gro -n 2kyv_DMPC_cg_ion.ndx -o 2kyv_DMPC_cg_monomer.gro

When prompted, choose *r_1-52*.

Create an index file for the monomer and add a group for the backbone *GO* and *GN* beads.

.. code-block:: bash
  
  echo -e "a GO GN \nq" | gmx make_ndx -f 2kyv_DMPC_cg_monomer.gro -o 2kyv_DMPC_cg_monomer.ndx

Create a position restraint file for the monomer backbone.

.. code-block:: bash
  
  gmx genrestr -f 2kyv_DMPC_cg_monomer.gro -n 2kyv_DMPC_cg_monomer.ndx -fc 100 100 100 -o posre_BB.itp

When prompted, choose *GO_GN*.

Edit each ``topol_Protein_chain_*.itp`` (A to E) to include the new position restraints:

.. list-table::
   :align: center
   :widths: 50 50
   :header-rows: 1

   * - Topology before editing
     - Topology after editing
   * - | ; Include Position restraint file 
       | #ifdef POSRES    
       | #include \"posre_Protein.itp\" 
       | #endif
       |   
       |    
       |   
       |   
      
     - | ; Include Position restraint file 
       | #ifdef POSRES    
       | #include \"posre_Protein.itp\" 
       | #endif
       |
       | #ifdef POSREBB    
       | #include \"posre_BB.itp\" 
       | #endif

Use a similar procedure to set the positional restraints on lipid's phosphates.

.. code-block:: bash
  
  gmx editconf -f 2kyv_DMPC_cg_ion.gro -n 2kyv_DMPC_cg_ion.ndx -o DMPC_cg.gro

When prompted, choose *Lipid*.

Create an index file of the membrane and add a group for phosphates (BFO beads).

.. code-block:: bash
  
  echo -e "a BFO \nq" | gmx make_ndx -f DMPC_cg.gro -o DMPC_cg.ndx

Create a position restraint file for phosphate groups in z coordinate.

.. code-block:: bash
  
  gmx genrestr -f DMPC_cg.gro -n DMPC_cg.ndx -fc 0 0 100 -o posre_Pz.itp

When prompted, choose *BFO*.

Edit ``topol_Lipid_chain_F.itp`` to include the new position restraints and define the flags *POSREZ* to switch on these restraints in the input file (.mdp).

.. list-table::
   :align: center
   :widths: 50 50
   :header-rows: 1

   * - Topology before editing
     - Topology after editing
   * - | ; Include Position restraint file 
       | #ifdef POSRES    
       | #include \"posre_Lipid_chain_F.itp\" 
       | #endif
       | 
       | 
       | 
       | 
              
     - | ; Include Position restraint file 
       | #ifdef POSRES    
       | #include \"posre_Lipid_chain_F.itp\" 
       | #endif
       |
       | #ifdef POSREZ    
       | #include \"posre_Pz.itp\" 
       | #endif

6.5. Run the simulation
________________________

.. important:: 

  By default, in this tutorial we will use input files for GROMACS on GPU (``sirah.ff/tutorial/6/GPU``). Example input files for using GROMACS on CPU can be found at: ``sirah.ff/tutorial/6/CPU``.

The folder ``sirah.ff/tutorial/6/GPU/`` contains typical input files for energy minimization
(``em1_CGLIPROT.mdp`` and ``em2_CGLIPROT.mdp``), equilibration (``eq1_CGLIPROT.mdp`` and ``eq2_CGLIPROT.mdp``) and production (``md_CGLIPROT.mdp``) runs. Please
check carefully the input flags therein.

Make a new folder for the run:

.. code-block:: bash

  mkdir run; cd run

**Energy Minimization of side chains by restraining the backbone**:

.. code-block:: bash

  gmx grompp -f ../sirah.ff/tutorial/6/GPU/em1_CGLIPROT.mdp -p ../topol.top -n ../2kyv_DMPC_cg_ion.ndx -c ../2kyv_DMPC_cg_ion.gro -r ../2kyv_DMPC_cg_ion.gro -o 2kyv_DMPC_cg_em1.tpr 

.. code-block:: bash

  gmx mdrun -deffnm 2kyv_DMPC_cg_em1 &> EM1.log &  

**Energy Minimization of the whole system**:

.. code-block:: bash

  gmx grompp -f ../sirah.ff/tutorial/6/GPU/em2_CGLIPROT.mdp -p ../topol.top -n ../2kyv_DMPC_cg_ion.ndx -c 2kyv_DMPC_cg_em1.gro -r 2kyv_DMPC_cg_em1.gro -o 2kyv_DMPC_cg_em2.tpr 

.. code-block:: bash

  gmx mdrun -deffnm 2kyv_DMPC_cg_em2 &> EM2.log & 

**Equilibration 1**:

Position restraints are defined in ``eq1_CGLIPROT.mdp`` file for protein backbone in xyz and phosphate groups (BFO beads) in z coordinate by setting keywords ``-DPOSREBB`` and ``-DPOSREZ``, respectively.

.. code-block:: bash

  gmx grompp -f ../sirah.ff/tutorial/6/GPU/eq1_CGLIPROT.mdp -p ../topol.top -n ../2kyv_DMPC_cg_ion.ndx -c 2kyv_DMPC_cg_em2.gro -r 2kyv_DMPC_cg_em2.gro -o 2kyv_DMPC_cg_eq1.tpr 

.. code-block:: bash

  gmx mdrun -deffnm 2kyv_DMPC_cg_eq1 &> EQ1.log & 

**Equilibration 2**:

.. code-block:: bash

  gmx grompp -f ../sirah.ff/tutorial/6/GPU/eq2_CGLIPROT.mdp -p ../topol.top -n ../2kyv_DMPC_cg_ion.ndx -c 2kyv_DMPC_cg_eq1.gro -r 2kyv_DMPC_cg_eq1.gro -o 2kyv_DMPC_cg_eq2.tpr 

.. code-block:: bash

  gmx mdrun -deffnm 2kyv_DMPC_cg_eq2 &> EQ2.log & 

**Production (1000ns)**:

.. code-block:: bash

  gmx grompp -f ../sirah.ff/tutorial/6/GPU/md_CGLIPROT.mdp -p ../topol.top -n ../2kyv_DMPC_cg_ion.ndx -c 2kyv_DMPC_cg_eq2.gro -r 2kyv_DMPC_cg_eq2.gro -o 2kyv_DMPC_cg_md.tpr 

.. code-block:: bash

  gmx mdrun -deffnm 2kyv_DMPC_cg_md &> MD.log & 

6.6. Visualizing the simulation
________________________________

That’s it! Now you can analyze the trajectory.

Process the output trajectory at folder ``run/`` to account for the Periodic Boundary Conditions (PBC):

.. code-block:: bash

  gmx trjconv -s 2kyv_DMPC_cg_em1.tpr -f 2kyv_DMPC_cg_md.xtc -o 2kyv_DMPC_cg_md_pbc.xtc -n ../2kyv_DMPC_cg_ion.ndx -pbc mol

When prompted, choose *System* for output.

Now you can check the simulation using VMD:

.. code-block:: bash

  vmd ../2kyv_DMPC_cg_ion.psf ../2kyv_DMPC_cg_ion.gro 2kyv_DMPC_cg_md_pbc.xtc -e ../sirah.ff/tools/sirah_vmdtk.tcl

.. note::
    
    The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. Use the command ``sirah-help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.
