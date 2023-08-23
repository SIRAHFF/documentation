This tutorial shows how to use the SIRAH force field to perform a coarse grained (CG) simulation of a protein embedded in a lipid bilayer using explicit solvent (called WatFour, WT4). The main references for
this tutorial are: `SIRAH 2.0 <https://doi.org/10.1021/acs.jctc.9b00006>`__, `SIRAH Lipids <https://doi.org/10.1021/acs.jctc.9b00435>`__, and `SIRAH Tools <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`__.
We strongly advise you to read these articles before starting the tutorial.

.. note::

	We strongly advise you to read and complete :ref:`Tutorial 5 <Tutorial 5>` and :ref:`Tutorial 6 <Tutorial 6>` before starting.
	

.. important::

    Check :ref:`download <download amber>` section for download and set up details before to start this tutorial.     
    Since this is **tutorial 7**, remember to replace ``X.X`` in your folder directory. The files corresponding to this tutorial can be found in: ``sirah_[version].amber/tutorial/7/``
	
7.1. Build CG representations
______________________________

.. caution::

	The mapping to CG requires the correct protonation state of each residue at a given pH. See :doc:`FAQs <../FAQ>` and :ref:`Tutorial 5 <Tutorial 5>` for cautions while preparing and mapping atomistic proteins to SIRAH.
	
Map the atomistic structure of protein `2KYV <https://www.rcsb.org/structure/2KYV>`__ to its CG representation:  

.. code-block:: bash

  ./sirah.amber/tools/CGCONV/cgconv.pl -i -i sirah.amber/tutorial/7/2kyv.pqr -o 2kyv_cg.pdb 
  
The input file ``-i`` 2kyv.pqr contains the atomistic representation of `2KYV <https://www.rcsb.org/structure/2KYV>`__ structure at pH 7.0, while the output ``-o`` 2kyv_cg.pdb is its SIRAH CG representation. 

.. important::

	If you already have an atomistic protein within a membrane, then you can simply map the entire system to SIRAH (this is highly recommended) and skip the step 2, however clipping the membrane patch may be required to set a correct solvation box (see bellow). By default no mapping is applied to lipids, as there is no standard naming convention for them. So users are requested to append a MAP file from the list in :ref:`Table 1 <table>`, by setting the flag ``-a`` in **cgconv.pl**. We recommend using `PACKMOL <https://m3g.github.io/packmol/>`__ for building the system. Reference building-block structures are provided at folder sirah.amber/PDB/, which agree with the mapping scheme ``sirah.amber/tools/CGCONV/maps/tieleman_lipid.map``. See :doc:`FAQs <../FAQ>` for cautions on mapping lipids to SIRAH and tips on using fragment-based topologies.  

.. tip::

  This an advanced usage of the script **cgconv.pl**, you can learn other capabilities from its help:
  ``./sirah.amber/tools/CGCONV/cgconv.pl -h``

7.1.1. Embed the protein in a lipid bilayer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We will show one possible way to do it by starting from a pre-equilibrated CG membrane patch.

In ``sirah_[version].amber/tutorial/7/`` you will find a pre-stabilized CG DMPC bilayer patch, concatenate it with the previously generated CG representation of the phospholamban (PLN) pentamer (PDB code: `2KYV <https://www.rcsb.org/structure/2KYV>`__):

.. code-block:: bash

	head -qn -1 2kyv_cg.pdb ./sirah.amber/tutorial/7/DMPC_cg.pdb > 2kyv_DMPC_cg_init.pdb

Luckily, we already oriented the protein inside the membrane. For setting up your own system you can go to `Orientations of Proteins in Membranes (OPM) database <https://opm.phar.umich.edu/>`__ and (if your structure is available) use the dummy atoms provided there to make them match with your membrane model (see **Figure 1**).

.. figure:: /../images/Tuto7.png
   :align: center
   :width: 100%

   **Figure 1.** Protein oriented inside the membrane from the OPM database with dummy atoms represented as orange spheres.
   
   

7.1.2. Delete close contact lipid molecules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We need to use VMD to to delete lipid molecules in close contact with the protein. For a proper treatment and visualization of the system in VMD you must first generate the molecular topology and initial coordinate files.
  
Use a text editor to create the file ``geninit.leap`` including the following lines:

.. code-block:: console

    # Load SIRAH force field
    addPath ./sirah.amber
    source leaprc.sirah

    # Load model
    ProtMem = loadpdb 2kyv_DMPC_cg_init.pdb

    # Save Parms
    saveAmberParmNetcdf ProtMem 2kyv_DMPC_cg_init.prmtop 2kyv_DMPC_cg_init.ncrst

    # EXIT
    quit

Run the LEaP application to generate the molecular topology and initial coordinate files:

.. code-block:: bash

    tleap -f geninit.leap

.. caution::

    Warning messages about long, triangular or square bonds in ``leap.log`` file are fine and expected due to the CG topology.

Now, open the files on VMD:

.. code-block:: bash

	vmd 2kyv_DMPC_cg_init.prmtop 2kyv_DMPC_cg_init.ncrst -e sirah.amber/tools/sirah_vmdtk.tcl

.. tip::

    VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right
    ones. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
    Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages.

In the VMD main window, select *Graphics* > *Representations*. In the *Selected Atoms* text entry type:

.. code-block:: text

	not (same residue as (sirah_membrane within 3.5 of sirah_protein) or (sirah_membrane and x < 5 or x > 142 or y < 2 or y > 140))

.. important::

	In the first part of the selection, lipid molecules in close contact with the protein are removed. The second one is made to “trim” the membrane patch, deleting lipids with acyl chains located in the periodic boundary images. This is frequent when using pre-equilibrated membrane patches and is necessary to avoid clashes in the following steps.

To save the refined protein-membrane system, in the VMD main window click on ``2kyv_DMPC_cg_init.prmtop``, then select *File* > *Save Coordinates*. In the *Selected atoms* option choose the selection you have just created and Save as ``2kyv_DMPC_cg.pdb``.
	
From now on it is just normal AMBER stuff!

7.2. Prepare LEaP input
________________________

Use a text editor to create the file ``gensystem.leap`` including the following lines:

.. code-block:: console

    # Load SIRAH force field
    addPath ./sirah.amber
    source leaprc.sirah

    # Load model
    Lipid = loadpdb DMPC64_cg.pdb

    # Add solvent, counterions and 0.15M NaCl
    # Tuned solute-solvent closeness for best hydration
    solvateBox Lipid WT4BOX {0 0 40} 0.7
    addIonsRand Lipid NaW 33 ClW 33

    # Save Parms
    saveAmberParmNetcdf protein DMPC64_cg.prmtop DMPC64_cg.ncrst

    # EXIT
    quit

.. seealso::

   The available ionic species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW). One ion pair (e.g. NaW-ClW) each 34 WT4 molecules renders a salt concentration of ~0.15M (see :ref:`Appendix <Appendix>` for details). Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`_.

7.3. Run LEaP
_______________

Run the LEaP application to generate the molecular topology and initial coordinate files:

.. code-block:: bash

    tleap -f gensystem.leap

.. caution::

    Warning messages about long, triangular or square bonds in ``leap.log`` file are fine and expected due to the CG topology.

This should create a topology file ``2kyv_DMPC_cg.prmtop`` and a coordinate file ``2kyv_DMPC_cg.ncrst``.

Use VMD to check how the CG model looks:

.. code-block:: bash

  vmd 2kyv_DMPC_cg.prmtop 2kyv_DMPC_cg.ncrst -e ./sirah.amber/tools/sirah_vmdtk.tcl

By selecting +X, +Y and +Z periodic images from the *Periodic* tab in the *Graphical Representations* window you will see unwanted water near the
hydrophobic region of the membrane and small vacuum slices at box boundaries. In the following step we will fix these issues by deleting those water molecules and reducing the box dimensions a few angstroms. See :doc:`FAQs <../FAQ>` for issues on membrane systems in AMBER. If you do not find your issue please start a discussion in our `github discussion page F&Q <https://github.com/SIRAHFF/documentation/discussions>`_.

.. tip::

    VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right
    ones. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
    Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages.

7.4. Resize the box with CPPTRAJ
_________________________________

.. note::

	As PACKMOL does not consider periodicity while building up the system, increasing the XY box sides a few Angstroms may be required to avoid bad contacts between images.
		
Use a text editor to create the file ``resize_box.cpptraj`` including the following lines:

.. code-block:: console

    # Set reference coordinate file
    reference 2kyv_DMPC_cg.ncrst
     
    # Remove water using distance-based masks
    strip :WT4&(@BCT1,BCT2,BC13,BC23<:12.0) parmout 2kyv_DMPC_cg.prmtop

    # New box dimensions
    box x 128 y 128 z 114

    # Amber NetCDF Restart generation
    trajout 2kyv_DMPC_cg_nb.ncrst

    # Do it!
     go

    # EXIT
     quit

.. caution::
	This is a critical step when preparing membrane systems to simulate with AMBER. In this case, the new box dimensions were set after some trial and error tests to allow for limited overlap between periodic box images. An excessive overlap may lead to important atom clashes an eventual system explosion during minimization/simulation, while insufficient overlap may impact the membrane cohesivity at PBC boundaries leading to pore formations or other issues.
	
Run the CPPTRAJ application application to adjust the size of the simulation box:

.. code-block:: bash

    cpptraj -p 2kyv_DMPC_cg.prmtop -y 2kyv_DMPC_cg.ncrst -i resize_box.cpptraj

Once again, use VMD to check the PBC images in the new box of the system:

.. code-block:: bash

  vmd 2kyv_DMPC_cg.prmtop 2kyv_DMPC_cg_nb.ncrst -e ./sirah.amber/tools/sirah_vmdtk.tcl
  
	
7.5. Run the simulation
________________________

Make a new folder for the run:

.. code-block:: bash

    mkdir -p run; cd run

The folder ``sirah.amber/tutorial/7/`` contains typical input files for energy minimization
(``em1_Prot-Lip.in`` and ``em2_Prot-Lip.in``), heating (``heat_Prot-Lip.in``), equilibration (``eq_Prot-Lip.in``) and production (``md_Prot-Lip.in``) runs. Please check carefully the input flags.

.. tip::

    **Some flags used in AMBER**

   - ``-i``: Input file.
   - ``-o``: Output file.
   - ``-p``: Parameter/topology file.
   - ``-c``: Coordinate file.
   - ``-r``: Restart file.
   - ``-x``: Trajectory file.
   - ``ref``: Reference file


.. warning::

    These input files are executed by the **GPU** implementation of *pmemd.cuda*, due to the system size we do not recommend the use of **CPU** implementations of AMBER.

.. note::

    Other available implementations that could be used: ``sander``  or ``pmemd``, both **CPU** implementations of AMBER. 

	
**Energy Minimization of side chains by restraining the backbone:**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/7/em1_Prot-Lip.in -p ../2kyv_DMPC_cg.prmtop -c ../2kyv_DMPC_cg_nb.ncrst -ref ../2kyv_DMPC_cg_nb.ncrst -o 2kyv_DMPC_cg_em_1.out -r 2kyv_DMPC_cg_em_1.ncrst &

**Energy Minimization of of the whole system:**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/7/em2_Prot-Lip.in -p ../2kyv_DMPC_cg.prmtop -c 2kyv_DMPC_cg_em_1.ncrst -ref 2kyv_DMPC_cg_em_1.ncrst -o 2kyv_DMPC_cg_em_2.out -r 2kyv_DMPC_cg_em_2.ncrst &
 
**Heating:**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/7/heat_Prot-Lip.in -p ../2kyv_DMPC_cg.prmtop -c 2kyv_DMPC_cg_em_2.ncrst -ref 2kyv_DMPC_cg_em_2.ncrst -o 2kyv_DMPC_cg_eq_0.out -r 2kyv_DMPC_cg_eq_0.ncrst -x 2kyv_DMPC_cg_eq_0.nc &

.. important::

	To avoid “*skinnb errors*” on GPU due to large box size fluctuations, the system must be equilibrated by several “short” runs using a large *skinnb* value. The number and length of the runs may vary according to the characteristic stabilization times of the system. For more information visit the `AMBER tutorial on lipids <http://ambermd.org/tutorials/advanced/tutorial16/>`__.
	
**Periodic box equilibration in GPU code (500 ps x 9):**

.. code-block:: bash

	for i in $(seq 1 9)
	do
		echo "running equilibration $i"
		pmemd.cuda -O \
		-i ../sirah.amber/tutorial/7/heat_Prot-Lip.in \
		-p ../2kyv_DMPC_cg.prmtop \
		-c 2kyv_DMPC_cg_eq_$(($i -1)).ncrst \
		-ref 2kyv_DMPC_cg_eq_$(($i -1)).ncrst \
		-o 2kyv_DMPC_cg_eq_$i.out \
		-r 2kyv_DMPC_cg_eq_$i.ncrst \
		-x 2kyv_DMPC_cg_eq_$i.nc
	done &
  
**Production (1000ns):**

.. code-block:: bash

   pmemd.cuda -O -i ../sirah.amber/tutorial/7/md_Prot-Lip.in -p ../2kyv_DMPC_cg.prmtop -c 2kyv_DMPC_cg_eq_9.ncrst -o 2kyv_DMPC_cg_md.out -r 2kyv_DMPC_cg_md.ncrst -x 2kyv_DMPC_cg_md.nc &


7.6. Visualizing the simulation
_______________________________

That’s it! Now you can analyze the trajectory.


Process the output trajectory to account for the Periodic Boundary Conditions (PBC):

.. code-block:: bash

    echo -e "autoimage\ngo\nquit\n" | cpptraj -p ../2kyv_DMPC_cg.prmtop -y 2kyv_DMPC_cg_md.nc -x 2kyv_DMPC_cg_md_pbc.nc --interactive


Now you can check the simulation using VMD:

.. code-block::

    vmd ../2kyv_DMPC_cg.prmtop 2kyv_DMPC_cg_md_pbc.nc -e ../sirah.amber/tools/sirah_vmdtk.tcl

.. note::

    The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD.
	
