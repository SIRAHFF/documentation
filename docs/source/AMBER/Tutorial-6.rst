This tutorial shows how to use the SIRAH force field to perform a coarse grained (CG) simulation of a
DMPC bilayer in explicit solvent (called WatFour, WT4). The main references for
this tutorial are: `Barrera et al. SIRAH Lipids <https://doi.org/10.1021/acs.jctc.9b00435>`_, and `Machado et al. SIRAH Tools <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_.
We strongly advise you to read these articles before starting the tutorial.

.. note::

	If you are not familiar with membrane stuff we strongly recommend you to first perform the `AMBER tutorial on lipids <http://ambermd.org/tutorials/advanced/tutorial16/>`__.
	

.. important::

    Check :ref:`install <Download amber>` section for download and set up details before to start this tutorial.     
    Since this is **tutorial 6**, remember to replace ``X.X`` in your folder directory. The files corresponding to this tutorial can be found in: ``sirah_[version].amber/tutorial/6/``
	
6.1. Build CG representations
______________________________

Map the atomistic structure of the preassembled DMPC bilayer to its CG representation:  

.. code-block:: bash

  ./sirah.amber/tools/CGCONV/cgconv.pl -i ./sirah.amber/tutorial/6/DMPC64.pdb -o DMPC64_cg.pdb -a sirah.amber/tools/CGCONV/maps/tieleman_lipid.map  
  
The input file ``-i`` DMPC64.pdb contains the atomistic representation of the DMPC bilayer, while the output ``-o`` DMPC64_cg.pdb is its SIRAH CG representation. The flag ``-a`` is a mapping scheme for lipids.

.. important::

	By default no mapping is applied to lipids, as there is no standard naming convention for them. So users are requested to append a MAP file from the list in :ref:`Table 1 <table>`, by setting the flag ``-a`` in **cgconv.pl**. We recommend using `PACKMOL <https://m3g.github.io/packmol/>`_ for building the system. Reference building-block structures are provided at folder ``sirah.amber/PDB/``, which agree with the mapping scheme in ``sirah.amber/tools/CGCONV/maps/tieleman_lipid.map``. The provided DMPC bilayer contains 64 lipid molecules per leaflet distributed in a 6.4 \* 6.4 nm surface, taking into account an approximate area per lipid of 0.64 nm\ :sup:`2` \ at 333 K  . The starting configuration was created with the input file ``sirah.amber/tutorial/6/DMPC_bilayer.pkm``. See :doc:`FAQ` for cautions on mapping lipids to SIRAH and tips on using fragment-based topologies.   

.. tip::

  This an advanced usage of the script **cgconv.pl**, you can learn other capabilities from its help:
  ``./sirah.amber/tools/CGCONV/cgconv.pl -h``

The input file ``DMPC64.pdb`` contains all the heavy atoms composing the lipids, while the output ``DMPC64_cg.pdb`` preserves a few of them. Please check both PDB and PQR structures using VMD:	

.. code-block:: bash

  vmd -m sirah.amber/tutorial/6/DMPC64.pdb DMPC64_cg.pdb


From now on it is just normal AMBER stuff!

6.2. Prepare LEaP
___________________

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

   The available ionic species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW). One ion pair (e.g. NaW-ClW) each 34 WT4 molecules renders a salt concentration of ~0.15M (see Appendix 1). Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`_.

6.3. Run LEaP
_______________

Run the LEaP application to generate the molecular topology and initial coordinate files:

.. code-block:: bash

    tleap -f gensystem.leap

.. caution::

    Warning messages about long, triangular or square bonds in ``leap.log`` file are fine and expected due to the CG topology.

This should create a topology file ``DMPC64_cg.prmtop`` and a coordinate file ``DMPC64_cg.ncrst``.

Use VMD to check how the CG model looks:

.. code-block:: bash

  vmd DMPC64_cg.prmtop DMPC64_cg.ncrst -e ./sirah.amber/tools/sirah_vmdtk.tcl

By selecting +X, +Y and +Z periodic images from the Periodic tab in the Graphical Representations window you will see small vacuum slices at box boundaries. In the following step we will fix this issue by reducing the box dimensions a few angstroms. See :doc:`FAQ` for issues on membrane systems in AMBER.

.. tip::

    VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right
    ones. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
    Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages.

6.4. Resize the box with CPPTRAJ
_________________________________

.. note::

	As PACKMOL does not consider periodicity while building up the system, increasing the XY box sides a few Angstroms may be required to avoid bad contacts between images.
		
Use a text editor to create the file ``resize_box.cpptraj`` including the following lines:

.. code-block:: console

    # New box dimensions
    box x 66 y 66 z 132
 
    # Amber NetCDF Restart generation
    trajout DMPC64_cg_nb.ncrst

    # Do it!
    go

    # EXIT
    quit

Run the CPPTRAJ application to generate the molecular topology and initial coordinate files:

.. code-block:: bash

    cpptraj -p DMPC64_cg.prmtop -y DMPC64_cg.ncrst -i resize_box.cpptraj

Once again, use VMD to check the PBC images in the new box of the system:

.. code-block:: bash

  vmd DMPC64_cg.prmtop DMPC64_cg_nb.ncrst -e ./sirah.amber/tools/sirah_vmdtk.tcl
  
	
6.5. Run the simulation
________________________

Make a new folder for the run:

.. code-block:: bash

    mkdir -p run; cd run

The folder ``sirah.amber/tutorial/6/`` contains typical input files for energy minimization
(``em_Lipid.in``), heating (``heat_Lipid.in``), equilibration (``eq_Lipid.in``) and production (``md_Lipid.in``) runs. Please check carefully the
input flags.

.. tip::

    **Some flags used in AMBER**

   - ``sander``: The AMBER program for molecular dynamics simulations.
   - ``pmemd.cuda``: The GPU implementation of AMBER program's Particle Mesh Ewald Molecular Dynamics for simulations.
   - ``-i``: Input file.
   - ``-o``: Output file.
   - ``-p``: Parameter/topology file.
   - ``-c``: Coordinate file.
   - ``-r``: Restart file.
   - ``-x``: Trajectory file.
   - ``ref``: Reference file


.. caution::

	These input files are executed by the GPU implementation of *pmemd*.
	
**Energy Minimization:**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/6/em_Lipid.in.in -p ../DMPC64_cg.prmtop -c ../DMPC64_cg_nb.ncrst -o DMPC64_cg_em.out -r DMPC64_cg_em.ncrst &
 
**Heating:**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/6/heat_Lipid.in.in -p ../DMPC64_cg.prmtop -c DMPC64_cg_em.ncrst -o DMPC64_cg_eq_0.out -r DMPC64_cg_eq_0.ncrst -x DMPC64_cg_eq_0.nc &

.. important::

	To avoid “skinnb errors” on GPU due to large box size fluctuations, the system must be equilibrated by several “short” runs using a large skinnb value. The number and length of the runs may vary according to the characteristic stabilization times of the system. For more information visit the `AMBER tutorial on lipids <http://ambermd.org/tutorials/advanced/tutorial16/>`__.
	
**Periodic box equilibration in GPU code (500 ps x 9):**

.. code-block:: bash

	for i in $(seq 1 9)
	do
		echo "running equilibration $i"
		pmemd.cuda -O \
		-i ../sirah.amber/tutorial/6/eq_Lipid.in \
		-p ../DMPC64_cg.prmtop \
		-c DMPC64_cg_eq_$(($i -1)).ncrst \
		-o DMPC64_cg_eq_$i.out \
		-r DMPC64_cg_eq_$i.ncrst \
		-x DMPC64_cg_eq_$i.nc
	done &
  
**Production (1000ns):**

.. code-block:: bash

   pmemd.cuda -O -i ../sirah.amber/tutorial/6/md_Lipid.in.in -p ../DMPC64_cg.prmtop -c DMPC64_cg_eq_9.ncrst -o DMPC64_cg_md.out -r DMPC64_cg_md.ncrst -x DMPC64_cg_md.nc &



6.6. Visualizing the simulation
_______________________________

That’s it! Now you can analyze the trajectory.


Process the output trajectory to account for the Periodic Boundary Conditions (PBC):

  .. code-block:: bash

      echo -e "autoimage\ngo\nquit\n" | cpptraj -p ../DMPC64_cg.prmtop -y DMPC64_cg_md.nc -x DMPC64_cg_md_pbc.nc --interactive


Now you can check the simulation using VMD:

.. code-block::

    vmd ../DMPC64_cg.prmtop DMPC64_cg_md_pbc.nc -e ../sirah.amber/tools/sirah_vmdtk.tcl

.. note::

    The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD.
	
You can also use CPPTRAJ to calculate the area per lipid:

.. code-block:: bash

     cpptraj -p ../DMPC64_cg.prmtop -i ../sirah.amber/tutorial/6/area_lipid.cpptraj

Use Grace to plot the results:

.. code-block:: bash

     xmgrace apl_DMPC64_310K.dat

.. note::

    To calculate the area per lipid, divide the membrane's area by the DMPC molecules per leaflet:   
	
	.. math::
		\frac{Area}{Lipid} = \frac{Box(x) * Box(y)}{64} 


And density profiles and bilayer thickness: 

.. code-block:: bash

     cpptraj -p ../DMPC64_cg.prmtop -i ../sirah.amber/tutorial/6/dens_profile.cpptraj

Use Grace to plot the results:

.. code-block:: bash

     xmgrace -nxy dens_profile_DMPC64_310K.dat

.. note::

    The thickness of the bilayer is the distance between the two peaks corresponding to the position of phosphate beads (BFO) along the z-axis.


.. _table: 

Table 1. Available mapping files (MAPs) at folder ``sirah.amber/tools/CGCONV/maps/`` for converting atomistic lipid structures to SIRAH models. **Important!** MAPs can not inter-convert different name conventions, e.g. amber_lipid.map won’t generate fragment-based residues from residue-based force fields. Due to possible nomenclature conflicts, users are advised to check and modify the MAPs as required.

+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
|      **Map**            | **Type**\* | **Compatibility**                        | **Source**                                                     |
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| amber_lipid.map         |   F        | AMBER Lipid11-17 force fields            | | `AMBER <http://ambermd.org/>`__                              |
|                         |            |                                          | | `HTMD <https://software.acellera.com/htmd/index.html>`__     | 
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| GAFF_lipid.map          |   R        | AMBER GAFF force field                   | `LipidBook <https://lipidbook.org/>`__                         |
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| charmm_lipid.map        |   R        | | CHARMM 27/36 force field, and “CHARMM  | | `CHARMM-GUI <https://charmm-gui.org/>`__                     |
|                         |            | | compatible” GAFF nomenclature          | | `GROMACS <https://www.gromacs.org/>`__                       |  
|                         |            |                                          | | `LipidBook <https://lipidbook.org/>`__                       | 
|                         |            |                                          | | `MemBuilder <http://bioinf.modares.ac.ir/software/mb/>`__    |     
|                         |            |                                          | | `HTMD <https://software.acellera.com/htmd/index.html>`__     |
|                         |            |                                          | | `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_               |   
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| slipids.map             |   R        | Stockholm lipids force field             | | `SLIPIDS <http://www.fos.su.se/~sasha/SLipids/About.html>`__ |
|                         |            |                                          | | `MemBuilder <http://bioinf.modares.ac.ir/software/mb/>`__    |   
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| OPLSA-AA_2014_lipid.map |   R        | All-atoms lipids for OPLS force field    | | `Maciejewski <https://doi.org/10.1021/jp5016627>`__          |
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| OPLSA-UA_lipid.map      |   R        | United-atom lipids for OPLS force field  | `LipidBook <https://lipidbook.org/>`__                         |
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| GROMOS43a1_lipid.map    |   R        | | United-atom lipids for GROMOS 43a1 and | | `LipidBook <https://lipidbook.org/>`__                       |
|                         |            | | CKP force fields                       | | `MemBuilder <http://bioinf.modares.ac.ir/software/mb/>`__    |      
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| GROMOS43a1-s3_lipid.map |   R        | | United-atom lipids for GROMOS 43a1-s3  | | `GROMACS <https://www.gromacs.org/>`__                       |
|                         |            | | force field                            | | `LipidBook <https://lipidbook.org/>`__                       |
|                         |            |                                          | | `MemBuilder <http://bioinf.modares.ac.ir/software/mb/>`__    |
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| GROMOS53a6_lipid.map    |   R        | | United-atom lipids for GROMOS 53a6     | | `GROMACS <https://www.gromacs.org/>`__                       |
|                         |            | | force field                            | | `LipidBook <https://lipidbook.org/>`__                       | 
|                         |            |                                          | | `MemBuilder <http://bioinf.modares.ac.ir/software/mb/>`__    | 
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| tieleman_lipid.map      |   R        | | Berger lipids as implemented by        | | `Tieleman <https://doi.org/10.1021/ja0624321>`__             |
|                         |            | | Tieleman et al. for GROMOS             | | `LipidBook <https://lipidbook.org/>`__                       |  
|                         |            | | force fields.                          |                                                                |
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+

\* Fragment-based (F) or Residue-based (R) topology.