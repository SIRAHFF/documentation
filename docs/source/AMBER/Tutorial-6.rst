.. note::

   Please report bugs, errors or enhancement requests through `Issue Tracker <https://github.com/SIRAHFF/documentation/issues>`_ or if you have a question about SIRAH open a `New Discussion <https://github.com/SIRAHFF/documentation/discussions>`_.
   
This tutorial shows how to use the SIRAH force field to perform a coarse grained (CG) simulation of a
DMPC bilayer in explicit solvent (called WatFour, WT4). The main references for
this tutorial are: `Barrera et al. <https://doi.org/10.1021/acs.jctc.9b00435>`_ and `Machado & Pantano <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_.
We strongly advise you to read these articles before starting the tutorial. You may also find interesting `this book chapter <https://pubs.aip.org/books/monograph/137/chapter-abstract/58880922/Simulating-Transmembrane-Proteins-with-the-Coarse?redirectedFrom=fulltext>`_.

.. note::

	If you are not familiar with membrane stuff we strongly recommend you to first perform the `Amber tutorial on lipids <http://ambermd.org/tutorials/advanced/tutorial16/>`__.
	

.. important::

    Check the :ref:`Setting up SIRAH <download amber>` section for download and set up details before starting this tutorial.
    Since this is **Tutorial 6**, remember to replace ``X.X`` and the files corresponding to this tutorial can be found in: ``sirah_[version].amber/tutorial/6/``
	
6.1. Build CG representations
______________________________

Map the atomistic structure of the preassembled DMPC bilayer to its CG representation:  

.. code-block:: bash

  ./sirah.amber/tools/CGCONV/cgconv.pl -i ./sirah.amber/tutorial/6/DMPC64.pdb -o DMPC64_cg.pdb -a sirah.amber/tools/CGCONV/maps/tieleman_lipid.map  
  
The input file ``-i`` DMPC64.pdb contains the atomistic representation of the DMPC bilayer, while the output ``-o`` DMPC64_cg.pdb is its SIRAH CG representation. The flag ``-a`` is a mapping scheme for lipids.

.. important::

	By default, no mapping is applied to lipids, as there is no standard naming convention for them. So users are requested to append a MAP file from the list in :ref:`Table 1 <table>`, by setting the flag ``-a`` in ``cgconv.pl``. We recommend using `PACKMOL <https://m3g.github.io/packmol/>`_ for building the system. Reference building-block structures are provided at folder ``sirah.amber/PDB/``, which agree with the mapping scheme in ``sirah.amber/tools/CGCONV/maps/tieleman_lipid.map``. The provided DMPC bilayer contains 64 lipid molecules per leaflet distributed in a 6.4 \* 6.4 nm surface, taking into account an approximate area per lipid of 0.64 nm\ :sup:`2` \ at 333 K. The starting configuration was created with the input file ``sirah.amber/tutorial/6/DMPC_bilayer.pkm``. See :doc:`FAQs <../FAQ>` for cautions on mapping lipids to SIRAH and tips on using fragment-based topologies. If you do not find your issue please start a discussion in our `github discussion page F&Q <https://github.com/SIRAHFF/documentation/discussions>`_.   

.. tip::

  This an advanced usage of the script **cgconv.pl**, you can learn other capabilities from its help by typing:

  .. code-block:: bash

    ./sirah.amber/tools/CGCONV/cgconv.pl -h

The input file ``DMPC64.pdb`` contains all the heavy atoms composing the lipids, while the output ``DMPC64_cg.pdb`` preserves a few of them. Please check both PDB structures using VMD:	

.. code-block:: bash

  vmd -m sirah.amber/tutorial/6/DMPC64.pdb DMPC64_cg.pdb


From now on it is just normal Amber stuff!

6.2. Prepare LEaP input
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

   The available electrolyte species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW) which represent solvated ions in solution. One ion pair (e.g., NaW-ClW) each 34 WT4 molecules results in a salt concentration of ~0.15M (see :ref:`Appendix <Appendix>` for details). Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`_.

6.3. Run LEaP
_______________

Run the LEaP application to generate the molecular topology and initial coordinate files:

.. code-block:: bash

    tleap -f gensystem.leap

.. note::

    Warning messages about long, triangular or square bonds in ``leap.log`` file are fine and expected due to the CG topology.

This should create a topology file ``DMPC64_cg.prmtop`` and a coordinate file ``DMPC64_cg.ncrst``.

Use VMD to check how the CG model looks:

.. code-block:: bash

  vmd DMPC64_cg.prmtop DMPC64_cg.ncrst -e ./sirah.amber/tools/sirah_vmdtk.tcl

By selecting +X, +Y and +Z periodic images from the *Periodic* tab in the *Graphical Representations* window you will see small vacuum slices at box boundaries. In the following step we will fix this issue by reducing the box dimensions a few angstroms. See :doc:`FAQs <../FAQ>` for issues on membrane systems in Amber.

.. tip::

    VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right
    ones, according to the CG representation. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
    Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.

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

Run the CPPTRAJ application to adjust the size of the simulation box:

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
(``em_Lipid.in``), heating (``heat_Lipid.in``), equilibration (``eq_Lipid.in``) and production (``md_Lipid.in``) runs. Please check carefully the input flags.

.. tip::

    **Some commonly used flags in AMBER**

   - ``-i``: Input file.
   - ``-o``: Output file.
   - ``-p``: Parameter/topology file.
   - ``-c``: Coordinate file.
   - ``-r``: Restart file.
   - ``-x``: Trajectory file.
   - ``-ref``: Reference file


.. warning::

	These input files are executed by the **GPU** implementation of ``pmemd.cuda``. Other available modules are ``sander`` or ``pmemd``, which are both **CPU** implementations of Amber.
	
	However, this simulation is time consuming owing to the system’s size. A parallel or CUDA implementation of Amber is advised.

	
**Energy Minimization:**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/6/em_Lipid.in -p ../DMPC64_cg.prmtop -c ../DMPC64_cg_nb.ncrst -o DMPC64_cg_em.out -r DMPC64_cg_em.ncrst &
 
**Heating:**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/6/heat_Lipid.in -p ../DMPC64_cg.prmtop -c DMPC64_cg_em.ncrst -o DMPC64_cg_eq_0.out -r DMPC64_cg_eq_0.ncrst -x DMPC64_cg_eq_0.nc &

.. important::

	To avoid “*skinnb errors*” on GPU due to large box size fluctuations, the system must be equilibrated by several “short” runs using a large *skinnb* value. The number and length of the runs may vary according to the characteristic stabilization times of the system. For more information visit the `Amber tutorial on lipids <http://ambermd.org/tutorials/advanced/tutorial16/>`__.
	
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

   pmemd.cuda -O -i ../sirah.amber/tutorial/6/md_Lipid.in -p ../DMPC64_cg.prmtop -c DMPC64_cg_eq_9.ncrst -o DMPC64_cg_md.out -r DMPC64_cg_md.ncrst -x DMPC64_cg_md.nc &



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

    The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. Use the command ``sirah-help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.
	
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


