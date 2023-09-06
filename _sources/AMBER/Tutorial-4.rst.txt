This tutorial shows how to use the SIRAH force field to perform a coarse grained (CG) simulation of a
closed circular DNA using the Generalized Born model (GB) for implicit solvent.

The main references for this tutorial are: `SIRAH DNA <https://pubs.acs.org/doi/abs/10.1021/ct900653p>`_ (latest parameters are those reported in: `WAT4 <https://pubs.acs.org/doi/abs/10.1021/ct100379f>`_) and `SIRAH Tools <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_. We strongly advise you to read these articles before starting the tutorial.

.. important::

    Check :ref:`download <download amber>` section for download and set up details before to start this tutorial.
    Since this is **tutorial 4**, remember to replace ``X.X``, the files corresponding to this tutorial can be found in: ``sirah_[version].amber/tutorial/4/``


4.1. Build CG representations
______________________________

Map the atomistic structure of the closed circular DNA to its CG representation:

.. code-block:: bash
   
   ./sirah.amber/tools/CGCONV/cgconv.pl -i ./sirah.amber/tutorial/4/ccdna.pdb | sed -e 's/DCX A 1/CW5 A 1/'  -e 's/DCX B 1/CW5 B 1/' -e 's/DGX A 100/GW3 A 100/' -e 's/DGX B 100/GW3 B 100/' > ccdna_cg.pdb

.. note::

	5' and 3' end residues mast be renamed to AW5, TW5, GW5 or CW5 and AW3, TW3, GW3 or CW3 to represent the corresponding Adenine, Thymine, Guanine or Cytosine extremes in a closed circular DNA.

.. tip::

    This is the basic usage of the script **cgconv.pl**, you can learn other capabilities from its help:
    ``./sirah.amber/tools/CGCONV/cgconv.pl -h``

The input file ``-i`` ccdna.pdb contains all the heavy atoms composing the DNA molecule, while the output ``-o`` ccdna_cg.pdb preserves a few of them.

Please check both PDB structures using VMD:

.. code-block:: bash

	vmd -m ./sirah.amber/tutorial/4/ccdna.pdb ccdna_cg.pdb

From now on it is just normal AMBER stuff!

4.2. Prepare LEaP input
________________________

Use a text editor to create the file ``gensystem.leap`` including the following lines:

.. code-block:: console

	# Load SIRAH force field
	addPath ./sirah.amber
	source leaprc.sirah

	# Load model
	dna = loadpdb ccdna_cg.pdb

	# Make a covalently closed circular DNA
	bond dna.1.PX dna.100.C1X
	bond dna.101.PX dna.200.C1X

	# Save Parms
	saveAmberParmNetcdf dna ccdna_cg.prmtop ccdna_cg.ncrst

	# EXIT
	quit

4.3. Run LEaP
______________

Run the LEAP application to generate the molecular topology and initial coordinate files:

.. code-block:: bash

    tleap -f gensystem.leap

.. note::

    Warning messages about long, triangular or square bonds in ``leap.log`` file are fine and
    expected due to the CG topology.

This should create a topology file ``ccdna_cg.prmtop`` and a coordinate file ``ccdna_cg.ncrst``.

Use VMD to check how the CG model looks like:

.. code-block::

	vmd ccdna_cg.prmtop ccdna_cg.ncrst -e ./sirah.amber/tools/sirah_vmdtk.tcl


.. tip::

    VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right
    ones. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
    Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages.

4.4. Run the simulation
________________________

Make a new folder for the run:

.. code-block:: bash

    mkdir -p run; cd run

The folder ``sirah.amber/tutorial/4/GPU`` contains typical input files for energy minimization (``em_GBSA.in``), equilibration (``eq_GBSA.in``) and production (``md_GBSA.in``) runs. Please check carefully the input flags therein.

.. tip::

    **Some flags used in AMBER**

   - ``-i``: Input file.
   - ``-o``: Output file.
   - ``-p``: Parameter/topology file.
   - ``-c``: Coordinate file.
   - ``-r``: Restart file.
   - ``-x``: Trajectory file.

.. caution::

    These input files are executed by the **GPU** implementation of ``pmemd.cuda``. Other available implementations that could be used: ``sander``  or ``pmemd``, both **CPU** implementations of AMBER.

.. note::

   You can find example input files for CPU and GPU, within ``sirah.amber/tutorial/4/``


**Energy Minimization:**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/4/GPU/em_GBSA.in -p ../ccdna_cg.prmtop -c ../ccdna_cg.ncrst -o ccdna_cg_em.out -r ccdna_cg_em.ncrst &

**Equilibration:**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/4/GPU/eq_GBSA.in -p ../ccdna_cg.prmtop -c ccdna_cg_em.ncrst -o ccdna_cg_eq.out -r ccdna_cg_eq.ncrst -x ccdna_cg_eq.nc &

**Production (100ns):**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/4/GPU/md_GBSA.in -p ../ccdna_cg.prmtop -c ccdna_cg_eq.ncrst -o ccdna_cg_md.out -r ccdna_cg_md.ncrst -x ccdna_cg_md.nc &


4.5. Visualizing the simulation
________________________________

That’s it! Now you can load, visualize and analize the trajectory file in VMD:

.. code-block::

	vmd ../ccdna_cg.prmtop ../ccdna_cg.ncrst ccdna_cg_md.nc -e ../sirah.amber/tools/sirah_vmdtk.tcl

.. note::

    The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. Use the command ``sirah-help`` in the Tcl/Tk console of VMD to access the manual pages.