This tutorial shows how to perform a coarse grained (CG) simulation of a double stranded DNA using the Generalized Born model for implicit solvent (GB) and the SIRAH force field. The main references
for this tutorial are: `Dans et al <https://pubs.acs.org/doi/abs/10.1021/ct900653p>`_ (latest parameters are those reported `here <https://pubs.acs.org/doi/abs/10.1021/ct100379f>`_) and `Machado and Pantano <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_. We strongly advise you to read these articles before starting the tutorial.

.. important::

    Check :ref:`setting up SIRAH <download amber>` section for download and set up details before to start this tutorial.
    Since this is **tutorial 1**, remember to replace ``X.X``, the files corresponding to this tutorial can be found in: ``sirah_[version].amber/tutorial/1/``

.. note::

    If you are not familiar with DNA stuff we strongly recommend you to first perform the `AMBER
    tutorial on DNA <http://ambermd.org/tutorials/basic/tutorial1>`_.

1.1. Build CG representations
______________________________

Map the atomistic structure of a 20-mer DNA to its CG representation:

.. code-block:: bash

  ./sirah.amber/tools/CGCONV/cgconv.pl -i ./sirah.amber/tutorial/1/dna.pdb  -o dna_cg.pdb

The input file ``-i`` dna.pdb contains all the heavy atoms composing the DNA molecule, while the  output ``-o`` dna_cg.pdb preserves a few of them.
Please check both PDB structures using VMD:

.. code-block:: bash

    vmd -m ./sirah.amber/tutorial/1/dna.pdb dna_cg.pdb

.. tip::

    This is the basic usage of the script **cgconv.pl**, you can learn other capabilities from its help:
    ``./sirah.amber/tools/CGCONV/cgconv.pl -h``

From now on it is just normal AMBER stuff!

1.2. Prepare LEaP input
________________________

Use a text editor to create the file ``gensystem.leap`` including the following lines:

.. code-block:: console

    # Load SIRAH force field
    addPath ./sirah.amber
    source leaprc.sirah

    # Load model
    dna = loadpdb dna_cg.pdb

    # Save Parms
    saveAmberParmNetcdf dna dna_cg.prmtop dna_cg.ncrst

    # EXIT
    quit


1.3. Run LEaP
______________

Run the LEAP application to generate the molecular topology and initial coordinate files:

.. code-block:: bash

    tleap -f gensystem.leap

.. caution::

    Warning messages about long, triangular or square bonds in ``leap.log`` file are fine and
    expected due to the CG topology.


This should create a topology file dna_cg.prmtop and a coordinate file dna_cg.ncrst.

Use VMD to check how the CG model looks like:

.. code-block:: bash

    vmd dna_cg.prmtop dna_cg.ncrst -e ./sirah.amber/tools/sirah_vmdtk.tcl

.. tip::

    VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right
    ones, according to the CG representation. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
    Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages.


1.4. Run the simulation
________________________

Make a new folder for the run:

.. code-block:: bash

    mkdir -p run; cd run


The folder ``sirah.amber/tutorial/1/PMEMD.GPU/`` contains typical input files for energy minimization
(``em_GB.in``), equilibration (``eq_GB.in``) and production (``md_GB.in``) runs. Please check carefully the input
flags therein.

.. tip::

    **Some commonly used flags in AMBER**

   - ``-i``: Input file.
   - ``-o``: Output file.
   - ``-p``: Parameter/topology file.
   - ``-c``: Coordinate file.
   - ``-r``: Restart file.
   - ``-x``: Trajectory file.


.. caution::

    These input files are executed by the **GPU** implementation of ``pmemd.cuda``. Other available modules are ``sander``  or ``pmemd``, which are both **CPU** implementations of AMBER.

.. note::

    You can find example input files for CPU versions of sander and pmemd at folders ``SANDER/`` and  ``PMEMD.CPU/``, within ``sirah.amber/tutorial/1/``

**Energy Minimization:**

.. code-block:: bash

  pmemd.cuda -O -i ../sirah.amber/tutorial/1/PMEMD.GPU/em_GB.in -p ../dna_cg.prmtop -c ../dna_cg.ncrst -o dna_cg_em.out -r dna_cg_em.ncrst &

.. important::

    In the course of long MD simulations the capping residues may eventually separate, this effect is
    called helix fraying. To avoid such behavior is necessary to set Watson-Crick restraints for the capping base pairs of this CG DNA at the end of ``eq_GB.in`` and ``md_GB.in`` files. Check the files lines that start with *&rst*.


.. warning:: 

    If you are using SANDER, to avoid the helix frying, you must create a symbolic link to the file ``dna_cg.RST``, which
    contains the definition of Watson-Crick restraints for the capping base pairs of this CG DNA:


    .. code-block:: bash

        ln -s ../sirah.amber/tutorial/1/SANDER/dna_cg.RST

    
    The file dna_cg.RST can only be read by SANDER, PMEMD reads a different restrain format.


**Equilibration:**

.. code-block:: bash

        pmemd.cuda -O -i ../sirah.amber/tutorial/1/PMEMD.GPU/eq_GB.in -p ../dna_cg.prmtop -c dna_cg_em.ncrst -o dna_cg_eq.out -r dna_cg_eq.ncrst -x dna_cg_eq.nc &

**Production (100ns):**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/1/PMEMD.GPU/md_GB.in -p ../dna_cg.prmtop -c dna_cg_eq.ncrst -o dna_cg_md.out -r dna_cg_md.ncrst  -x dna_cg_md.nc &


1.5. Visualizing the simulation
________________________________

Now you can load, visualize and analize the trajectory in VMD:

.. code-block::

    vmd ../dna_cg.prmtop ../dna_cg.ncrst dna_cg_md.nc -e ../sirah.amber/tools/sirah_vmdtk.tcl

.. note::

    The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. Use the command ``sirah-help`` in the Tcl/Tk console of VMD to access the manual pages.

