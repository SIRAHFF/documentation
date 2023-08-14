This tutorial shows how to perform a coarse grained (CG) simulation of a double stranded DNA using the Generalized Born model for implicit solvent (GB) and the SIRAH force field. The main references
for this tutorial are: `Dans et al. SIRAH DNA <https://pubs.acs.org/doi/abs/10.1021/ct900653p>`_ (latest parameters are those reported in: `Darré et al. WAT4?) <https://pubs.acs.org/doi/abs/10.1021/ct100379f>`_, `Machado et al. SIRAH Tools <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_. We strongly advise you to read these articles before starting the tutorial.

.. important::

    Check :ref:`install <Download amber>` section for download and set up details before to start this tutorial.
    **This is tutorial 1**, remember to replace ``X.X``, the files corresponding to this tutorial can be found in: ``sirah_[version].amber/tutorial/1/``

1. Build CG representations
____________________________

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

2. Prepare leap
_______________

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


3. Run LEAP
____________

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
    ones. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
    Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages.


4. Run the simulation
_______________________

Make a new folder for the run:

.. code-block:: bash

    mkdir -p run; cd run

In the course of long MD simulations the capping residues may eventually separate, this effect is
called helix fraying. To avoid such behavior create a symbolic link to the file ``dna_cg.RST``, which
contains the definition of Watson-Crick restraints for the capping base pairs of this CG DNA:

.. code-block:: bash

    ln -s ../sirah.amber/tutorial/1/SANDER/dna_cg.RST

.. important::

    The file dna_cg.RST can only be read by SANDER, PMEMD reads a different restrain format.

The folder ``sirah.amber/tutorial/1/SANDER/`` contains typical input files for energy minimization
(em_GB.in), equilibration (eq_GB.in) and production (md_GB.in) runs. Please check carefully the input
flags therein.

.. tip::

    **Some flags used in AMBER**

   - ``sander``: The AMBER program for molecular dynamics simulations.
   - ``-i``: Input file.
   - ``-o``: Output file.
   - ``-p``: Parameter/topology file.
   - ``-c``: Coordinate file.
   - ``-r``: Restart file.
   - ``-x``: Trajectory file.

**Energy Minimization:**

   .. code-block:: bash

      $ sander -O -i ../sirah.amber/tutorial/1/SANDER/em_GB.in -p ../dna_cg.prmtop -c ../dna_cg.ncrst -o dna_cg_em.out -r dna_cg_em.ncrst &

**Equilibration:**

    .. code-block:: bash

        $ sander -O -i ../sirah.amber/tutorial/1/SANDER/md_GB.in -p ../dna_cg.prmtop -c dna_cg_eq.ncrst -o dna_cg_md.out -r dna_cg_md.ncrst

**Production (100ns):**

    .. code-block:: bash

        sander -O -i ../sirah.amber/tutorial/1/SANDER/md_GB.in -p ../dna_cg.prmtop -c dna_cg_eq.ncrst -o dna_cg_md.out -r dna_cg_md.ncrst  -x dna_cg_md.nc &

.. note::

    You can find example input files for CPU and GPU versions of pmemd at folders PMEMD.CPU/ and PMEMD.GPU/ within sirah.amber/tutorial/1/

5. Visualising the simulation
______________________________

Now you can load, visualize and analize the trajectory file in VMD:

.. code-block::

    vmd ../dna_cg.prmtop ../dna_cg.ncrst dna_cg_md.nc -e ../sirah.amber/tools/sirah_vmdtk.tcl

.. hint::

    The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD.
