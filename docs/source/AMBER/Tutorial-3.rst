This tutorial shows how to apply a multiscale representation of DNA, which is available in the SIRAH
force field, to study protein-DNA interactions. In this approach, the protein and its binding region are
treated with an atomistic force field, while the contextual DNA is represented at coarse-grained (CG)
level. The simulated system is composed of the human TATA binding protein (TBP PDB: `1C9B <https://www.rcsb.org/structure/1c9b>`_) bounded to the TATA box at the promoter region (-64 to +13) of the Human Immunodeficiency Virus type 1 (HIV-1, GenBank: `K03455 <https://www.ncbi.nlm.nih.gov/nuccore/K03455.1>`_, **Figure 1**). Solvent effects are considered implicitly using the Generalized Born model (GB).

The main references for this tutorial are: `All-atoms/CG DNA <https://pubs.rsc.org/en/Content/ArticleLanding/2011/CP/c1cp21248f>`_ and `SIRAH Tools <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_. We strongly advise you to read these articles before starting the tutorial.


.. figure:: /../images/Tuto3.png
   :alt: map to buried treasure

   **Figure 1.** Promoter region of HIV-1. The TATA box is highlighted in orange, while base pairs at 0.7nm of TBP are colored in yellow.


.. important::

    Check :ref:`download <download amber>` section for download and set up details before to start this tutorial.
    Since this is **tutorial 3**, remember to replace ``X.X``, the files corresponding to this tutorial can be found in: ``sirah_[version].amber/tutorial/3/``


3.1. Build CG representations
______________________________

Map only the atomistic structure of nucleotides outside the TATA box of the HIV-1 promoter to CG:

.. code-block:: bash
   
   ./sirah.amber/tools/CGCONV/cgconv.pl -i ./sirah.amber/tutorial/3/tbphiv.pdb -R 1-32,45-110,123-154 -o dna_hyb.pdb

	
The input file ``-i`` tbphiv.pdb contains all the heavy atoms composing the protein-DNA system. Mapped
residues are selected through option ``-R``. The selection considers a buffer of two base pairs at each
side of the TATA box (**Figure 1**). The atomistic/CG interface must always be a step of B form DNA. The
resulting coordinates are saved in the output ``-o`` dna_hyb.pdb.

Please check both PDB structures in VMD:

.. code-block:: bash

   vmd -m ./sirah.amber/tutorial/3/tbphiv.pdb dna_hyb.pdb

.. tip::

  This is the basic usage of the script **cgconv.pl**, you can learn other capabilities from its help:
  ``./sirah.amber/tools/CGCONV/cgconv.pl -h``

From now on it is just normal AMBER stuff!


3.2. Prepare LEaP
__________________

Use a text editor to create the file ``gensystem.leap`` including the following lines:

.. code-block:: console

   # Load AMBER force field (parm14SB/parmbsc0)
   source oldff/leaprc.ff14SB

   # Load SIRAH force field
   addPath ./sirah.amber
   source leaprc.sirah

   # Load model
   dna = loadpdb dna_hyb.pdb
   
   # Save Parms
   saveAmberParmNetcdf dna dna_hyb.prmtop dna_hyb.ncrst
   
   # EXIT
   quit

.. note::

   Notice: According to AMBER version 10, 11 or 12 (14) the source file for parm99SB/bsc0 is leaprc.ff99bsc0, leaprc.ff10 or leaprc.ff12SB respectively.

3.3. Run LEaP
______________

Run the LEaP application to generate the molecular topology and initial coordinate files:

.. code-block:: bash

    tleap -f gensystem.leap

.. warning::

    Warning messages about long, triangular or square bonds in ``leap.log`` file are fine and
    expected due to the CG topology of some residues.

This should create a topology file ``dna_hyb.prmtop`` and a coordinate file ``dna_hyb.ncrst``.

Use VMD to check how the multiscale model looks like:

.. code-block:: bash

   vmd dna_hyb.prmtop dna_hyb.ncrst -e ./sirah.amber/tools/sirah_vmdtk.tcl

.. tip::

    VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right
    ones. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
    Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages.

3.4. Run the simulation
________________________

Make a new folder for the run:

.. code-block:: bash

    mkdir -p run; cd run

The folder ``sirah.amber/tutorial/3/PMEMD.GPU/`` contains typical input files for energy minimization
(``em_HYB.in``), equilibration (``eq_HYB.in``) and production (``md_HYB.in``) runs. Please check carefully the input flags therein.

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

   You can find example input files for CPU versions of sander and pmemd at folders ``SANDER/`` and  ``PMEMD.CPU/``, within ``sirah.amber/tutorial/3/``

   This simulation is very time consuming owing to the system's size, so try a parallel or CUDA implementation of AMBER.

**Energy Minimization:**

.. code-block:: bash

   sander -O -i ../sirah.amber/tutorial/3/SANDER/em_HYB.in -p ../dna_hyb.prmtop -c ../dna_hyb.ncrst -o dna_hyb_em.out -r dna_hyb_em.ncrst &

.. important::

    In the course of long MD simulations the capping residues may eventually separate, this effect is
    called helix fraying. To avoid such behavior is necessary to set Watson-Crick restraints for the capping base pairs of this CG DNA at the end of ``eq_GB.in`` and ``md_GB.in`` files. Check the files lines that start with *&rst*.


.. warning:: 

    If you are using SANDER to avoid such behavior create a symbolic link to the file ``dna_cg.RST``, which
    contains the definition of Watson-Crick restraints for the capping base pairs of this CG DNA:


    .. code-block:: bash

        ln -s ../sirah.amber/tutorial/3/SANDER/dna_cg.RST

    
    The file dna_cg.RST can only be read by SANDER, PMEMD reads a different restrain format.


**Equilibration:**

.. code-block:: bash

   sander -O -i ../sirah.amber/tutorial/3/SANDER/eq_HYB.in -p ../dna_hyb.prmtop -c dna_hyb_em.ncrst -o dna_hyb_eq.out -r dna_hyb_eq.ncrst -x dna_hyb_eq.nc &

**Production (10ns):**

.. code-block:: bash

   sander -O -i ../sirah.amber/tutorial/3/SANDER/md_HYB.in -p ../dna_hyb.prmtop -c dna_hyb_eq.ncrst -o dna_hyb_md.out -r dna_hyb_md.ncrst -x dna_hyb_md.nc &

.. note::

    You can find example input files for CPU and GPU versions of pmemd at folders PMEMD.CPU/ and PMEMD.GPU/ within sirah.amber/tutorial/3/

3.5. Visualizing the simulation
________________________________

That’s it! Now you can load, visualize and analize the trajectory file in VMD:

.. code-block::

   vmd ../dna_hyb.prmtop ../dna_hyb.ncrst dna_hyb_md.nc -e ../sirah.amber/tools/sirah_vmdtk.tcl

.. note::

    The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. Use the command ``sirah-help`` in the Tcl/Tk console of VMD to access the manual pages.
