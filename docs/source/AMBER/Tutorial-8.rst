.. note::

   Please report bugs, errors or enhancement requests through `Issue Tracker <https://github.com/SIRAHFF/documentation/issues>`_ or if you have a question about SIRAH open a `New Discussion <https://github.com/SIRAHFF/documentation/discussions>`_.
   
This tutorial shows how to use the SIRAH force field to perform a coarse grained (CG) simulation of a
glycoprotein in explicit solvent (called WatFour, WT4). The main references for
this tutorial are: `Darré et al. <https://pubs.acs.org/doi/abs/10.1021/ct100379f>`_, `Machado et al. <https://doi.org/10.1021/acs.jctc.9b00006>`__ and `Machado & Pantano  <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_. We strongly advise you to read these articles before starting the tutorial.

.. important::

    Check the :ref:`Setting up SIRAH <download amber>` section for download and set up details before starting this tutorial.
    Since this is **Tutorial 8**, remember to replace ``X.X`` and the files corresponding to this tutorial can be found in: ``sirah_[version].amber/tutorial/8/``


8.1. Build CG representations
_____________________________

.. caution::

	The mapping to CG requires the correct protonation state of each residue at a given pH. We recommend using the `CHARMM-GUI server <https://www.charmm-gui.org/>`_ with the **Glycan Reader & Modeler** to prepare your system, choosing the output naming scheme of AMBER for best compatibility. An account is required to access any of the CHARMM-GUI Input Generator modules, and it can take up to 24 hours to obtain one. 

This example use the structure of glycoprotein `1GYA <https://www.rcsb.org/structure/1GYA>`_, using CHARMM-GUI server are obtained the files parm7 and rst7, these files are converted to pdb format with the naming scheme of AMBER/GLYCAM using AmberTool as follow:

.. code-block:: bash

	ambpdb -p ./sirah.amber/tutorial/8/step3_input.parm7 -c ./sirah.amber/tutorial/8/step3_input.rst7 > 1GYA_glycam.pdb

From the file 1GYA_glycam.pdb generated delete the solvent, rename to 1GYA_glycam_NoW.pdb and then map the protonated atomistic structure to its CG representation:   

.. code-block:: bash

	./sirah.amber/tools/CGCONV/cgconv.pl -i ./sirah.amber/tutorial/8/1GYA_glycam_NoW.pdb -o 1GYA_cg.pdb  
  
The input file ``-i`` 1GYA_glycam_NoW.pdb contains the atomistic representation of `1GYA <https://www.rcsb.org/structure/1GYA>`_ structure at pH **7.0**, while the output ``-o`` 1GYA_cg.pdb is its SIRAH CG representation.

.. tip::

	This is the basic usage of the script **cgconv.pl**, you can learn other capabilities from its help by typing:

	.. code-block:: bash

		./sirah.amber/tools/CGCONV/cgconv.pl -h	
		
.. note::

	**Pay attention to residue names when mapping structures from other atomistic force fields or experimental structures.** Although we provide compatibility for naming schemes in PDB, GMX, GROMOS, CHARMM and OPLS, there might always be some ambiguity in the residue naming, specially regarding protonation states, that may lead to a wrong mapping. For example, SIRAH Tools always maps the residue name “HIS” to a Histidine protonated at the epsilon nitrogen (:math:`N_{\epsilon}`) regardless the actual proton placement. Similarly, protonated Glutamic and Aspartic acid residues must be named “GLH” and “ASH”, otherwise they will be treated as negative charged residues. In addition, protonated and disulfide bonded Cysteines must be named “CYS” and “CYX” respectively. These kind of situations need to be carefully checked by the users. In all cases the residues preserve their identity when mapping and back-mapping the structures. Hence, the total charge of the protein should be the same at atomistic and SIRAH levels. You can check the following mapping file to be sure of the compatibility: ``sirah.amber/tools/CGCONV/maps/sirah_prot.map`` and ``sirah_glycans.map``.    

  
.. important::

	By default, charged termini are used. However, it is possible to set them neutral by renaming the residues from **s**\[code\] to **a**\[code\] (Nt-acetylated) or **m**\[code\] (Ct-amidated) after mapping to CG, where \[code\] is the root residue name in SIRAH. For example, to set a neutral N-terminal Histidine protonated at epsilon nitrogen (:math:`N_{\epsilon}`) rename it from “sHe” to “aHe”.


Please check both PDBs structures using VMD:	

.. code-block:: bash

  vmd -m sirah.amber/tutorial/8/1GYA_glycam_NoW.pdb 1GYA_cg.pdb


From now on it is just normal Amber stuff!


8.2. Prepare LEaP input
_________________________

Use a text editor to create the file ``gensystem.leap`` including the following lines:

.. code-block:: console

    # Load SIRAH force field
    addPath ./sirah.amber
    source leaprc.sirah

    # Load model
    glycoprot = loadpdb 1GYA_cg.pdb

    charge glycoprot

    # N-Glycosilation
    bond glycoprot.111.GO2  glycoprot.110.GO2
    bond glycoprot.110.GO2  glycoprot.109.GO6

    bond glycoprot.112.GO2  glycoprot.109.GO3
    bond glycoprot.109.GO2  glycoprot.108.GO6

    bond glycoprot.114.GO2 glycoprot.113.GO2
    bond glycoprot.113.GO2 glycoprot.108.GO3

    bond glycoprot.108.GO2  glycoprot.107.GO4
    bond glycoprot.107.GNac glycoprot.106.GO4
    # ASN 65
    bond glycoprot.106.GNac glycoprot.65.BND


    # Add solvent, counterions and 0.15M NaCl
    # Tuned solute-solvent closeness for best hydration
    solvateOct glycoprot WT4BOX 20 0.7
    addIonsRand glycoprot NaW 40 ClW 41

    # Save topology
    saveAmberParmNetcdf glycoprot 1GYA_cg.prmtop 1GYA_cg.ncrst

    # EXIT
    quit

.. caution::

    Each glycosidic bond must be defined explicitly in LEaP using the command bond, e.g.: “*bond unit.ri.beadi unit.rj.beadj*”. Where *ri* and *rj* correspond to the residue index in the topology file starting from 1, which may differ from the biological sequence in the PDB file. And *beadi* and *beadj* are the names of the beads involved in the glycosidic bond of the corresponding residues.

    The above also applies to each disulfide bond, e.g.: “*bond unit.ri.BSG unit.rj.BSG*”. You can try the command *pdb4amber* to get those indexes from the atomistic structure, but be aware that it may not work if the Cysteine residues are too far away (in this case result in an empty file):

    .. code-block:: bash

		pdb4amber -i sirah.amber/tutorial/8/1GYA_glycam_NoW.pdb -o 1GYA_aa.pdb && cat 1GYA_aa_sslink


	
.. seealso::

       The available electrolyte species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW) which represent solvated ions in solution. One ion pair (e.g., NaW-ClW) each 34 WT4 molecules results in a salt concentration of ~0.15M (see :ref:`Appendix <Appendix>` for details). Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`__.
	   

8.3. Run LEaP 
____________________

Run the LEaP application to generate the molecular topology and initial coordinate files:

.. code-block:: bash

    tleap -f gensystem.leap

.. note::

    Warning messages about long, triangular or square bonds in ``leap.log`` file are fine and expected due to the CG topology of some residues.


This should create a topology file ``1GYA_cg.prmtop`` and a coordinate file ``1GYA_cg.ncrst``.

Use VMD to check how the CG model looks like and particularly the presence of glycosidic bonds:

.. code-block:: bash

  vmd 1GYA_cg.prmtop 1GYA_cg.ncrst -e ./sirah.amber/tools/sirah_vmdtk.tcl


.. tip::

    VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right
    ones, according to the CG representation. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
    Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.

8.4. Run the simulation
_______________________

Make a new folder for the run:

.. code-block:: bash

    mkdir -p run; cd run

The folder ``sirah.amber/tutorial/8/`` contains typical input files for energy minimization
(``em1_WT4.in`` and ``em2_WT4.in``), relaxation (or equilibration) (``eq1_WT4.in`` and ``eq2_WT4.in``) and production (``md_WT4.in``) runs. Please check carefully the
input flags therein, in particular the definition of flag *chngmask=0* at *&ewald* section is **mandatory**.

.. tip::

    **Some commonly used flags in Amber**

   - ``-i``: Input file.
   - ``-o``: Output file.
   - ``-p``: Parameter/topology file.
   - ``-c``: Coordinate file.
   - ``-r``: Restart file.
   - ``-x``: Trajectory file.
   - ``-ref``: Reference file

.. caution::

	These input files are executed by the **GPU** implementation of ``pmemd.cuda``. Other available modules are ``sander`` or ``pmemd``, which are both **CPU** implementations of Amber.

.. note::

	The same input files can be used to run on CPU with the modules ``pmemd`` or ``sander``.
	
	
**Energy Minimization of side chains and solvent by restraining the protein backbone and glycan rings:**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/8/em1_WT4.in -p ../1GYA_cg.prmtop -c ../1GYA_cg.ncrst -ref ../1GYA_cg.ncrst -o 1GYA_cg_em1.out -r 1GYA_cg_em1.ncrst &
 
**Energy Minimization of whole system:**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/8/em2_WT4.in -p ../1GYA_cg.prmtop -c ../1GYA_cg_em1.ncrst -o 1GYA_cg_em2.out -r 1GYA_cg_em2.ncrst &

**Solvent Relaxation (or equlibration) in NPT:**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/8/eq1_WT4.in -p ../1GYA_cg.prmtop -c 1GYA_cg_em2.ncrst -ref 1GYA_cg_em2.ncrst -o 1GYA_cg_eq1.out -r 1GYA_cg_eq1.ncrst -x 1GYA_cg_eq1.nc &
  
.. caution::

	Option **restraintmask=:'1-114'** in input file ``eq1_WT4.in`` must be set specifically for each system to restrain all glycoprotein’s residues.

**Soft Relaxation to improve side chain and glycan solvation (NPT):**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/8/eq2_WT4.in -p ../1GYA_cg.prmtop -c 1GYA_cg_eq1.ncrst -ref 1GYA_cg_eq1.ncrst -o 1GYA_cg_eq2.out -r 1GYA_cg_eq2.ncrst -x 1GYA_cg_eq2.nc &
  

**Production (1000ns):**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/8/md_WT4.in -p ../1GYA_cg.prmtop -c 1GYA_cg_eq2.ncrst -o 1GYA_cg_md.out -r 1GYA_cg_md.ncrst -x 1GYA_cg_md.nc &



8.5. Visualizing the simulation
________________________________

That’s it! Now you can analyze the trajectory.
Process the output trajectory to account for the Periodic Boundary Conditions (PBC):

  .. code-block:: bash

      echo -e "autoimage\ngo\nquit\n" | cpptraj -p ../1GYA_cg.prmtop -y 1GYA_cg_md.nc -x 1GYA_cg_md_pbc.nc --interactive

Load the processed trajectory in VMD:

.. code-block::

    vmd ../1GYA_cg.prmtop ../1GYA_cg.ncrst 1GYA_cg_md.nc -e ../sirah.amber/tools/sirah_vmdtk.tcl

.. note::

     The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. Use the command ``sirah-help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.
