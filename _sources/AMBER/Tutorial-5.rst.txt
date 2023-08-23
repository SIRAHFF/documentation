This tutorial shows how to use the SIRAH force field to perform a coarse grained (CG) simulation of a
protein in explicit solvent (called WatFour, WT4). The main references for
this tutorial are: `WAT4 <https://pubs.acs.org/doi/abs/10.1021/ct100379f>`_, `SIRAH 2.0 <https://doi.org/10.1021/acs.jctc.9b00006>`_ and `SIRAH Tools <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_. We strongly advise you to read these articles before starting the tutorial.

.. important::

    Check :ref:`download <download amber>` section for download and set up details before to start this tutorial.     
    Since this is **tutorial 5**, remember to replace ``X.X`` in your folder directory. The files corresponding to this tutorial can be found in: ``sirah_[version].amber/tutorial/5/``


5.1. Build CG representations
_____________________________

.. caution::

	The mapping to CG requires the correct protonation state of each residue at a given pH. We recommend using the `PDB2PQR server <https://server.poissonboltzmann.org/pdb2pqr>`_ and choosing the output naming scheme of AMBER for best compatibility. Be aware that modified residues lacking parameters such as: MSE (seleno MET), TPO (phosphorylated THY), SEP (phosphorylated SER) or others are deleted from the PQR file by the server. In that case, mutate the residues to their unmodified form before submitting the structure to the server.    

Map the protonated atomistic structure of protein `1CRN <https://www.rcsb.org/structure/1CRN>`_ to its CG representation:   

.. code-block:: bash

  ./sirah.amber/tools/CGCONV/cgconv.pl -i ./sirah.amber/tutorial/5/1CRN.pqr -o 1CRN_cg.pdb  
  
The input file ``-i`` 1CRN.pqr contains the atomistic representation of `1CRN <https://www.rcsb.org/structure/1CRN>`_ structure at pH **7.0**, while the output ``-o`` 1CRN_cg.pdb is its SIRAH CG representation.


.. note::

	**Pay attention to residue names when mapping structures from other atomistic force fields or experimental structures.** Although we provide compatibility for naming schemes in PDB, GMX, GROMOS, CHARMM and OPLS, there always may be some ambiguity in the residue naming, specially regarding protonation states, that may lead to a wrong mapping. For example, SIRAH Tools always maps the residue name “HIS” to a Histidine protonated at epsilon nitrogen (:math:`N_{\epsilon}`) regardless the actual proton placement. Similarly, protonated Glutamic and Aspartic acid residues must be named “GLH” and “ASH”, otherwise they will be treated as negative charged residues. In addition, protonated and disulfide bonded Cysteines must be named “CYS” and “CYX” respectively. These kind of situations need to be carefully checked by the users. In all cases the residues preserve their identity when mapping and back-mapping the structures. Hence, the total charge of the protein should be the same at atomistic and SIRAH level. You can check the following mapping file to be sure of the compatibility: ``sirah.amber/tools/CGCONV/maps/sirah_prot.map``.    

  
.. important::

	By default charged termini are used, but it is possible to set then neutral by renaming the residues from **s**\[code\] to **a**\[code\] (Nt-acetylated) or **m**\[code\] (Ct-amidated) after mapping to CG, where \[code\] is the root residue name in SIRAH. For example, to set a neutral N-terminal Histidine protonated at epsilon nitrogen (:math:`N_{\epsilon}`) rename it from “sHe” to “aHe”.

.. tip::

  This is the basic usage of the script **cgconv.pl**, you can learn other capabilities from its help:
  ``./sirah.amber/tools/CGCONV/cgconv.pl -h``	

Please check both PDB and PQR structures using VMD:	

.. code-block:: bash

  vmd -m sirah.amber/tutorial/5/1CRN.pqr 1CRN_cg.pdb


From now on it is just normal AMBER stuff!


5.2. Prepare LEaP
_________________

Use a text editor to create the file ``gensystem.leap`` including the following lines:

.. code-block:: console

    # Load SIRAH force field
    addPath ./sirah.amber
    source leaprc.sirah

    # Load model
    protein = loadpdb 1CRN_cg.pdb

    # Info on system charge
    charge protein  
	
    # Set S-S bridges
    bond protein.3.BSG protein.40.BSG
    bond protein.4.BSG protein.32.BSG
    bond protein.16.BSG protein.26.BSG

    # Add solvent, counterions and 0.15M NaCl
    # Tuned solute-solvent closeness for best hydration
    solvateOct protein WT4BOX 20 0.7
    addIonsRand protein NaW 22 ClW 22

    # Save Parms
    saveAmberParmNetcdf protein 1CRN_cg.prmtop 1CRN_cg.ncrst

    # EXIT
    quit

.. caution::

    Each disulfide bond must be defined explicitly in LEaP using the command bond, e.g.: “*bond unit.ri.BSG unit.rj.BSG*”. Where *ri* and *rj* correspond to the residue index in the topology file starting from 1, which may differ from the biological sequence in the PDB file. You can try the command *pdb4amber* to get those indexes from the atomistic structure, but be aware that it may not work if the Cysteine residues are too far away:	

    .. code-block:: bash

	   pdb4amber -i sirah.amber/tutorial/5/1CRN.pqr -o 1CRN_aa.pdb && cat 1CRN_aa_sslink

	
.. seealso::

   The available ionic species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW). One ion pair (e.g. NaW-ClW) each 34 WT4 molecules renders a salt concentration of ~0.15M (see :ref:`Appendix <Appendix>` for details).
   Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`_.

5.3. Run LEaP
______________

Run the LEaP application to generate the molecular topology and initial coordinate files:

.. code-block:: bash

    tleap -f gensystem.leap

.. warning::

    Warning messages about long, triangular or square bonds in ``leap.log`` file are fine and expected due to the CG topology.

This should create a topology file ``1CRN_cg.prmtop`` and a coordinate file ``1CRN_cg.ncrst``.

Use VMD to check how the CG model looks like and particularly the presence of disulfide bonds:

.. code-block:: bash

  vmd 1CRN_cg.prmtop 1CRN_cg.ncrst -e ./sirah.amber/tools/sirah_vmdtk.tcl


.. tip::

    VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right
    ones. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
    Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages.

5.4. Run the simulation
_______________________

Make a new folder for the run:

.. code-block:: bash

    mkdir -p run; cd run

The folder ``sirah.amber/tutorial/5/`` contains typical input files for energy minimization
(``em1_WT4.in`` and ``em2_WT4.in``), equilibration (``eq1_WT4.in`` and ``eq2_WT4.in``) and production (``md_WT4.in``) runs. Please check carefully the
input flags therein, in particular the definition of flag *chngmask=0* at *&ewald* section is **mandatory**.

.. tip::

    **Some flags used in AMBER**

   - ``-i``: Input file.
   - ``-o``: Output file.
   - ``-p``: Parameter/topology file.
   - ``-c``: Coordinate file.
   - ``-r``: Restart file.
   - ``-x``: Trajectory file.
   - ``-ref``: Reference file

.. caution::

    These input files are executed by the **GPU** implementation of ``pmemd.cuda``. Other available implementations that could be used: ``sander``  or ``pmemd``, both **CPU** implementations of AMBER.


**Energy Minimization of side chains and solvent by restraining the backbone:**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/5/em1_WT4.in -p ../1CRN_cg.prmtop -c ../1CRN_cg.ncrst -ref ../1CRN_cg.ncrst -o 1CRN_cg_em1.out -r 1CRN_cg_em1.ncrst &
 
**Energy Minimization of whole system:**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/5/em2_WT4.in -p ../1CRN_cg.prmtop -c ../1CRN_cg_em1.ncrst -o 1CRN_cg_em2.out -r 1CRN_cg_em2.ncrst &

**Solvent Equilibration (NPT):**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/5/eq1_WT4.in -p ../1CRN_cg.prmtop -c 1CRN_cg_em2.ncrst -ref 1CRN_cg_em2.ncrst -o 1CRN_cg_eq1.out -r 1CRN_cg_eq1.ncrst -x 1CRN_cg_eq1.nc &
  
.. caution::

	Option **restraintmask=:'1-46'** in input file ``eq1_WT4.in`` must be set specifically for each system to embrace all protein’s residues.

**Soft equilibration to improve side chain solvation (NPT):**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/5/eq2_WT4.in -p ../1CRN_cg.prmtop -c 1CRN_cg_eq1.ncrst -ref 1CRN_cg_eq1.ncrst -o 1CRN_cg_eq2.out -r 1CRN_cg_eq2.ncrst -x 1CRN_cg_eq2.nc &
  

**Production (1000ns):**

.. code-block:: bash

	pmemd.cuda -O -i ../sirah.amber/tutorial/5/md_WT4.in -p ../1CRN_cg.prmtop -c 1CRN_cg_eq2.ncrst -o 1CRN_cg_md.out -r 1CRN_cg_md.ncrst -x 1CRN_cg_md.nc &


.. important::

	The same input files can be used to run on CPU *pmemd* or *sander*.
	

5.5. Visualizing the simulation
________________________________

That’s it! Now you can analyze the trajectory.
Process the output trajectory to account for the Periodic Boundary Conditions (PBC):

  .. code-block:: bash

      echo -e "autoimage\ngo\nquit\n" | cpptraj -p ../1CRN_cg.prmtop -y 1CRN_cg_md.nc -x 1CRN_cg_md_pbc.nc --interactive

**Load the processed trajectory in VMD:**

.. code-block::

    vmd ../1CRN_cg.prmtop ../1CRN_cg.ncrst 1CRN_cg_md.nc -e ../sirah.amber/tools/sirah_vmdtk.tcl

.. note::

    The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages.
