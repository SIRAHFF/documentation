.. caution::
    
    This tutorial was developed for the `Workshop and Course on Molecular, Physical, and Computational Virology <https://www.sirahff.com/events.html#Curso_virologia>`_. Please be aware that the scripts and commands given here might not work right on your system due to software versions diferences.

.. 
   _Please report bugs, errors or enhancement requests through `Issue Tracker <https://github.com/SIRAHFF/documentation/issues>`_ or if you have a question about SIRAH open a `New Discussion <https://github.com/SIRAHFF/documentation/discussions>`_.
   
This tutorial shows how to use the SIRAH force field to perform a coarse grained (CG) simulation of a
non-enveloped Virus-like Particle (VLP) in explicit solvent (called WatFour, WT4). The main references for
this tutorial are: `Darré et al. <https://pubs.acs.org/doi/abs/10.1021/ct100379f>`_, `Machado et al. <https://doi.org/10.1021/acs.jctc.9b00006>`_, `Machado et al. <https://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b00659>`__ and `Soñora et al. <https://pubs.acs.org/doi/10.1021/acs.jcim.0c01205>`_. We strongly advise you to read these articles before starting the tutorial.


.. important::

    Check the :ref:`Setting up SIRAH <download amber>` section for download and set up details before starting this tutorial.
    Since this is **Tutorial 9**, remember to replace ``X.X`` and the files corresponding to this tutorial can be found in: ``sirah_[version].amber/tutorial/9/``

.. warning::
    
    This tutorial uses tools like Packmol and PDB2PQR that are not part of the AmberTools suite. You can find the latest version of Packmol `here <https://m3g.github.io/packmol/>`__. You can install PDB2PQR locally by:

    .. code-block:: bash
        
        conda install conda-forge::pdb2pqr









9.1. Structure preprocessing
_____________________________

.. note::

   This tutorial provides the necessary steps to assemble viral particles for SIRAH CG simulation. However, it is crucial to address any potential errors in the PDB files before proceeding. These preprocessing steps can prevent issues related to incorrect atom labeling, missing residues, clashes, or other common artifacts that are present in raw PDB data. View Prof. Jodi Hadden-Perilla's fantastic video on `Tips and Tricks for MD Simulation of Viruses <https://www.youtube.com/watch?v=OHf1FqOLxTI>`_ for more details.
   
  

9.1.1 Working with the Biological Assembly of an entire VLP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. important::

   In this tutorial, we will use a Biological Assembly from the PDB to work with an already-assembled Porcine Circovirus 2 (PCV2) VLP. **In** :ref:`9.1.3 <assalt>`, **we demonstrate how to generate a VLP of a non-enveloped virus from its asymmetric unit if your system is not already assembled**.


The PCV2 `6OLA <https://www.rcsb.org/structure/6OLA>`_ structure, which includes the VLP and genomic DNA fragments in its inner core, is used in this tutorial. Download the Biological Assembly from the PDB, making sure you are in the ``sirah.amber/tutorial/9/`` folder:

.. code-block:: bash

  wget https://files.rcsb.org/download/6OLA-assembly1.cif.gz  

Decompress the Biological Assembly:

.. code-block:: bash

  gzip -d 6OLA-assembly1.cif.gz
  
Open the cif file with VMD:

.. code-block:: bash

  vmd 6OLA-assembly1.cif

After loading the structure, we must saved it in PDB format by going to **File** > **Save Coordinates** > **Selected Atoms** > **all** option. Enter **File type** as pdb and hit **Save**. Assign filename as ``6ola_assembly.pdb`` and click **OK**. 

You can close VMD now.

Some structures may have structural **issues** or **artifacts** that affect their preparation or simulation. Thus, it is crucial to check the details of the structure and its integrity prior to any step. 

Regarding the **6OLA** structure, there are three critical issues that require correction: 

* **Issue 1**. The structure contains 60 DNA fragments and each DNA fragment includes a phosphate atom from the subsequent nucleotide. However, this nucleotide was not modeled; only the phosphate was included. Thus, when we attempt to generate the topology using LEaP, the process is stopped since it fails to reconstruct a nucleotide that is absent. This should be corrected by removing the extra phosphate and oxigen atoms. We can correct this using the following command:

..
    _sed '/\(P\|OP1\|OP2\)[[:space:]]\+DC[[:space:]]\+[A-Z][[:space:]]\+1/d' 6ola_assembly.pdb > 6ola_edited.pdb

.. code-block:: bash

    sed '/\(P[[:space:]]\+1DC[[:space:]]\+[A-Z][[:space:]]\+1\|OP11DC[[:space:]]\+[A-Z][[:space:]]\+1\|OP21DC[[:space:]]\+[A-Z][[:space:]]\+1\)/d' 6ola_assembly.pdb > 6ola_modificated_1.pdb

This line of code creates a new PDB file by removing the lines containing the **P**, **OP1**, and **OP2** atoms from the original structure.

* **Issue 2**. Another issue is the incorrect assignment of residue names according to their protonation state. Following AMBER force field rules, histidine residue name can vary depending on whether the protonation occurs at the delta nitrogen (HID), epsilon nitrogen (HIE), or if the residue is doubly protonated (HIP). Different residue names can also be given to ASP (ASH), GLU (GLH), and CYS (CYM or CYX if they are part of a disulfide bond).

In order to comply with AMBER force field rules, the PDB file must be edited to update residues protonations. To achieve this, we will first inspect the PDB structure in VMD to identify the relevant protonations.

.. note::

    In this part of the tutorial, we established the protonation states using the proton locations as provided in the original Biological Assembly of `6OLA <https://www.rcsb.org/structure/6OLA>`_ structure. 

Open VMD:

.. code-block:: bash

    vmd 6ola_modificated_1.pdb  

Once in VMD, let's examine one of the fragments (or protomers) of the PCV2 capsid. To do this, go to the **Graphics** tab and select **Representations**. In the **Graphics Representations** menu, go to **Selected Atoms** and make the selection of histidines from fragment 0 using the following selection command:

.. code-block:: bash

    resname HIS and fragment 0

You should run this command again and check the protons position for all the other protonable residues, such as ASP, GLU, and CYS. 

.. important::

    The only residues in this structure that appear to have distinct protonation states are the three histidines. Use your mouse to navigate and identify which nitrogens are protonated to correctly adjust their residue's names. Following proton location, HIS 113 is protonated at the epsilon nitrogen (:math:`N_{\epsilon}`), HIS 122 is protonated at the delta nitrogen (:math:`N_{\delta}`), and HIS 160 is protonated at the epsilon nitrogen (:math:`N_{\epsilon}`). 

Close VMD. Use the command below to modify the names of the histidine's residues:

.. code-block:: bash
    
    sed -e '/HIS[[:space:]]\+[A-Z][[:space:]]113/s/HIS/HIE/' -e '/HIS[[:space:]]\+[A-Z][[:space:]]122/s/HIS/HID/' -e '/HIS[[:space:]]\+[A-Z][[:space:]]160/s/HIS/HIE/' 6ola_modificated_1.pdb > 6ola_modificated_2.pdb

..
    _sed -e '/ HIS . 113 /s/HIS/HIE/' -e '/ HIS . 122 /s/HIS/HID/' -e '/ HIS . 160 /s/HIS/HIE/' 6ola_edited.pdb > 6ola_edited_modificated.pdb


* **Issue 3**. Another important step for this structure is to correct the atom index numbers, ensuring they range from 0 to 99999. Additionally, a `TER` record should be added at the end of each chain, and an `END` record should be included at the end of the structure.

To do this, we can use the ``fixpdb.tcl`` VMD script found in the ``sirah.amber/tutorial/9/`` folder: 

.. code-block:: bash

    vmd -dispdev text 6ola_modificated_2.pdb -e ./sirah.amber/tutorial/9/fixpdb.tcl  

Then

.. code-block:: bash

    mv 6ola_modificated_2_final.pdb 6ola_modificated_final.pdb 

If all stages were completed successfully, the final PDB file ``6ola_modificated_final.pdb`` is generated.

* **Issue 4 (Issue 2 - Alternative)**. Instead of assigning protonation states by hand, you can use tools that take into account the surrounding residues and the pH conditions to do it automatically. You can use servers, like the `CHARMM-GUI server <https://www.charmm-gui.org/>`_ and `PDB2PQR server <https://server.poissonboltzmann.org/pdb2pqr>`_. The ``6ola_modificated_final.pdb`` should work without any issues, though it may take a while to execute.

.. tip::

   An additional approach to circumvent potential server execution issues is to submit a minimal fragment of your symmetrized VLP that includes all of the interfaces between subunits as input.

.. caution::

    The mapping to CG requires the correct protonation state of each residue at a given pH. We recommend using the `CHARMM-GUI server <https://www.charmm-gui.org/>`_ and use the **PDB Reader & Manipulator** to prepare your system. An account is required to access any of the CHARMM-GUI Input Generator modules, and it can take up to 24 hours to obtain one. 
    
    Other option is the `PDB2PQR server <https://server.poissonboltzmann.org/pdb2pqr>`_ and choosing the output naming scheme of AMBER for best compatibility. This server was utilized to generate the *PQR* file featured in this tutorial. Be aware that modified residues lacking parameters such as: MSE (seleno MET), TPO (phosphorylated TYR), SEP (phosphorylated SER) or others are deleted from the PQR file by the server. In that case, mutate the residues to their unmodified form before submitting the structure to the server.

Given the scale and complexity of the system, using the ``pdb2pqr`` script to perform the calculations locally is one method to get around potential server problems. Therefore, we can assign protonation states using pH ~ 7.0 with the following command:

 
.. code-block:: bash

    pdb2pqr --ff=AMBER --keep-chain --ffout=AMBER --titration-state-method=propka --with-ph=7.0 --include-header ./6ola_modificated_final.pdb ./6ola_modificated_final.pqr  

A PQR file ``6ola_modificated_final.pqr`` is generated containing protonation state of each residue at the given pH.

9.1.2 Calculate system charge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

	As a result of the system's size and the various methods available for assigning hydrogen atoms, the most straightforward method of verifying the current charge of the VLP is to calculate the number of charged amino acids and phosphate atoms within the system.   

To do this, we can use the ``CountCharge.sh`` script found in the ``sirah.amber/tutorial/9/`` folder by:     

.. code-block:: bash

    ./sirah.amber/tutorial/9/CountCharge.sh 6ola_modificated_final.pdb   

This script evaluates the number of positively and negatively charged residues in the system and provides the net charge that can be used in later steps of preparing the VLP for simulations.

.. tip::

	The file ``CountCharge.sh`` was made for this system, we suggest that you inspect the  file to gain a more comprehensive understanding of the script and then change it to work with your own system. 

.. code-block:: console

    +--------------------+----------+
    |      Molecule      | Count    |    
    +====================+==========+
    | LYS (CA)           |  540     |  
    | ARG (CA)           | 1020     |
    | GLU (CA)           |  240     |         
    | ASP (CA)           |  720     |   
    | Phosphate (P)      |  180     |   
    |--------------------|----------|
    |Total Positv Charge | 1560     |          
    |Total Negatv Charge | 1140     |      
    |Net Charge          |  420     |           
    +--------------------|----------+

.. _assalt:

9.1.3 Packing the assimetric unit to build an entire VLP (Alternative)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. caution::

    This step should be performed only when the input file (.pdb or .cif) is the asymmetric unit. For more details about asymetric unit and biological assemblies, check the `RCSB PDB website <https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/biological-assemblies>`_. 
    
In complex systems like viruses, some biological assemblies are constructed from the asymmetric unit, which in virology may be represented by a single protein or a small number of them, referred to as a capsomer or protomer. In these cases, it is essential to assemble the whole VLP by employing symmetry operations (rotations and translations) with the symmetry matrices included in the original PDB file. Thus, based on the available structure, it may be necessary to initially assemble the entire VLP.

In ``sirah.amber/tutorial/9/example_3JCI`` you will find the asymmetric unit of the capsid protein of the PCV2 (PDB code: `3JCI <https://www.rcsb.org/structure/3JCI>`__) along with the necessary files to assemble its VLP. 

The first step to assemble a VLP is to minimize the protomer. To do this, we need to create a ``3JCI.leap`` file to generate the molecular topology and initial coordinate files:

.. code-block:: console

    # Load AMBER force field
    source leaprc.protein.ff14SB

    # Load model
    pcv2 = loadPdb 3jci.pdb

    # Save Parms
    saveamberparm pcv2 3JCI_aa.prmtop 3JCI_aa.inpcrd

    # EXIT
    quit

Run the LEaP application with:

.. code-block:: bash

    tleap -f 3JCI.leap

Make a new folder to run the minimization:

.. code-block:: bash

    mkdir -p run; cd run
  
Run the minimization with sander:

.. code-block:: bash

    mpirun -np 4 sander.MPI -O -i ../em_cpu.in -p ../3JCI_aa.prmtop -c ../3JCI_aa.inpcrd -ref ../3JCI_aa.inpcrd -inf mdinfo -o 3JCI_aa_em.out -r 3JCI_aa_em.rst7 &
  
.. note::

    The same input file can be used to run on CPU with the modules ``pmemd``, ``sander`` or ``sander.MPI``. Or the GPU implementation, ``pmemd.cuda``.
    
Once the protomer has been minimized, return to the original folder:

.. code-block:: bash

  cd ..  
 
And build the VLP using the ``capsid.tcl`` VMD script: 

.. code-block:: bash

  vmd -dispdev text -e ./sirah.amber/tutorial/9/capsid.tcl  

.. note::

    Due to its complexity, we suggest that you inspect the ``capsid.tcl`` file to gain a more comprehensive understanding of the procedure.
  
That’s it! You can now use VMD to see what the VLP looks like:

.. code-block:: bash

  vmd PCV2_capsid_OK.pdb


9.2. Build CG representations
_____________________________

.. caution::
	
	We are using the ``6ola_modificated_final.pdb`` file for this section of the tutorial, which makes use of the assigned protonation states and hydrogen atoms as given in the original Biological Assembly VLP. But you may also use the automatically calculated protonation states file ``6ola_modificated_final.pqr``.

Map the atomistic structure of VLP that was prepared in the previous step to its CG representation:   

.. code-block:: bash

  ./sirah.amber/tools/CGCONV/cgconv.pl -i 6ola_modificated_final.pdb -o 6ola_modificated_final_cg.pdb
  

The input file ``-i`` 6ola_modificated_final.pdb has the atomistic model of the VLP that has had its structure issues fixed,, while the output ``-o`` 6ola_modificated_final_cg.pdb is its SIRAH CG representation.

.. tip::

    This is the basic usage of the script **cgconv.pl**, you can learn other capabilities from its help by typing:

    .. code-block:: bash

        ./sirah.amber/tools/CGCONV/cgconv.pl -h 
        
.. note::

    **Pay attention to residue names when mapping structures from other atomistic force fields or experimental structures.** Although we provide compatibility for naming schemes in PDB, GMX, GROMOS, CHARMM and OPLS, there might always be some ambiguity in the residue naming, specially regarding protonation states, that may lead to a wrong mapping. For example, SIRAH Tools always maps the residue name “HIS” to a Histidine protonated at the epsilon nitrogen (:math:`N_{\epsilon}`) regardless the actual proton placement. Similarly, protonated Glutamic and Aspartic acid residues must be named “GLH” and “ASH”, otherwise they will be treated as negative charged residues. In addition, protonated and disulfide bonded Cysteines must be named “CYS” and “CYX” respectively. These kind of situations need to be carefully checked by the users. In all cases the residues preserve their identity when mapping and back-mapping the structures. Hence, the total charge of the protein should be the same at atomistic and SIRAH levels. You can check the following mapping file to be sure of the compatibility: ``sirah.amber/tools/CGCONV/maps/sirah_prot.map``.    

  
.. important::

    By default, charged termini are used, but it is possible to set them neutral by renaming the residues from **s**\[code\] to **a**\[code\] (Nt-acetylated) or **m**\[code\] (Ct-amidated) after mapping to CG, where \[code\] is the root residue name in SIRAH. For example, to set a neutral N-terminal Histidine protonated at epsilon nitrogen (:math:`N_{\epsilon}`) rename it from “sHe” to “aHe”.


Please check both structures using VMD: 

.. code-block:: bash

  vmd -m 6ola_modificated_final.pdb 6ola_modificated_final_cg.pdb



9.3. Wrapping up VLP system with Packmol
_________________________________________________

.. warning::

    Before packing the VLP system, it is necessary to estimate the number of coarse-grained water molecules (WT4) per layer (see **Figure 1**), i.e., the inner-virus and outer-virus regions. In this tutorial, we will use the following radius to pack the system with Packmol (left) and LEaP (right).

    .. figure:: /../images/water_box_virus.png
     :align: center
     :width: 100%

    **Figure 1.** Dimensions of the solvation box used in this tutorial.

To calculate the approximate number of solvent molecules in each layer, you can use a function in the ``calc_n.awk`` file that estimates the number of molecules in a given volume. 

.. code-block:: bash

  awk -f ./sirah.amber/tutorial/9/calc_n.awk ./sirah.amber/tutorial/9/layers_radius.dat  

The two lines in the ``layers_radius.dat`` are read by this code. One line calculates the number of WT4 outside the VLP, and the other line calculates the number of WT4 inside the VLP (see **Figure 1**).

.. note::

    Each line of the ``layers_radius.dat`` has five parameters. The first and second ones are WT4's molecular weight and density, respectively. The third and fourth parameters are the starting and end point of the compartment and defines its radius in Angstroms (Å). The fourth value needs to be greater than the third. The fifth and final parameter gives a label to each compartment. Depending on your system, you may need to adjust the radius values in the ``layers_radius.dat`` file. 

The output of this command is: 

.. code-block:: bash

	WT4_out = 18235
	WT4_in = 6393

The next step is to use the program Packmol to put together the system's parts: VLP, solvent (both inside and outside), and ions (inside) (see **Figure 1**).

.. 
    _This tutorial step has been simplified to eliminate the need for numerous Packmol and tLEaP rounds. As a result, we provide the precise quantity of ion and water molecules required to complete one Packmol step and one tLEaP step. It is important to keep in mind, though, that the whole process usually involves making the whole water box to guess how many water molecules are there and calculate the right salt concentration before the final packed system is reached. 	

To accurately add the right number of solvent and ion molecules using Packmol, it is important to consider the charge of the VLP. As previously determined, the net charge of the VLP is **+420**. 

To balance this, we need to introduce at least **420** negative ions (Cl⁻) into the inner core of the VLP. However, previous simulations conducted by us and other researchers demonstrated through trial and error that the inner core of the virion must have a minimum of **620 Cl⁻** ions to not only neutralize the VLP’s charge but also compensate for the lack of a complete genome within the particle. As a result, the inner solvent layer contains::

	620 Cl⁻ beads (ClW in its CG naming).

	6393 - 620 = 5773 WT4 molecules at the inner layer.

.. caution::

     It is important to be aware that the charges may vary when using the ``6ola_modificated_final.pqr``, and as a result, the number of ions necessary to neutralize the charge imbalance may differ.

This information is provided to the ``PCV2.pkm`` file, which is then utilized by Packmol to pack the system. Use the following command to execute Packmol: 

.. code-block:: bash

  packmol < ./sirah.amber/tutorial/9/PCV2.pkm >> PCV2_packmol.log &  

.. important::

    Edit the ``PCV2.pkm`` file to fix the number of solvent and ions according to your system.

.. 
    _After running the first iteration with Packmol, we used the Packmol output as input for TLEaP to add the necessary water molecules to create a truncated octahedral box. The output from TLEaP results in a system with a total of approximately 23,450 WT4 water molecules, though this number may vary slightly. Initially, the system contained 14,010 WT4 molecules, meaning that 9,440 WT4 molecules were added in this step using TLEaP. It’s important to remember that biological systems under physiological conditions are typically simulated at a salt concentration of 150 mM. Therefore, to reach this concentration, we had added the corresponding number of ion pairs for the 9,440 newly added WT4 molecules which corresponds to 278 ion pairs. At the end we have 556 ions, which reduces the updated number of water molecules to 8,884.

.. 
    _In this new iteration with Packmol, we will add 278 NaW and 278 ClW ions to the region outside the viral particle. Then, we will use the Packmol output as input for TLEaP once again to add the 8,884 WT4 water molecules and create the truncated octahedral box. This setup will allow us to begin the energy minimization and equilibration process.

In order for LEaP to recognize proteins and DNA as distinct molecules, we also need to add `TER` records at the end of each chain and a `END` record at the end of the PDB file.

Thus, we can use the ``fixpdb.tcl`` VMD script found in the ``sirah.amber/tutorial/9/`` folder again: 

.. code-block:: bash

    vmd -dispdev text PCV2_6OLA_CG_solv.pdb -e ./sirah.amber/tutorial/9/fixpdb.tcl  


A PDB file ``PCV2_6OLA_CG_solv_final.pdb`` is generated and ready for LEaP.

From now on it is just normal Amber stuff!


9.4 Prepare LEaP input
_________________________________________________

In this step, we will set the truncated octahedral box and the salt concentration to 150 mM with LEaP. There are 35,720 WT4 molecules in the system after LEaP (an increase of 11,712 WT4). To reach 150 mM NaCl concentration on the outer layer of WT4, we need to add::

    35720/34 = 1050 ionic pairs

.. seealso::

       The available electrolyte species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW) which represent solvated ions in solution. One ion pair (e.g., NaW-ClW) each 34 WT4 molecules results in a salt concentration of ~0.15M (see :ref:`Appendix <Appendix>` for details). Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`__.

However, keep in mind that the VLP currently has a 200 negative charge::
    
    420 - 620 = -200

So, in order to achieve a more balanced system, we can split this charge. In this case, the number of counterions to add is::

    1,050 + 100 = 1,150 NaW

    1,050 - 100 = 950 ClW

.. seealso::

       Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`__.


Use a text editor to create the file ``6OLA_gensystem.leap`` including the following lines:

.. code-block:: console

    # Load SIRAH force field
    addPath ./sirah.amber
    source leaprc.sirah

    # Load model
    pcv2 = loadpdb PCV2_6OLA_CG_solv_final.pdb

    # Info on system charge
    charge pcv2  

    # Add solvent and build a truncated octahedral box
    # Tuned solute-solvent closeness for best hydration
    solvateOct pcv2 WT4BOX 17 0.9

    # Add counterions and 0.15M NaCl
    addIonsRand pcv2 NaW 1150 ClW 950

    # Save Parms
    saveAmberParmNetcdf pcv2 PCV2_6OLA_CG.prmtop PCV2_6OLA_CG.ncrst

    # Save pdb file
    savepdb pcv2 PCV2_6OLA_CG.pdb

    # EXIT
    quit


.. caution::

    There is only one CYS in this VLP, and there are no disulfide bridges in it. However, if you are working with a system that includes these bonds, please take note of the following:

    Each disulfide bond must be defined explicitly in LEaP using the command bond, e.g.: “*bond unit.ri.BSG unit.rj.BSG*”. Where *ri* and *rj* correspond to the residue index in the topology file starting from 1, which may differ from the biological sequence in the PDB file. You can try the command *pdb4amber* to get those indexes from the atomistic structure, but be aware that it may not work if the Cysteine residues are too far away: 

    .. code-block:: bash

       pdb4amber -i YOUR_VLP.pqr -o YOUR_VLP_aa.pdb && cat YOUR_VLP_aa_sslink


9.5 Run LEaP 
_________________________________________________

Run LEaP to generate the molecular topology and initial coordinate files:

.. code-block:: bash

    tleap -f 6OLA_gensystem.leap

.. note::

    Warning messages about long, triangular or square bonds in ``leap.log`` file are fine and expected due to the CG topology of some residues.


This should create a topology file ``PCV2_6OLA_CG.prmtop``, a coordinate file ``PCV2_6OLA_CG.ncrst`` and a pdb file ``PCV2_6OLA_CG.pdb``.

Use VMD to check how the CG model looks like and particularly the presence of disulfide bonds whether your system has them:

.. code-block:: bash

  vmd PCV2_6OLA_CG.prmtop PCV2_6OLA_CG.ncrst -e ./sirah.amber/tools/sirah_vmdtk.tcl


.. tip::

    VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right
    ones, according to the CG representation. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
    Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.


9.6. Run the simulation
_______________________

Make a new folder for the run:

.. code-block:: bash

    mkdir -p run; cd run

The folder ``sirah.amber/tutorial/9/`` contains typical input files for energy minimization
(``VLP_cg_em1.in`` to ``VLP_cg_em5.in``), equilibration (``VLP_cg_eq1.in`` to ``VLP_cg_eq4.in``) and production (``VLP_cg_md1.in``) runs. Please check carefully the
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


**Energy Minimization of solvent and ions by restraining the whole protein and DNA fragments:**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/9/VLP_cg_em1.in -p PCV2_6OLA_CG.prmtop -c PCV2_6OLA_CG.ncrst -ref PCV2_6OLA_CG.ncrst -o PCV2_6OLA_CG_em1.out -r PCV2_6OLA_CG_em1.ncrst &
 

**Energy Minimization of side chains and solvent by restraining the backbone of proteins and DNA:**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/9/VLP_cg_em2.in -p PCV2_6OLA_CG.prmtop -c PCV2_6OLA_CG_em1.ncrst -ref PCV2_6OLA_CG_em1.ncrst -o PCV2_6OLA_CG_em2.out -r PCV2_6OLA_CG_em2.ncrst &
 

**Energy Minimization of whole system:**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/9/VLP_cg_em3.in -p PCV2_6OLA_CG.prmtop -c PCV2_6OLA_CG_em2.ncrst -ref PCV2_6OLA_CG_em2.ncrst -o PCV2_6OLA_CG_em3.out -r PCV2_6OLA_CG_em3.ncrst &


**Energy Minimization of side chains and solvent by restraining the backbone:**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/9/VLP_cg_em4.in -p PCV2_6OLA_CG.prmtop -c PCV2_6OLA_CG_em3.ncrst -ref PCV2_6OLA_CG_em3.ncrst -o PCV2_6OLA_CG_em4.out -r PCV2_6OLA_CG_em4.ncrst &
 

**Energy Minimization of whole system:**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/9/VLP_cg_em5.in -p PCV2_6OLA_CG.prmtop -c PCV2_6OLA_CG_em4.ncrst -ref PCV2_6OLA_CG_em4.ncrst -o PCV2_6OLA_CG_em5.out -r PCV2_6OLA_CG_em5.ncrst &

**Solvent Equilibration (NPT) by restraining the whole protein and DNA fragments:**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/9/VLP_cg_eq1.in -p PCV2_6OLA_CG.prmtop -c PCV2_6OLA_CG_em5.ncrst -ref PCV2_6OLA_CG_em5.ncrst -o PCV2_6OLA_CG_eq1.out -r PCV2_6OLA_CG_eq1.ncrst -x PCV2_6OLA_CG_eq1.nc &

**Soft equilibration to improve side chain solvation (NPT) by restraining the backbone of proteins and DNA:**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/9/VLP_cg_eq2.in -p PCV2_6OLA_CG.prmtop -c PCV2_6OLA_CG_eq1.ncrst -ref PCV2_6OLA_CG_eq1.ncrst -o PCV2_6OLA_CG_eq2.out -r PCV2_6OLA_CG_eq2.ncrst -x PCV2_6OLA_CG_eq2.nc &

**Soft solvent Equilibration (NPT) by restraining the proteins' backbone:**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/9/VLP_cg_eq3.in -p PPCV2_6OLA_CG.prmtop -c PCV2_6OLA_CG_eq2.ncrst -ref PCV2_6OLA_CG_eq2.ncrst -o PCV2_6OLA_CG_eq3.out -r PCV2_6OLA_CG_eq3.ncrst -x PCV2_6OLA_CG_eq3.nc &


**Whole System Equilibration (NPT) without restraints:**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/9/VLP_cg_eq4.in -p PCV2_6OLA_CG.prmtop -c PCV2_6OLA_CG_eq3.ncrst -ref PCV2_6OLA_CG_eq3.ncrst -o PCV2_6OLA_CG_eq4.out -r PCV2_6OLA_CG_eq4.ncrst -x PCV2_6OLA_CG_eq4.nc &
  

**Production (1000ns):**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/9/VLP_cg_md1.in -p PCV2_6OLA_CG.prmtop -c PCV2_6OLA_CG_eq4.ncrst -ref PCV2_6OLA_CG_eq4.ncrst -o PCV2_6OLA_CG_md1.out -r PCV2_6OLA_CG_md1.ncrst -x PCV2_6OLA_CG_md1.nc &



9.7. Visualizing the simulation
________________________________

That’s it! Now you can analyze the trajectory.
Process the output trajectory to account for the Periodic Boundary Conditions (PBC):

.. code-block:: bash

      echo -e "center (@GN,GC,GO)\n image familiar\n go\n quit\n" | cpptraj -p PCV2_6OLA_CG.prmtop -y PCV2_6OLA_CG_md1.nc -x PCV2_6OLA_CG_md1_pbc.nc --interactive

Load the processed trajectory in VMD:

.. code-block::

    vmd PCV2_6OLA_CG.prmtop PCV2_6OLA_CG.ncrst PCV2_6OLA_CG_md1_pbc.nc -e ../sirah.amber/tools/sirah_vmdtk.tcl

.. note::

     The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. Use the command ``sirah-help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.
