.. note::

   Please report bugs, errors or enhancement requests through `Issue Tracker <https://github.com/SIRAHFF/documentation/issues>`_ or if you have a question about SIRAH open a `New Discussion <https://github.com/SIRAHFF/documentation/discussions>`_.
   
This tutorial shows how to use the SIRAH force field to perform a coarse grained (CG) simulation of a protein in explicit solvent (called WatFour, WT4). The main references for
this tutorial are: `Cantero et al. <https://doi.org/10.1021/acs.jpcb.4c03278>`_, `Darré et al. <https://pubs.acs.org/doi/abs/10.1021/ct100379f>`_, `Machado et al. <https://doi.org/10.1021/acs.jctc.9b00006>`__ and `Machado & Pantano  <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_. We strongly advise you to read these articles before starting the tutorial.

.. important::

    Check the :ref:`Setting up SIRAH <download amber>` section for download and set up details before starting this tutorial.
    Since this is **Tutorial 1**, remember to replace ``X.X`` and the files corresponding to this tutorial can be found in: ``sirah_[version].amber/tutorial/NAMD/1/``




1.1. Build CG representations
_____________________________

.. note::

    In this section, we use the inherent functionality of NAMD to parse AMBER files. The tLEaP module of AmberTools is used to generate parameters and topology information.  


.. caution::

	The mapping to CG requires the correct protonation state of each residue at a given pH. We recommend using the `CHARMM-GUI server <https://www.charmm-gui.org/>`_ and use the **PDB Reader & Manipulator** to prepare your system. An account is required to access any of the CHARMM-GUI Input Generator modules, and it can take up to 24 hours to obtain one. 
	
	Other option is the `PDB2PQR server <https://server.poissonboltzmann.org/pdb2pqr>`_ and choosing the output naming scheme of AMBER for best compatibility. This server was utilized to generate the *PQR* file featured in this tutorial. Be aware that modified residues lacking parameters such as: MSE (seleno MET), TPO (phosphorylated TYR), SEP (phosphorylated SER) or others are deleted from the PQR file by the server. In that case, mutate the residues to their unmodified form before submitting the structure to the server.

Map the protonated atomistic structure of protein `1CRN <https://www.rcsb.org/structure/1CRN>`_ to its CG representation:   

.. code-block:: bash

  ./sirah.amber/tools/CGCONV/cgconv.pl -i ./sirah.amber/tutorial/NAMD/1/1CRN.pqr -o 1CRN_cg.pdb  
  
The input file ``-i`` 1CRN.pqr contains the atomistic representation of `1CRN <https://www.rcsb.org/structure/1CRN>`_ structure at pH **7.0**, while the output ``-o`` 1CRN_cg.pdb is its SIRAH CG representation.

.. tip::

	This is the basic usage of the script **cgconv.pl**, you can learn other capabilities from its help by typing:

	.. code-block:: bash

		./sirah.amber/tools/CGCONV/cgconv.pl -h	
		
.. note::

	**Pay attention to residue names when mapping structures from other atomistic force fields or experimental structures.** Although we provide compatibility for naming schemes in PDB, GMX, GROMOS, CHARMM and OPLS, there might always be some ambiguity in the residue naming, specially regarding protonation states, that may lead to a wrong mapping. For example, SIRAH Tools always maps the residue name “HIS” to a Histidine protonated at the epsilon nitrogen (:math:`N_{\epsilon}`) regardless the actual proton placement. Similarly, protonated Glutamic and Aspartic acid residues must be named “GLH” and “ASH”, otherwise they will be treated as negative charged residues. In addition, protonated and disulfide bonded Cysteines must be named “CYS” and “CYX” respectively. These kind of situations need to be carefully checked by the users. In all cases the residues preserve their identity when mapping and back-mapping the structures. Hence, the total charge of the protein should be the same at atomistic and SIRAH levels. You can check the following mapping file to be sure of the compatibility: ``sirah.amber/tools/CGCONV/maps/sirah_prot.map``.    

  
.. important::

	By default, charged termini are used, but it is possible to set them neutral by renaming the residues from **s**\[code\] to **a**\[code\] (Nt-acetylated) or **m**\[code\] (Ct-amidated) after mapping to CG, where \[code\] is the root residue name in SIRAH. For example, to set a neutral N-terminal Histidine protonated at epsilon nitrogen (:math:`N_{\epsilon}`) rename it from “sHe” to “aHe”.


Please check both PDB and PQR structures using VMD:	

.. code-block:: bash

  vmd -m sirah.amber/tutorial/NAMD/1/1CRN.pqr 1CRN_cg.pdb


From now on it is just normal Amber stuff!


1.2. Prepare with LEaP
_____________________________

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
    solvatebox protein WT4BOX 20 0.7
    addIonsRand protein NaW 22 ClW 22

    # Save Parms
    saveAmberParm protein 1CRN_cg_ionized.prmtop 1CRN_cg_ionized.rst
    savepdb protein 1CRN_cg_ionized.pdb

    # EXIT
    quit


.. caution::

    Each disulfide bond must be defined explicitly in LEaP using the command bond, e.g.: “*bond unit.ri.BSG unit.rj.BSG*”. Where *ri* and *rj* correspond to the residue index in the topology file starting from 1, which may differ from the biological sequence in the PDB file. You can try the command *pdb4amber* to get those indexes from the atomistic structure, but be aware that it may not work if the Cysteine residues are too far away:	

    .. code-block:: bash

	   pdb4amber -i sirah.amber/tutorial/NAMD/1/1CRN.pqr -o 1CRN_aa.pdb && cat 1CRN_aa_sslink

	
.. seealso::

       The available electrolyte species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW) which represent solvated ions in solution. One ion pair (e.g., NaW-ClW) each 34 WT4 molecules results in a salt concentration of ~0.15M (see :ref:`Appendix <Appendix>` for details). Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`__.
	   

1.3. Run LEaP 
____________________

Run the tLEaP application to generate the molecular topology and initial coordinate files:

.. code-block:: bash

    tleap -f gensystem.leap

.. note::

    Warning messages about long, triangular or square bonds in ``leap.log`` file are fine and expected due to the CG topology of some residues.


This should create a topology file ``1CRN_cg_ionized.prmtop`` and a coordinate file ``1CRN_cg_ionized.rst``. The last line of ``1CRN_cg_ionized.rst`` file contains the cell dimension information needed in the NAMD configuration file, for additional information, please refer to the `Using the AMBER force field in NAMD documentation <https://ambermd.org/namd/namd_amber.html>`_.

.. note::
    
    To check the last line of the 1CRN_cg_ionized.rst you can use:

    .. code-block:: bash

        tail -n 1 1CRN_cg_ionized.rst
    
    Save this information so it can be used in NAMD input files.


For this tutorial, the cell dimensions are:

.. code-block:: bash

       73.3223400  70.2433400  72.8663400  90.0000000  90.0000000  90.0000000

The first three values represent the x, y, and z dimensions. The remaining three values define an orthorhombic box.

Use VMD to check how the CG model looks like and particularly the presence of disulfide bonds:

.. code-block:: bash

  vmd 1CRN_cg_ionized.prmtop 1CRN_cg_ionized.rst -e ./sirah.amber/tools/sirah_vmdtk.tcl


.. tip::

    VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right
    ones, according to the CG representation. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
    Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.


1.4. Create Backbone and Protein restraints 
________________________________________________

In NAMD, it is necessary to assign a PDB file that contains the system's restraints. Generally, this is done by assigning values in the last column of the PDB file that corresponds to the B-factor. These files can be generated using a VMD script. 

Use a text editor to create the file ``restraints.tcl`` including the following lines:

.. code-block:: console

    # Load SIRAH tools
    source sirah.amber/tools/sirah_vmdtk.tcl
    
    # Upload the system
    mol new 1CRN_cg_ionized.prmtop type parm7 waitfor all
    mol addfile 1CRN_cg_ionized.pdb type pdb waitfor all

    # Clean B column
    set all [atomselect top all]
    $all set beta 0

    # Setup the backbone restraints
    set bb [atomselect top "sirah_protein and sirah_backbone"]
    $bb set beta 2.4
    $all writepdb bb_restraints.pdb

    # Setup the protein restraints
    $all set beta 0
    set prot [atomselect top "sirah_protein"]
    $prot set beta 2.4
    $all writepdb prot_restraints.pdb

    # Exit
    quit

Run the VMD script to generate the pdb restriction file:

.. code-block:: bash

    vmd -dispdev text –e restraints.tcl 


1.5. Run the simulation
_______________________

1.5.1 NAMD2
~~~~~~~~~~~~~

Make a new folder for the run:

.. code-block:: bash

    mkdir -p run; cd run
  
Copy the input files to the folder:

.. code-block:: bash

    cp sirah.amber/tutorial/NAMD/1/NAMD2/*.conf .

The folder ``sirah.amber/tutorial/NAMD/1/NAMD2`` contains typical input files for energy minimization (``em1.conf`` and ``em2.conf``), heating (``heat.conf``), equilibration (``eq1.conf`` and ``eq2.conf``) and production (``md.conf``) runs.  Please carefully review the input files, paying especially attention to the cell dimension values, names, and restrictions.


.. note::

    The same input files can be used to run on CPU or GPU. However, in NAMD2 (CUDA), the number of processors used (+p option) significantly affects performance. By contrast, in NAMD3 (CUDA), this value does not directly correlate to higher performance.
    
.. tip::

    If you have more than one GPU card, be sure you set the GPU number properly. For example, in order to utilize GPU 0, it is necessary to execute this command prior to running:

    .. code-block:: bash
        
        export CUDA_VISIBLE_DEVICES=0

    You might also specify which and how many CPUs you need to used. For example, if you require 24 CPUs:

    .. code-block:: bash

        namd2 +p24 namd_input.conf > namd_input.log &

    To indicate which cores to use:

    .. code-block:: bash

        namd2 +setcpuaffinity +pemap 0-23 +p24 namd_input.conf > namd_input.log &

    
**Energy Minimization of side chains and solvent by restraining the backbone:**

.. code-block:: bash

    namd2 +p8 em1.conf > em1.log &

.. note::

    In this stage, the restriction file ``bb_restraints.pdb`` is assigned to the consref and conskfile flags.
 
**Energy Minimization of whole system:**

.. code-block:: bash

    namd2 +p8 em2.conf > em2.log &

**Solvent Equilibration (NPT):**

.. code-block:: bash

    namd2 +p8 heat.conf > heat.log &
  
.. note::

    In this stage, the restriction file ``prot_restraints.pdb`` is assigned to the consref and conskfile flags.

**Solvent equilibration (NPT):**

.. code-block:: bash

    namd2 +p8 eq1.conf > eq1.log &
	
**Soft equilibration to improve side chain solvation (NPT):**

.. code-block:: bash

    namd2 +p8 eq2.conf > eq2.log &

.. note::

    In this stage, the restriction file ``bb_restraints.pdb`` is assigned to the consref and conskfile flags. The constraintScaling is set to 0.1 to soften the restraints.

**Production (1000ns):**

.. code-block:: bash

    namd2 +p8 md.conf > md.log &


1.5.2 NAMD3
~~~~~~~~~~~~~

Make a new folder for the run:

.. code-block:: bash

    mkdir -p run; cd run

Copy the input files to the folder:

.. code-block:: bash

    cp sirah.amber/tutorial/NAMD/1/NAMD3/*.conf .

The folder ``sirah.amber/tutorial/NAMD/1/NAMD3`` contains typical input files for energy minimization (``em1.conf`` and ``em2.conf``), heating (``heat.conf``), equilibration (``eq1.conf`` and ``eq2.conf``) and production (``md.conf``) runs.  Please carefully review the input files, paying especially attention to the cell dimension values, names, and restrictions.


.. note::

    The same input files can be used to run on on CPU or GPU. However, in NAMD2 (CUDA), the number of processors used (+p option) significantly affects performance. By contrast, in NAMD3 (CUDA), this value does not directly correlate to higher performance.

.. warning::
    
    To use GPU cards in NAMD3, you need to enable the GPU-resident mode with the CUDASOAintegrate option. In the input files for this tutorial the option is enabled.

    .. code-block:: console

        #NAMD3 parameters 
        if {1} {                                ;# Enable the block
        CUDASOAintegrate     on                 ;# Enable GPU-resident mode.
        }


.. tip::

    If you have more than one GPU card, be sure you set the GPU number properly. For example, in order to utilize GPU 0, it is necessary to execute this command prior to running:

    .. code-block:: bash
        
        export CUDA_VISIBLE_DEVICES=0

    You might also specify which and how many CPUs you need to used. For example, if you require 4 CPUs:

    .. code-block:: bash

        namd3 +p4 namd_input.conf > namd_input.log &

    To indicate which cores to use:

    .. code-block:: bash

        namd3 +setcpuaffinity +pemap 0-3 +p4 namd_input.conf > namd_input.log &

    
**Energy Minimization of side chains and solvent by restraining the backbone:**

.. code-block:: bash

    namd3 +p4 em1.conf > em1.log &

.. note::

    In this stage, the restriction file ``bb_restraints.pdb`` is assigned to the consref and conskfile flags.
 
**Energy Minimization of whole system:**

.. code-block:: bash

    namd3 +p4 em2.conf > em2.log &

**Solvent Equilibration (NPT):**

.. code-block:: bash

    namd3 +p4 heat.conf > heat.log &
  
.. note::

    In this stage, the restriction file ``prot_restraints.pdb`` is assigned to the consref and conskfile flags.

**Solvent equilibration (NPT):**

.. code-block:: bash

    namd3 +p4 eq1.conf > eq1.log &


**Soft equilibration to improve side chain solvation (NPT):**

.. code-block:: bash

    namd3 +p4 eq2.conf > eq2.log &

.. note::

    In this stage, the restriction file ``bb_restraints.pdb`` is assigned to the consref and conskfile flags. The constraintScaling is set to 0.1 to soften the former restraints.

**Production (1000ns):**

.. code-block:: bash

    namd3 +p4 md.conf > md.log &


1.6. Visualizing the simulation
________________________________

That’s it! Now you can analyze the trajectory.

Load the processed trajectory in VMD:

.. code-block::

    vmd ../1CRN_cg_ionized.prmtop MD.dcd -e sirah.amber/tools/sirah_vmdtk.tcl

.. note::

     The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. Use the command ``sirah-help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.


1.7 How to modify input files
________________________________

The provided NAMD input files were created using the information from this tutorial; for other systems, check and carefully adjust the inputs in accordance with the tips bellow.

1.7.1 All files
~~~~~~~~~~~~~~~~~

For all NAMD input (``*.conf``) files:

1. Change input file names:

.. code-block:: console

    amber         on                           ;# Turns on the AMBER Format inputs
    parmfile      your_system.prmtop           ;# File containing the force field parameters
    coordinates   your_system.pdb              ;# File containing initial coordinates

2. Change output file names:

.. code-block:: console

    set outputName     your_output_name        ;# Base name of output simulation files

1.7.2 Energy minimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For ``em1.conf``:

1. Assign the file name of the restriction file in consref and conskfile line:

.. code-block:: console

    #Constraints on the Backbone 
    if {1} {                                  ;# If 1 read the block 
    constraints          on                   ;# Turns on constraints 
    consref              your_restraints.pdb  ;# Reference PDB file for constraint positions 
    conskfile            your_restraints.pdb  ;# File containing constraint force constants 
    constraintScaling    1                    ;# Scaling factor for constraint forces 
    consexp              2                    ;# Exponent for constraint potential 
    conskcol             B                    ;# Column for constraint forces in the pdb file
    }

2. Set periodic cell information with the values from your rst file. NAMD will not read the box information from the crd/rst file generated from tLEaP. Instead you will have to specify the box information using the cellBasisVector flag. Rewrite the information from the last line of the rst file to this format. To check the last line you can use:

.. code-block:: bash

    tail -n 1 your_system.rst

.. code-block:: bash

   dx    dy    dz    90.0000000  90.0000000  90.0000000

.. code-block:: console

    # Periodic Cell 
    if {1} { 
    cellBasisVector1    dx         00.0000    00.0000 
    cellBasisVector2    00.0000    dy         00.0000 
    cellBasisVector3    00.0000    00.0000    dz 
    cellOrigin          00.0000    00.0000    00.0000
    }

    
.. caution::
    
    The method of setting periodic cell information for NAMD by directly utilizing the final line values from the rst file is only applicable to cubic periodic boxes. For non-cubic periodic box, such as truncated octahedron or rhombic dodecahedron, the values are determined by specific calculations using the length of the edge of the box for each kind. Refer to the section *using non-cubic periodic boxes* in the tutorial `Using the Amber force field in NAMD <https://ambermd.org/namd/namd_amber.html>`_.     


For ``em2.conf``:

1. Now, in the **input parameters** section, include the block that sets the end of the previous minimization stage as the starting point:

.. code-block:: console

    if {1} {                                    ;# This block checks for a restart files from a previous simulation (Minimization_01) 
    set inputname      yourMinimization_01      ;# Root name for restart files                                   
    binCoordinates     $inputname.restart.coor  ;# Read coordinate restart file 
    #binVelocities     $inputname.restart.vel   ;# Read velocity restart file 
    extendedSystem     $inputname.restart.xsc   ;# Read extended system restart file 
    }

The flag `binVelocities` is commented because the previous minimization stage does not assign velocity values. 

2. From now on, the `Periodic Cell` block at **simulation parameters** section will remain disabled (if {0}) because the box size will be read from the previous stage.

3. At this stage, the restraint block in the **additional parameters** section is disabled (if {0}) so that the minimization is carried out without any restraints. However, if you need to use restraints you can enable this section (if {1}).

    
1.7.3 Heating
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For ``heat.conf``, we will gradually bring the system up to the target temperature while maintaining control of the thermostat and barostat. This stage is performed with restraints on the backbone ``bb_restraints.pdb``.

1. In the **simulation parameters** section, the temperature and pressure control blocks are enabled (if {1}). Particularly the langevin and langevinPiston options, which control the thermostat and barostat, respectively. Change it as needed.

.. code-block:: console

    # Temperature control 
    if {1} {                     ;# If 0 don't read the block 
    langevin           on        ;# Turns off Langevin thermostat 
    langevinDamping    50        ;# Damping coefficient  gamma of 50/ps 
    langevinHydrogen   off       ;# Turns off Langevin dynamics for hydrogen atoms (if on) 
    langevinTemp       60        ;# Temperature target for Langevin dynamics (if on) 
    } 
          
    # Pressure  control 
    if {1} {                     ;# If 0 don't read the block                                                    
    useGroupPressure     no      ;# Disables group pressure control 
    useFlexibleCell      no      ;# Disables flexible cell (not for water boxes) 
    useConstantArea      no      ;# Disables constant area control (not for water boxes) 
    langevinPiston       on      ;# Enable Langevin piston for pressure control (off by default) 
    langevinPistonTarget 1.01325 ;# Target pressure for Langevin piston (in bar) (if on). 
    langevinPistonPeriod 200.    ;# Period of the Langevin piston (if on) 
    langevinPistonDecay  100.    ;# Decay time of the Langevin piston (if on)
    langevinPistonTemp   60      ;# Temperature target for Langevin piston (if on) 
    }

2. In the **execution instructions** section, there is a script that in `nSteps` progressively increases the temperature by adjusting the parameters of the barostat and thermostat. You can modify it to suit your desired temperature:

.. code-block:: console

    set Temp    300                               ;# Set temperature  target   
    set barostat 1                                ;# Set pressure target 
    set nSteps    600                             ;# Defines the number of simulation steps to run per temperature increment
    # for loop iterates through a temperature range 
    for {set t 60} {$t <= $Temp} {incr t} {run $nSteps;langevintemp $t;if {$barostat} {langevinpistontemp $t}}


1.7.4 Equilibration
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For ``eq1.conf``:

1. Enable the **additional parameters** section (if {1}) so that the equilibration is carried out with restraints to the whole protein. In order to accomplish this, the restriction file ``prot_restraints.pdb`` should be assigned to the consref and conskfile lines:

.. code-block:: console

    # Constraints of protein 
    if {1} {                                   ;# If 1 read the block 
    constraints           on                   ;# Turns on constraints 
    consref               your_restraints.pdb  ;# Reference PDB file for constraint positions 
    conskfile             your_restraints.pdb  ;# File containing constraint force constants 
    constraintScaling     1                    ;# Scaling factor for constraint forces 
    consexp               2                    ;# Exponent for constraint potential 
    conskcol              B                    ;# Column for constraint forces in the pdb file
    }



For ``eq2.conf``:

1. Enable the **additional parameters** section (if {1}) so that the equilibration is carried out with a small restraint on the backbone. In order to accomplish this, the restriction file ``bb_restraints.pdb`` should be assigned to the consref and conskfile lines and the constraintScaling should be set to 0.1:

.. code-block:: console

    # Constraints of protein 
    if {1} {                                   ;# If 1 read the block 
    constraints           on                   ;# Turns on constraints 
    consref               your_restraints.pdb  ;# Reference PDB file for constraint positions 
    conskfile             your_restraints.pdb  ;# File containing constraint force constants 
    constraintScaling     0.1                  ;# Scaling factor for constraint forces 
    consexp               2                    ;# Exponent for constraint potential 
    conskcol              B                    ;# Column for constraint forces in the pdb file
    }

1.7.5 Production
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For ``md.conf``:

1. At this stage, the restraint block in the **additional parameters** section is disabled (if {0}) so that the production stage is carried out without any restraints.






