.. note::

   Please report bugs, errors or enhancement requests through `Issue Tracker <https://github.com/SIRAHFF/documentation/issues>`_ or if you have a question about SIRAH open a `New Discussion <https://github.com/SIRAHFF/documentation/discussions>`_.
   
This tutorial shows how to apply the multiscale solvation approach of SIRAH force field in Steered Molecular Dynamics (SMD) simulations. In this tutorial, an immunoglobulin domain known as `I10 <https://royalsocietypublishing.org/doi/10.1098/rsob.160114>`_ derived from Titin is enclosed by CG waters referred to as WatFour (WT4). These waters are further embedded within supra-coarse-grained (SCG) molecules called WatElse (WLS), which serve as a representation of bulk water. The main references for this tutorial are: `Darré et al. <https://pubs.acs.org/doi/abs/10.1021/ct100379f>`_ and `Machado et al. <https://doi.org/10.1021/acs.jctc.7b00659>`__.

.. note::

	We strongly advise you to read and complete :ref:`Tutorial 2 <Tutorial 2 G>` and :ref:`Tutorial 3 <Tutorial 3 G>` before starting. We also recommend you to perform the `Umbrella Sampling <htthttp://www.mdtutorials.com/gmx/umbrella/index.html>`__ tutorial of GROMACS to get familiar with pulling simulations parameters.
	

.. important::

    Check the :ref:`Setting up SIRAH <download gromacs>` section for download and set up details before starting this tutorial.
    Since this is **Tutorial 7**, remember to replace ``X.X``, the files corresponding to this tutorial can be found in: ``sirah.ff/tutorial/7/``
	

7.1. Setting pulling direction
________________________________

.. caution::
  
  The basic idea behind SMD is to apply an external force or mechanical perturbation to a specific part of the biomolecular system under study and observe its response. This external force is applied to one atom or groups of atoms as a "steering" force, guiding or pulling the atoms along a predefined path or direction. Thus, prior to proceeding, we need to remember to prepare the PDB file to specify the pulling direction. In this example, the **Z-axis** was employed to orient the protein terminus and to determine the box size.

In ``sirah_[version].ff/tutorial/7/`` you will find the I10 domain (PDB code: `4QEG <https://www.rcsb.org/structure/4QEG>`__) already aligned along the Z-axis, ``I10_z.pdb`` file. 

.. figure:: /../images/TutorialSMD1.png
   :align: center
   :width: 100%

   **Figure 1.** I10 domain N-terminal (blue) and C-terminal (red) aligned to the Z-axis from ``I10_z.pdb``.

To get the protein aligned along the Z-axis, we used the following commands in VMD's *Tk Console* (*Extensions* > *Tk Console*):

.. code-block:: console
   
  #Select the protein 
  set all [atomselect top "protein"]

  #Change protein coordinates to align to its center 
  $all moveby  [vecinvert [measure center $all]]

  #Rotate around x-axis 50 degrees 
  $all move [transaxis x 50]

  #Rotate around y-axis -11 degrees
  $all move [transaxis y -11]

  #Write PDB file with the protein aligned to z-axis
  $all writepdb I10_z.pdb

.. note::

  For setting up your own system you can open your own PDB in VMD and probe different alternatives to get the correct direction.



7.2. Build CG representations
______________________________

.. caution::

  The mapping to CG requires the correct protonation state of each residue at a given pH. We recommend using the `CHARMM-GUI server <https://www.charmm-gui.org/>`_ and use the **PDB Reader & Manipulator** to prepare your system. An account is required to access any of the CHARMM-GUI Input Generator modules, and it can take up to 24 hours to obtain one. 
  
  Other option is the `PDB2PQR server <https://server.poissonboltzmann.org/pdb2pqr>`_ and choosing the output naming scheme of AMBER for best compatibility. This server was utilized to generate the *PQR* file featured in this tutorial. Be aware that modified residues lacking parameters such as: MSE (seleno MET), TPO (phosphorylated THY), SEP (phosphorylated SER) or others are deleted from the PQR file by the server. In that case, mutate the residues to their unmodified form before submitting the structure to the server.

  See :ref:`Tutorial 3 <Tutorial 3 G>` for cautions while preparing and mapping atomistic proteins to SIRAH.
	
Map the atomistic structure of the I10 domain to its CG representation:  

.. code-block:: bash

  ./sirah.ff/tools/CGCONV/cgconv.pl -i sirah.ff/tutorial/7/I10_z.pqr -o I10_cg.pdb 
  
The input file ``-i`` I10_z.pqr contains the atomistic representation of `4QEG <https://www.rcsb.org/structure/4QEG>`__ structure at pH **7.0** and aligned to the **Z-axis**, while the output ``-o`` I10_cg.pdb is its SIRAH CG representation. 

.. tip::

  This is the basic usage of the script **cgconv.pl**, you can learn other capabilities from its help by typing:

  .. code-block:: bash

    ./sirah.ff/tools/CGCONV/cgconv.pl -h


Please check both PDB and PQR structures using VMD:

.. code-block:: bash 
  
  vmd -m sirah.ff/tutorial/7/I10_z.pqr I10_cg.pdb

From now on it is just normal GROMACS stuff!


7.3. PDB to GROMACS format
__________________________

Use ``pdb2gmx`` to convert your PDB file into GROMACS format: 

.. code-block:: bash

  gmx pdb2gmx -f I10_cg.pdb -o I10_cg.gro

When prompted, choose *SIRAH force field* and then *SIRAH solvent models*.
In this specific case, the charge of the system is -5.000 e; this will be addressed later.

.. note:: 

  By default charged terminal are used but it is possible to set them neutral with option ``-ter``

.. note::

  Warning messages about long, triangular or square bonds are fine and expected due to the CG topology of some residues.

.. caution::

  However, missing atom messages are errors which probably trace back to the
  mapping step. In that case, check your atomistic and mapped structures and do not carry on the
  simulation until the problem is solved.


7.4. Solvate the system
_______________________

.. danger::

  Since we are doing a SMD simulation, we need to carefully set the box dimensions to provide enough solvent to the "stretched" protein and keep a good separation between the surfaces of periodic images.

  In this system, I10 has 88 amino acids and the maximum extension size of each amino acid is 0.34 nm. Thus, "stretched" protein will be 88 * 0.34 = 29.92 nm (~ 30.0 nm) long. In order to accommodate the pulling, GROMACS stipulates a minimum box size double this value, i.e. 60 nm for the Z-axis. However, for optimal results, it is recommended that the dimensions of the box be 2.5 to 3 times greater than the maximum length of the protein when in its extended conformation. Therefore, for this tutorial the box used is 10 10 90 nm.

  .. figure:: /../images/TutorialSMD2.png
   :align: center
   :width: 100%

   **Figure 2.** Dimensions of the multiscale solvation box used in this tutorial.

In order to have a multiscale solvent approach using WT4 and WLS, two steps are needed to solvate the systems. 

First, define the simulation region of the system to be enclosed by WT4 (pink in **Figure 2**)

.. code-block:: bash 
  
  gmx editconf -f I10_cg.gro -o I10_cg_box.gro -box 10 10 20 -bt triclinic -c

.. note::

  At this step, if you don't want to use a multiscale solvent method, the whole box dimension (10 10 90 nm) can be used to add only WT4 molecules. 

Add WT4 molecules:

.. code-block:: bash 

  gmx solvate -cp I10_cg_box.gro -cs sirah.ff/wt416.gro -o I10_cg_solv1.gro

.. note:: 

  Before GROMACS version 5.x, the command *gmx solvate* was called *genbox*.

Edit the [ molecules ] section in ``topol.top`` to include the number of added WT4 molecules:

.. list-table::
   :align: center
   :widths: 50 50
   :header-rows: 1

   * - Topology before editing
     - Topology after editing
   * - | [ molecules ] 
       | ; Compound #mols 
       | Protein_chain_A    1    
       |  

    
     - | [ molecules ] 
       | ; Compound #mols
       | Protein_chain_A    1    
       | WT4             6281

.. hint::
  
  If you forget to read the number of added WT4 molecules from the output of *solvate*, then use the following command line to get it 

  .. code-block:: console

    grep -c WP1 I10_cg_solv1.gro

.. caution::
  
  The number of added WT4 molecules, **6281**, may change according to the software version.

Remove misplaced WT4 molecules within 0.3 nm of protein:

.. code-block:: bash
  
  echo q | gmx make_ndx -f I10_cg_sol1.gro -o I10_cg_sol1.ndx

.. code-block:: bash 

  gmx grompp -f sirah.ff/tutorial/7/GPU/em1_CGPROT.mdp -p topol.top -po delete1.mdp -c I10_cg_sol1.gro -o I10_cg_sol1.tpr -maxwarn 2

.. caution::
  
  New GROMACS versions may complain about the non-neutral charge of the system, aborting the generation of the TPR file by command grompp. We will neutralize the system later, so to overcame this issue, just allow warning messages by adding the following keyword to the grompp command line: ``-maxwarn 2``

.. code-block:: bash 

  gmx select -f I10_cg_sol1.gro -s I10_cg_sol1.tpr -n I10_cg_sol1.ndx -on rm_close_wt4.ndx -select 'not (same residue as (resname WT4 and within 0.3 of group Protein))'

.. code-block:: bash 

  gmx editconf -f I10_cg_sol1.gro -o I10_cg_sol2.gro -n rm_close_wt4.ndx


Edit the [ molecules ] section in ``topol.top`` to correct the number of WT4 molecules, **6261**.

.. hint::
  
  If you forget to read the number of added WT4 molecules from the output of *solvate*, then use the following command line to get it 

  .. code-block:: console

    grep -c WP1 I10_cg_solv2.gro


Now, we include the second solvent layer of solvent with WLS molecules (green in **Figure 2**):

.. code-block:: bash 
  
  gmx editconf -f I10_cg_sol2.gro -o I10_cg_box2.gro -box 10 10 90 -bt triclinic -c

.. hint::

  We can check the final box dimensions with VMD:

  .. code-block:: bash

    vmd I10_cg_box2.gro

  In the *Tk Console* (*Extensions* > *Tk Console*) use the command:

  .. code-block:: bash

   pbc box


Add WLS molecules:

.. code-block:: bash 

  gmx solvate -cp I10_cg_box2.gro -cs sirah.ff/wlsbox.gro -o I10_cg_sol3.gro

.. note:: 

  Before GROMACS version 5.x, the command *gmx solvate* was called *genbox*.


Edit the [ molecules ] section in ``topol.top`` to include the number of added WLS molecules:

.. list-table::
   :align: center
   :widths: 50 50
   :header-rows: 1

   * - Topology before editing
     - Topology after editing
   * - | [ molecules ] 
       | ; Compound #mols 
       | Protein_chain_A    1    
       | WT4             6261
       |   

    
     - | [ molecules ] 
       | ; Compound #mols
       | Protein_chain_A    1    
       | WT4             6261
       | WLS             4697

.. hint::
  
  If you forget to read the number of added WLS molecules from the output of *solvate*, then use the following command line to get it 

  .. code-block:: console

    grep -c LN1 I10_cg_solv3.gro

.. caution::
  
  The number of added WLS molecules, **4697**, may change according to the software version.

Remove misplaced WLS molecules within 7.8 nm of protein:

.. code-block:: bash
  
  echo q | gmx make_ndx -f I10_cg_sol3.gro -o I10_cg_sol3.ndx

.. code-block:: bash 

  gmx grompp -f sirah.ff/tutorial/7/GPU/em1_CGPROT.mdp -p topol.top -po delete3.mdp -c I10_cg_sol3.gro -o I10_cg_sol3.tpr -maxwarn 2

.. caution::
  
  New GROMACS versions may complain about the non-neutral charge of the system, aborting the generation of the TPR file by command grompp. We will neutralize the system later, so to overcame this issue, just allow warning messages by adding the following keyword to the grompp command line: ``-maxwarn 2``

.. code-block:: bash 

  gmx select -f I10_cg_sol3.gro -s I10_cg_sol3.tpr -n I10_cg_sol3.ndx -on rm_close_wls.ndx -select 'not (same residue as (resname WLS and within 7.8 of group Protein))'

.. code-block:: bash 

  gmx editconf -f I10_cg_sol3.gro -o I10_cg_sol4.gro -n rm_close_wls.ndx

Edit the [ molecules ] section in ``topol.top`` to correct the number of WLS molecules, **4582**.

.. hint::
  
  If you forget to read the number of added WLS molecules from the output of *solvate*, then use the following command line to get it 

  .. code-block:: console

    grep -c LN1 I10_cg_solv4.gro


.. note::
  
  Consult ``sirah.ff/0ISSUES`` and :doc:`FAQs <../FAQ>` for information on known solvation issues.


Add CG counterions and 0.15M NaCl:

.. code-block:: bash

  gmx grompp -f sirah.ff/tutorial/7/GPU/em1_CGPROT.mdp -p topol.top -po delete4.mdp -c I10_cg_sol4.gro -o I10_cg_sol4.tpr -maxwarn 3

.. caution::
  
  New GROMACS versions may complain about the non-neutral charge of the system, aborting the generation of the TPR file by command grompp. We are about to neutralize the system, so to overcame this issue, just allow warning messages by adding the following keyword to the grompp command line: ``-maxwarn 3``

.. code-block:: bash

  gmx genion -s I10_cg_sol4.tpr -o I10_cg_ion.gro -np 187 -pname NaW -nn 182 -nname ClW

When prompted, choose to substitute *WT4* molecules by *ions*.

.. note:: 

  The available electrolyte species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW) which represent solvated ions in solution. One ion pair (e.g., NaW-ClW) each 34 WT4 molecules results in a salt concentration of ~0.15M (see :ref:`Appendix <Appendix>` for details). Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`_.

Edit the [ molecules ] section in ``topol.top`` to include the CG ions and the correct number of WT4, WLS, and ions.

.. list-table::
   :align: center
   :widths: 50 50
   :header-rows: 1

   * - Topology before editing
     - Topology after editing
   * - | [ molecules ] 
       | ; Compound #mols 
       | Protein_chain_A    1    
       | WT4             6261
       | WLS             4582
       |
       |  

    
     - | [ molecules ] 
       | ; Compound #mols
       | Protein_chain_A    1    
       | WT4             5892
       | NaW              187
       | ClW              182
       | WLS             4582

.. caution::

  Following the above order is important: the number of WT4 comes first, then the number of ions, and finally WLS.

Before running the simulation it may be a good idea to visualize your molecular system. CG molecules are not recognized by molecular visualizers and will not display correctly. To fix this problem you may
generate a PSF file of the system using the script ``g_top2psf.pl``:

.. code-block:: bash

  ./sirah.ff/tools/g_top2psf.pl -i topol.top -o I10_cg_ion.psf

.. note::

  This is the basic usage of the script ``g_top2psf.pl``, you can learn other capabilities from its help:
  
  .. code-block:: bash

    ./sirah.ff/tools/g_top2psf.pl -h


Use VMD to check how the CG system looks like:

.. code-block::

  vmd I10_cg_ion.psf I10_cg_ion.gro -e sirah.ff/tools/sirah_vmdtk.tcl

.. tip::

  VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right ones, according to the CG representation. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
  Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.


To achive a proper interaction between the protein and solvent, we will perform a equilibration step applying restraints over the protein backbone.

Create an index file including a group for the backbone GN and GO beads:

.. code-block:: bash 

  echo -e "a GN GO\n\nq" | gmx make_ndx -f I10_cg_ion.gro -o I10_cg_ion.ndx

.. note::

  WT4 and CG ions (NaW, ClW) are automatically set to the group *SIRAH-Solvent*.

Generate restraint files for the backbone *GN* and *GO* beads:

.. code-block:: bash

 gmx genrestr -f I10_cg.gro -n I10_cg_ion.ndx -o bkbres.itp

.. code-block:: bash

  gmx genrestr -f I10_cg.gro -n I10_cg_ion.ndx -o bkbres_soft.itp -fc 100 100 100

When prompted, choose the group *GN_GO*

Add the restraints to ``topol.top``:

.. list-table:: 
   :align: center
   :widths: 50 50
   :header-rows: 1

   * - Topology before editing
     - Topology after editing
   * - | ; Include Position restraint file
       | #ifdef POSRES
       | #include "posre.itp"
       | #endif
       
       | ; Include water topology
       | #include"../sirah.ff/solv.itp" 
       |
       
       | #ifdef POSRES_WATER       
       | ; Position restraint for each water oxygen
       | [ position_restraints ]
       | ;  i funct       fcx        fcy        fcz
       |    1    1       1000       1000       1000
       | #endif

       |
       |

       |
       |
       |
       |
       |
       |

       |
       |
       |
       |
       
                  
     - | ; Include Position restraint file
       | #ifdef POSRES
       | #include "posre.itp"
       | #endif
        
       | ; Backbone restraints
       | #ifdef GN_GO
       | #include "bkbres.itp"
       | #endif
     
       | ; Backbone soft restrains
       | #ifdef GN_GO_SOFT
       | #include "bkbres_soft.itp"
       | #endif

       | ; Include water topology
       | #include"../sirah.ff/solv.itp"

       | #ifdef POSRES_WATER       
       | ; Position restraint for each water oxygen
       | [ position_restraints ]
       | ;  i funct       fcx        fcy        fcz
       |    1    1       1000       1000       1000
       | #endif

       | ; Solvent restrains
       | #ifdef Sirah_solvent
       | #include "posre_sirah_solvent.itp"
       | #endif
       |


7.5. Run the simulation
________________________

.. important:: 

  By default in this tutorial we will use input files for GROMACS on GPU (``sirah.ff/tutorial/7/GPU``). 

The folder ``sirah.ff/tutorial/7/GPU/`` contains typical input files for energy minimization (``em1_CGPROT.mdp``, ``em2_CGPROT.mdp``, and ``em3_CGPROT.mdp``), equilibration (``eq1_CGPROT.mdp`` and ``eq2_CGPROT.mdp``), production (``md_CGPROT.mdp``) and SMD (``SMD_Force_CGPROT.mdp`` and ``SMD_Velocity_CGPROT.mdp``) runs. Please check carefully the input flags therein.

Make a new folder for the run:

.. code-block:: bash

  mkdir run; cd run

**Energy Minimization of side chains by restraining the backbone and Sirah-solvent**:

.. code-block:: bash

  gmx grompp -f ../sirah.ff/tutorial/7/GPU/em1_CGPROT.mdp -p ../topol.top -po em1.mdp -n ../I10_cg_ion.ndx -c ../I10_cg_ion.gro -r ../I10_cg_ion.gro -o I10_cg_em1.tpr

.. code-block:: bash

  gmx mdrun -deffnm I10_cg_em1 &> EM1.log &

**Energy Minimization of side chains by restraining the backbone**:

.. code-block:: bash

  gmx grompp -f ../sirah.ff/tutorial/7/GPU/em2_CGPROT.mdp -p ../topol.top -po em2.mdp -n ../I10_cg_ion.ndx -c I10_cg_em1.gro -o I10_cg_em2.tpr 

.. code-block:: bash

  gmx mdrun -deffnm I10_cg_em2 &> EM2.log & 

**Energy Minimization of whole system**:

.. code-block:: bash

  gmx grompp -f ../sirah.ff/tutorial/7/GPU/em3_CGPROT.mdp -p ../topol.top -po em3.mdp -n ../I10_cg_ion.ndx -c I10_cg_em2.gro -o I10_cg_em3.tpr 

.. code-block:: bash

  gmx mdrun -deffnm I10_cg_em3 &> EM3.log & 

**Solvent equilibration**:

.. code-block:: bash

  gmx grompp -f ../sirah.ff/tutorial/7/GPU/eq1_CGPROT.mdp -p ../topol.top -po eq1.mdp -n ./I10_cg_ion.ndx -c I10_cg_em2.gro -r I10_cg_em2.gro -o I10_cg_eq1.tpr 

.. code-block:: bash

  gmx mdrun -deffnm I10_cg_eq1 &> EQ1.log & 

**Soft equilibration to improve side chain solvation**:

.. code-block:: bash

  gmx grompp -f ../sirah.ff/tutorial/7/GPU/eq2_CGPROT.mdp -p ../topol.top -po eq2.mdp -n ../I10_cg_ion.ndx -c I10_cg_eq1.gro -r I10_cg_eq1.gro -o I10_cg_eq2.tpr

.. code-block:: bash

  gmx mdrun -deffnm I10_cg_eq2 &> EQ2.log & 


**SMD Force or velocity**: 

Here, we need to modify the index file to add the "Pull" groups, before running SMD force or velocity simulations. 

Copy the ``I10_cg_ion.ndx`` with a new name:

.. code-block:: bash

  cp ../I10_cg_ion.ndx ../I10_cg_ion_pull.ndx

Open the ``I10_cg_ion_pull.ndx`` in any text editor and manually add to the end of the ``I10_cg_ion_pull.ndx`` file these two new groups:

+-----------------+
| Groups to add   |
+=================+
| | [ pull1 ]     | 
| | 4             |
| | [ pull2 ]     | 
| | 437           |
+-----------------+

In this tutorial we are going to run only the **SMD Force** simulation:

.. code-block:: bash

  gmx grompp -f ../sirah.ff/tutorial/7/GPU/SMD_Force_CGPROT.mdp -p ../topol.top -po md.mdp -n ../I10_cg_ion_pull.ndx -c I10_cg_eq2.gro -o I10_cg_SMD_F.tpr

.. code-block:: bash
  
  gmx mdrun -deffnm I10_cg_SMD_F &> SMD_F.log &

However, you can also run a **SMD Velocity** simulation:

.. code-block:: bash

  gmx grompp -f ../sirah.ff/tutorial/7/GPU/SMD_Velocity_CGPROT.mdp -p ../topol.top -po md.mdp -n ../I10_cg_ion_pull.ndx -c I10_cg_eq2.gro -o I10_cg_SMD_V.tpr

.. code-block:: bash

  gmx mdrun -deffnm I10_cg_SMD_V &> SMD_V.log &

.. note::

  It is also possible to use the files you made to run a MD simulation of your system: 

  **Production (1000ns)**:

  .. code-block:: bash

    gmx grompp -f ../sirah.ff/tutorial/7/GPU/md_CGPROT.mdp -p ../topol.top -po md.mdp -n ../I10_cg_ion.ndx -c I10_cg_eq2.gro -o I10_cg_md.tpr

  .. code-block:: bash

    gmx mdrun -deffnm I10_cg_md &> MD.log &


7.6. Visualizing the simulation
________________________________

That’s it! Now you can analyze the trajectory.

GROMACS automatically creates plot results to the SMD simulations. The position versus simulation time is saved as ``I10_cg_md_pullx.xvg``, and the force versus simulation time is saved as ``I10_cg_md_pullf.xvg``. 

You can plot the results using Xmgrace.

.. code-block:: bash
  
  xmgrace I10_cg_md_pullx.xvg

For ``I10_cg_md_pullx.xvg``, a plot similar to **Figure 3** will appear:

.. figure:: /../images/Fuerza_100pN.png
   :align: center
   :width: 100%

   **Figure 3.** Position vs Time plot created by GROMACS from the SMD simulation.


In addition, you can process the output trajectory at folder ``run/`` to account for the Periodic Boundary Conditions (PBC).

For **SMD Force**:

.. code-block:: bash

  gmx trjconv -s I10_cg_em1.tpr -f I10_cg_SMD_F.xtc -o I10_cg_SMD_F_pbc.xtc -n ../I10_cg_ion_pull.ndx -ur compact -center -pbc mol

For **SMD Velocity**:

.. code-block:: bash
  
  gmx trjconv -s I10_cg_em1.tpr -f I10_cg_SMD_V.xtc -o I10_cg_SMD_V_pbc.xtc -n ../I10_cg_ion_pull.ndx -ur compact -center -pbc mol

When prompted, choose *Protein* for centering and *System* for output.

.. note::

  If you had also run a MD simulation, you could use the following commands to account for PBC:

  .. code-block:: bash

    gmx trjconv -s I10_cg_em1.tpr -f I10_cg_md.xtc -o I10_cg_md_pbc.xtc -n ../I10_cg_ion.ndx -ur compact -center -pbc mol 


Now you can check the simulation using VMD:

.. code-block:: bash

  vmd ../I10_cg_ion.psf ../I10_cg_ion.gro I10_cg_SMD_F_pbc.xtc -e ../sirah.ff/tools/sirah_vmdtk.tcl

.. note::
    
    The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. Use the command ``sirah-help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.



