This tutorial shows how to visualize and analize trajectory files using the **SIRAH Tools Plugin** for VMD. The main reference
for this tutorial is `Machado & Pantano <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_. The protein−DNA complex used here is the 3′ Repair Exonuclease 1 (TREX 1, PDB id: `5YWS <https://www.rcsb.org/structure/5YWS>`_) and its simulation was previously dicussed by `Klein et al. 2020 <https://doi.org/10.1021/acs.jcim.0c00160>`_. We strongly advise you to read these articles before starting the tutorial.

.. warning::

    This tutorial requires the use of **VMD** and **AmberTools**. Please install these applications if you do not already have them. If you need to install any of these applications, please visit the following websites for installation instructions and additional information: `Install AmberTools <https://ambermd.org/GetAmber.php#ambertools>`_ and `Installation VMD guide <https://www.ks.uiuc.edu/Research/vmd/current/ig/node6.html>`_.


.. important::

    In this tutorial, we will focus on the **sirah_vmdtk.tcl** plugin. The utilization of the SIRAH Tools scripts ``cgconv.pl`` to map AA to CG, and ``g_top2psf.pl`` to convert GROMACS' topologies to PSF files, is covered elsewhere in the documentation. We recommend that you review the other tutorial sections.



Visualization of CG systems
____________________________

One of the features of the SIRAH Tools plugin for VMD, **sirah_vmdtk.tcl**, is that it improves the visualization and analysis of CG trajectories. This plugin sets the correct van der Waals radii, sets the correct coloring code by atom and residue types, and has macros for selecting molecular components on SIRAH trajectories.

After processing the output trajectory to account for the Periodic Boundary Conditions (PBC) (see :ref:`Amber <AMBER>`, :ref:`GROMACS <GROMACS>` or :ref:`NAMD <NAMD>` tutorials for examples of how you can do this), load the processed trajectory and the ``sirah_vmdtk.tcl`` file in VMD:

.. code-block:: bash

    vmd ../5YWS_cg.psf ../5YWS_cg.gro 5YWS_cg_md_pbc.nc -e ../sirah.amber/tools/sirah_vmdtk.tcl

.. tip::
	
	You can also load ``sirah_vmdtk.tcl`` inside VMD. Go to **Extensions** > **Tk Console** and enter:
	
	.. code-block:: bash

		source ../sirah.amber/tools/sirah_vmdtk.tcl
	
If ``sirah_vmdtk.tcl`` is not loaded in VMD, the trajectory will still be loaded, but the correct molecular connectivity, bead size, charges, etc. will not be generated (see Figure 1). 

.. figure:: /../images/sirah_tools_1.png
   :align: center
   :width: 100%

   **Figure 1.** SIRAH CG simulation loaded in VMD. A) The trajectory imported without the **sirah_vmdtk.tcl** plugin. B) Trajectory imported using the **sirah_vmdtk.tcl** plugin. 
	

SIRAH macros
~~~~~~~~~~~~~~~~

In VMD, a macro is a text that represents a selection. Macros are a useful feature of VMD when you use certain selections often. In the SIRAH Tools plugin for VMD, **sirah_vmdtk.tcl**, 10 macros are available:

* ``sirah_membrane`` - select only lipid residues;
* ``sirah_nucleic`` - select only acid nucleic residues;
* ``sirah_protein`` - select only protein residues;
* ``sirah_basic`` - select only basic amino acid residues;
* ``sirah_acidic`` - select only acidic amino acid residues;
* ``sirah_polar`` - select only polar amino acid residues;
* ``sirah_neutral`` - select only neutral amino acid residues;
* ``sirah_backbone`` - select the backbone beads of protein (GN, GC, and GO) or acid nucleic (PX, C5X, and O3');
* ``sirah_water`` - select all solvent available in the force field (WT4 and WLS);
* ``sirah_ions`` - select all ions available in the force field (KW, NaW, ClW, MgX, CaX, and Znx). 

To use any SIRAH macro, go to **Graphics** > **Representations** (See Figure 2) and click on the **Create Rep** button. In the **Selected Atoms** window, erase the text that appears there and type any of the available *Singlewords* above.

.. figure:: /../images/sirah_tools_2.png
   :align: center
   :width: 100%

   **Figure 2.** Our trajectory loaded in VMD with the default representation and without any selection.
   
.. tip::
	 To see all available VMD's *Singlewords*, click on the **Selections** tab in the **Graphics** > **Representations** window and look through the *Singlewords* panel. You can also double-click on a singleword, and then, click on the **Apply** button to use a VMD Macro.
	
Since our system is a protein−DNA complex with coordinating magnesium ions, we can use four macros: ``sirah_protein`` (see Figure 3A), ``sirah_nucleic`` (see Figure 3B), ``sirah_water`` (see Figure 3C), and ``sirah_ions`` (see Figure 3D).

.. figure:: /../images/sirah_tools_3.png
   :align: center
   :width: 100%

   **Figure 3.** Each of the macros of our systems shown as a single image. A) The ``sirah_protein`` macro using the **Licorice** as **Drawing Method**. The backbone beads of the protein are pink, and the sidechains ones are blue when using the **Coloring Method** Name. B) The ``sirah_nucleic`` macro using the *Licorice* as **Drawing Method**. For the DNA, the backbone beads PX, C5X, and O3' are colored dark yellow, cyan, and red, respectively, and the sidechain ones are colored red and blue using the **Coloring Method** Name. C) The ``sirah_water`` macro using the *CPK* as **Drawing Method**. In this system, only WT4 was used and is colored green by the **Coloring Method** Name. D) The ``sirah_ions`` macro using the *VDW* as **Drawing Method**. The ClW, NaW, and MgW ions are colored cyan, blue, and purple, respectively, when using the **Coloring Method** Name.


To use all representations together, just create different representations for each one (clicking on the **Create Rep** button in the **Graphics** > **Representations** window). Your protein-DNA complex should now look similar to the one in Figure 4.

.. figure:: /../images/sirah_tools_4.png
   :align: center
   :width: 100%

   **Figure 4.** Our protein-DNA complex after using the SIRAH macros ``sirah_backbone``, ``sirah_water``, and ``sirah_ions``. To improve visualization, we changed the **Material** of the water beads to Ghost and we hid some elements of the system.
   
.. tip:: 
	As you can see, macros can be very useful and when saving your work in a saved state (**File** > **Save Visualization State**), macros are included in the saved state file to be used later (**File** > **Load Visualization State**).


Structural analysis of CG systems
__________________________________


Secondary structural analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Backmapping analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
