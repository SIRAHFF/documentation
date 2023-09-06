This tutorial shows how to visualize and analize trajectory files using the **SIRAH Tools Plugin** for VMD. The main reference
for this tutorial is `SIRAH Tools <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_. The protein−DNA complex used here is the 3′ Repair Exonuclease 1 (TREX 1, PDB id: `5YWS <https://www.rcsb.org/structure/5YWS>`_) and its simulation was previously dicussed by `Klein et al. 2020 <https://doi.org/10.1021/acs.jcim.0c00160>`_. We strongly advise you to read these articles before starting the tutorial.

.. warning::

    This tutorial requires the use of **VMD**, **AmberTools**, and **R**. Please install these applications if you do not already have them. If you need to install any of these applications, please visit the following websites for installation instructions and additional information: `Install AmberTools <https://ambermd.org/GetAmber.php#ambertools>`_, `Installation VMD guide <https://www.ks.uiuc.edu/Research/vmd/current/ig/node6.html>`_ and the `The R Project for Statistical Computing <https://www.r-project.org/>`_ page.


.. important::

    In this tutorial, we will focus on the **sirah_vmdtk.tcl** plugin. The utilization of the SIRAH Tools scripts ``cgconv.pl`` to map AA to CG, and ``g_top2psf.pl`` to convert GROMACS' topologies to PSF files, is covered elsewhere in the documentation. We recommend that you review the other tutorial sections.



Visualization of CG systems
____________________________

One of the features of the SIRAH Tools plugin for VMD, **sirah_vmdtk.tcl**, is that it improves the visualization and analysis of CG trajectories. This plugin sets correct van der Waals radii, coloring code by atom and residue types and macros for selecting molecular components on SIRAH trajectories.

After processing the output trajectory to account for the Periodic Boundary Conditions (PBC) (see :ref:`AMBER <AMBER>`, :ref:`GROMACS <GROMACS>` or :ref:`NAMD <NAMD>` tutorials for examples of how you can do this), load the processed trajectory and the ``sirah_vmdtk.tcl`` file in VMD:

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

   **Figure 1.** SIRAH CG simulation loaded in VMD. A) The trajectory was imported without the **sirah_vmdtk.tcl** plugin. B) Trajectory imported using the **sirah_vmdtk.tcl** plugin. 
	

Structural analysis of CG systems
_________________________________