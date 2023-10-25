This tutorial shows how to perform *Root Mean Square Deviation* (RMSD), *Root Mean Square Fluctuation* (RMSF), *Radius of gyration* (Rg), and *Solvent-Accessible Aurface Area* (SASA) analyses on a SIRAH CG simulation trajectory file using VMD. The protein−DNA complex used here is the 3′ Repair Exonuclease 1 (TREX 1, PDB id: `5YWS <https://www.rcsb.org/structure/5YWS>`_) and its simulation was previously dicussed by `Klein et al. <https://doi.org/10.1021/acs.jcim.0c00160>`__ (check the :ref:`SIRAH Tools tutorial <SIRAH tools>` for additional analyses on this system). We strongly advise you to read the article before starting the tutorial.


Since VMD is a powerful molecular visualization software, you may find it convenient and familiar to perform trajectory analysis directly in VMD. In the following sections, we will concentrate on the practical side of these analyses in VMD and will not delve into the theory. 

.. warning::

    This tutorial requires the use of **VMD** and **Xmgrace**. Install these applications if you do not already have them. If you need to install these application, visit the following websites for installation instructions and additional information: `VMD installation  guide <https://www.ks.uiuc.edu/Research/vmd/current/ig/node6.html>`_ and `Xmgrace installation guide <https://plasma-gate.weizmann.ac.il/Grace/doc/UsersGuide.html#s2>`__.


.. important::

    In this tutorial, we used **VMD version 1.9.4**. A few details may vary depending on the software version and operating system.


Loading the Trajectory
-----------------------

After processing the output trajectory to account for the Periodic Boundary Conditions (PBC) (see :ref:`Amber <AMBER>`, :ref:`GROMACS <GROMACS>` or :ref:`NAMD <NAMD>` tutorials for examples of how you can do this), load the processed trajectory and the ``sirah_vmdtk.tcl`` file in VMD:

.. code-block:: bash

    vmd ../5YWS_cg.psf 5YWS_cg_md_pbc.xtc -e ../sirah.ff/tools/sirah_vmdtk.tcl

.. tip::
   
   You can also load ``sirah_vmdtk.tcl`` inside VMD. Go to *Extensions* > *Tk Console* and enter:
   
   **For AMBER**

   .. code-block:: bash

      source ../sirah.amber/tools/sirah_vmdtk.tcl
   
   **For GROMACS**
   
   .. code-block:: bash

      source ../sirah.ff/tools/sirah_vmdtk.tcl

.. attention::

    The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. These macros will facilite our atom selection within VMD. You can go to the :ref:`SIRAH Tools tutorial <SIRAH tools>` to learn more about the macros.


That’s it! Now you can analyze the trajectory.

RMSD
--------

The *Root Mean Square Deviation* (RMSD) is one of the most common analyses performed when analyzing a MD trajectory. Without delving into the theoretical aspects, it can be stated that the *RMSD* serves as a metric for assessing the degree of structural resemblance across different molecular conformations within a given system along a trajectory. 

To perform this analysis you can take advantage of the **RMSD Visualizer Tool**, a built-in plugin in VMD designed to compute *RMSD*, or using a **Tcl script**.


1.1 RMSD Visualizer Tool
__________________________

Once your trajectory is loaded, go to *Extensions* > *Analysis* > *RMSD Visualizer Tool*. It will open the interface dedicated to *RMSD* analysis.

In the *RMSD Visualizer Tool* interface, select **Top** in the *Molecule* selection. To analyze the protein backbone, enter **sirah_backbone and sirah_protein** in the *Atom selection* box (see **Figure 1A**).

.. tip::

   If you wish to focus on a specific bead, you can enter their name in the *Atom selection* box. For example, if we want something similar to a "carbon alpha" *RMSD*, we can type **name GC**.

Click on the ``ALIGN`` button. This will superimpose each frame of the trajectory to the reference frame (in this case the first frame, frame 0) based on the selected groups of atoms to minimize *RMSD*. This step is not required but is recommended to display only the differences that arise from structure fluctuations and not from the displacements and rotations of the molecule as a whole.  

.. tip::

   Every selection need to be aligned before calculating RMSD.

.. note::
   
   Since the VMD default *Reference* selection (*Molecule ID*: **self** and *Frame*: **0**) was used, all the atoms of the selected molecule will be rotated and translated to fit the structure of the trajectory first frame (frame number 0). But you can modify this by setting a reference molecule. It is important to remember, when doing an alignment using another molecule as the reference, the selections for both molecules need to have the exact same number of atoms.

With all the necessary settings in place, click on the ``RMSD`` button. The tool will perform the analysis, computing the *RMSD* values for the selected region over the trajectory. After the calculation is complete, click on the ``Plot result`` button and a new window will open with a plot showcasing the computed *RMSD* versus frames (see **Figure 1B**). 

.. figure:: /../images/sirah_analysis_1.png
   :align: center
   :width: 100%

   **Figure 1.** Protein *RMSD* of a SIRAH CG simulation using the *RMSD Visualizer Tool* of VMD. A) *RMSD Visualizer Tool* options used to perfom *RMSD*. B) The result *RMSD vs Frame* of the **sirah_backbone and sirah_protein** selection. 

.. tip::

   If your trajectory is too big, you can change the *Step size* parameter in the *Trajectory* section to skip frames. However, alignment does not accept a *Step size*.


Now, let's calculate the *RMSD* for the DNA molecule. In the *Atom selection* box, type in **sirah_backbone and sirah_nucleic**. Then, click on the ``ALIGN`` button and then on the ``RMSD`` button. A new *RMSD* line will appear in the results box (see **Figure 2A**). Select it and click on the ``Plot result`` button to open the *RMSD* versus frame plot window (see **Figure 2B**).

.. figure:: /../images/sirah_analysis_2.png
   :align: center
   :width: 100%

   **Figure 2.** DNA *RMSD* of a SIRAH CG simulation using the *RMSD Visualizer Tool* of VMD. A) *RMSD Visualizer Tool* options used to perfom *RMSD*. B) The result *RMSD vs Frame* of the **sirah_backbone and sirah_nucleic** selection. 

With the *RMSD Visualizer Tool*, you can also choose to plot the *RMSD* of both protein and DNA in the same window by selecting both calculations and clicking on the ``Plot result`` button (see **Figure 3**).

.. figure:: /../images/sirah_analysis_3.png
   :align: center
   :width: 100%

   **Figure 3.** Protein (blue) and DNA (red) *RMSD* of a SIRAH CG simulation using the *RMSD Visualizer Tool* of VMD. 

.. tip::
 
   If you prefer to use external tools for visualization plots, you can save the data in a text file (on the plot window go to *File* > *Export to Xmgrace* or *File* > *Export to ASCII vectors*) and then import the data into your preferred plotting software.


1.2. Tcl script
________________

The same actions can be taken on the scripting level using the *Tk Console*. Thus, you can create a ``rmsd_protein.tcl`` file to calculate *RMSD* of the protein and output a text file ``rmsd_prot.dat``:

.. code-block:: console
   
   # set output file name 
   set outfile [open rmsd_prot.dat w];

   # set reference as the first frame using protein backbone as selection
   set reference [atomselect top "sirah_backbone and sirah_protein" frame 0]
   # set trajectory selection also as the protein backbone
   set compare [atomselect top "sirah_backbone and sirah_protein"]

   # get the number of frames
   set N [molinfo top get numframes]

   # calculate RMSD for all frames
   for {set i 0} {$i < $N} {incr i} {
      # get the correct frame
      $compare frame $i
      # do the alignment
      $compare move [measure fit $compare $reference]
      # compute the RMSD
      set rmsd [measure rmsd $compare $reference]

      # print the RMSD in the output file
      puts $outfile "$i \t $rmsd"
   }
   close $outfile

With the ``rmsd_protein.tcl`` file in your work directory, go to *Extensions* > *Tk Console* and enter:
   
.. code-block:: bash

   source rmsd_protein.tcl
   
You can use Xmgrace to plot the result:

.. code-block:: bash

   xmgrace rmsd_prot.dat

You can create a ``rmsd_nucleic.tcl`` file by changing the following lines:
   
* ``set outfile [open rmsd_prot.dat w]`` to ``set outfile [open rmsd_nucl.dat w]``;
* ``set reference [atomselect top "sirah_backbone and sirah_protein" frame 0]`` to ``set reference [atomselect top "sirah_backbone and sirah_nucleic" frame 0]``;
* ``set compare [atomselect top "sirah_backbone and sirah_protein"]`` to ``set compare [atomselect top "sirah_backbone and sirah_nucleic"]``.

You can use Xmgrace to plot the results:

.. code-block:: bash

   xmgrace rmsd_prot.dat rmsd_nucl.dat


RMSF
---------

The *Root Mean Square Fluctuation* (RMSF) of a structure is the time average of the *RMSD* per residue. In contrast to the *RMSD*, which quantifies how much a structure deviates from a reference over time, the *RMSF* can disclose which system components are the most mobile. To perform this analysis, we are using a **Tcl script** directly using the *Tk Console*. 

Thus, you can create a ``rmsf_protein.tcl`` file to calculate *RMSF* of the protein and output a text file ``rmsf_prot.dat``:

.. code-block:: console
   
   # set output file name
   set outfile [open rmsf_prot.dat w];

   # set reference and selection of protein
   set reference [atomselect top "sirah_protein and name GC" frame 0]
   set sel [atomselect top "sirah_protein and name GC"]

   # get the number of frames
   set N [molinfo top get numframes]

   #do the alignment
   for {set i 0} {$i < $N} {incr i} {
      # get the correct frame
      $sel frame $i
      # do the alignment
      $sel move [measure fit $sel $reference]
   }

   # calculate rmsf for all trajectory frames
   set rmsf [measure rmsf $sel first 0 last -1 step 1]

   # print to file the rmsf by residue
   for {set i 0} {$i < [$sel num]} {incr i} {
      puts $outfile "[expr {$i+1}] [lindex $rmsf $i]"
   } 
   close $outfile

.. tip::

   If your trajectory is too big, you can change the ``step`` parameter in the ``set rmsf [measure rmsf $sel first 0 last -1 step 1]`` line to skip frames. You can also change ``first`` and ``last`` parameters if you have a frame range.    

With the ``rmsf_protein.tcl`` file in your work directory, go to *Extensions* > *Tk Console* and enter:
   
.. code-block:: bash

   source rmsf_protein.tcl
   
You can use Xmgrace to plot the result and a plot similar to **Figure 4** will appear:

.. code-block:: bash

   xmgrace rmsf_prot.dat

.. figure:: /../images/sirah_analysis_4.png
   :align: center
   :width: 100%

   **Figure 4.** *RMSF* plot of the GC bead of all the residues from the SIRAH CG simulation using the Xmgrace program. 

.. tip::
 
   The file ``rmsf_prot.dat`` can be plotted in any external tools for visualization plots.