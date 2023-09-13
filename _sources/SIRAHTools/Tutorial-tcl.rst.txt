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

The utility ``sirah_ss`` assigns secondary structures to CG proteins in SIRAH classifying aminoacids in *a-helix (H)*, *extended b-sheet (E)* or, otherwise, *coil (C)* conformations,based on the instantaneous values of the backbone’s torsional angles and Hydrogen bond-like (HB) interactions (`Darre et al <https://pubs.acs.org/doi/10.1021/ct5007746>`_). Function sirah_ss produces ASCII files of average and by-frame results, which can be visualized as a color matrix using the python script ``plot_ss.py``.

To use the ``sirah_ss`` tool it is necessary to load the ``sirah_vmdtk.tcl`` script. In the vmd Tk console we could see the options that the script gives us: 

.. code-block:: console

   sirah_help 

or instead view the help for the sirah_ss functionality

.. code-block:: console

   sirah_ss help

we will get the following output:

.. code-block:: console

      >>>> sirah_ss <<<<

      version: 1.1 [Dec 2015]

      Description:

        Command to calculate the secondary structure of SIRAH proteins.
        
        If mol is omitted then all options are applied to top molecule.
        By default the secondary structure data is stored in memory
        as the array sscache_data(mol,frame) and the coloring method
        'Secondary Structure' will display correctly. By default the
        script reads the sscache_data if available in memory and does
        not recalculate the conformation, unless flag ramach is set or
        the sscache_data is cleaned. The outname and noprint keywords
        control the generated output files.

        The current command version does not support PBC options.

      Usage: sirah_ss [options]

      Options: These are the optional arguments

        mol       Set the molecule ID, default top
        
        first     First frame to analyze
        
        last      Last frame to analyze
        
        sel       "VMD selection", default "all"
        
        outname   List {} of keywords and names for output files.
                  Keywords and default names:
                  
                    byframe   ss_by_frame.xvg
                    byres     ss_by_res.xvg
                    global    ss_global.xvg
                    mtx       ss.mtx
                    phi       phi.mtx
                    psi       psi.mtx
        
        load      Load secondary structure data from file (e.g ss.mtx) into molecule (mol)
        
        noprint   Flag to avoid printing standard results to out files
        
        ramach    Flag to write phi and psi angles of residue selection (sel) to out files
        
        nocahe    Flag to avoid saving data to sscache_data
        
        now       Flag to calculate the secondary structure at current frame.
                  No sscache_data or output file is saved.
        
        clean     Flag to clean the sscache_data of the selected molecule (mol) and exit
        
        help      Display help and exit


      Examples:

        1.        Load secondary structure data into molecules 2 and extend the calculation
                  sirah_ss mol 2 load ss.mtx
                  sirah_ss mol 2

        2.        Set custom output file names
                  sirah_ss outname {mtx myss.mtx byframe myss_by_frame.xvg}

In the help you can see the available flags for sirah_ss. By default the outputs of sirah_ss are the by frame file ``s_by_frame.xvg``, by residue ``ss_by_res.xvg``, average file ``ss_global.xvg`` and a matrix of secondary structure vs time ``ss.mtx``. 

We can also select on which mol within VMD to perform the analysis and on which frames, e.g., in the standard case if we would like to calculate it for the whole structure (loaded on top) on all frames: 

1. Make a selection, by defalut sirah_ss selection is *all*, but in the case that need change selection, for instance select only the protein in the Tk console we can typing: 

.. code-block:: console

   set protein [atomselect top sirah_protein]

Here we select only the protein beads.

2. Calculate the secondary structure in the selection

.. code-block:: console

   sirah_ss

This will by default calculate the secondary structure for all frames of the trajectory of the top molecule in our selection and by default generate four files (byres, byframe, global and mtx matrix).

We have prepared a python script ``plot_ss.py`` that can be used to plot any of the output files (byres, by frame and mtx matrix). 

.. important:: 

   The ``plot_ss.py`` script works properly with **python 3.9**. You can for example create a conda environment using: 

   .. code-block:: console

      conda create --name plot_ss python=3.9

   and activate it using:

   .. code-block:: console

      conda activate plot_ss

   Finally check de python version:

   .. code-block:: console

      python --version

Once we are sure we are working on a correct python version, we can use the ``plot_ss.py`` script.

.. code-block:: python

   python plot_ss.py -h

This will show us the options we have within the script, most of the flags are to change the appearance of the graph (colors, font size, label size). 

.. code-block:: console

   Create a PNG image from an input file (ss.mtx, ss_by_frame.xvg, or ss_by_res.xvg) generated by SIRAH.

   optional arguments:
     -h, --help            Show this help message and exit
     -i [input]            Input file name (ss.mtx, ss_by_frame.xvg, ss_by_res.xvg)
     -d [dpi]              DPI (dots per inch) for saving the figure (default: 300)
     -tu [tu]              Time unit (us or ns) (default: us)
     -dt [dt]              Time between consecutive frames (default: 1e-04)
     -H [helix-color]      Alpha-helix color (default: darkviolet)
     -E [beta-sheet color] Beta-sheet color (default: yellow)
     -C [coil color]       Coil color (default: aqua)
     -o [out name]         Image output name (default: input name)
     -wt [width]           Image width (inches) (default: 10)
     -ht [height]          Image height (inches) (default: 8)
     -xfs [xlab fontsize]  Fontsize for x-axis labels (default: 14)
     -yfs [ylab fontsize]  Fontsize for y-axis labels (default: 14)
     -yticks [y # ticks]   Number of major tick locators on the y-axis (default: 10)
     -xticks [x # ticks]   Number of major tick locators on the x-axis (default: 10)
     -xtsize [xtsize]      Size of ticks on the x-axis (default: 12)
     -ytsize [ytsize]      Size of ticks on the y-axis (default: 12)
     -title [title]        Title of the plot (default: None)
     -ttsize [title_size]  Size of the plot title (default: 16)
     --version             Print version and exit

the basic use of the script is: 

.. code-block:: console

   python plot_ss.py -i filename 

``filename`` is any of the files ``ss_by_res.xvg``, ``ss_by_frame.xvg`` or ``ss.mtx``.

In the case of the ``ss.mtx`` matrix it may take some time (no more than a couple of minutes), so be patient. 

.. tip::

   To use the flags that modify the colors any of the following entries is valid for example if we would like to use red as the color for the a-helix:

   ``-H red``
   ``-H r````
   ``-H "red``
   ``-H "r"``
   ``-H "#FF0000"``

   In the case of HEX code in mandatory use quotes.

The script generates the graphs shown in figure 5, remember that you must plot one by one, selecting the filename according to your interest.





Backmapping analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
