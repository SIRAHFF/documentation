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

   **Figure 4.** Our protein-DNA complex after using the SIRAH macros ``sirah_backbone``, ``sirah_water``, and ``sirah_ions``. To improve visualization, we used Licorice as **Drawing Method** for ``sirah_backbone`` and CPK for ``sirah_water`` and ``sirah_ions``. In addition, the **Material** of the water beads was changed to Ghost and we hid some elements of the system.
   
.. tip:: 

	As you can see, macros can be very useful and when saving your work in a saved state (**File** > **Save Visualization State**), macros are included in the saved state file to be used later (**File** > **Load Visualization State**).


Structural analysis of CG systems
__________________________________

Besides the features that enhance the visualization of SIRAH CG simulations, two additional SIRAH Tools features can be used to analyze the trajectory: ``sirah_ss`` and ``sirah_backmap``. They will be discussed below.


Secondary structural analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The utility ``sirah_ss`` assigns secondary structures to CG proteins in SIRAH, classifying residues into *α-helix (H)*, *extended β-sheet (E)* or, otherwise, *coil (C)* conformations, based on the instantaneous values of the backbone’s torsional angles and Hydrogen bond-like (HB) interactions (see `Darre et al <https://pubs.acs.org/doi/10.1021/ct5007746>`_ for more information). The ``sirah_ss`` feature produces ASCII files of average and by-frame results, which can be visualized as color plots using the python script ``plot_ss.py``.

With the ``sirah_vmdtk.tcl`` file loaded, we can access the ``sirah_help`` feature by going to **Extensions** > **Tk Console** and entering 

.. code-block:: console

   sirah_help 

to see the all the SIRAH tools available options. Or, we can go straight to the help for the ``sirah_ss`` feature by typing

.. code-block:: console

   sirah_ss help

The output will be the following:

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


All the available options for ``sirah_ss`` are documented in this help. By default, four outputs are created by ``sirah_ss``: ``ss_by_frame.xvg``, ``ss_by_res.xvg``, ``ss_global.xvg`` and ``ss.mtx``. The ``ss_by_frame.xvg`` file gives the secondary structure percentages of the three classifications (H, E, or C) by each frame of the simulations. The ``ss_by_res.xvg`` file gives the secondary structure percentages of H, E, and C by each residue. The ``ss_global.xvg`` file shows the overall percentages and standard deviation of H, E, and C of the entire simulation. And the ``ss.mtx`` file gives a matrix of the secondary structure transitions, between H, E, or C, of the residues versus the simulation time.

We can also choose on which *mol* and on which frames to do the analysis within VMD. For example, if we want to calculate the secondary structure for the whole structure (loaded on *top*) on all frames, we can:

1. Make a VMD selection. By defalut sirah_ss selection is *all*, but if we need to change the selection, for example to choose only the protein, we can type in the **Tk Console**:

	.. code-block:: console

		set protein [atomselect top sirah_protein]

	Here we selected only the protein beads.

2. Calculate the secondary structure information of the selection by typing:

	.. code-block:: console

		sirah_ss

This will, by default, calculate the secondary structure for all frames of the trajectory of our selection and generate the four files.

Furthermore, we have prepared a python script ``plot_ss.py`` that can be used to plot three of the four files (``ss_by_frame.xvg``, ``ss_by_res.xvg`` and ``ss.mtx``). 

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

Once we know we are using the right version of Python, we can use the ``plot_ss.py`` script by typing:

.. code-block:: python

   python plot_ss.py -h

This will show us the options we have within the script, most of the flags are for changing the appearance of the graph (colors, font size, label size, etc). 

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

The basic use of the script is: 

.. code-block:: console

   python plot_ss.py -i filename 

where ``filename`` is one of the files ``ss_by_res.xvg``, ``ss_by_frame.xvg`` or ``ss.mtx``.

In the case of the ``ss.mtx`` matrix it may take some time (no more than a couple of minutes), so be patient. 

.. tip::

   To use the flags that modify the colors, any of the following entries is valid. For example, if we would like to use red as the color for the α-helix, we can type:

   ``-H red``
   ``-H r````
   ``-H "red``
   ``-H "r"``
   ``-H "#FF0000"``

   In the case of HEX code usage, quotation marks are required.

The script generates the plots shown in Figure 5, remember that you must plot one by one, selecting the filename according to your interest.

.. figure:: /../images/sirah_tools_5.png
   :align: center
   :width: 100%

   **Figure 5.** Outcomes of ``plot_ss.py`` where the secondary structural elements α-helix (H), extended β-sheet (E), and coil (C) are colored purple, yellow, and blue, respectively. A) The ``ss.mtx`` matrix plot shows the protein's secondary structure transitions during the 3.0 μs MD simulation. B)  The ``ss_by_frame.xvg`` plot shows the three secondary structure elements (H, E, and C) percentages by each frame of the 3.0 μs MD simulation. C) The ``ss_by_res.xvg`` plot shows the three secondary structure elements (H, E, and C) percentages by each residue during the 3.0 μs MD simulation. Here, we put the three plots together and edited the figure to use only one legend. However, the script makes three separate plots, each with its own legend.

.. note::
	
	The files ``ss_by_res.xvg``, ``ss_by_frame.xvg``, and ``ss.mtx`` can also be plotted using other plotting programs such as `Grace <https://plasma-gate.weizmann.ac.il/Grace/>`_ or `R <https://www.r-project.org/>`_. In addition, the ``plot_ss.py`` script can also be modified or improved according to the necessity of the user.



Backmapping analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
