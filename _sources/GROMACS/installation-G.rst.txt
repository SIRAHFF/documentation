Required Software
__________________

**GROMACS**

- GROMACS 4.5.5 or later version properly installed in your computer. 

.. tip::

	`Download <https://www.gromacs.org/Downloads>`_ lastest version of GROMACS and check `installation guide <https://manual.gromacs.org/documentation/current/install-guide/index.html>`_. 


**VMD**

The molecular visualization program `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_, version 1.9.3 or later is (`freely available download <https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD>`_). Check the `installation guide <https://www.ks.uiuc.edu/Research/vmd/current/ig/node6.html>`_ for more instructions.


Prior knowledge
_______________

How to perform a standard atomistic molecular dynamics simulations with GROMACS and basic usage of
VMD. 


.. tip::

    **GROMACS** 

    If you are not familiar with atomistic molecular dynamics simulations we strongly recommend you to perform the usage tutorials of GROMACS: `Introduction to Molecular Dynamics <https://tutorials.gromacs.org/md-intro-tutorial.html>`_, `Introduction to Membrane-Protein Simulation <https://tutorials.gromacs.org/membrane-protein.html>`_ and `others <http://www.mdtutorials.com/gmx/>`_. For more details and documentation, check `GROMACS Manual <https://manual.gromacs.org/>`_. 


    **VMD**

    If you are not familiar with VMD we strongly recommend you to perform the basic usage tutorial of VMD (`using VMD <https://www.ks.uiuc.edu/Training/Tutorials/vmd/tutorial-html/index.html>`_). For more details and documentation, check `VMD User Guide <https://www.ks.uiuc.edu/Research/vmd/current/ug/ug.html>`_ and the `official VMD documentation <https://www.ks.uiuc.edu/Research/vmd/current/docs.html#tutorials>`_.


Download and setting SIRAH
___________________________

Download `SIRAH for GROMACS <https://github.com/SIRAHFF/documentation/releases/download/GROMACS/sirah_x2.2_20-07.tgz>`_ and uncompress it into your working directory. 

.. code-block:: bash

    tar -xzvf sirah_x2.2_20-07.tgz

In case of using another version of sirah:

.. code-block:: bash

	tar -xzvf ``sirah_[version].tgz``
	

You will get a folder ``sirah_[version].ff/`` containing the force field definition, the SIRAH Tools in
``sirah_[version].ff/tools/``, molecular structures to build up systems in ``sirah_[version].ff/PDB``, the required material to perform the tutorial in ``sirah_[version].ff/tutorial/X.X/``.

Make a new folder for this tutorial:

.. caution::

  Remember to change **X.X** to the number that corresponds to the tutorial you are going to do.

.. code-block:: bash

	mkdir tutorialX.X; cd tutorialX.X

Create the following symbolic links in the folder tutorialX.X:

.. code-block:: bash

	ln -s ../sirah_[version].ff sirah.ff

.. code-block:: bash

	ln -s sirah.ff/residuetypes.dat

.. code-block:: bash
	
	ln -s sirah.ff/specbond.dat

.. warning::

	Files residuetypes.dat and specbond.dat are essential for the correct definition of molecular
	groups and auto-detection of disulfide bonds and cyclic DNA polymers.