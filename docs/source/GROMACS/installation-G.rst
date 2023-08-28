Required Software
__________________


- GROMACS 4.5.5 or later version properly installed in your computer. 
- The molecular visualization program VMD (freely available at www.ks.uiuc.edu/Research/vmd).

Prior knowledge
________________

How to perform a standard atomistic molecular dynamic simulation with GROMACS.

Set Up SIRAH
_____________

Download the file ``sirah_[version].gmx.tgz`` from <www.sirahff.com> and uncompress it into your
working directory. 

.. note:: 
	
	[version] should be replaced with the actual package version e.g.: x2_18-09
	tar -xzvf ``sirah_[version].gmx.tgz``

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