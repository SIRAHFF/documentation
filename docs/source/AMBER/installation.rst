Required Software
__________________

.. _vmd: www.ks.uiuc.edu/Research/vmd

AMBER 16 and AMBER Tools 16 or later versions properly installed and running in your computer. The molecular visualization program VMD 1.9.3 or later version (`freely available download <vmd>`_).


Prior knowledge
_______________

How to perform a standard atomistic molecular dynamic simulation with AMBER and basic usage of
VMD. If you are not familiar with DNA stuff we strongly recommend you to first perform the `AMBER
tutorial on DNA <http://ambermd.org/tutorials/basic/tutorial1>`_.

Download and set up SIRAH Force Field
______________________________________

Download the file ``sirah_[version].amber.tgz`` from www.sirahff.com and uncompress it into your
working directory. Notice: ``[version]`` should be replaced with the actual package version e.g.: x2_18-09

.. code-block:: bash

    $ tar -xzvf sirah_[version].amber.tgz

You will get a folder ``sirah_[version].amber/`` containing the force field definition, the SIRAH Tools in
``sirah_[version].amber/tools/``, molecular structures to build up systems in ``sirah_[version].amber/PDB/``,
frequently asked questions in ``sirah_[version].amber/tutorial/SIRAH_FAQs.pdf`` and the required
material to perform the tutorial in ``sirah_[version].amber/tutorial/1/``

Make a new folder for this tutorial in your working directory:

.. code-block:: bash

    $ mkdir tutorial1; cd tutorial1

Create the following symbolic link in the folder tutorial1:

.. code-block:: bash

    ln -s ../sirah_[version].amber sirah.amber
