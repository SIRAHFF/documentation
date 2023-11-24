.. note:: 

    SIRAH Force Field is distributed with `AMBER <https://ambermd.org/index.php>`_ and `AmberTools <https://ambermd.org/AmberTools.php>`_ offical releases.

Required Software
__________________

**AMBER**

Amber and AmberTools (version 16 or later) versions properly installed and running in your computer. 

.. tip::

    `Download <https://ambermd.org/GetAmber.php#ambertools>`_ the lastest version of Amber and AmberTools and refer to the *How to obtain Amber23* section for more details on getting Amber23. In addition, go to the `Install Amber <https://ambermd.org/Installation.php>`_ page for specific instructions based on your operating system. If you have Amber already installed, it is unnecessary to install AmberTools separately.


**VMD**

The molecular visualization program `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_, version 1.9.3 or later (`freely available download <https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD>`_). Check VMD's `Installation guide <https://www.ks.uiuc.edu/Research/vmd/current/ig/node6.html>`_ for more instructions.


Prior knowledge
_______________

How to perform a standard atomistic molecular dynamics simulations with AMBER and basic usage of
VMD. 


.. tip::

    **AMBER** 

    If you are not familiar with atomistic molecular dynamics simulations, we strongly recommend you to perform the basic usage tutorials of Amber and AmberTools: `1. Building Systems <https://ambermd.org/tutorials/BuildingSystems.php>`_, `3. Creating Stable Systems and Running MD <https://ambermd.org/tutorials/Relaxation.php>`_ and `5. Case Studies <https://ambermd.org/tutorials/Introductory.php>`_. For more details and documentation, check `Amber23 Manual <https://ambermd.org/doc12/Amber23.pdf>`_. 


    **VMD**

    If you are not familiar with VMD, we strongly recommend you to perform the basic usage tutorial of VMD (`Using VMD <https://www.ks.uiuc.edu/Training/Tutorials/vmd/tutorial-html/index.html>`_). For more details and documentation, check `VMD User Guide <https://www.ks.uiuc.edu/Research/vmd/current/ug/ug.html>`_ and the official `VMD documentation <https://www.ks.uiuc.edu/Research/vmd/current/docs.html#tutorials>`_.

Download and Setting up SIRAH
_______________________________

.. note::

   Please report bugs, errors or enhancement requests through `Issue Tracker <https://github.com/SIRAHFF/documentation/issues>`_ or if you have a question about SIRAH open a `New Discussion <https://github.com/SIRAHFF/documentation/discussions>`_.

Download `SIRAH for Amber <https://github.com/SIRAHFF/documentation/releases/download/AMBER/sirah_x2.3_23-11.amber.tar.gz>`_ and uncompress it into your working directory.

.. code-block:: bash

    tar -xzvf sirah_x2.3_23-11.amber.tar.gz

If you are using a different version of SIRAH, type:

.. code-block:: bash

    tar -xzvf sirah_[version].amber.tgz


You will get the folder ``sirah_[version].amber/`` containing the force field definition, SIRAH Tools in the 
``sirah_[version].amber/tools/`` folder, molecular structures to build up systems in ``sirah_[version].amber/PDB/`` and the required material to perform the tutorials in ``sirah_[version].amber/tutorial/X.X/``

.. caution::

	Remember to change **X.X** to the number that corresponds to the tutorial you are going to do.

Make a new folder for this tutorial in your working directory:

.. code-block:: bash

    mkdir tutorialX.X; cd tutorialX.X

Create the following symbolic link in the folder tutorialX.X:

.. code-block:: bash

    ln -s ../sirah_x2.3_23-11.amber sirah.amber

.. note:: 

    SIRAH Force Field is also distributed with `AMBER <https://ambermd.org/index.php>`_ and `AmberTools <https://ambermd.org/AmberTools.php>`_ official releases. To use the native AMBER version of SIRAH, create a symbolic link located in $AMBERHOME:

    .. code-block:: bash
	
		ln -s $AMBERHOME/dat/SIRAH sirah.amber

    Check the `AMBER Manual <https://ambermd.org/doc12/Amber23.pdf>`_ section **3.11.2** for more details. 

    However, if you want the lastest parameters and implementations strongly advise you to use the developers version of SIRAH from GitHub. 