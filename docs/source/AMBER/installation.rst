Required Software
__________________

**AMBER**

AMBER (version 16) and AMBER Tools (version 16) or later versions properly installed and running in your computer. 

.. tip::

    `Download <https://ambermd.org/GetAmber.php#ambertools>`_ lastest version of amber. See *How to obtain Amber23* section. 
    `Install amber <https://ambermd.org/Installation.php>`_. Select according to the operating system you are using 
    `Install amberTools <https://ambermd.org/GetAmber.php#ambertools>`_.  If you have amber installed, it is not necessary to install ambertools separately.


**VMD**

The molecular visualization program `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_, version 1.9.3 or later is (`freely available download <https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD>`_). Check the `installation guide <https://www.ks.uiuc.edu/Research/vmd/current/ig/node6.html>`_ for more instructions.


Prior knowledge
_______________

How to perform a standard atomistic molecular dynamics simulations with AMBER and basic usage of
VMD. 


.. tip::

    **AMBER** 

    If you are not familiar with atomistic molecular dynamics simulations we strongly recommend you to perform the basic usage tutorials of Amber and AmberTools: `1.Building Systems <https://ambermd.org/tutorials/BuildingSystems.php>`_, `3.Creating Stable Systems and Running MD <https://ambermd.org/tutorials/Relaxation.php>`_ and `5.Case Studies <https://ambermd.org/tutorials/Introductory.php>`_. For more details and documentation, check `Amber23 Manual <https://ambermd.org/doc12/Amber23.pdf>`_. 


    **VMD**

    If you are not familiar with VMD we strongly recommend you to perform the basic usage tutorial of VMD (`using VMD <https://www.ks.uiuc.edu/Training/Tutorials/vmd/tutorial-html/index.html>`_). For more details and documentation, check `VMD User Guide <https://www.ks.uiuc.edu/Research/vmd/current/ug/ug.html>`_ and the `official VMD documentation <https://www.ks.uiuc.edu/Research/vmd/current/docs.html#tutorials>`_.

Download and setting SIRAH
___________________________

Download `SIRAH for amber <https://github.com/SIRAHFF/documentation/releases/download/AMBER/sirah_x2.2_20-08.amber.tgz>`_ and uncompress it into your working directory.

.. code-block:: bash

    tar -xzvf sirah_x2.2_20-08.amber.tgz

In case of using another version of sirah:

.. code-block:: bash

    tar -xzvf sirah_[version].amber.tgz


You will get a folder ``sirah_[version].amber/`` containing the force field definition, the SIRAH Tools in
``sirah_[version].amber/tools/``, molecular structures to build up systems in ``sirah_[version].amber/PDB/`` and the required material to perform tutorials in ``sirah_[version].amber/tutorial/X.X/``

.. caution::

  Remember to change **X.X** to the number that corresponds to the tutorial you are going to do.

Make a new folder for this tutorial in your working directory:

.. code-block:: bash

    mkdir tutorialX.X; cd tutorialX.X

Create the following symbolic link in the folder tutorialX.X:

.. code-block:: bash

    ln -s ../sirah_x2.2_20-08.amber sirah.amber
