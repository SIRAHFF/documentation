.. SIRAH documentation master file, created by
   sphinx-quickstart on Wed Aug  9 14:07:26 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SIRAH FF documentation
==================================

.. caution::

   **This website is currently undergoing construction, which means that some of its components are constantly being updated and expanded upon.**

.. image:: ../images/pagina_principal.gif
  :width: 100%


`SIRAH <http://www.sirahff.com/>`_ is a simple and intuitive coarse-grained force field for molecular dynamics simulations. It offers plug-and-play topologies and parameters for aqueous solvent, phospholipids, DNA, metal ions, and proteins, including post-translational modifications and glycosylation. SIRAH current version has the ability to leverage GPU acceleration in popular MD codes. It provides a collection of scripts, SIRAH tools, to facilitate setting up a SIRAH CG simulation. 


The SIRAH force field offers several advantages compared to atomistic models in simulating biomolecular systems:

.. image:: ../images/pros.png
  :width: 30
  :align: left

**Speed and Efficiency**: Coarse-grained (CG) models, such as SIRAH, provide a significant speedup in simulations compared to atomistic models. This allows for the study of larger and more complex biological systems over longer time scales.

.. image:: ../images/pros.png
  :width: 30
  :align: left

**Multiscale Capabilities**: SIRAH offers a multiscale approach, allowing for the combination of different levels of resolution within the same simulation setup. This includes the ability to covalently link all-atom and CG models in nucleic acid chains and the incorporation of a multi-resolution model for the solvent.

.. image:: ../images/pros.png
  :width: 30
  :align: left

**Self-Consistent Force Field**: SIRAH aims to provide a self-consistent force field that describes the intra and intermolecular interactions among moieties of diverse chemical nature in a holistic view of biological ensembles' structural and dynamic features. This ensures that the force field is capable of accurately representing the interactions within biomolecular systems.

The force field is accessible to the entire scientific community without a license. Check out the :doc:`About SIRAH` section for further information.

.. note::

   Please report bugs, errors or enhancement requests through `Issue Tracker <https://github.com/SIRAHFF/documentation/issues>`_ or if you have a question about SIRAH open a `New Discussion <https://github.com/SIRAHFF/documentation/discussions>`_.


Follow us in our social media profiles: |google-sirah| |youtube-sirah| |twitter-sirah| |github-sirah| 
 

.. |google-sirah| raw:: html

   <a href="https://scholar.google.com.uy/citations?hl=es&user=Rqr6Jw4AAAAJ" class="author-social" target="_blank"><i class="fa fa-google" style="font-size:30px;"></i></a>
   
.. |youtube-sirah| raw:: html

   <a href="https://www.youtube.com/channel/UCLK25yNsJPlej_fEKdLc1qw" class="author-social" target="_blank"><i class="fa fa-youtube-play" style="font-size:30px;"></i></a>
   
.. |twitter-sirah| raw:: html

   <a href="https://twitter.com/SIRAHForcefield" class="author-social" target="_blank"><i class="fa fa-twitter" style="font-size:30px;"></i></a>

.. |github-sirah| raw:: html

   <a href="https://github.com/SIRAHFF" class="author-social" target="_blank"><i class="fa fa-github" style="font-size:30px;"></i></a>  

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Manual

   About SIRAH
   Citation
   Further reading
   FAQ
   Developers

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Tutorials

   Tutorials amber
   Tutorials gromacs
   Tutorials namd
   Tutorials sirahtools
   Tutorials analysis
