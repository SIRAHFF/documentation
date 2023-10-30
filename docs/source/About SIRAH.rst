About SIRAH
===========

`SIRAH <http://www.sirahff.com/>`_ is a coarse-grained force field for molecular dynamics simulations that offers plug-and-play topologies and parameters for aqueous solvent, phospholipids, DNA, metal ions, and proteins, including post-translational modifications and glycosylation. It was designed, implemented, and developed by the `Biomolecular Simulations Group <https://pasteur.uy/laboratorios/simulaciones-biomoleculares/)>`_ at the `Institut Pasteur de Montevideo <https://pasteur.uy/>`_, Uruguay, which is also responsible for its regular maintenance.  

The current version, **SIRAH 2.2**, includes:   

* A description of different protonation states and most common post-translational modifications for protein residues;       

* Compatibility for mapping various force field atom types and experimental structures;      

* The ability to leverage GPU acceleration in Amber and GROMACS codes;     

* A collection of scripts, known as **SIRAH Tools**, that facilitate the process of mapping all-atom files to CG representations, backmapping, visualizing, and analyzing SIRAH trajectories directly on the popular VMD program for molecular visualization.     

Advantages
--------------

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

Limitations
--------------

However, it is important to note that the SIRAH force field also has limitations compared to atomistic models:

.. image:: ../images/cons.png
  :width: 30
  :align: left

**Atomic-Level Precision**: CG models, including SIRAH, sacrifice atomic-level precision in order to achieve computational efficiency. This means that certain details, such as the explicit representation of individual water molecules or specific side chain conformations, may not be accurately captured.

.. image:: ../images/cons.png
  :width: 30
  :align: left

**Limited Molecular Diversity**: While SIRAH incorporates parameters for the most prevalent biological molecules, the molecular diversity within biological systems is extensive, making it challenging to encompass all relevant biomolecules. Creating arbitrary molecular topologies may require additional data from experimental databases and physicochemical intuition.

.. image:: ../images/cons.png
  :width: 30
  :align: left

**Protein Folding Simulations**: Although SIRAH has been successful in reproducing spontaneous aggregation, folding of small peptides, and specific recognition of proteins, the unbiased formation of large helical segments is still challenging for the force field.

In summary, the SIRAH force field offers speed, efficiency, and multiscale capabilities for simulating biomolecular systems. However, it sacrifices atomic-level precision and may have limitations in handling molecular diversity and protein folding simulations compared to atomistic models.

License 
-------
The SIRAH force field is in the public domain and is therefore accessible to the entire scientific community without a license. However, if SIRAH was beneficial to your work, please remember to :doc:`cite us <Citation>`.

SIRAH Tools is free software. It can be redistributed and/or modified under the terms of the `GNU General Public License <https://www.gnu.org/licenses/gpl-3.0.html>`_ as published by the `Free Software Foundation <https://www.fsf.org/>`_.


Acknowledgments  
---------------

The SIRAH force field is developed by the `Biomolecular Simulations Group <https://pasteur.uy/laboratorios/simulaciones-biomoleculares/)>`_ and supported partially by the Institut Pasteur de Montevideo and FOCEM (MERCOSUR Structural Convergence Fund - COF 03/11).   


Follow us  
---------------

Follow us in our social media profiles: |google-sirah| |youtube-sirah| |twitter-sirah| |github-sirah| 
 

.. |google-sirah| raw:: html

   <a href="https://scholar.google.com.uy/citations?hl=es&user=Rqr6Jw4AAAAJ" class="author-social" target="_blank"><i class="fa fa-google" style="font-size:30px;"></i></a>
   
.. |youtube-sirah| raw:: html

   <a href="https://www.youtube.com/channel/UCLK25yNsJPlej_fEKdLc1qw" class="author-social" target="_blank"><i class="fa fa-youtube-play" style="font-size:30px;"></i></a>
   
.. |twitter-sirah| raw:: html

   <a href="https://twitter.com/SIRAHForcefield" class="author-social" target="_blank"><i class="fa fa-twitter" style="font-size:30px;"></i></a>

.. |github-sirah| raw:: html

   <a href="https://github.com/SIRAHFF" class="author-social" target="_blank"><i class="fa fa-github" style="font-size:30px;"></i></a>  
