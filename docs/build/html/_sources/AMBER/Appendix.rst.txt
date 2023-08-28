Calculating ionic concentrations
________________________________

.. math::   

.. role:: raw-latex(raw)
            :format: latex html

.. raw:: html

   <script type="text/javascript" src="http://localhost/mathjax/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>

.. raw:: latex html

       
           \[\rho_{\text{WT4}} = \rho_{\text{H}_2\text{O}} = 1000 \, \text{g/L} \]

            \[MW_{\text{H}_2\text{O}} = 18 \, \text{g/mol} \]

            \[1 \text{WT4} \sim 11 \text{H}_2\text{O} \]

            \[\text{M}=\frac{\text{mol}}{V} \quad ; \quad \text{n}={\text{N}_A} \quad ; \quad \rho = \frac{\text{m}}{V} \quad ; \quad \text{m} = \text{mol} \cdot \text{M} \cdot \text{W}\]

            \[\text{V} = \frac{\text{m}}{\rho} = \frac{\text{mol} \cdot \text{MW}_{\text{H}_2\text{O}}}{\rho} = \frac{\text{n}_{\text{H}_2\text{O}} \cdot \text{MW}_{\text{H}_2\text{O}}}{{\text{N}_A} \cdot \rho}\]


            \[\text{M}=\frac{\text{mol}}{V} = \frac{\text{n}_{\text{ion}}}{\text{N}_{\text{A}} \cdot \text{V}} = \frac{\text{n}_{\text{ion}}}{\text{N}_{\text{A}}} \frac{\text{N}_{\text{A}} \cdot \rho}{\text{n}_{\text{H}_2\text{O}} \cdot \text{MW}_{\text{H}_2\text{O}}} = \frac{\text{n}_{\text{ion}} \cdot 1000}{\text{n}_{\text{WT}_4} \cdot \text{(11)} \cdot \text{(18)}} \sim 5 \frac{\text{n}_{\text{ion}}}{\text{n}_{\text{WT}_4}}  \]


Number of *WT4* molecules per ion at 0.15M: :raw-latex:`\(\text{n}_{WT4} = 5 {{\text{n}_{\text{ion}}} \over \text{M}} = {5(1) \over 0.15} \sim 34 \)`


.. _table: 

Available mapping files
_______________________________


Table 1. Available mapping files (MAPs) at folder ``sirah.amber/tools/CGCONV/maps/`` for converting atomistic lipid structures to SIRAH models. **Important!** MAPs can not inter-convert different name conventions, e.g. amber_lipid.map won’t generate fragment-based residues from residue-based force fields. Due to possible nomenclature conflicts, users are advised to check and modify the MAPs as required.

+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
|      **Map**            | **Type**\* | **Compatibility**                        | **Source**                                                     |
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| amber_lipid.map         |   F        | AMBER Lipid11-17 force fields            | | `AMBER <http://ambermd.org/>`__                              |
|                         |            |                                          | | `HTMD <https://software.acellera.com/htmd/index.html>`__     | 
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| GAFF_lipid.map          |   R        | AMBER GAFF force field                   | `LipidBook <https://lipidbook.org/>`__                         |
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| charmm_lipid.map        |   R        | | CHARMM 27/36 force field, and “CHARMM  | | `CHARMM-GUI <https://charmm-gui.org/>`__                     |
|                         |            | | compatible” GAFF nomenclature          | | `GROMACS <https://www.gromacs.org/>`__                       |  
|                         |            |                                          | | `LipidBook <https://lipidbook.org/>`__                       | 
|                         |            |                                          | | `MemBuilder <http://bioinf.modares.ac.ir/software/mb/>`__    |     
|                         |            |                                          | | `HTMD <https://software.acellera.com/htmd/index.html>`__     |
|                         |            |                                          | | `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_               |   
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| slipids.map             |   R        | Stockholm lipids force field             | | `SLIPIDS <http://www.fos.su.se/~sasha/SLipids/About.html>`__ |
|                         |            |                                          | | `MemBuilder <http://bioinf.modares.ac.ir/software/mb/>`__    |   
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| OPLSA-AA_2014_lipid.map |   R        | All-atoms lipids for OPLS force field    | | `Maciejewski <https://doi.org/10.1021/jp5016627>`__          |
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| OPLSA-UA_lipid.map      |   R        | United-atom lipids for OPLS force field  | `LipidBook <https://lipidbook.org/>`__                         |
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| GROMOS43a1_lipid.map    |   R        | | United-atom lipids for GROMOS 43a1 and | | `LipidBook <https://lipidbook.org/>`__                       |
|                         |            | | CKP force fields                       | | `MemBuilder <http://bioinf.modares.ac.ir/software/mb/>`__    |      
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| GROMOS43a1-s3_lipid.map |   R        | | United-atom lipids for GROMOS 43a1-s3  | | `GROMACS <https://www.gromacs.org/>`__                       |
|                         |            | | force field                            | | `LipidBook <https://lipidbook.org/>`__                       |
|                         |            |                                          | | `MemBuilder <http://bioinf.modares.ac.ir/software/mb/>`__    |
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| GROMOS53a6_lipid.map    |   R        | | United-atom lipids for GROMOS 53a6     | | `GROMACS <https://www.gromacs.org/>`__                       |
|                         |            | | force field                            | | `LipidBook <https://lipidbook.org/>`__                       | 
|                         |            |                                          | | `MemBuilder <http://bioinf.modares.ac.ir/software/mb/>`__    | 
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+
| tieleman_lipid.map      |   R        | | Berger lipids as implemented by        | | `Tieleman <https://doi.org/10.1021/ja0624321>`__             |
|                         |            | | Tieleman et al. for GROMOS             | | `LipidBook <https://lipidbook.org/>`__                       |  
|                         |            | | force fields.                          |                                                                |
+-------------------------+------------+------------------------------------------+----------------------------------------------------------------+

\* Fragment-based (F) or Residue-based (R) topology.



