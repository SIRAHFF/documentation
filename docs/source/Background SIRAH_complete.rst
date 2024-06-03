Background
==================

In recent decades, the field of molecular dynamics (MD) simulations has undergone significant progress and development, allowing for the investigation of biological systems at nanoseconds to microsecond scales. However, complex systems involving millions of atoms or large-size protein assemblies are still too computationally expensive to be studied in atomistic detail without the aid of specialized supercomputers. Thus, the computational cost associated with MD simulations has driven the development of cost-effective approximations, particularly Coarse-Grained (CG) models, which attempt to increase system complexity and spatiotemporal sampling. The resolution of CG models has become essential for investigating molecular processes at the nano and meso dimensions and addressing mechanisms that have been previously inaccessible to traditional modeling methods. CG models consist of effective particles, referred to as beads, that represent groups of corresponding atoms. CG models can vary in resolution, ranging from supra-CG levels where a single bead represents an entire protein, to near-atomistic models that preserve most of the chemical characteristics [:ref:`1 <ref1>`]. 

Currently, one of the most prominent CG models is the `SIRAH <http://www.sirahff.com/>`_ (Southamerican Initiative for a Rapid and Accurate Hamiltonian) force field, which was developed by the `Biomolecular Simulations Group <https://pasteur.uy/laboratorios/simulaciones-biomoleculares/)>`_ `Institut Pasteur de Montevideo <https://pasteur.uy/>`_.  It covers parameters and topologies for aqueous solvent, phospholipids, DNA, metal ions, and proteins. A recent update introduced modifications to bonded and non-bonded parameters, protonation states, post-translational modifications, and compatibility improvements for different force fields [:ref:`2 <ref2>`]. 


The SIRAH force field for CG simulations
-----------------------------------------

The distinguishing features of the SIRAH force field for CG simulations lie in its approach to mapping from atomistic to CG representations and its strategic selection of interaction potentials. The mapping procedure entails strategically placing effective interactive beads at pivotal atoms involved in the structure or at atoms that form crucial intermolecular interactions (**Figure 1**). The distribution of these beads corresponds to intended interactions of functional groups based on size and charge, resulting in a heterogeneous distribution with higher bead density in regions that establish more diverse interactions.

.. figure:: ../images/mousepad-old.png
   :align: center
   :width: 100%
   
   **Figure 1.** SIRAH force field CG representation.   
           
SIRAH employs a classical two-body Hamiltonian, facilitating its use in various MD engines without the need for extensive learning or format changes. This choice enables anyone familiar with standard all-atoms MD simulations in engines like AMBER or GROMACS to run CG simulations using SIRAH seamlessly. The classical Hamiltonian requires the determination of numerous parameters, but SIRAH's mapping strategy significantly reduces this burden. Equilibrium distances are derived directly from statistical data, quantum-level calculations, or canonical conformations, minimizing the number of parameters to be determined.

The initial CG model for DNA served as the foundation for SIRAH, with force constants, partial charges, and Lennard-Jones (LJ) parameters derived through trial and error simulations on DNA segments. The approach of transferring and adapting parameters based on similar functional groups ensures analogous interaction parameters for diverse molecular moieties. SIRAH's versatility is exemplified through its CG models for different biomolecular families. In the following sections, we provide a synopsis of the CG models developed by SIRAH; however, for a more comprehensive material, please refer to [:ref:`3 <ref3>`].


The CG DNA model
-----------------

The SIRAH's DNA model involves six effective beads representing each of the four CG nucleotides (**Figure 2**) [:ref:`4 <ref4>`]. The mapping strategy considers the 5' - 3' prime polarity and electrostatic complementarity between A-T and G-C base pairs. The backbone is represented by two beads at the phosphate and C5' Carbon positions, while three beads on the Watson-Crick edge ensure base pair recognition. The five-membered sugar ring is depicted by a single bead situated at the C1' position, linking the backbone to the Watson-Crick edge.

.. figure:: ../images/mousepad-old-dna.png
   :align: center
   :width: 80%
   
   **Figure 2.** SIRAH force field DNA CG representation.   

By selecting this mapping option, the specific base-pair recognition of the B-form DNA is maintained, and the distortion effects of mismatches are captured precisely. However, it has limitations, since it excludes less frequent inter-nucleotide interactions, such as sugar edge or Hoogsteen base pairs.

In SIRAH's DNA model, bead sizes determined by LJ parameters are heterogeneous, maintaining correct stacking distances in a double-stranded configuration. The effective beads representing bases adopt LJ sizes from the Barcelona force field, while those representing the backbone have larger sizes. The partial charges, assigned to ensure electrostatic recognition, are determined to reproduce the electrostatic potential of the force field in the grooves of a double-stranded structure. The initial mass distribution allows MD simulations with a timestep of 5 fs. 

This CG DNA model reproduces the structure and dynamics of double-stranded DNA comparable to atomistic force fields and demonstrates spontaneous formation of large "bubbles" within DNA, fraying, rehybridization, and matches experimentally determined persistence lengths of single-stranded filaments [:ref:`4 <ref4>`].


The WatFour model for CG explicit solvent
------------------------------------------------

In tandem with the DNA model development, a CG aqueous solvent was created, featuring CG water and monovalent electrolytic ions (sodium, potassium, chloride) (**Figure 3**) [:ref:`5 <ref5>`]. 

.. figure:: ../images/mousepad-old-solvent.png
   :align: center
   :width: 60%
   
   **Figure 3.** SIRAH force field Solvent CG representation.   

Unlike typical CG water, SIRAH's WatFour (WT4) model aimed to replicate the structure of an elementary water cluster, including a central water surrounded by four identical molecules forming a tetrahedron. Hydrogen atoms were removed, and only the oxygen atoms at the tetrahedron's vertices were retained, connected by harmonic bonds. This flexible tetrahedral structure generated its own dielectric permittivity and electrostatic screening by adding partial charges to the four beads, creating a quadrupole with two positively and two negatively charged beads. The partial charges were adopted from the SPC water model to ensure compatibility with fully atomistic water models for multiscale simulations [:ref:`5 <ref5>`].

Iterative fitting was performed on the LJ energy well depth, which corresponded to the experimental diffusion coefficient of pure water at 300 K, and the bead size, which mirrored the second solvation peak of water. The mass of the beads was set to achieve a density of 1 kg/dl. The WT4 model, which resembled a bulkier "water molecule," corresponded to the second apex of the radial distribution function for atomistic water.

Monovalent ions in SIRAH, represented by single beads with a net charge of +/- 1e, were developed based on neutron scattering data, reflecting the chemical identity of sodium, potassium, and chloride ions [:ref:`5 <ref5>`]. The ions' depth of the LJ well matched that of the WT4 beads, offering the flexibility to adjust ionic strength by modifying added salt in the simulation box. The incorporation of electrolytic ions and accurate electrostatic description using the Particle Mesh Ewald summation methods contribute to the relevant features of SIRAH.


The CG protein model
---------------------

The CG protein model in SIRAH employs varying bead sizes to reflect different amino acid interactions. The latest version [:ref:`2 <ref2>`], refined in 2019, has significantly improved the ability to reproduce protein structures. The atomistic to CG mapping of protein side chains follows the DNA model philosophy, with effective beads placed at selected atoms along side chains, representing hydrophobic, aromatic, and polar interactions (**Figure 4**).

.. figure:: ../images/mousepad-old-amino.png
   :align: center
   :width: 90%
   
   **Figure 4.** SIRAH force field amino acids CG representation.  

Hydrophobic amino acids are neutral beads at specific positions with an LJ diameter of 0.42 nm. Aromatic amino acids use smaller beads, 0.35 nm, for stacking-like interactions, with partial charges on certain residues to preserve Hydrogen bond possibilities. Polar amino acids retain beads in functional groups, while acidic and basic amino acids have partial charges which add up to a net charge of +/- 1e.

The aminoacidic backbone is represented with three beads for Nitrogen, Cα Carbon, and carboxylic Oxygen positions, facilitating easy transformation between all-atoms and CG. Bonded parameters for amino acids follow the rules outlined for DNA, with force constants for bond and angular stretching adapted from the same set of parameters. This approach has been found to be effective and time-efficient. 

In version 2.2 [:ref:`2 <ref2>`], all bead masses are set to 50 a.u., and common post-translational modifications, including phosphorylation and acetylation, and different protonation states are available. In addition, divalent ion parameters for Zinc, Magnesium, and Calcium, derived from statistical analyses and validated through multiple CG simulations, enable SIRAH simulations of a wide range of metal-bound macromolecules.


CG models for phospholipids
---------------------------------------

Following the completion of DNA, aqueous solvent, and protein models, the SIRAH force field aimed to incorporate a suitable CG lipid representation for simulating membrane proteins. Focusing on prototypical phospholipids, including phosphatidyl-choline (PC), -ethanolamine (PE), and –serine (PS) heads, along with myristoyl (M), palmitic (P), and oleic (O) acyl chains, SIRAH enabled simulations of diverse eukaryotic membrane components [:ref:`6 <ref6>`]. Utilizing the existing functional groups in the force field, parameterization of these lipids required minimal modifications, ensuring compatibility and accurate replication of lipid bilayer mechanical properties such as thickness, areas per lipid, order parameter, etc., and their dependence with the temperature.

During protein simulations embedded in lipid bilayers, spurious insertions of acyclic tails into the protein core were observed. To address this, specific interactions between hydrophobic protein side chains and acyl chains were set outside Lorentz-Berthelot combination rules, yielding accurate representations of the SarcoEndoplasmic Reticulum Calcium (SERCA) pump's tilted orientation in a DMPC bilayer [:ref:`6 <ref6>`]. This modification facilitated simulations of electrostatics-driven opening of Connexin 26 channels, demonstrating predictive power in identifying mutations inhibiting channel opening [:ref:`7 <ref7>`]. The approach was also employed for cost-effective simulations of entire viral capsids and envelopes, allowing construction and simulation of a Zika Virus-Like Particles on a multi-microsecond time scale [:ref:`8 <ref8>`].


Multiscale simulations
-----------------------

The development of the SIRAH force field in a classical two-body Hamiltonian framework has facilitated multiscale simulations, eliminating the need for non-Hamiltonian interaction terms and ensuring efficiency without communication delays between software modules. 

Two multiscale implementations in SIRAH are emphasized: first, an all-atoms/CG model covalently linking both resolutions within a nucleic acid chain [:ref:`9 <ref9>`]; second, a multiresolution solvent model allowing the mixture of fully atomistic solutes with a shell of atomistic solvent surrounded by CG water, applicable to highly solvated systems like viral capsids [:ref:`10 <ref10>`]. 

A triple solvation scheme, treating water at all-atoms, CG, and supraCG levels, is also available. This is particularly useful for complex cellular systems and has been applied to assemble and simulate VLPs systems in an onion-shaped configuration using CG water (WT4) and supra-CG solvent (WLS) [:ref:`10 <ref10>`]. MD simulations of entire VLPs, such as those studying Flaviviruses with membranes and proteinaceous envelopes, offer crucial insights into their dynamics and are vital for understanding biological systems at a level accessible only through computer simulations [:ref:`8 <ref8>`].


Overwriting combination rules
--------------------------------

The SIRAH force field introduces a modification in the calculation of LJ interactions to address issues with electrolytic ions in proteins and DNA. Unlike traditional MD packages using Lorentz-Berthelot (LB) combination rules, SIRAH employs an "outside-of-LB trick" that allows specific LJ parameters for certain bead pairs, enabling the fine-tuning of interactions. This approach provides adaptability to regulate interactions applying only to specific bead pairs, in accordance with various physicochemical settings [:ref:`3 <ref3>`].

SIRAH comprises 56 different bead types, with 197 interactions defined outside LB combination rules among 1540 possible pair combinations [:ref:`3 <ref3>`]. The modifications include cation-π interactions between aromatic residues and Lysine, methylated Lysine, and zwitterionic N-terminal beads. The force field corrects the size of backbone beads, crucial for forming α helices and Hydrogen bonds, ensuring compatibility with compact structures. It facilitates the formation of secondary structure elements and enhances interactions with other force field components.


Performance
------------

The latest version of the SIRAH force field leverages GPU implementations in GROMACS and AMBER, enabling CG simulations on desktop computers at a rate of a few microseconds per day for medium-sized systems. Larger systems of around a million particles can achieve speeds of hundreds of nanoseconds per day [:ref:`2 <ref2>`]. 

Recently, to illustrate SIRAH's performance, a comparison was made between a SIRAH CG simulation and an atomistic simulation (Amber's FF14SB) of the SARS-CoV-2 Spike protein's receptor binding domain (RBD) with human ACE2 and the amino acid transporter B0AT1 (see [:ref:`3 <ref3>`]). The CG model exhibited a 60-fold speedup, simulating approximately 660 ns per day with a 20 fs time step, compared to the atomistic model's 11 ns per day with a 2 fs time step, using the same system. 

Nevertheless, it is essential to take into account the constraints of the force field beyond its speed implications, as various force fields may possess distinct capabilities. Thus, exercise caution when making direct comparisons between CG force fields, considering their distinct strengths and drawbacks.


Limitations
------------

Although the SIRAH force field offers speed, efficiency, and multiscale capabilities for simulating biomolecular systems, it has some limitations such as: 

* It potentially compromises precision in both structural and energetic aspects. SIRAH, similar to other CG force fields, faces limitations in scenarios that demand atomic-level precision, such as interactions mediated by single water molecules or ligands with specific binding sites. Examples like potassium channels or aquaporins, where individual water molecules play a crucial role, may be challenging for CG models that combine multiple water molecules into a single effective bead.

* Protein folding simulations are not extensively explored. Although SIRAH is successful in reproducing spontaneous aggregation and small peptide folding, the unbiased formation of large helical segments remains challenging. 

* The molecular diversity in biological systems is vast, making it nearly impossible to encompass all relevant biomolecules. Establishing a generally valid methodology for creating arbitrary molecular topologies involves converting new topologies from all-atom to CG, relying on experimental data, organic chemistry knowledge, and physicochemical intuition.


Perspectives
-------------

The rapid advancement of computer power has established MD as a valuable tool in biomedical sciences for understanding intricate processes and vast biological systems. Developing force fields that are universally applicable to all biological molecular families and enable communication at different levels of molecular resolution is a crucial and complex task. Recently, the SIRAH force field expanded its scope by incorporating glycans to simulate polysaccharide chains and protein glycosylation (see [:ref:`11 <ref11>`]). In addition, the lipid diversity will be enhanced by including sphingomyelins, ceramides, and cholesterol, crucial components of endoplasmic reticulum membranes and flaviviral envelopes. Additionally, testing parameters for POPG, a lipid found in bacterial membranes, is underway to improve the realism of antibiotic peptide mode of action descriptions. In the medium term, there are plans to introduce a coarse-grained model for RNA, which is crucial for the description of viral particles and a major area of focus for the group's ongoing research.


References
-------------

.. _ref1:

[1] Borges-Araújo, L.; Patmanidis, I.; Singh, A. P.; Santos, L. H. S.; Sieradzan, A. K.; Vanni, S.; Czaplewski, C.; Pantano, S.; Wataru Shinoda, W.; Monticelli, L.; Liwo, A.; Marrink, S. J.; Souza, P. C. T. Pragmatic Coarse-Graining of Proteins: Models and Applications. Journal of Chemical Theory and Computation. 2023. |Review-2| |Review2-cit| 

.. |Review-2| image:: https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.3c00733-blue?color=blue
   :alt: Access the paper
   :target: https://doi.org/10.1021/acs.jctc.3c00733
   
.. |Review2-cit| image:: https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1021%2Facs.jctc.3c00733
   :alt: Citation
   :target: https://scholar.google.com.uy/scholar?cites=14982031192725054357

.. _ref2:

[2] Machado, M. R.; Barrera, E. E.; Klein, F.; Soñora, M.; Silva, S.; Pantano, S. The SIRAH 2.0 Force Field: Altius, Fortius, Citius. Journal of Chemical Theory and Computation 2019, 15, 2719–2733. |SIRAH2.0|  |SIRAH2.0-cit|

.. |SIRAH2.0| image:: https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.9b00006-blue
   :alt: Access the paper
   :target: https://doi.org/10.1021/acs.jctc.9b00006

.. |SIRAH2.0-cit| image:: https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1021%2Facs.jctc.9b00006
   :alt: Citation
   :target: https://scholar.google.com/scholar?oi=bibs&hl=es&cites=5136612330374064800

.. _ref3:

[3] Klein, F.; Soñora, M.; Santos, L. H.; Frigini, E. N.; Ballesteros-Casallas, A.; Machado, M. R.; Pantano, S. The SIRAH force field: a suite for simulations of complex biological systems at the coarse-grained and multiscale levels. Journal of Structural Biology 2023, 107985. |Review| |Review-cit|

.. |Review| image:: https://img.shields.io/badge/DOI-10.1016%2Fj.jsb.2023.107985-blue
   :alt: Access the paper
   :target: https://doi.org/10.1016/j.jsb.2023.107985
   
.. |Review-cit| image:: https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1016%2Fj.jsb.2023.107985
   :alt: Citation
   :target: https://scholar.google.com/scholar?cites=11014340861876399425

.. _ref4:

[4] Dans, P. D.; Zeida, A.; Machado, M. R.; Pantano, S. A Coarse Grained Model for Atomic-Detailed DNA Simulations with Explicit Electrostatics. Journal of Chemical Theory and Computation 2010, 6, 1711–1725. |DNA| |DNA-cit|

.. |DNA| image:: https://img.shields.io/badge/DOI-10.1021%2Fct900653p-blue
   :alt: Access the paper
   :target: https://doi.org/10.1021/ct900653p

.. |DNA-cit| image:: https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1021%2Fct900653p
   :alt: Citation
   :target: https://scholar.google.com/scholar?oi=bibs&hl=es&cites=12499613729973955498

.. _ref5:

[5] Darré, L.; Machado, M. R.; Dans, P. D.; Herrera, F. E.; Pantano, S. Another Coarse Grain Model for Aqueous Solvation: WAT FOUR? Journal of Chemical Theory and Computation 2010, 6, 3793–3807. |Solvent| |Solvent-cit|

.. |Solvent| image:: https://img.shields.io/badge/DOI-10.1021%2Fct100379f-blue
   :alt: Access the paper
   :target: https://doi.org/10.1021/ct100379f

.. |Solvent-cit| image:: https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1021%2Fct100379f
   :alt: Citation
   :target: https://scholar.google.com/scholar?oi=bibs&hl=es&cites=11533073503238221292

.. _ref6:

[6] Barrera, E. E.; Machado, M. R.; Pantano, S. Fat SIRAH: Coarse-Grained Phospholipids To Explore Membrane–Protein Dynamics. Journal of Chemical Theory and Computation 2019, 15, 5674–5688. |FatSirah| |FatSirah-cit|

.. |FatSirah| image:: https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.9b00435-blue
   :alt: Access the paper
   :target: https://doi.org/10.1021/acs.jctc.9b00435
   
.. |FatSirah-cit| image:: https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1021%2Facs.jctc.9b00435
   :alt: Citation
   :target: https://scholar.google.com/scholar?oi=bibs&hl=es&cites=13191972720970339574

.. _ref7:

[7] Zonta, F.; Buratto, D.; Crispino, G.; Carrer, A.; Bruno, F.; Yang, G.; Mammano, F.; Pantano, S. Cues to Opening Mechanisms From in Silico Electric Field Excitation of Cx26 Hemichannel and in Vitro Mutagenesis Studies in HeLa Transfectans. Frontiers in Molecular Neuroscience 2018, 11, 170. |MemProt-1| |MemProt-cit|     
   
.. |MemProt-1| image:: https://img.shields.io/badge/DOI-10.3389%2Ffnmol.2018.00170-blue
   :alt: Access the paper
   :target: https://doi.org/10.3389/fnmol.2018.00170
   
.. |MemProt-cit| image:: https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.3389%2Ffnmol.2018.00170
   :alt: Citation
   :target: https://scholar.google.com/scholar?cites=7027056542531206464&as_sdt=2005&sciodt=0,5&hl


.. _ref8:

[8] Soñora, M.; Martínez, L.; Pantano, S.; Machado, M. R. Wrapping Up Viruses at Multiscale Resolution: Optimizing PACKMOL and SIRAH Execution for Simulating the Zika Virus. Journal of Chemical Information and Modeling 2021, 61, 408–422. |VLP2| |VLP2-cit|     

.. |VLP2| image:: https://img.shields.io/badge/DOI-10.1021%2Facs.jcim.0c01205-blue
   :alt: Access the paper
   :target: https://doi.org/10.1021/acs.jcim.0c01205
      
.. |VLP2-cit| image:: https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1021%2Facs.jcim.0c01205
   :alt: Access the paper
   :target: https://scholar.google.com/scholar?cites=8645160591236740149

.. _ref9:

[9] Machado, M. R.; Zeida, A.; Darré, L.; Pantano, S. From quantum to subcellular scales: multi-scale simulation approaches and the SIRAH force field. Interface Focus 2019, 9, 20180085. |MC2|  |MC2-cit| 

.. |MC2| image:: https://img.shields.io/badge/DOI-10.1098%2Frsfs.2018.0085-blue?label=DOI
   :alt: Access the paper
   :target: https://doi.org/10.1098/rsfs.2018.0085

.. |MC2-cit| image:: https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1098%2Frsfs.2018.0085
   :alt: Citation
   :target: https://scholar.google.com/scholar?cites=5473055142318037579

.. _ref10:

[10] Machado, M. R.; González, H. C.; Pantano, S. MD Simulations of Virus like Particles with Supra CG Solvation Affordable to Desktop Computers. Journal of Chemical Theory and Computation 2017, 13, 5106–5116. |MC1| |MC1-cit|  

.. |MC1| image:: https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.7b00659-blue
   :alt: Access the paper
   :target: https://doi.org/10.1021/acs.jctc.7b00659

.. |MC1-cit| image:: https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1021%2Facs.jctc.7b00659
   :alt: Citation
   :target: https://scholar.google.com/scholar?cites=16637391138490147245

.. _ref11:

[11] Garay, P. G.; Machado, M. R.; Verli, H.; Pantano, S. SIRAH late harvest: coarse-grained models for protein glycosylation. Journal of Chemical Theory and Computation 2024. |GLY| |GLY-cit|

.. |GLY| image:: https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.3c00783-blue
   :alt: Access the paper
   :target: https://pubs.acs.org/doi/10.1021/acs.jctc.3c00783

.. |GLY-cit| image:: https://img.shields.io/endpoint?url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1021%2Facs.jctc.3c00783
   :alt: Citation
