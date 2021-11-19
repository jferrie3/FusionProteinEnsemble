# FusionProteinEnsemble
 
## Algorithm Descriptions
### FusionProteinModeler
The FusionProteinModeler script takes sequences of fusion protein constituents and linkers as inputs along with the structures of each constituent to generate a structure of the full-length fusion protein along with additional files useful for subsequent sampling.
### FastFloppyLinker
The FastFloppyLinker script take structural outputs from FusionProteinModeler along with information and fragment libraries for the linker segments and information about multimerization and symmetry, when applicable, to sample the structural dynamics of fusion proteins.

## Directory Descriptions
### FusionProteinModeler_Input_Output_Files
Contains the inputs for running FusionProteinModeler for all Halo-SNAPf fusion proteins and the resulting outputs
### FastFloppyLinker_Output_Files
Contains the resultant structural ensembles and inter-fluorophores distances from FastFloppyLinker simulations
### Cluster_Submit_Scripts
Contains the submit scripts used for running FastFloppyLinker on the Savio cluster at UC Berkeley
### Accessory_Scripts
Contains extra helpful scripts related to this project

## Installation Guide
__Operating System:__ Linux (64-bit) or Windows Ubuntu Subsystem

__Programming Language__: Python
This code was specifically written and tested in Python3.7 (python3.7.5)

__Required Python Packages:__
- Pyrosetta
	- This was specifically written and tested with PyRosetta 2020.03+release.f1c7e97=py37_0
- biopython
- pyrosetta_help (https://github.com/matteoferla/pyrosetta_help)

__Anaconda Environment:__
An anaconda environment containing all necessary packages can be found in the anaconda folder. Build time for this Anaconda environment takes on the order of mintues to hours depending on the number of processors used and is largely dependent on the PyRosetta build time. On a normal computer this is expected to take ~30 minutes. With this file you can generate a local version of this environment using the command:

```conda env create -f lion.yml```

__Additional Reccommended Software Packages:___
- Blast (blast-2.2.26) (https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/)

## Running FusionProteinModeler
FusionProteinModeler is run using the command-line, where inputs are specified using an arguement parser. An example run command is as follows:
```
run ../FusionProteinModeler.py -p Halo-SNAPf.pdbs -f Halo-SNAPf.fastas -frags_rapx Halo-SNAPf_Reweighted.rapx.ss2 -frags_ppred Halo-SNAPf_Reweighted.ppred.ss2 -noblast
```
Each of the parser flags are described below:
```
-f  --Sequence_FASTAs  Input the sequences files of each of the input domains/linkers in order form N to C
-p  --Input_PDBs PDB structure each of each non-linker domain
-md  --Multimer_Domain  Enter the PDB file name of the monomeric PDB used in the Input_PDBs arguement
-mp  --Multimer_PDB  PDB structure of the multimerizing domain in the multimeric form for symmetry construction
-floppy  --FastFloppyLinker  Running with FastFloppyLinker will automatically perform the FastFloppyLinker simulation after creating the fusion protein
-noblast  --NoBlast  Turn off the PsiBlast search for fragment selection
-legacyblast  --LegacyBlast  Perform the Blast search using the legacy blastgpg
-nofrags  --NoFrags  Turn off the Fragment Picker and use of Fragments in FastFloppyLinker
-frags_ppred  --PPred_Fragment_File  PsiPred file Required for the Fragment Picker
-frags_rapx  --RapX_Fragment_File  RaptorX file Required for the Fragment Picker
-o  --Output_Name  Name of Output PDB structure NOTE: Do not add .pdb
```
### Acquiring the Necessary Inputs
#### Generating FASTA files
FASTA files should directly correpond to cloned proteins. When determining linker sequences, in addition to including the cloned linker sequences in each linker FASTA, N/C-terminal tails that are not included within the structured regions of supplied PDBs should be included as linkers.
#### Fragment Library Construction
Construction of fragment libraries has been automated in the script, but requires secondary structure prediction files which can be obtained from [RaptorX Property](http://raptorx.uchicago.edu/StructurePropertyPred/predict/), [PsiPred PSIPRED 4.0](http://bioinf.cs.ucl.ac.uk/psipred/) or other servers of your choice.
After obtaining these prediction files, the script _Diso_SS2_Reweight_Opt_FFT.py_ on the [AbInitiVO and FastFloppyTail] (https://github.com/jferrie3/AbInitioVO-and-FastFloppyTail) GitHub page can be used to generate reweighted version of these predictions that accounts for predicted disordered probabilities.
Additionally, follow directions for employing various Reweighting methods using disordered probability predictions to ensure correct fragment generation.
#### Sequence BLAST
BLAST is automated within the script to generate useful inputs for fragment library generation. This can be turned off with the -noblast flag. Historically, Legacy BLAST has been used to generate inputs for fragment library curation and can be run by adding -legacyblast to the command-line and is the reccommended approach. This requires extra scripts to be configured within the LegacyBLAST folder.
Regardless of the BLAST strategy used, all BLAST will require download of the non-redundant database.
#### Multimeric Protein Simulation using Symmetry
To simulate a protein that is known to have a multimeric interaction, all FASTAs (-f inputs) and PDBs (-p inputs) are input as monomers. By providing the multimer containing PDB of the multimerizing domain as the -mp input and the monomeric -p supplied PDB as the -md input, symmetry files and an symmetrized output fusion protein will be automatically generated.
It should be noted, that this method is designed to only have a single multimerizing domain as the possible basis for the simulation.
### Outputs of FusionProteinModeler
__PDB file containing the fusion protein structure__
__FASTA file containing the fusion protein sequence__
__Linker file containing regions of the fusion protein that are linkers to sample with FastFloppyLinker__
__Loops file containing regions of the fusion protein that are loops within folded domains to sample with Backrub motions during FastFloppyLinker__
__Fragment library for sampling linker regions with FastFloppyLinker__

## Running FastFloppyLinker
FastFloppyLinker is run using the command-line, where inputs are specified using an arguement parser.
```
run FastFloppyTail.py -inpdb Halo-SNAPf_WorkingFusion.pdb -linkers Halo-SNAPf_WorkingFusion.linkers -loops Halo-SNAPf_WorkingFusion.loops -t_frag Halo-SNAPf_WorkingFusion.200.3mers
```
Each of the parser flags are described below:
```
-in  --Input_FASTA_File  Name of the text file containing the FASTA sequence of the protein of interest. Carot should not be in same line as sequence, UniProt format preferred.
-ftnstruct  --Number_of_FloppyTail_Structures  Number of structures to sample during FloppyTail portion. Default = 400
-t_frag  --Three_Mer_Frag_Library  Name of the file containing the three-mer fragment library generated by disorder corrected method
-linkers  --Linker_File  File containing the starting and ending residue number for the linkers in the fusion protein
-m_num  --Multimer_Number  Number of multimers whose linkers need to be sampled. Only supply for multimers not being run under symmetry
-m_anchor  --Multimer_Anchor  Side of the fusion which forms the multimer. Options N or C for N or C-terminal fusion
-symm  --Symmetry_File  Supply the symmetry file if you want the protein to be run under symmetry (all moves perfectly mirrored)
-cycles  --FloppyTail_Cycles  Number of sampling cycles within each stage of FloppyTail, identical to increase_cycles flag in C++ ClassicAbInitio and AbRelax. Default 1
-refinesubset  --Refine_Subset  Only subjects the lowest X% of structures to Relax refinement where X is specified as input following flag. Default 0
-relnstruct  --Number_of_Relax_Structures  Number of independent full-atom Relax sampling trajectories from a single AbInitio structure. Default 0
-diso  --Disorder_Probability_Prediction_File  File containing per residue disorder probability prediction in RaptorX format. Generally acquired from RaptorX prediciton.
-inpdb  --Input_PDB_File  Name of the text file containing the PDB structure of the protein of interest. All residues are required, missing residues are not constructed
-pymol  --PYMOL  Running with this flag will utilize PyMOL Mover for visualization
-rg  --RG  Ability to provide radius of gyration via Ferrie et. al. JPCB 2020 rg score term
-clusterid --Cluster_ID Batch submit script ID for output record keeping when running on a cluster
-verbose --Verbose Will unmute information coming from PyRosetta
```
### Acquiring the Necessary Inputs
Input structure (-inpdb), 3mer fragment files (-t_frag) and linker segments (-linkers) are generated from the FusionProteinModeler script. The linker file directly informs the FastFloppyLinker script which regions to sample while the fragment library directly provides the associated fragments.

#### Symmetry Specifications
Required symmetry file (-symm) is generated by FusionProteinModeler while the number of monomers within the multimer (-m_num) and the location of the multimerizing domain (-m_anchor). Currently, there is only support for the multimerizing domain to be located at the N or C terminus by modulating the FoldTree for sampling.


