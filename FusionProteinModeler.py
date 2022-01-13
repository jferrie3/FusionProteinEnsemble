## SPECIFIC LOCATIONS ##
VALL_LOCATION = '/mnt/d/VirtualMachines/VBoxShared/vall.jul19.2011.gz'
BLAST_DATABASE_LOCATION = '/mnt/e/Blast/nr'
ATOM_PAIR_CONSTRAINT_LOCATION = '/mnt/d/TG_Protein_Modeling/generate_atom_pair_constraint_file.py'
RENUMBER_PDB_LOCATION= '/mnt/d/TG_Protein_Modeling/renumber_pdb.py'

## IMPORTS ##
### PYROSETTA IMPORTS ###
from pyrosetta import * 
init(extra_options = '-extra_patch_fa /mnt/d/TG_Protein_Modeling/OGQ_SidechainConjugation.txt /mnt/d/TG_Protein_Modeling/TMR_SidechainConjugation.txt')
from pyrosetta_help import Blueprinter
### BIOPYTHON IMPORTS ###
from Bio.Seq import Seq
from Bio.PDB import *
from Bio import pairwise2
from Bio import Align
from Bio.SeqUtils import seq1
from Bio.Blast.Applications import NcbipsiblastCommandline as psiblast
### GENERAL PYTHON IMPORTS ###
import argparse
import os
from os.path import exists
import fileinput
import numpy as np
import multiprocessing
import sys

## ARGUEMENT PARSER ##
parser = argparse.ArgumentParser(description='Program')
parser.add_argument('-f', '--Sequence_FASTAs', action='store', type=str, required=True,
	help='Input the sequences files of each of the input domains/linkers in order form N to C')
parser.add_argument('-p', '--Input_PDBs', action='store', type=str, required=True,
	help='PDB structure each of each non-linker domain')
parser.add_argument('-md', '--Multimer_Domain', action='store', type=str, required=False,
	help='Enter the PDB file name of the monomeric PDB used in the Input_PDBs arguement ')
parser.add_argument('-mp', '--Multimer_PDB', action='store', type=str, required=False,
	help='PDB structure of the multimerizing domain in the multimeric form for symmetry construction')
parser.add_argument('-floppy', '--FastFloppyLinker', action='store_true', required=False,
	help='Running with FastFloppyLinker will automatically perform the FastFloppyLinker simulation after creating the fusion protein')
parser.add_argument('-noblast', '--NoBlast', action='store_true', required=False,
	help='Turn off the PsiBlast search for fragment selection')
parser.add_argument('-legacyblast', '--LegacyBlast', action='store_true', required=False,
	help='Perform the Blast search using the legacy blastgpg')
parser.add_argument('-nofrags', '--NoFrags', action='store_true', required=False,
	help='Turn off the Fragment Picker and use of Fragments in FastFloppyLinker')
parser.add_argument('-frags_ppred', '--PPred_Fragment_File', action='store', required=False,
	help='PsiPred file Required for the Fragment Picker')
parser.add_argument('-frags_rapx', '--RapX_Fragment_File', action='store', required=False,
	help='RaptorX file Required for the Fragment Picker')
parser.add_argument('-o', '--Output_Name', action='store', type=str, required=False,
	help='Name of Output PDB structure NOTE: Do not add .pdb', default='WorkingFusion')
args = parser.parse_args()

if not args.NoFrags:
	if not args.PPred_Fragment_File and not args.RapX_Fragment_File:
		sys.exit()
		
## DEFINED OBJECTS PYTHON ##
### GetSequence: For importing SINGLE sequences from FASTA files ###
def GetSequence(fasta):
	fasta_file = open(fasta, 'r')
	fasta_lines = fasta_file.readlines()
	fasta_counter = 0
	fasta_sequence = ' '
	for fasta_line in fasta_lines:
		if '>' not in fasta_line:
			if fasta_counter == 0:
				if '\n' in fasta_line:
					fasta_sequence = fasta_line.split('\n')[0]
				else:
					fasta_sequence = fasta_line
				fasta_counter = 1
			else:
				if '\n' in fasta_line:
					fasta_sequence = fasta_sequence + fasta_line.split('\n')[0]
				else:
					fasta_sequence = fasta_sequence + fasta_line
	return fasta_sequence	
	
### GetSequences: For importing MULTIPLE sequences from a SINGLE FASTA files ###
def GetSequences(fasta):
	fasta_file = open(fasta, 'r')
	fasta_lines = fasta_file.readlines()
	fasta_counter = 0
	fasta_sequence = ''
	fasta_sequences = []
	for fasta_line in fasta_lines:
		if '>' in fasta_line and len(fasta_sequence) > 1:
			fasta_sequences.append(fasta_sequence)
			fasta_sequence = ''
			fasta_counter = 0
		if '>' not in fasta_line:
			if fasta_counter == 0:
				if '\n' in fasta_line:
					fasta_sequence = fasta_line.split('\n')[0]
				else:
					fasta_sequence = fasta_line
				fasta_counter = 1
			else:
				if '\n' in fasta_line:
					fasta_sequence = fasta_sequence + fasta_line.split('\n')[0]
				else:
					fasta_sequence = fasta_sequence + fasta_line
	fasta_sequences.append(fasta_sequence)				
	return fasta_sequences
	
### THREE LETTER - ONE LETTER DICTIONARY ###
three_one_dictionary = {'C':'CYS', 'D':'ASP', 'S':'SER', 'Q':'GLN', 'K':'LYS',
     'I':'ILE', 'P':'PRO', 'T':'THR', 'F':'PHE', 'N':'ASN', 
     'G':'GLY', 'H':'HIS', 'L':'LEU', 'R':'ARG', 'W':'TRP', 
     'A':'ALA', 'V':'VAL', 'E':'GLU', 'Y':'TYR', 'M':'MET'}
     
### ParseAlignments: Performs alignments between target/scaffold and identifies isertions, deletions, and mutations ###
def ParseAlignments(target_seq,scaffold_seq):
	aligner = Align.PairwiseAligner()
	aligner.open_gap_score = -10
	alignments = aligner.align(target_seq,scaffold_seq)
	formatted_alignment = str(alignments[0])
	print(formatted_alignment)
	formatted_alignment = formatted_alignment.split('\n')
	parsed_target_seq = formatted_alignment[0]
	parsed_scaffold_seq = formatted_alignment[2]
	parsed_encoding = formatted_alignment[1]
	indels = []
	mutations = []
	start_res = 0
	temp_insert = ''
	insert_start = 0
	insert_end = 0
	deletion_start = 0
	deletion_end = 0
	scaffold_seq_numbering = []
	current_seq_numbering = 0
	#### Identifying the starting residue ####
	for encoding_idx,encoding in enumerate(parsed_encoding):
		if encoding != '-' and start_res == 0:
			start_res = encoding_idx
		if parsed_scaffold_seq[encoding_idx] != '-':
			current_seq_numbering += 1
			scaffold_seq_numbering.append(current_seq_numbering)
		else:
			scaffold_seq_numbering.append(0)
	#### Parsing the Alignment Encodings ####	
	for encoding_idx,encoding in enumerate(parsed_encoding):
		if encoding == '-':
			if parsed_target_seq[encoding_idx] == '-': ## Deletion
				if insert_start != 0:
					insert_end = scaffold_seq_numbering[encoding_idx]
					indel_code = 'I'
					indels.append([temp_insert, insert_start, insert_end, indel_code])
					temp_insert = ''
					insert_start = 0
					insert_end = 0 
				if deletion_start == 0:
					deletion_start = scaffold_seq_numbering[encoding_idx]
					continue
				else:
					continue
			else: ## Insertion
				temp_insert += parsed_target_seq[encoding_idx]
				if deletion_start != 0:
					deletion_end = scaffold_seq_numbering[encoding_idx] - 1
					indel_code = 'D'
					indels.append(['X', deletion_start, deletion_end, indel_code])
					deletion_start = 0
					deletion_end = 0
				if insert_start == 0:
					if encoding_idx == 0:
						insert_start = 1
					else:
						insert_start = scaffold_seq_numbering[encoding_idx-1]
					continue
		else:
			if insert_start != 0:
				insert_end = scaffold_seq_numbering[encoding_idx]
				indel_code = 'I'
				indels.append([temp_insert, insert_start, insert_end, indel_code])
				temp_insert = ''
				insert_start = 0
				insert_end = 0
			elif deletion_start != 0:
				deletion_end = encoding_idx
				indel_code = 'D'
				indels.append(['X', deletion_start, deletion_end, indel_code])
				deletion_start = 0
				deletion_end = 0
			if encoding == '.':
				mutations.append([scaffold_seq_numbering[encoding_idx], parsed_scaffold_seq[encoding_idx], parsed_target_seq[encoding_idx]]) 
			if encoding == '|':
				continue
	if deletion_start != 0:
		deletion_end = max(scaffold_seq_numbering)
		indel_code = 'D'
		indels.append(['X', deletion_start, deletion_end, indel_code])			
	if insert_start != 0:
		insert_end = max(scaffold_seq_numbering)
		indel_code = 'I'
		indels.append([temp_insert, insert_start, insert_end, indel_code])
	return indels, mutations	

### LinkerInfo: Extract the sequence, start and end position of the linker in the fusion protein ###
def LinkerInfo(fusionpose, linkerseq):
	aligner = Align.PairwiseAligner()
	aligner.open_gap_score = -5
	fusionpose_seq = Seq(fusionpose.sequence())
	alignments = aligner.align(fusionpose_seq,linkerseq)
	formatted_alignment = str(alignments[0])
	formatted_alignment = formatted_alignment.split('\n')
	parsed_fusion_seq = formatted_alignment[0]
	parsed_linker_seq = formatted_alignment[2]
	parsed_encoding = formatted_alignment[1]
	linker_start_idx = 0
	linker_end_idx = 0
	linker_seq = ''
	for encoding_idx,encoding in enumerate(parsed_encoding):
		if linker_start_idx != 0 and linker_end_idx != 0:
			break
		elif encoding != '-' and linker_start_idx == 0:
			linker_start_idx = encoding_idx - 2
		elif encoding == '-' and linker_start_idx != 0:
			linker_end_idx = encoding_idx + 1
	linker_start = linker_start_idx + 1		
	linker_end = linker_end_idx + 1
	linker_seq = parsed_fusion_seq[linker_start_idx:linker_end_idx]
	return linker_start, linker_end, linker_seq	

### Remodeler: Performs mutations, insertions, deletions using RosettaRemodel
def Remodeler(pose, mutations, indels):
	blueprint = Blueprinter.from_pose(pose)
	## Make mutations
	for mutation in mutations:
		blueprint.mutate(mutation[0], mutation[2])
	## Perform indels
	for indel in indels:
		if indel[3] == 'I':
			if indel[1] == pose.total_residue():
				for indel_res in indel[0][::-1]:
					blueprint.insert_after(indel[2], 'PIKAA '+indel_res)
			else:		
				for indel_res in indel[0]:
					blueprint.insert_before(indel[2], 'PIKAA '+indel_res)
		if indel[3] == 'D':
			blueprint.del_span(indel[1], indel[2])
	print(blueprint)		
	for resrow in blueprint.rows:
		if len(resrow) > 0:
			if resrow[0] == 0:
				resrow[2] = 'L'
			if '[' in resrow[len(resrow)-1]:
				resrow[3] = 'PIKAA '+resrow[1]	
		print(resrow)
	blueprint.set('mutant.blu')
	#remodel_mover = blueprint.get_remodelmover()
	pyrosetta.rosetta.basic.options.set_boolean_option('remodel:quick_and_dirty', True)
	pyrosetta.rosetta.basic.options.set_string_option('remodel:generic_aa', 'A')
	remodel_mover = pyrosetta.rosetta.protocols.forge.remodel.RemodelMover()
	remodel_mover.register_options()
	remodel_mover.redesign_loop_neighborhood(True)
	SF1 = create_score_function("ref2015")
	#SF1.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint, 10.0)
	remodel_mover.fullatom_scorefunction(SF1)
	remodel_mover.apply(pose)
	return pose, blueprint
	
### Fuser ###
def MultiFuser(fasta_sequences, pdb_structures):
	SF2 = create_score_function("ref2015_cart")
	FLpose = Pose()
	### Track the linker regions so they can be passed to FastFloppyTail
	linker_segments = []
	### Append the Linker to the N-term protein and idealize 
	for Npose_idx,Npose in enumerate(pdb_structures):
		if Npose_idx == len(pdb_structures)-1:
			break
		if Npose_idx > 0:
			Npose.assign(FLpose)
		Lpose = pose_from_sequence(fasta_sequences[Npose_idx*2 + 1])
		Cpose = pdb_structures[Npose_idx + 1]
		start_res = Npose.total_residue()
		end_res = start_res + Lpose.total_residue()
		idealize_vector = pyrosetta.rosetta.utility.vector1_unsigned_long()
		linker_segments.append([start_res, end_res])
		for res_i in range(start_res, end_res+1):
			idealize_vector.append(res_i)
		#idealize_vector.append(end_res)
		rm_lower(Lpose.conformation(), 1)
		rm_upper(Npose.conformation(), start_res)
		task = standard_packer_task(Lpose)
		task.restrict_to_repacking()
		SF1 = create_score_function("ref2015")
		packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(SF1, task)
		packer.apply(Lpose)
		movemap = MoveMap()
		movemap.set_bb(True)
		movemap.set_chi(True)
		#movemap.set(pyrosetta.rosetta.core.id.DOF_Type.THETA, True)
		minmover= pyrosetta.rosetta.protocols.minimization_packing.MinMover()
		minmover.movemap(movemap)
		minmover.score_function(SF1)
		minmover.min_type('lbfgs_armijo_nonmonotone') #Call particular minimization type (dfpmin_armijo_nonmonotone) originally linmin
		minmover.tolerance(0.0001)
		minmover.apply(Lpose)
		Npose = pyrosetta.rosetta.protocols.grafting.insert_pose_into_pose(Npose, Lpose, Npose.total_residue(), Npose.total_residue())					
		ft = FoldTree()
		ft.simple_tree(Npose.total_residue())
		Npose.fold_tree(ft)
		pyrosetta.rosetta.protocols.idealize.basic_idealize(Npose, idealize_vector, SF2, 1)
		### Append the C-term protein to the Linker and idealize 
		start_res = Npose.total_residue()
		end_res = start_res + 1
		idealize_vector = pyrosetta.rosetta.utility.vector1_unsigned_long()
		idealize_vector.append(start_res)
		idealize_vector.append(end_res)
		rm_lower(Cpose.conformation(), 1)
		rm_upper(Npose.conformation(), start_res)
		#Npose.append_pose_by_jump(Cpose, start_res)
		Npose = pyrosetta.rosetta.protocols.grafting.insert_pose_into_pose(Npose, Cpose, Npose.total_residue(), Npose.total_residue())
		ft = FoldTree()
		ft.simple_tree(Npose.total_residue())
		Npose.fold_tree(ft)
		pyrosetta.rosetta.protocols.idealize.basic_idealize(Npose, idealize_vector, SF2, 0)
		FLpose.assign(Npose)
	'''	
	## After creating the Fusions Protein	
	SF2 = create_score_function("ref2015_cart")
	movemap_cart = MoveMap()
	movemap_cart.set_bb(False)
	movemap_cart.set_chi(False)
	for linker in linker_segments:
		for res_i in range(linker[0], linker[1] + 1):
			movemap_cart.set_bb(res_i, True)
			movemap_cart.set_chi(res_i, True)
	movemap_cart.set(pyrosetta.rosetta.core.id.DOF_Type.THETA, True)
	minmover_cart = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
	minmover_cart.movemap(movemap_cart)
	minmover_cart.score_function(SF2)
	minmover_cart.min_type('lbfgs_armijo_nonmonotone') #Call particular minimization type (dfpmin_armijo_nonmonotone) originally linmin
	minmover_cart.tolerance(0.0001)
	minmover_cart.cartesian(False)
	minmover_cart.omega(True)
	minmover_cart.apply(FLpose)
	pmm.apply(FLpose)
	'''
	### Dump the fusion protein ###
	FLseq = FLpose.sequence()
	with open(args.Output_Name + '.fasta', 'a') as FLfasta:
		if args.LegacyBlast:
			FLseq2 = '>FusionOutput\n' + FLseq + '\n'
			FLfasta.write(FLseq2)
		else:	
			FLfasta.write(FLseq)
	return linker_segments, FLseq, FLpose

### Fusion3DDomainAlign: Aligns the fusion protein with the multimerizing domain ###
def Fusion3DDomainAlign(FLpdb, Mpdb):
	pdbparser = PDBParser()		
	FLstruct = pdbparser.get_structure('fusion', FLpdb)
	Mstruct = pdbparser.get_structure('multi', Mpdb)
	FLchains = {chain.id:seq1(''.join(residue.resname for residue in chain))for chain in FLstruct.get_chains()}
	Mchains = {chain.id:seq1(''.join(residue.resname for residue in chain))for chain in Mstruct.get_chains()}
	FLseq = FLchains['A']
	Mseq = Mchains['A']
	aligner = Align.PairwiseAligner()
	aligner.open_gap_score = -10 
	alignments = aligner.align(FLseq, Mseq)
	formatted_alignment = str(alignments[0])
	formatted_alignment = formatted_alignment.split('\n')
	parsed_FLpose_seq = formatted_alignment[0]
	parsed_Mpose_seq = formatted_alignment[2]
	parsed_encoding = formatted_alignment[1]
	FL_res = 0
	FL_res_list = []
	M_res = 0
	M_res_list = []
	for encoding_idx, encoding in enumerate(parsed_encoding):
		if parsed_FLpose_seq[encoding_idx] != '-':
			FL_res += 1
		if parsed_Mpose_seq[encoding_idx] != '-':
			M_res += 1
		if encoding == '|':
			if M_res != 1 and M_res != len(Mseq):
				FL_res_list.append(FL_res)
				M_res_list.append(M_res)
			else:
				continue
	FLatoms = []
	FL_sup_atoms = []
	M_sup_atoms = []
	sup = Superimposer()
	for atom in FLstruct[0]['A'].get_atoms():
		FLatoms.append(atom)
		if atom.full_id[3][1] in FL_res_list:
			FL_sup_atoms.append(atom)
	for atom in Mstruct[0]['A'].get_atoms():
		if atom.full_id[3][1] in M_res_list:
			M_sup_atoms.append(atom)
	sup.set_atoms(M_sup_atoms, FL_sup_atoms)
	sup.apply(FLstruct)
	io = PDBIO()
	io.set_structure(FLstruct)
	out_pdb_name = 'AlignedFusionProtein.pdb'
	io.save(out_pdb_name)
	return out_pdb_name

### GenerateSymmFile: Generate symmetry file for running sampling with symmetry ###
#### Note: This only works if the fusion protein has been aligned to the monomeric protein ####
#### Note: Additionally the chain specifier (-a) for the symm generation should match the chain used as the monomer ####
def GenerateSymmFile(multimer_pdb_file, monomer_pdb_file, fusion_pdb_file):
	### Create the symmetry file ###
	cmd = './make_symmdef_file.pl -a A -r 10 -p ' + str(multimer_pdb_file) + ' > ' + str(args.Output_Name + '.symm')
	os.system(cmd)
	### Modify the center of mass to reflect the correct residue for the fusion protein ###
	fusion_pose = pose_from_pdb(fusion_pdb_file)
	monomer_pose = pose_from_pdb(monomer_pdb_file)
	multimer_pose = pose_from_pdb(multimer_pdb_file)
	multimer_number = round(multimer_pose.total_residue()/monomer_pose.total_residue())
	xyz_com = np.array(pyrosetta.rosetta.protocols.geometry.center_of_mass(monomer_pose, 1, monomer_pose.total_residue()))
	closest_dist = 10000.0
	closest_res = 0
	closest_res_xyz = xyz_com
	for res in range(1, fusion_pose.total_residue()):
			xyz_res = np.array(fusion_pose.residue(res).xyz('C'))
			distance = np.linalg.norm(xyz_com - xyz_res)
			if distance < closest_dist:
				closest_dist = distance
				closest_res = res
	with fileinput.FileInput(args.Output_Name + '.symm', inplace=True) as symmfile:
		for line in symmfile:
			if 'COM' in line:
				print(line.replace('COM', str(closest_res)))
			else:
				print(line)
	print('Made Symm File for ' + str(multimer_number) + '-mer Multimer')
	return multimer_number
	
### PerformBlast: Performs a BLAST search on the linker sequence ###
def PerformBlast():
	psiblast_cline = psiblast(cmd = 'psiblast',\
	db = BLAST_DATABASE_LOCATION,\
	query = args.Output_Name + '.fasta',\
	evalue = 0.000001,\
	inclusion_ethresh = 0.000001,\
	max_hsps = 100000,\
	num_alignments = 100000,\
	num_iterations = 2,\
	num_threads = int(multiprocessing.cpu_count()/2),\
	outfmt = 6,\
	out = args.Output_Name+'.psi',\
	out_pssm = args.Output_Name+'.pssm',\
	out_ascii_pssm = args.Output_Name+'.ascii.pssm',\
	save_pssm_after_last_round = True,\
	show_gis = True)
	stdout, stderr = psiblast_cline()
	return 'Blast Complete'

def PickFragments(FLseq, linker_segments):
	## Doing everything through PyRosetta Initialization
	linker_query = ''
	for linker in linker_segments:
		for res_i in range(linker[0], linker[1] + 1):
			linker_query += str(res_i) + ' '
	if args.LegacyBlast:
		init(extra_options='\
		-in::file::vall ' + VALL_LOCATION + ' \
		-in::file::fasta ' + args.Output_Name + '.fasta \
		-frags::picking::query_pos ' + linker_query + '\
		-frags::ss_pred ' + str(args.PPred_Fragment_File) + ' psipred ' + str(args.RapX_Fragment_File) + ' raptorx \
		-out::file::frag_prefix ' + args.Output_Name + ' \
		-frags::describe_fragments ' + args.Output_Name + '.fsc \
		-frags::scoring::config /mnt/d/TG_Protein_Modeling/quota-protocol.wghts \
		-frags::frag_sizes 3 \
		-frags::n_candidates 1000 \
		-frags::n_frags 200 \
		-frags::picking::quota_config_file /mnt/d/TG_Protein_Modeling/quota.def \
		-in::file::checkpoint OUTPUT.chk')
	elif args.NoBlast:
		init(extra_options='\
		-in::file::vall ' + VALL_LOCATION + ' \
		-in::file::fasta ' + args.Output_Name + '.fasta \
		-frags::picking::query_pos ' + linker_query + '\
		-frags::ss_pred ' + str(args.PPred_Fragment_File) + ' psipred ' + str(args.RapX_Fragment_File) + ' raptorx \
		-out::file::frag_prefix ' + args.Output_Name + ' \
		-frags::describe_fragments ' + args.Output_Name + '.fsc \
		-frags::scoring::config /mnt/d/TG_Protein_Modeling/quota_noblast-protocol.wghts \
		-frags::frag_sizes 3 \
		-frags::n_candidates 1000 \
		-frags::n_frags 200 \
		-frags::picking::quota_config_file /mnt/d/TG_Protein_Modeling/quota.def')
	else:	
		init(extra_options='\
		-in::file::vall ' + VALL_LOCATION + ' \
		-in::file::fasta ' + args.Output_Name + '.fasta \
		-frags::picking::query_pos ' + linker_query + '\
		-frags::ss_pred ' + str(args.PPred_Fragment_File) + ' psipred ' + str(args.RapX_Fragment_File) + ' raptorx \
		-out::file::frag_prefix ' + args.Output_Name + ' \
		-frags::describe_fragments ' + args.Output_Name + '.fsc \
		-frags::scoring::config quota-protocol.wghts \
		-frags::frag_sizes 3 \
		-frags::n_candidates 1000 \
		-frags::n_frags 200 \
		-frags::picking::quota_config_file quota.def \
		-in::file::pssm ' + args.Output_Name + '.ascii.pssm')
	
	fragpicker = pyrosetta.rosetta.protocols.frag_picker.FragmentPicker()
	fragpicker.parse_quota_command_line()
	fragpicker.parse_command_line()
	fragment_score_scheme = pyrosetta.rosetta.std.ostringstream()
	fragpicker.show_scoring_methods(fragment_score_scheme)
	fragpicker.quota_protocol()
	fragpicker.save_fragments()
	
	return 'Fragments Picked'

## USERFUL PYROSETTA THINGS ##
pmm = PyMOLMover()
rm_lower = pyrosetta.rosetta.core.conformation.remove_lower_terminus_type_from_conformation_residue
rm_upper = pyrosetta.rosetta.core.conformation.remove_upper_terminus_type_from_conformation_residue

## IMPORT THE PROTEIN STRUCTURES AND SEQUENCES ##
fasta_sequences = []
pdb_structures = []
scaffold_sequences = []

### IMPORT SEQUENCES ###
if len(args.Sequence_FASTAs.split(',')) == 1:
	fasta_sequences = GetSequences(args.Sequence_FASTAs.split(',')[0])
else:	
	for sequence in args.Sequence_FASTAs.split(','):
		fasta_sequences.append(GetSequence(sequence))

### IMPORT STRUCTURES ###
if len(args.Input_PDBs.split(',')) == 1:
	pdb_list_file = open(args.Input_PDBs, 'r')
	pdb_lines = pdb_list_file.readlines()
	for pdb_line in pdb_lines:
		pdb_structures.append(pdb_line.split('\n')[0])
else:
	for pdb in args.Input_PDBs.split(','):
		pdb_structures.append(pdb)

## MODIFY THE IMPORTED PROTEIN STRUCTURES ##
for pdb_idx, pdb in enumerate(pdb_structures):
	### Renumber and perform constraied relax
	renum_pdb = pdb.split('.pdb')[0] + '_Renumbered.pdb'
	dump_pdb_name = pdb.split('.pdb')[0] + '_Remodeled.pdb'
	if exists(dump_pdb_name):
		pdb_structures[pdb_idx] = pose_from_pdb(dump_pdb_name)
		scaffold_sequences.append(str(pdb_structures[pdb_idx].sequence()))
	else:
		os.system('python3 ' + RENUMBER_PDB_LOCATION + ' -p ' + str(pdb) + ' -o ' + str(renum_pdb))
		constraints = str(pdb.split('.pdb')[0]) + '_AtomPairConstraints.cst'
		os.system('python3 '+ ATOM_PAIR_CONSTRAINT_LOCATION + ' ' + str(renum_pdb) + ' ' + str(constraints))
		pdb_structures[pdb_idx] = pose_from_pdb(renum_pdb)
		scaffold_sequences.append(str(pdb_structures[pdb_idx].sequence()))
		SF1 = create_score_function("ref2015")
		SF1.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint, 10.0)
		movemap = MoveMap()
		movemap.set_bb(True)
		movemap.set_chi(True)
		movemap.set_jump(False)
		#movemap.set(pyrosetta.rosetta.core.id.DOF_Type.THETA, True)
		minmover= pyrosetta.rosetta.protocols.minimization_packing.MinMover()
		minmover.movemap(movemap)
		minmover.score_function(SF1)
		minmover.min_type('lbfgs_armijo_nonmonotone') #Call particular minimization type (dfpmin_armijo_nonmonotone) originally linmin
		minmover.tolerance(0.0001)
		fastrelax = pyrosetta.rosetta.protocols.relax.FastRelax()
		fastrelax.min_type('lbfgs_armijo_nonmonotone')
		fastrelax.set_scorefxn(SF1)
		fastrelax.set_movemap(movemap)
		### Adding constraints and performing the relax ###
		set_constraints = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
		set_constraints.constraint_file(str(constraints))
		set_constraints.apply(pdb_structures[pdb_idx])
		fastrelax.apply(pdb_structures[pdb_idx])
		### Extract sequences
		seq = fasta_sequences[pdb_idx*2]
		pose_seq = Seq(pdb_structures[pdb_idx].sequence())
		### Idenify and perform isertions, deletions, mutations.
		indels, mutations = ParseAlignments(seq, pose_seq)
		pdb_structures[pdb_idx], blue = Remodeler(pdb_structures[pdb_idx], mutations, indels)
		### Minimize the structures
		movemap = MoveMap()
		movemap.set_bb(True)
		movemap.set_chi(True)
		movemap.set(pyrosetta.rosetta.core.id.DOF_Type.THETA, False)
		minmover= pyrosetta.rosetta.protocols.minimization_packing.MinMover()
		minmover.movemap(movemap)
		minmover.score_function(SF1)
		minmover.min_type('lbfgs_armijo_nonmonotone') #Call particular minimization type (dfpmin_armijo_nonmonotone) originally linmin
		minmover.tolerance(0.0001)
		#minmover.cartesian(True)
		minmover.omega(False)
		#### Applying the MinMover ####
		minmover.apply(pdb_structures[pdb_idx])
		pdb_structures[pdb_idx].dump_pdb(dump_pdb_name)

## MAKE THE PROTEIN FUSION ##
#FLpose = Pose()
linker_segments = []
loop_segments = []
FLseq = ''
#Flpose, linker_segments = MultiFuser(fasta_sequences, pdb_structures)
linker_segments, FLseq, FLpose = MultiFuser(fasta_sequences, pdb_structures)


## IDENTIFY THE LOOPS
scaffold_FLseq = ''
for seq_idx, seq_item in enumerate(fasta_sequences):
	if seq_idx % 2 == 0:
		scaffold_FLseq += scaffold_sequences[int(seq_idx/2)]
	else:
		scaffold_FLseq += seq_item
		
FLindels, FLmutations = ParseAlignments(scaffold_FLseq, FLseq)
for indel in FLindels:
	if indel[3] == 'D':
		if int(indel[2] - indel[1]) > 4:
			loop_segments.append([indel[1], indel[2]])

### SAVE FINAL OUTPUT POSE ###
FLpose.dump_pdb(args.Output_Name + '.pdb')

## PREPARE FOR SUBSEQUENT SAMPLING ##
### SYMMETRY ###
### Set up the fusion for symmetry based off of template PDB ### 
if args.Multimer_Domain:
	if args.Multimer_Domain == 'N':
		aligned_fusion_pdb = Fusion3DDomainAlign(args.Output_Name + '.pdb', renum_N_pdb)
		multimer_number = GenerateSymmFile(args.Multimer_PDB, renum_N_pdb, aligned_fusion_pdb)
	elif args.Multimer_Domain == 'C':
		aligned_fusion_pdb = Fusion3DDomainAlign(args.Output_Name + '.pdb', renum_C_pdb)
		multimer_number = GenerateSymmFile(args.Multimer_PDB, renum_C_pdb, aligned_fusion_pdb)

### LINKER FLEXIBILITY/FRAGMENTS ###
linker_output = ''
for linker_line in linker_segments:
	for linker_item in linker_line:
		linker_output = linker_output + str(linker_item)
		linker_output = linker_output + ' '
	linker_output = linker_output + '\n'
with open(args.Output_Name + '.linkers', 'a') as linkerfile:
	linkerfile.write(linker_output)
	
### LINKER FLEXIBILITY/FRAGMENTS ###
loops_output = ''
for loops_line in loop_segments:
	for loops_item in loops_line:
		loops_output = loops_output + str(loops_item)
		loops_output = loops_output + ' '
	loops_output = loops_output + '\n'
with open(args.Output_Name + '.loops', 'a') as loopsfile:
	loopsfile.write(loops_output)	

#### GENERATING FRAGMENTS FOR FASTFLOPPYLINKER ####
##### RUNNING BLAST #####
if not args.NoBlast: 
	if not args.NoFrags:
		if args.LegacyBlast:
			run_legacy_blast = '/mnt/e/Blast_for_Checkpoint_File.pl ' + args.Output_Name + '.fasta'
			os.system(run_legacy_blast)
			print("Blast Completed")
		else:	
			bl_com_stat = PerformBlast()
			print(bl_com_stat)
	
##### RUNNING FRAGMENT PICKER #####
if not args.NoFrags:
	pf_com_stat = PickFragments(FLseq, linker_segments)
	print(pf_com_stat)
	
#### RUNNING FASTFLOPPYLINKER ####
if args.FastFloppyLinker:
	run_cmd = 'python ../FastFloppyLinker2.py -inpdb ' + str(args.Output_Name + '.pdb') + ' -linkers ' + str(args.Output_Name + '.linkers') 
	if not args.NoFrags:
		run_cmd += ' -t_frag ' + str(args.Output_Name + '.200.3mers')
	if args.Multimer_Domain:
		run_cmd += ' -symm ' + str(args.Output_Name + '.symm')
		run_cmd += ' -m_anchor ' + args.Multimer_Domain
	os.system(run_cmd)

