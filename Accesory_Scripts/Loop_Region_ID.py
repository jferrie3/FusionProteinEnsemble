### BIOPYTHON IMPORTS ###
from Bio.Seq import Seq
from Bio.PDB import *
from Bio import pairwise2
from Bio import Align
from Bio.SeqUtils import seq1
from Bio.Blast.Applications import NcbipsiblastCommandline as psiblast
import sys

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
	return formatted_alignment
	#### Identifying the starting residue ####

loop_in_pdb = GetSequence('../Halo-1xIg-SNAPf/Halo-1xIg-SNAPf_WorkingFusion.fasta')
loop_id_pdb = GetSequence(str(sys.argv[1]))

loop_residues = []
loop_file = ('../Halo-1xIg-SNAPf/Halo-1xIg-SNAPf_WorkingFusion.loops')
loops_file = open(loop_file, 'r')
loops_lines = loops_file.readlines()
for loops_line in loops_lines:
	start_res = int(loops_line.split(' ')[0])
	end_res = int(loops_line.split(' ')[1])
	### Need to add the length of a monomer to each of the linker residue numbers to get it to sample all linkers
	for residue in range(start_res, end_res + 1):
		loop_residues.append(residue)

alignments = ParseAlignments(loop_in_pdb, loop_id_pdb)
loop_segments = []
segment_start = 0
segment_length = loop_residues[len(loop_residues)-1] - loop_residues[0]
if '-' in alignments[0]:
	resi_count = 0
	for resi_idx,resi in enumerate(alignments[0]):
		if resi != '-':
			resi_count += 1
		if resi_count == loop_residues[0]:
			segment_start = resi_idx + 1
			break		
	loop_segments.append([segment_start, segment_start + segment_length])
else:
	resi_count = 0
	temp_count = 0
	for resi_idx,resi in enumerate(alignments[2]):
		temp_count += 1
		if resi != '-':
			resi_count += 1
		if temp_count == loop_residues[0] - 1:
			segment_start = resi_count + 1
			break		
	loop_segments.append([segment_start, segment_start + segment_length])

	
loops_output = ''
for loops_line in loop_segments:
	for loops_item in loops_line:
		loops_output = loops_output + str(loops_item)
		loops_output = loops_output + ' '
	loops_output = loops_output + '\n'
with open(sys.argv[2], 'a') as loopsfile:
	loopsfile.write(loops_output)