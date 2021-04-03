'''
@author Issac Chaudry

A utility to check sequencing results.

'''

from Bio.pairwise2 import *
from Bio import SeqIO
from Bio import Align
from Bio import pairwise2
from colorama import Fore, Back, Style, init
from tqdm import tqdm
import math
import os
import sys
import csv
import difflib
import json
import argparse
import os


###############
#Configuration
###############

#Set for the colorama module
init(autoreset=True)

root_path = '/'.join(str(os.path.realpath(__file__)).split('/')[:-1])

#Read configuration from config.json file
with open(root_path + '/config.json', 'r') as in_file:
	file_reader = json.load(in_file)

	#Alignment parameters
	match_score = file_reader['alignment_parameters'][0]['match_score']
	mismatch_score = file_reader['alignment_parameters'][0]['mismatch_score']
	open_gap_score = file_reader['alignment_parameters'][0]['open_gap_score']
	extend_gap_score = file_reader['alignment_parameters'][0]['extend_gap_score']
	target_end_gap_score = file_reader['alignment_parameters'][0]['target_end_gap_score']
	query_end_gap_score = file_reader['alignment_parameters'][0]['query_end_gap_score']
	max_alignments = file_reader['alignment_parameters'][0]['max_alignments']

	#Path to construct library (as a CSV)
	construct_lib_path = file_reader['construct_lib']

	construct_search_threshold = file_reader['construct_search_threshold']


#Path to the folder that was downloaded from Eurofins and holds all the sequencing results.
sequencing_results_path = '1855463/' 


####################
#Accessory Functions
####################

def add_construct():
	'''
	Appends a construct to the construct library designated in the config.json
	'''

	#Path to the csv file to append to the construct_lib
	path = 'none'
	while not os.path.exists(path):
		path = input('Enter path to csv with constructs to add:')

		if not os.path.exists(path):
			print('Path not found. Ctrl+C to exit.')

	#Parse from inputted path
	to_write = {}
	with open(path, 'r') as file:
		file_reader = csv.reader(file)

		for row in file_reader:
			if len(row) != 2:
				continue
			to_write[row[0].strip()] = row[1].strip()

	#Write to the construct lib
	with open(construct_lib_path, 'a') as lib:
		lib.write('\n')
		for cname in to_write.keys():
			lib.write(cname + ',' + to_write[cname])

	print('Done appending.')

def clean_name(name):
	'''
	Cleans the construct name for the search function.
	'''
	return name.lower().replace('-','').replace('_', '').replace(' ', '')

def ask_to_continue():
	print('Press enter to continue...')
	Style.RESET_ALL
	input()


'''
NOT NEEDED! Switched to using .seq files.
'''
def remove_phd_comments(file_name):
	'''
	Removes the comment section from the PHD file because and NULL values will break the biopython parser.

	Parameters
	----------
	file_name: str
		The PHD file to remove the comments for. 

	Returns
	--------
	file_name: str
		The same filename that was passed in.
	'''

	with open(file_name, 'r') as f_reader:

		#Pull all the lines from the PHD file
		lines = f_reader.readlines()

		#A list of all the lines that belong in the new output PHD file
		to_write = []

		#True if writing the current line in the file, False otherwise
		writing = True

		#Go through all of the lines and pull the ones to write in the new file
		for line in lines:
			if writing:
				to_write.append(line)

			#Stop writing once at the start of the comments
			if line.strip() == "BEGIN_COMMENT":
				writing = False

			#Start writing again after the comments section
			if line.strip() == "END_COMMENT":
				to_write.append(line)
				writing = True

		#Write the new lines to same file
		with open(file_name, 'w') as f_writer:
			f_writer.writelines(to_write)

		return file_name

def check_sequencing(sequencing_result, construct_sequence, a_id):
	'''
	This function takes care of most of the pipeline. It checks the sequenceing and prints the alingment (and statistics).

	Parameters
	----------
	sequencing_result: str
		The sequence obtained from Europhins (or other sequencing service)
	construct_sequence: str
		The sequence of the target construct.
	a_id: str
		An identifier for this alignment
	'''

	Style.RESET_ALL

	aligner = Align.PairwiseAligner()
	aligner.mode = 'local'
	aligner.match_score = match_score
	aligner.mismatch_score = mismatch_score
	aligner.open_gap_score = open_gap_score
	aligner.extend_gap_score = extend_gap_score
	aligner.target_end_gap_score = target_end_gap_score
	aligner.query_end_gap_score = query_end_gap_score

	alignments = aligner.align(sequencing_result, construct_sequence)

	if len(alignments) > max_alignments:
		tqdm.write(Fore.RED + 'Too many alignments for ' + a_id + '. Trimming...')
		Style.RESET_ALL
		alignments = [alignments[i] for i in range(max_alignments)]

	for alignment in alignments:
		#Parse out alignement components
		#aligned_sequencing_result = list(str(alignment.format()).split('\n'))[0]
		#aligned_construct_sequence = list(str(alignment.format()).split('\n'))[2]
		#align_chars = list(str(alignment.format()).split('\n'))[1]

		aligned_sequencing_result = list(str(format(alignment)).split('\n'))[0]
		aligned_construct_sequence = list(str(format(alignment)).split('\n'))[2]
		align_chars = list(str(format(alignment)).split('\n'))[1]

		#Calculate and print alignment statistics
		percent_id = 100 * (align_chars.count('|')/len(aligned_construct_sequence.strip()))
		
		percent_gaps = 100 * (align_chars.count('-')/len(aligned_construct_sequence.strip())) 

		percent_mismatches = 100 * (align_chars.count('.')/len(aligned_construct_sequence.strip()))

		query_cover = 100 * (alignment.aligned[1][-1][-1] - alignment.aligned[1][0][0]) / len(construct_sequence)

		print('#'*os.get_terminal_size().columns)
		print('Alignment ID: ', a_id)
		print('Alignment score: ', alignment.score)

		print('Percent ID: ', percent_id, '\tPercent Gap: ', percent_gaps, '\tPercent Mismatches: ', percent_mismatches, '\tQuery Coverage: ', query_cover)
		print('#'*os.get_terminal_size().columns)

		########################################
		#Print out the alignment in clean format
		########################################

		#Determine the start positions for each sequence.
		sequencing_result_aligned_regions = alignment.aligned[0]
		construct_sequence_aligned_regions = alignment.aligned[1]
		sequencing_result_start_pos = sequencing_result_aligned_regions[0][0]
		construct_sequence_start_pos = construct_sequence_aligned_regions[0][0]

		#Strip white space from the sequences and alignment characters
		aligned_sequencing_result_stripped = aligned_sequencing_result.strip()
		aligned_construct_sequence_stripped = aligned_construct_sequence.strip()
		align_chars_stripped = align_chars.strip()


		char_per_line = os.get_terminal_size().columns - 20
		num_lines = math.floor(len(align_chars_stripped)/char_per_line)
		char_last_line = len(align_chars_stripped) % char_per_line

		pos_sequencing_result = sequencing_result_start_pos
		pos_construct_sequence = construct_sequence_start_pos
		pos_alignment_char = 0

		for l in range(num_lines):
			print(Fore.YELLOW + str(pos_sequencing_result) + '->' + str(pos_sequencing_result + char_per_line - 1))
			print(Fore.YELLOW + aligned_sequencing_result_stripped[pos_sequencing_result:pos_sequencing_result+char_per_line])
			print(align_chars_stripped[pos_alignment_char:pos_alignment_char+char_per_line])
			print(Fore.GREEN + aligned_construct_sequence_stripped[pos_construct_sequence:pos_construct_sequence+char_per_line])
			print(Fore.GREEN + str(pos_construct_sequence) + '->' + str(pos_construct_sequence + char_per_line - 1))
			print()
			pos_sequencing_result += char_per_line
			pos_construct_sequence += char_per_line
			pos_alignment_char += char_per_line

		print(Fore.YELLOW + str(pos_sequencing_result) + '->' + str(sequencing_result_aligned_regions[-1][-1]))
		print(Fore.YELLOW + aligned_sequencing_result_stripped[pos_sequencing_result:pos_sequencing_result+char_last_line])
		print(align_chars_stripped[pos_alignment_char:pos_alignment_char+char_last_line])
		print(Fore.GREEN + aligned_construct_sequence_stripped[pos_construct_sequence:pos_construct_sequence+char_last_line])
		print(Fore.GREEN + str(pos_construct_sequence) + '->' + str(construct_sequence_aligned_regions[-1][-1]))

		Style.RESET_ALL
		print('\n')

##################
#Start of pipeline
##################
def main():

	#Store all constructs in a dictionary from the CSV file
	construct_lib = {}

	with open(construct_lib_path, 'r') as file:
		reader = csv.reader(file)
		for row in reader:
			if len(row) != 2:
				continue
			if len(row[0].strip()) != 0:
				construct_lib[row[0].strip()] = row[1].strip()

	#Remove the comments section from each of the PHD files so the Biopython module can parse it
	#seq_results_files = [remove_phd_comments(sequencing_results_path + file) for file in os.listdir(sequencing_results_path) if 'phd' in file]
	seq_results_files = sorted([(sequencing_results_path + file) for file in os.listdir(sequencing_results_path) if 'seq' in file])

	#Iterate through all of the PHD files
	for file in seq_results_files:
		print('Checking ', file.split('/')[-1])

		#Pull construct sequence
		construct_name = input('Enter construct to check against: ')
		construct_sequence = ''

		matched_constructs = []

		for c in construct_lib.keys():
			if construct_name == c:
				construct_sequence = construct_lib[construct_name]

			if clean_name(construct_name) in clean_name(c):
				matched_constructs.append(c)

			if (len([d[0] for d in difflib.ndiff(clean_name(c), clean_name(construct_name)) if d[0].strip() != '']) < 10):
				matched_constructs.append(c)

		if len(construct_sequence) == 0:
			if len(matched_constructs) == 0:
				print(Fore.RED + 'No valid construct found for ' + str(file))
				Style.RESET_ALL
				ask_to_continue()
				continue

			if(len(matched_constructs) > 1):
				got_name = False
				reponse = ''
				while not got_name:
					print(Fore.RED + 'Multiple possible constructs for ' + file + ', ' + construct_name + '. Select correct one or enter "q" to skip.')
					Style.RESET_ALL

					for i in range(len(matched_constructs)):
						print(str(i + 1) + '. ' + matched_constructs[i])

					response = input()

					if response == 'q':
						got_name = True
						continue

					try:
						response = int(response)
					except ValueError:
						print(Fore.RED + 'Invalid selection. Try again or enter "q" to skip.')
						continue

					if int(response) - 1 in range(len(matched_constructs)):
						construct_name = matched_constructs[int(response) - 1]
						got_name = True
						continue
					else:
						print(Fore.RED + 'Invalid selection. Try again or enter "q" to skip.')

				if response == 'q':
					continue


			construct_sequence = construct_lib[construct_name]

			if construct_sequence == '':
				print(Fore.RED + 'No construct found. Skipping...')
				continue

			#Pull the sequence from the phd file
			sequencing_result = SeqIO.read(file, 'fasta').seq.upper()

		#Perform the alingment and print the results
		Style.RESET_ALL
		check_sequencing(sequencing_result, construct_sequence, file + ', ' + construct_name)

	Style.RESET_ALL

if __name__ == '__main__':

	#Set up parser for input
	parser = argparse.ArgumentParser()
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument('-p', '--path', help='the path to the directory holding the sequencing results.')
	group.add_argument('--add_construct', help='add this flag to add a construct to the construct library.', action='store_true')

	args = parser.parse_args()

	if args.add_construct:
		add_construct()

	if args.path:
		sequencing_results_path = args.path
		if sequencing_results_path[-1] is not '/':
			sequencing_results_path += '/'
		main()



	


		






