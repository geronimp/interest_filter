#!/usr/bin/env python
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Joel Boyd"
__copyright__ = "Copyright 2017"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"
__version__ = "0.0.1"

###############################################################################

import logging
import subprocess
import argparse
import yaml
import json
import shutil
import os
from string import Template
from annotate_wrapper import AnnotateWrapper

###############################################################################

debug={1:logging.CRITICAL,
       2:logging.ERROR,
       3:logging.WARNING,
       4:logging.INFO,
       5:logging.DEBUG}

###############################################################################
# Functions

def check_params(args):
	'''
	Runs checks on the to ensure the inputs are valid.
	
	Parameters
	----------
	args 	- Obj. Argparse object. 
	
	Raises
	------
	Exception if any of the tests are failed
	'''


	if not(args.bin or args.fasta):
		raise Exception("No inputs found. Please provide something to either the --bin or the --fasta flags.")

###############################################################################
# Classes
class Checks:

	'''
	Parent class for the Check obejct. Contains template commands for each 
	check type (currently only one).
	'''

	TSV = 'tsv'
	TAB	= '\t'
	CSV = ','
	COLON = ':'
	TAX_SEP = '; '
	
	####### PARAMETER FIELDS #######
	# General
	TYPE 				= 'type'
	NAME 				= 'name'
	SUMMARY 			= 'summary'
	PARAMS 				= 'params'
	INTERPRETATION 		= 'interpretation'
	
	# Check types
	GRAFTM 				= 'graftm'

	# self.GRAFTM
	GPKG   				= 'gpkg'
	EVALUE 				= 'evalue'
	TARGET_ANNOTATION 	= 'target_annotation'

	# Downstream analysis
	HEATMAP_KOS			= 'heatmap_kos'

	param_templates \
			= { GRAFTM: [GPKG, EVALUE] }
	interpretation_templates \
			= { GRAFTM: [TARGET_ANNOTATION] }
	command_templates \
			= { GRAFTM: Template('graftM graft --threads $threads --force --evalue $evalue --forward $genome --output_directory $output --graftm_package $gpkg') }

class Check(Checks):

	def __init__(self, check_file):
		'''
		Parse the input check file.

		Parameters
		----------
		check_file 	- String. Path to .json file containing check information		
		'''
		check_json = json.load(open(check_file))

		for field, entry in check_json.items():
			field = str(field)
			entry = str(entry)

			if field == self.TYPE:
				self.type = entry
			
			elif field == self.NAME:
				self.name = entry
			
			elif field == self.PARAMS:
				self.parameter_dictionary \
					= {str(parameter_field): str(parameter_entry) for 
					   parameter_field, parameter_entry in 
					   check_json[field].items()}
			
			elif field == self.INTERPRETATION:
				self.interpretation_dictionary \
					= {str(interpretation_field): str(interpretation_entry) for 
					   interpretation_field, interpretation_entry in 
					   check_json[field].items()}
			
			elif field == self.SUMMARY:
				self.summary_dictionary \
					= {str(summary_field): str(summary_entry) for 
					   summary_field, summary_entry in 
					   check_json[field].items()}
	
	def cmd(self, genome_path, output, threads):
		'''
		Generate a command to be run by subprocess.

		Parameters
		----------
		genome_path	- Path to a genome bin to be used in the command.
		output 		- Path to the output directory to which the results of
					  the command will be written.
		threads 	- String. Number of threads to use in the command
		Output
		------
		A string of a command that is specific to a particular type of check
		'''

		self.bin_base = os.path.join(output, os.path.basename(genome_path))
		parameters 	  = dict(genome  = genome_path, 
							 output  = self.bin_base, 
							 threads = threads)

		for field in self.param_templates[self.type]:
			parameters[field] = self.parameter_dictionary[field]

		cmd = self.command_templates[self.type].substitute(parameters)
		
		print cmd
		return cmd

	def _interpret_target_annotation(self, results, target_annotations):
		'''
		Read in GraftM read annotations and determine whether they match the 
		target annotations


		Parameters
		----------
		results 			= String. Path to file containing annotation results.
		target_annotations 	= String. A comma separated list of the annotations
							  to look for in the results file
		
		Output
		------
		A dictionary containing the count of observed target annotations for the 
		result file.
		'''

		target_annotations = target_annotations.split(self.CSV)
		result_annotations = None
		counts 			   = {x:0 for x in target_annotations}
		
		for file in os.listdir(results):
			subfile = os.path.join(results, file)
			if os.path.isdir(subfile):
				for file in os.listdir(subfile):
					if file.endswith(self.TSV):
						result_annotations = os.path.join(subfile, file)
	
		if result_annotations:
			for line in open( result_annotations ):
				protein_name, classification = line.strip().split(self.TAB)
				for rank in classification.split(self.TAX_SEP):
					if rank in target_annotations:
						counts[rank]+=1
		return counts

	def interpret(self, results_file):
		'''
		Run through interpretation types for the check, and pass them to specific functions. 
		Currently only target_annotation is supported.

		Parameters
		----------
		results_file - String. Path to file containing the results to interpret.
		
		Output
		------
		A list of interpretation results.
		'''
		results_list = []
		for interpretation in self.interpretation_templates[self.type]:
			
			if interpretation == self.TARGET_ANNOTATION:
				output_result = self._interpret_target_annotation(results_file, self.interpretation_dictionary[interpretation])
				# Filter out genomes that did not have any of the target annotation.
				if max(output_result.values())>0: 
					results_list.append(output_result)
			
			else:
				raise Exception("Interpretation type not recognised.")
		return results_list

class RunInterestingChecks:

	'''
	Runs 'Checks' on genome bins. A check is specified using a json file in
	the following format:

	{
		"name": "", 			A string description of the check.
		"type": "", 			The type of check. Currently only "graftm" is supported
		"params": {}, 			Dictionary containing parameters associated with the check. 
		 						For the GraftM check, only "gpkg" and "evalue" is required.
		"interpretation": {}, 	Dictionary specifying how the output should be interpreted.
		 						For the GraftM check, only "target_annotation" is required.		
		"summary": {}, 			Dictionary specifying how the output of Enrichm should be 
		 						summarised, should a genome pass the check. For the GraftM
		 						check, only 'heatmap_kos' is required.
		"version": 1 			Version of the check 
	 }	

	So a real example, where we want to check for genomes with McrA proteins that fall within 
	the novel clade:

	{
		"name": "novel_mcrA_check", 
		"type": "graftm",
		"params": {"gpkg": "/srv/db/graftm/7//7.57.methylcoenzymem_reductase_novel.mcrA.gpkg",
				   "threads": "1e-10"}, 
		"interpretation": {"target_annotation": "99_Novel"}, 
		"summary": {"heatmap_kos": "data/novel_mcrA_kos.tsv"}, 
		"version": 1
	}

	This format is fairly flexible and should allow for different checks, 
	interpretations, and summary methods be added over time.
	'''

	OUTPUT_YAML     = 'results.yaml'
	OUTPUT_ANNTOATE = 'enrichm_annotate'
	CHECKS 			= 'checks'

	###########################################################################
	###########################################################################

	def __init__(self, output_directory, threads, force):
		'''		
		Parameters
		----------
		output_directory 	- String. Output directory to create
		threads				- String. Number of threads to use for dependency 
							  software
		force 				- Bool. Overwrite a directory with the same name 
							  as the specified output
		'''
		self.output_directory 		 = output_directory
		self.threads 		  		 = threads

		self.output_yaml 	  		 = os.path.join(self.output_directory, self.OUTPUT_YAML)
		self.annotate_wrapper_output = os.path.join(self.output_directory, self.OUTPUT_ANNTOATE)
		self.checks_subdirectory 	 = os.path.join(self.output_directory, self.CHECKS)

		if force:
			if os.path.isdir(self.output_directory):
				logging.warning('Removing directory with same name as designated output directory: %s' \
											 % (self.output_directory))
				shutil.rmtree(self.output_directory)

		# Generate output directory structure
		os.mkdir(self.output_directory) 		# Main output directory
		os.mkdir(self.checks_subdirectory)		# Subdirectory to store the output of files of each check
		os.mkdir(self.annotate_wrapper_output) 	# Subdirectory of enrichM annotations for files that pass checks.

	###########################################################################
	###########################################################################

	def _compile_yaml(self, check_results, output_yaml):
		'''
		Write results dictionary to file in yaml format

		Parameters
		----------
		check_results 	- Dictionary. To be written in yaml format to file
		output_yaml 	- String. File to which the yaml formatted dictionary
						  will be written
		'''

		with open(output_yaml, 'w') as output_yaml_io:
			yaml.dump(check_results, output_yaml_io, default_flow_style=False)

	###########################################################################
	###########################################################################

	def _run_check(self, check, bins_list):
		'''
		Excecute jobs associated with check (e.g. run graftM)

		Parameters
		----------
		check 		- Object. Check object.
		bins_list	- List. List of strings to files contianing genome or 
					  assembly bins.

		Output
		------
		This will probably change with time, but currently it is the path to 
		the key results file output by the commmand run by the check object.
		'''

		check_results = {}
		check_subdir = os.path.join(self.checks_subdirectory, check.name)
		counter = 1

		while os.path.exists(check_subdir):
			new_name = check_subdir + '.' + str(counter)
			logging.warning('check file with duplicate name: %s renaming to %s' % (check.name, new_name))
			check_subdir = new_name
			counter += 1
		
		os.mkdir(check_subdir)

		for bin_path in bins_list:

			logging.info('Working on genome %s' % bin_path)
			subprocess.call(check.cmd(bin_path, check_subdir, self.threads), shell=True)
			
			logging.info('Interpreting results')
			result = check.interpret(check.bin_base)
			
			if result:
				check_results[bin_path] = result

		return check_results

	###########################################################################
	###########################################################################

	def _summary(self, check, results):
		'''
		Run through summarising methods within the check object (e.g. Build a
		summary heatmap)

		Parameters
		----------
		check 		- Object. Check object.
		results 	- String. This will probably change with time, but currently
					  it is the path to the key results file output by the 
					  commmand run by the check object.
		'''
		
		seen = {}

		if Check.HEATMAP_KOS in check.summary_dictionary: 
			heatmap_kos = check.summary_dictionary[Check.HEATMAP_KOS]
			logging.info('Annotating the genomes that passed the check')
			for genome in results.keys():
				ar=AnnotateWrapper(genome, heatmap_kos, self.threads, self.annotate_wrapper_output, check.name)
				if genome in seen:
					logging.info('Skipping annotation for a genome that has already been processed by enrichm: %s' % (genome)) 
					ar.heatmap(ar.self.ko_matrx_output)
				else:
					seen[genome] = [ar.heatmap_output]
					ar.main()

	###########################################################################
	###########################################################################

	def _gather_bins(self, bin, fasta):
		
		'''
		Builds full paths for each genome/assembly in the bin file, and adds
		in any extra 'fastas' that were included.
		Parameters
		----------
		bin 			- String. Path to folder containing genomes to process
		fasta 			- List. List of paths to files containing genomes

		Output
		------
		List of strings to files contianing genome or assembly bins.
		'''

		bin_list = []	
		if bin:
			for bin_file in os.listdir(bin):
				bin_path = os.path.join(bin, bin_file)
				bin_list.append(bin_path)
		if fasta:
			bin_list += fasta

		return bin_list

	###########################################################################
	###########################################################################

	def do(self, bin, fasta, check_file_list, dont_annotate):
		'''
		Iterate through specified 'check files' 

		Parameters
		----------
		bin 			- String. Path to folder containing genomes to process
		fasta 			- List. List of paths to files containing genomes
		check_file_list	- List. List of check files in json format.		
		dont_annotate 	- Bool. Whether to annotate (True) or not (False)
		'''
		
		check_results_dict = {}
		
		bin_list = self._gather_bins(bin, fasta)
		if any(bin_list):
			for check_file in check_file_list:
				check = Check(check_file) # Load check

				logging.info('Running check: %s' % check.name)
				results = self._run_check(check, bin_list)
				check_results_dict[check.name] = results

				if dont_annotate:
					logging.info('Skipping annotation') 
					shutil.rmtree(self.annotate_wrapper_output)
				else:
					logging.info('Summarising results')
					self._summary(check, results)		
			
		logging.info('Writing check results to an output YAML file: %s' % self.output_yaml)
		self._compile_yaml( check_results_dict, self.output_yaml )

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='''RunInterestingChecks = Run Interesting Checks''')

	parser.add_argument('--bin', type=str, help='Directory containing genomes')
	parser.add_argument('--fasta', type=str, nargs='+', help='Fasta file', default=[])	
	parser.add_argument('--checks', type=str, help='Checks', nargs = '+', required=True)
	parser.add_argument('--dont_annotate', action='store_true', help='Dont run EnrichM to annotate genomes that pass the check filters')
	parser.add_argument('--force', action='store_true', help='Overwrite previous run with the same name')
	parser.add_argument('--output_directory', type=str, help='Output alignment', required=True)
	parser.add_argument('--threads', help='Number of threads to use from PROKKA and EnrichM. Default = 5.', default = '5')
	parser.add_argument('--log', help='Output logging information to file', type=str, default=False)
	parser.add_argument('--verbosity', help='1 - 5, 1 being silent, 5 being noisy indeed. Default = 4', type=int, default=4)
	
	args = parser.parse_args()

	if args.log:
		if os.path.isfile(args.log):
			raise Exception("File %s exists" % args.log)
		logging.basicConfig(filename=args.log, level=debug[args.verbosity], format='%(asctime)s %(levelname)s: %(message)s', datefmt='[%Y-%m-%d %I:%M:%S %p]')
	else:
		logging.basicConfig(level=debug[args.verbosity], format='%(asctime)s %(levelname)s: %(message)s', datefmt='[%Y-%m-%d %I:%M:%S %p]')
	
	check_params(args)

	ric = RunInterestingChecks(args.output_directory, args.threads, args.force)
	
	ric.do(args.bin, args.fasta, args.checks, args.dont_annotate)
	