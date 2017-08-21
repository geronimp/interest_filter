#!/usr/bin/env python
###############################################################################
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program. If not, see <http://www.gnu.org/licenses/>.        #
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
import tempfile
import os 

###############################################################################

class AnnotateWrapper:

    FAA = '.faa'
    ENRICHM_KO_MATRIX = 'ko_frequency_table.tsv'

    def __init__(self, genome, kos_file, threads, output_directory, heatmap_output):
        self.dir_path = os.path.dirname(os.path.realpath(__file__))

        self.genome  = genome
        self.threads = threads
        basename     = os.path.splitext(os.path.basename(self.genome))[0]
        
        self.ko_set  = set()
        for ko in open(kos_file):
            self.ko_set.add(ko.strip())

        if not os.path.isdir(output_directory):
            os.mkdir(output_directory)

        if output_directory:
            base_dir = os.path.join(output_directory, basename)  
        else:
            base_dir = basename

        os.mkdir(base_dir)

        self.enrichm_output     = os.path.join(base_dir, '%s.enrichm.output' % basename)
        self.heatmap_output     = os.path.join(base_dir, '%s.heatmap.pdf' % heatmap_output)
        self.ko_matrix_output   = os.path.join(self.enrichm_output, self.ENRICHM_KO_MATRIX)

    def do_enrichm(self):
        '''
        Call enrichM
        '''
        cmd = 'enrichm annotate --ko --threads %s --output %s --genome_files %s' \
                    % (self.threads, self.enrichm_output, self.genome)
        subprocess.call(cmd, shell=True)
        return self.ko_matrix_output
    
    def do_heatmap(self, ko_matrix_path):
        '''
        Generate a heatmap for kos in the self.ko_set from the enrichM output ko matrix

        Parameters
        ----------
        ko_matrix_path  - String. Path to file containing the full KO matrix to subset 
                          for the heatmap
        '''

        with tempfile.NamedTemporaryFile(suffix='.tsv') as ko_matrix_subset_io:
            original_ko_matrix = open(ko_matrix_path)
            ko_matrix_subset_io.write(original_ko_matrix.readline()) # Write header
            for line in original_ko_matrix:
                ko = line.split('\t')[0]
                if ko in self.ko_set:
                    ko_matrix_subset_io.write(line)
            ko_matrix_subset_io.flush()
            cmd = 'Rscript %s/basic_heatmap.R  -f %s -o %s' \
                        % (self.dir_path, ko_matrix_subset_io.name, self.heatmap_output)
            ### ~ TODO: Remove hard-coded r script. 
            subprocess.call(cmd, shell=True)

    def main(self):
        '''
        Run enrichM on the specified genome bin.
        '''
        logging.info('Running enrichm on genome: %s' % (self.genome))
        ko_matrix_path = self.do_enrichm()

        logging.info('Generating heatmap for %i kos' % len(self.ko_set))
        self.do_heatmap(ko_matrix_path)
