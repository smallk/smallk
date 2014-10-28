# -*- coding: utf-8 -*-
# Copyright 2013,2014 Georgia Institute of Technology.
#
# Licensed under the Apache License, Version 2.0 (the â€œLicenseâ€); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by
# applicable law or agreed to in writing, software distributed under the License
# is distributed on as â€œAS ISâ€ BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the specific language
# governing permissions and limitations under the License.

#-----------------------------------------------------------------------------

import numpy as np #must be imported prior to use of pysmallk library
from pysmallk import smallkapi
import argparse

# define a parser for the purpose of dynamic placement of the data_dir variable
parser = argparse.ArgumentParser(description="Run SmallK via python binding")
parser.add_argument('--data_dir', action='store', required=False, metavar='data_dir', default="../xdata_data")
args = parser.parse_args()

# define input file names
filepath_w = args.data_dir + "/nmf_init_w.csv"
filepath_h = args.data_dir + "/nmf_init_h.csv"
filepath_matrix = args.data_dir + "/reuters.mtx"
dictpath_dict = args.data_dir + "/reuters_dictionary.txt"

# print metadata
print 'Smallk major version:', smallkapi.get_major_version()
print 'Smallk minor version:', smallkapi.get_minor_version()
print 'Smallk patch level:', smallkapi.get_patch_level()
print 'Smallk version string:', smallkapi.get_version_string()

try:
	# verify input values
	smallkapi.load_matrix(filepath=filepath_matrix)
	assert(smallkapi.is_matrix_loaded())
	smallkapi.nmf(8, 'BPP', precision=2, min_iter=10,
		max_iter=10000, tol=0.000001, max_threads=3)

	inputs = smallkapi.get_inputs()

	if (2 != inputs['precision']):
		raise Exception('SetOutputPrecision failed')
	if (10 != inputs['min_iter']):
		raise Exception('SetMinIter failed')
	
	if (10000 != inputs['max_iter']):
		raise Exception('SetMaxIter failed')
	if (0.000001 != inputs['tol']):
		raise Exception('SetNmfTolerance failed')
	if (3 != inputs['max_threads']):
		raise Exception('SetMaxThreads failed')


    # factor the Reuters matrix using BPP and k == 8
    # use initializer matrices for W and H
    
	smallkapi.nmf(8, 'BPP', min_iter=1, precision=6, 
		infile_W=filepath_w, infile_H=filepath_h,
		max_threads=8, tol=0.005)


	# run a hierarchical clustering problem and generate 5 clusters
    # use XML format for the clustering result file
	smallkapi.set_outputformat("XML")
	smallkapi.load_dictionary(dictpath_dict)
	smallkapi.hiernmf2(5)

except Exception as e:
	print e.args

# always call finalize() when done with analysis; this helps clean up system variables
smallkapi.finalize()
