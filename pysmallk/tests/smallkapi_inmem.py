# -*- coding: utf-8 -*-
# Copyright 2013,2014 Georgia Institute of Technology.
#
# Licensed under the Apache License, Version 2.0 (the “License”); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by
# applicable law or agreed to in writing, software distributed under the License
# is distributed on as “AS IS” BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the specific language
# governing permissions and limitations under the License.

try:
	import numpy as np # must import numpy before importing the .so files
	import pysmallk
	import sys
	from scipy.sparse import csc_matrix
	import argparse
	from scipy.io import mmread, mmwrite

except ImportError:
	raise
	print 'ImportError: smallkapi test failed'
	raise

sk = pysmallk.SmallkAPI()

# define a parser for the purpose of dynamic placement of the data_dir variable
parser = argparse.ArgumentParser(description="Run SmallK via python binding")
parser.add_argument('--indir', action='store', required=False, metavar='data_dir', default="../xdata_data/")
args = parser.parse_args()

# print metadata
print 'Smallk major version:', sk.get_major_version()
print 'Smallk minor version:', sk.get_minor_version()
print 'Smallk patch level:', sk.get_patch_level()
print 'Smallk version string:', sk.get_version_string()


# make output path end in "/" if necessary
if args.indir[-1] != "/":
    args.indir += "/"

infile = 'reuters.mtx'
w_file = 'nmf_init_w.csv'
h_file = 'nmf_init_h.csv'
w_truth = '/test/nmf_result_w.csv'
h_truth = '/test/nmf_result_h.csv'
indict = 'reuters_dictionary.txt'

pathtofile = args.indir+infile
pathtoW = args.indir + w_file
pathtoH = args.indir + h_file
pathtoWtruth = args.indir + w_truth
pathtoHtruth = args.indir + h_truth
pathtodict = args.indir + indict

matrix = csc_matrix(mmread(pathtofile))
data = matrix.data
row_indices = matrix.indices
col_offsets = matrix.indptr


try:
	# verify input values
	sk.load_matrix(filepath=args.matrixfile)
	assert(sk.is_matrix_loaded())
	sk.nmf(8, 'BPP', precision=2, min_iter=10,
		max_iter=10000, tol=0.000001, max_threads=3)

	inputs = sk.get_inputs()

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

except Exception as e:
	print e.args

#set up for nmf routine
k = 8
algorithm = 'BPP'

sk.load_matrix(buffer=data, row_indices=row_indices, col_offsets=col_offsets,
	height=matrix.shape[0], width=matrix.shape[1], nz=len(data))
sk.nmf(k, algorithm, infile_W=pathtoW, infile_H=pathtoH,
	min_iter=1, precision=6)
h = sk.get_H()
w = sk.get_W()

h_test = np.genfromtxt(pathtoHtruth,delimiter=',')
w_test = np.genfromtxt(pathtoWtruth,delimiter=',')

if np.allclose(w, w_test, rtol=.01, atol=.01):
	print 'W matrix test passed'
else:
	print 'W matrix test failed'

if np.allclose(h, h_test, rtol=.01, atol=.01):
	print 'H matrix test passed'
else:
	print 'H matrix test failed'


with open(pathtodict) as dictionary:
    terms = dictionary.read().split("\n")
    terms.pop()

# sk.hiernmf2(16, dict_filepath=pathtodict)
sk.hiernmf2(5, dictionary=terms)

sk.finalize()

