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
	import numpy as np
	import pysmallk
	from scipy.sparse import csc_matrix
	from scipy.io import mmread, mmwrite
	import argparse
except ImportError:
	print 'ImportError: hierclust test failed'
	raise

h = pysmallk.Hierclust()

# define a parser for the purpose of dynamic placement of the data_dir variable
parser = argparse.ArgumentParser(description="Run SmallK via python binding")
parser.add_argument('--indir', action='store', required=False, metavar='data_dir', default="../xdata_data/")
args = parser.parse_args()

infile = 'reuters.mtx'
indict = 'reuters_dictionary.txt'
assign_truth = '/test/reuters_assignments_12.csv'
initdir = 'test/matrices.reuters/'

# make output path end in "/" if necessary
if args.indir[-1] != "/":
    args.indir += "/"

pathtofile = args.indir+infile
pathtodict = args.indir + indict
pathtoassign = args.indir + assign_truth
pathtoinit = args.indir + initdir

matrix = csc_matrix(mmread(pathtofile))
data = matrix.data
row_indices = matrix.indices
col_offsets = matrix.indptr


h.load_matrix(buffer=data, row_indices=row_indices, col_offsets=col_offsets,
	height=matrix.shape[0], width=matrix.shape[1], nz=len(data))
h.load_dictionary(filepath=pathtodict)

h.cluster(12, initdir=pathtoinit)

assign = h.get_assignments()

# if flat==1:

# 	terms_indices = h.get_flat_top_terms()
# 	print terms_indices

# 	with open(pathtodict) as dictionary:
# 	    dictionary = dictionary.read().split("\n")
# 	    dictionary.pop()

# 	terms = [[x for i, x in enumerate(dictionary) if i == idx] for idx in terms_indices]

assign_test = np.genfromtxt(pathtoassign, delimiter=',', skip_footer=1)

string = 'assignment file test passed'
for i in range(len(assign)):
	if assign_test[i] != assign[i]:
		string = 'assignment file test failed'
		break
print string


h.finalize()