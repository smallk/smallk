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
	from pysmallk import flatclust as f
	import argparse
except ImportError:
	print 'ImportError: flatclust test failed'
	raise


# define a parser for the purpose of dynamic placement of the data_dir variable
parser = argparse.ArgumentParser(description="Run SmallK via python binding")
parser.add_argument('--indir', action='store', required=False, metavar='data_dir', default="../xdata_data/")
args = parser.parse_args()

infile = 'rnd_256_256.csv'
w_file = 'flatclust_init_w.csv'
h_file = 'flatclust_init_h.csv'
indict = 'reuters_dictionary.txt'
assign_truth = '/test/flatclust_rnd_assignments_16.csv'

# make output path end in "/" if necessary
if args.indir[-1] != "/":
    args.indir += "/"

pathtofile = args.indir+infile
pathtoW = args.indir + w_file
pathtoH = args.indir + h_file
pathtodict = args.indir + indict
pathtoassign = args.indir + assign_truth

matrix = np.loadtxt(pathtofile, delimiter=',', dtype=float)

f.load_matrix(matrix=matrix, column_major=False)
f.load_dictionary(dictfile=pathtodict)


f.cluster(16, infile_W=pathtoW, infile_H=pathtoH, algorithm='HALS', min_iter=1, max_iter=5000)

assign = f.get_assignments()

# terms_indices = f.get_flat_top_terms()

# with open(pathtodict) as dictionary:
#     dictionary = dictionary.read().split("\n")
#     dictionary.pop()

# terms = [[x for i, x in enumerate(dictionary) if i == idx] for idx in terms_indices]

assign_test = np.genfromtxt(pathtoassign, delimiter=',')

string = 'assignment file test passed'
for i in range(len(assign)):
	if assign_test[i] != assign[i]:
		string = 'assignment file test failed'
		break
print string


f.finalize()

