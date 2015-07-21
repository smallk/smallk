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
	import os
	import sys
	import argparse
	import numpy as np
	from scipy.sparse import csc_matrix
	from scipy.io import mmread, mmwrite

	import pysmallk
except ImportError:
	print 'ImportError: preprocessor test failed'
	raise


# define a parser for the purpose of dynamic placement of the data_dir variable
parser = argparse.ArgumentParser(description="Run SmallK via python binding")
parser.add_argument('--indir', action='store', required=False, metavar='data_dir', default="../xdata_data/")
args = parser.parse_args()

p = pysmallk.Preprocessor()

# make output path end in "/" if necessary
if args.indir[-1] != "/":
    args.indir += "/"

infile = "matrix.mtx"
indict = "dictionary.txt"
indocs = "documents.txt"
filetest = "test/reduced_matrix_20news.mtx"
dicttest = "test/reduced_dictionary_20news.txt"
doctest = "test/reduced_documents_20news.txt"

pathtofile = args.indir+infile
pathtodict = args.indir+indict
pathtodocs = args.indir+indocs
pathtofiletest = args.indir+filetest
pathtodicttest = args.indir+dicttest
pathtodoctest = args.indir+doctest

# The constituent parts of a sparse matrix are the nonzero elements, contained in 'data',
# the row indices for the nonzero elements, contained in 'row_indices', and the column
# offsets indicating where each new column begins (should be of size width+1), contained
# in 'col_offsets'.

# Here, we read in the mtx file and convert it to a column-major sparse matrix and extract
# the data, row_indices, and col_offsets values.
matrix = csc_matrix(mmread(pathtofile))
data = matrix.data
row_indices = matrix.indices
col_offsets = matrix.indptr

with open(pathtodict) as dictionary:
    terms = dictionary.read().split("\n")
    terms.pop()

with open(pathtodocs) as docs:
    docids = docs.read().split("\n")
    docids.pop()

height = len(terms)
width = len(docids)
nz = len(data)


#set input values
maxiter = 8
docsperterm = 3
termsperdoc = 5
precision = 4
boolean_mode = 0


# load datasets
p.load_matrix(height=height, width=width, nz=nz, buffer=data, row_indices=row_indices, col_offsets=col_offsets)
p.load_dictionary(dictionary=terms)
p.load_documents(documents=docids)


#preprocess the dataset
p.preprocess(maxiter=maxiter, docsperterm=docsperterm, termsperdoc=termsperdoc, boolean_mode=boolean_mode)

#return the reduced dataset to python
reduced_docs = p.get_reduced_documents()
reduced_dict = p.get_reduced_dictionary()
reduced_scores = p.get_reduced_scores()
reduced_row_indices = p.get_reduced_row_indices()
reduced_col_offsets = p.get_reduced_col_offsets()
reduced_height = len(reduced_dict)
reduced_width = len(reduced_docs)

#print out the reduced dataset to confirm results
outdir = './'
outfile = outdir+"reduced_matrix.mtx"
outdict = outdir+"reduced_dictionary.txt"
outdocs = outdir+"reduced_documents.txt"

with open(outdict, 'w') as out:
	out.writelines("%s\n" % item  for item in reduced_dict)

with open(outdocs, 'w') as out:
	out.writelines("%s\n" % item  for item in reduced_docs)

 
reduced_matrix = csc_matrix((reduced_scores,reduced_row_indices, reduced_col_offsets), shape=(reduced_height, reduced_width), dtype=float)

test_matrix = csc_matrix(mmread(pathtofiletest))

with open(pathtodicttest) as dictionary:
    test_dict = dictionary.read().split("\n")
    test_dict.pop()

with open(pathtodoctest) as docs:
    test_docids = docs.read().split("\n")
    test_docids.pop()

#check dictionary
# if reduced_dict != test_dict:
# 	print 'preprocessor dictionary test failed'
# else:
# 	print 'preprocessor dictionary test passed'

# #check documents
# if reduced_docs != test_docids:
# 	print 'preprocessor documents test failed'
# else:
# 	print 'preprocessor documents test passed'

#check dictionary
if np.allclose(reduced_matrix.data, test_matrix.data, rtol=.01, atol=.01):
	print 'preprocessor matrix test passed'
else:
	print 'preprocessor matrix test failed'


