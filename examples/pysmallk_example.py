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


import sys
import numpy as np # must import numpy before importing the .so files
from scipy.sparse import csc_matrix
import argparse
from scipy.io import mmread, mmwrite
import pysmallk

# define a parser for the purpose of dynamic placement of the data_dir variable
parser = argparse.ArgumentParser(description="Run SmallK via python binding")
parser.add_argument('--indir', action='store', required=False, metavar='data_dir', default="../../../xdata_data/")
args = parser.parse_args()

# make output path end in "/" if necessary
if args.indir[-1] != "/":
    args.indir += "/"


print "*****************************************************"
print "*                                                   *"
print "*           Demonstrating use of the                *"
print "*           Preprocessor and SmallKAPI.             *"
print "*                                                   *"
print "*****************************************************"

infile = "matrix.mtx"
indict = "dictionary.txt"
indocs = "documents.txt"

pathtofile = args.indir+infile
pathtodict = args.indir+indict
pathtodocs = args.indir+indocs

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


p = pysmallk.Preprocessor()
# load datasets to be pruned by the preprocessor
p.load_inputmatrix(height=height, width=width, nz=nz, buffer=data, row_indices=row_indices, col_offsets=col_offsets)
p.load_dictionary(dictionary=terms)
p.load_documents(documents=docids)


#preprocess the dataset
p.preprocess()

#return the reduced dataset to python
reduced_docs = p.get_reduced_documents()
reduced_dict = p.get_reduced_dictionary()
reduced_scores = p.get_reduced_scores()
reduced_row_indices = p.get_reduced_row_indices()
reduced_col_offsets = p.get_reduced_col_offsets()
reduced_height = len(reduced_dict)
reduced_width = len(reduced_docs)


sk = pysmallk.SmallkAPI()
#load the preprocessed dataset into smallk
sk.load_matrix(buffer=reduced_scores, row_indices=reduced_row_indices, col_offsets=reduced_col_offsets,
	height=reduced_height, width=reduced_width, nz=len(reduced_scores))

#run hierarchical nmf on the preprocessed dataset
sk.load_dictionary(dictionary=reduced_dict)
sk.hiernmf2(5)

sk.finalize()


print "*****************************************************"
print "*                                                   *"
print "*        Demonstrating use of Flatclust             *"
print "*                                                   *"
print "*****************************************************"


#generate a random matrix 256 x 256
a = np.random.random((256, 256))

pathtodict = args.indir + 'reuters_dictionary.txt'
with open(pathtodict) as dictionary:
    terms = dictionary.read().split("\n")
    terms.pop()

f = pysmallk.Flatclust()

#load in the random matrix and run flatclusting
f.load_matrix(matrix=a)
f.load_dictionary(dictionary=terms)
f.cluster(16, algorithm='HALS')
f.write_output('assignments', 'fuzzy', 'tree')

f.finalize()

