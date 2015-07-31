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

# must import numpy before importing the .so files
import numpy as np
import argparse
import pysmallk

# define a parser for the purpose of dynamic placement of the data_dir variable
parser = argparse.ArgumentParser(description="Run SmallK via python binding")
parser.add_argument('--indir', action='store', required=False, metavar='data_dir', default="../../smallk_data/")
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

# paths to matrix, dictionary, and documents files that have not
# already been preprocessed for use in SmallK
infile = "articles_matrix.mtx"
indict = "articles_dictionary.txt"
indocs = "articles_documents.txt"

pathtofile = args.indir + infile
pathtodict = args.indir + indict
pathtodocs = args.indir + indocs

# instantiate a Preprocessor object
p = pysmallk.Preprocessor()

# load the inputs files
p.load_matrix(filepath=pathtofile)
p.load_dictionary(filepath=pathtodict)
p.load_documents(filepath=pathtodocs)

# preprocess the dataset
p.preprocess()

# return the reduced dataset to python
reduced_docs = p.get_reduced_documents()
reduced_dict = p.get_reduced_dictionary()
reduced_scores = p.get_reduced_scores()
reduced_row_indices = p.get_reduced_row_indices()
reduced_col_offsets = p.get_reduced_col_offsets()
reduced_height = len(reduced_dict)
reduced_width = len(reduced_docs)

# Now let's use the ouptuts of the preprocessor as the inputs to SmallkAPI,
# without first writing them to the filesystem.

# instantiate a SmallkAPI object
sk = pysmallk.SmallkAPI()

# load the preprocessed dataset into smallk
sk.load_matrix(buffer=reduced_scores, row_indices=reduced_row_indices, col_offsets=reduced_col_offsets,
	height=reduced_height, width=reduced_width, nz=len(reduced_scores))

# load the dictionary into smallk
sk.load_dictionary(dictionary=reduced_dict)

# run nmf and write the results to the filesystem
sk.nmf(5, 'BPP')

# using the dataset already loaded, run hierarchical nmf and write the results to the filesystem
sk.hiernmf2(5)


print "*****************************************************"
print "*                                                   *"
print "*        Demonstrating use of Flatclust             *"
print "*                                                   *"
print "*****************************************************"

# generate a random matrix 256 x 256
a = np.random.random((256, 256))

# load in an associated dictionary
pathtodict = args.indir + 'reuters_dictionary.txt'
with open(pathtodict) as dictionary:
    terms = dictionary.read().split("\n")

# instantiate a Flatclust object
f = pysmallk.Flatclust()

# load in the random matrix and run flatclusting
f.load_matrix(matrix=a)
f.load_dictionary(dictionary=terms)
f.cluster(16, algorithm='HALS')
a = f.get_assignments()

print 'ASSIGNMENTS:'
print a

# clean up the environment
sk.finalize()
f.finalize()
