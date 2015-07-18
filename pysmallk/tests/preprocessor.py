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

#-----------------------------------------------------------------------------
try:
	import sys
	import numpy as np #must be imported prior to use of pysmallk library
	import pysmallk
except ImportError:
	print 'ImportError: preprocessor test failed'
	raise


# This program can be used as a command line tool to run flat NMF clustering.
# It also demonstrates how to use each of the flatclust functions, which can be 
# integrated into other code.

# identify file names
infile = "matrix.mtx"
indict = "dictionary.txt"
indocs = "documents.txt"

p = pysmallk.Preprocessor()

# use the parser function to parse the command line arguments from the user
args = p.parser()

# make output path end in "/" if necessary
if args.indir[-1] != "/":
    args.indir += "/"

pathtofile = args.indir+infile
pathtodict = args.indir+indict
pathtodocs = args.indir+indocs

if args.outdir=="./":
    print "Default output directory: current directory"

outfile = args.outdir+"reduced_matrix.mtx"
outdict = args.outdir+"reduced_dictionary.txt"
outdocs = args.outdir+"reduced_documents.txt"


# load in the user-specified matrix, dictionary, and documents files
p.load_inputmatrix(filepath=pathtofile)
p.load_dictionary(pathtodict)
p.load_documents(pathtodocs)


# use the user-provided inputs to run the preprocessor 
# all keyword arguments are optional 
p.preprocess(maxiter=args.maxiter, docsperterm=args.docs_per_term,
	termsperdoc=args.terms_per_doc, boolean_mode=args.boolean_mode)

# write files to system using user-provided filenames
p.write_output(outfile, outdict, outdocs, precision=args.precision)
