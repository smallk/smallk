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
	from pysmallk import flatclust as f
except ImportError:
	print 'ImportError: flatclust test failed'
	raise

# This program can be used as a command line tool to run flat NMF clustering.
# It also demonstrates how to use each of the flatclust functions, which can be 
# integrated into other code.

# use the parser function to parse the command line arguments from the user
args = f.parser()

# load in the user-specified matrix and dictionary
f.load_matrix(filepath=args.matrixfile)
f.load_dictionary(dictfile=args.dictfile)

# use the user-provided inputs to run flat clustering
# all keyword arguments are optional 
f.cluster(args.clusters, infile_W=args.infile_W, infile_H=args.infile_H,
    algorithm=args.algorithm, min_iter=args.miniter, max_iter=args.maxiter, tol=args.tol,
    verbose=args.verbose, maxterms=args.maxterms, max_threads=args.maxthreads)

# write files to system using user-provided filenames
# assignfile and treefile should not include file extensions, e.g. 'tree_output'
# the library will append the appropriate ending to the filename based on the 'format' parameter
f.write_output(args.assignfile, args.treefile, outdir=args.outdir, format=args.format)

# always call finalize() when done with analysis; this helps clean up system variables
f.finalize()
