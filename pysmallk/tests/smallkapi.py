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
try:
	import sys
	import numpy as np #must be imported prior to use of pysmallk library
	from pysmallk import smallkapi as sk
	import argparse
except ImportError:
	print 'ImportError: smallkapi test failed'
	raise


args = sk.parser()

sk.load_matrix(filepath=args.matrixfile)

if args.hiernmf2:
	sk.hiernmf2(args.k, dict_filepath=args.dictfile, format=args.format, 
		maxterms=args.maxterms, hiernmf2tolerance=args.tol)
else:
	sk.nmf(args.k, args.algorithm, min_iter=args.miniter, precision=args.outprecision, 
		infile_W=args.infile_W, infile_H=args.infile_H,
		max_threads=args.maxthreads, tol=args.tol)


# always call finalize() when done with analysis; this helps clean up system variables
sk.finalize()
