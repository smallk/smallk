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

import sys
import numpy as np #must be imported prior to use of pysmallk library
import pysmallk 

# This program can be used as a command line tool to generate matrices.
# It also demonstrates how to use each of the matrixgen functions, which
# can be integrated into other code.

m = pysmallk.Matrixgen()

# use the parser function to parse the command line arguments from the user
args = m.parser()

if args.type == 'UNIFORM':
    A = m.uniform(args.height, args.width, center=args.rng_center, radius=args.rng_radius)

elif args.type == 'DENSE_DIAG':
    A = m.densediag(args.height, args.width, center=args.rng_center, radius=args.rng_radius)

elif args.type == 'SPARSE_DIAG':
    A = m.sparsediag(args.width, center=args.rng_center, radius=args.rng_radius)

elif args.type == 'IDENTITY':
    A = m.identity(args.height, args.width)

elif args.type == 'ONES':
    A = m.ones(args.height, args.width)

elif args.type == 'ZEROS':
    A = m.zeros(args.height, args.width)

elif args.type == 'SPARSE':
    A = m.sparse(args.height, args.width, args.nz_per_col)

# Demonstrates how to convert to numpy array for use in existing code:
# matrix = np.array(A)
# matrix.shape = (height,width)

# write the matrix to the filepath specified by the user
m.write_output(args.filename, precision=args.precision)