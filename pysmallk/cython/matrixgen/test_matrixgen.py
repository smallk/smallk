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

import argparse
import time
import libmatrixgen


rng = libmatrixgen.PyRandom()
rng.SeedFromTime()


parser = argparse.ArgumentParser()
parser.add_argument("--height",         action="store", required=True,  metavar="height")
parser.add_argument("--width",          action="store", required=True, metavar="width")
parser.add_argument("--filename",       action="store", required=True, metavar="filename")
parser.add_argument("--type",           action="store", required=False, metavar="type", default='UNIFORM',
                    choices=['UNIFORM', 'DENSE_DIAG', 'SPARSE_DIAG','IDENTITY', 'ONES', 'ZEROS', 'SPARSE'])
parser.add_argument("--rng_center",     action="store", required=False, metavar="rng_center", default=0.5)
parser.add_argument("--rng_radius",     action="store", required=False, metavar="rng_radius", default=0.5)
parser.add_argument("--precision",      action="store", required=False, metavar="precision", default=6)
parser.add_argument("--nz_per_col",     action="store", required=False, metavar="nz_per_col", default=1)

args = parser.parse_args()

#get values from user input
height = int(args.height)
width = int(args.width)
filename = str(args.filename)
center = float(args.rng_center)
radius = float(args.rng_radius)
precision = int(args.precision)
nz = int(args.nz_per_col)
typeInput = str(args.type)

#set time for tracking
t0 = time.time()

is_sparse = False

if typeInput == 'UNIFORM':
    A = libmatrixgen.PyUniform(width, height, center, radius, rng)

elif typeInput == 'DENSE_DIAG':
    A = libmatrixgen.PyDenseDiag(width, height, center, radius, rng)

elif typeInput == 'SPARSE_DIAG':
    is_sparse = True
    S = libmatrixgen.PySparseDiag(width, center, radius, rng)

elif typeInput == 'IDENTITY':
    A = libmatrixgen.PyIdentity(width, height)

elif typeInput == 'ONES':
    A = libmatrixgen.PyOnes(width, height)

elif typeInput == 'ZEROS':
    A = libmatrixgen.PyZeros(width, height)

elif typeInput == 'SPARSE':
    is_sparse = True
    S = libmatrixgen.PySparse(width, height, nz, rng)

if is_sparse:
    if not libmatrixgen.PyWriteMtxFile(filename, S, precision):
        print 'matrixgen error - sparse matrix file write failed'
    else:
        print '{} matrix succesfully written'.format(typeInput)
else:
    if not libmatrixgen.PyWriteDelimitedFile(A, height, height, width, filename, precision):
        print 'matrixgen error - file write failed'
    else: 
        print '{} matrix succesfully written'.format(typeInput)

t1 = time.time()
print "Matrix generation time: ", t1-t0, "s"
