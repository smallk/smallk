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

import cython 
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.string cimport string 
from cython.operator import dereference
import numpy as np
cimport numpy as np


cdef extern from "sparse_matrix_decl.hpp":
    cdef cppclass SparseMatrix[T]:
        SparseMatrix()
        unsigned int Width()
        unsigned int Height()
        void Reserve(const unsigned int height, const unsigned int width, const unsigned int nzmax)
        void BeginLoad()
        void EndLoad()
        void Load(const unsigned int row, const unsigned int col, const T& value)

cdef extern from "random.hpp":
    cdef cppclass Random:
        Random()
        void SeedFromTime()
        double RandomDouble(const double& center, const double& radius)


cdef extern from "sparse_matrix_ops.hpp":
    cdef cppclass RandomSparseMatrix:
        void RandomSparseMatrix(Random& rng, SparseMatrix[double]& A, 
                        const unsigned int nonzeros_per_column,
                        const unsigned int min_height, const unsigned int max_height,
                        const unsigned int min_width, const unsigned int max_width)

cdef extern from "sparse_matrix_io.hpp":
    cdef bool WriteMatrixMarketFile(const string& file_path,
                           const SparseMatrix[double]& S,
                           const unsigned int precision)

cdef extern from "delimited_file.hpp":
    cdef bool WriteDelimitedFile(const double* buffer, 
                        const unsigned int ldim,
                        const unsigned int height, 
                        const unsigned int width,
                        const string& filename,
                        const unsigned int precision,
                        const char DELIM)

#for whatever reason, without the following extern, an ImportError is
#generated complaining about finding __ZN12SparseMatrixIdE5ClearEv
#(SparseMatrix<double>::Clear())
cdef extern from "term_frequency_matrix.hpp":
    pass

@cython.boundscheck(False)
@cython.wraparound(False)

cdef class PyDoubleSparseMatrix:
    #original code uses SparseMatrix[T], but cython doesn't like templates
    cdef SparseMatrix[double] *thisptr
    def __cinit__(self):
        self.thisptr = new SparseMatrix[double]()
    def Width(self):
        return self.thisptr.Width()
    def Height(self):
        return self.thisptr.Height()
    cdef SparseMatrix[double]* getThis(self):
        return self.thisptr
    def Reserve(self, height, width, nzmax):
        self.thisptr.Reserve(height, width, nzmax)
    def BeginLoad(self):
        self.thisptr.BeginLoad()
    def EndLoad(self):
        self.thisptr.EndLoad()
    def Load(self, row, col, value):
        self.thisptr.Load(row, col, value)
        
cdef class PyRandom:
    cdef Random* thisptr
    def __cinit__(self):
        self.thisptr = new Random()
    def SeedFromTime(self):
        self.thisptr.SeedFromTime()
    cdef Random* getThis(self):
        return self.thisptr
    def RandomDouble(self, center, radius):
        return self.thisptr.RandomDouble(center, radius)

def PyUniform(unsigned int n, unsigned int m, double cent, double rad, PyRandom rng):
    cdef unsigned int c = 0
    cdef vector[double] A
    cdef unsigned int r
    while c != n:
        r = 0
        while r != m:
            A.push_back(rng.RandomDouble(cent, rad))
            r += 1
        c += 1
    return A

def PyDenseDiag(unsigned int n, unsigned int m, double cent, double rad, PyRandom rng):
    cdef unsigned int c = 0
    cdef vector[double] A
    cdef unsigned int r
    while c != n:
        r = 0
        while r != m:
            A.push_back(0.0)
            r += 1
        c += 1
    
    c = 0
    while c != n:
        A[c*m + c] = rng.RandomDouble(cent, rad)
        c += 1
    return A

def PyIdentity(unsigned int n, unsigned int m):
    cdef unsigned int c = 0
    cdef vector[double] A
    cdef unsigned int r
    while c != n:
        r = 0
        while r != m:
            A.push_back(0.0)
            r += 1
        c += 1
    
    c = 0
    while c != n:
        A[c*m + c] = 1.0
        c += 1
    return A
    
def PySparseDiag(unsigned int n, double cent, double rad, PyRandom rng):
    S = PyDoubleSparseMatrix()
    S.Reserve(n, n, n)
    S.BeginLoad()
    cdef unsigned int c = 0
    while c != n:
        S.Load(c, c, rng.RandomDouble(cent, rad))
        c += 1
    S.EndLoad()
    return S

def PyOnes(unsigned int n, unsigned int m):
    cdef unsigned int c = 0
    cdef unsigned int r
    cdef vector[double] A
    while c != n:
        r = 0
        while r != m:
            A.push_back(1.0)
            r += 1
        c += 1
    return A

def PyZeros(unsigned int n, unsigned int m):
    cdef unsigned int c = 0
    cdef unsigned int r
    cdef vector[double] A
    while c != n:
        r = 0
        while r != m:
            A.push_back(0.0)
            r += 1
        c += 1
    return A
    
def PySparse(unsigned int n, unsigned int m, unsigned int nz, PyRandom rng):
    S = PyDoubleSparseMatrix()
    RandomSparseMatrix(dereference(rng.getThis()), dereference(S.getThis()), nz, m, m, n, n)
    return S
        
def PyWriteMtxFile(string& filepath, PyDoubleSparseMatrix S, unsigned int precision):
    res = WriteMatrixMarketFile(filepath, dereference(S.getThis()), precision)
    return res

def PyWriteDelimitedFile(vector[double] A, unsigned int ldim, unsigned int height, unsigned int width, string& filename, unsigned int precision):
    res = WriteDelimitedFile(&A[0], ldim, height, width, filename, precision, ',')
    return res
