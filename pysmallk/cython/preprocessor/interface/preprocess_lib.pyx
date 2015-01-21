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

cdef extern from "sparse_matrix_io.hpp":
    cdef bool LoadMatrixMarketFile(const string& file_path,
                                   SparseMatrix[double]& A, 
                                   unsigned int& height,
                                   unsigned int& width,
                                   unsigned int& nnz)

    #cdef bool WriteMatrixMarketFile(const std::string& file_path,
    cdef bool WriteMatrixMarketFile(const string& file_path,
                           const SparseMatrix[double]& A,
                           const unsigned int precision)

cdef extern from "term_frequency_matrix.hpp":
    cdef cppclass TermFrequencyMatrix:
        TermFrequencyMatrix(const SparseMatrix[double]& S, bool boolean_mode)
        unsigned int Width()
        unsigned int Height()
        bool WriteMtxFile(const string& file_path, const double* scores, 
                          const unsigned int precision)

cdef extern from "preprocess.hpp":
    cdef unsigned int TT = 0xFFFFFFFF
    cdef unsigned int FF = 0

    bool preprocess_tf(TermFrequencyMatrix& A,
                       vector[unsigned int]& term_indices,
                       vector[unsigned int]& doc_indices,
                       vector[double]& scores,
                       const unsigned int MAX_ITER,
                       const unsigned int DOCS_PER_TERM,
                       const unsigned int TERM_PER_DOC)

@cython.boundscheck(False)
@cython.wraparound(False)

cdef class PyDoubleSparseMatrix:
    cdef SparseMatrix[double] *thisptr
    def __cinit__(self):
        self.thisptr = new SparseMatrix[double]()
    def Width(self):
        return self.thisptr.Width()
    def Height(self):
        return self.thisptr.Height()
    cdef SparseMatrix[double]* getThis(self):
        return self.thisptr

cdef class PyDoubleTermFrequencyMatrix:
    cdef TermFrequencyMatrix *thisptr
    def __cinit__(self, PyDoubleSparseMatrix A, bool b_mode):
        self.thisptr = new TermFrequencyMatrix(dereference(A.getThis()), b_mode)
    def Width(self):
        return self.thisptr.Width()
    def Height(self):
        return self.thisptr.Height()
    cdef TermFrequencyMatrix* getThis(self):
        return self.thisptr
    
def PyLoadMatrixMarketFile(string infile, PyDoubleSparseMatrix A): 
    cdef unsigned int height = 0, width = 0, nnz = 0 # assign values to get rid of warning
    cdef bool res = LoadMatrixMarketFile(infile, dereference(A.getThis()), height, width, nnz) 
    if not res:
        return res
    else:
        return (height, width, nnz)


def PyPreprocess_tf(PyDoubleTermFrequencyMatrix A, 
                    unsigned int height, unsigned int width, 
                    const unsigned int MAX_ITER, 
                    const unsigned int DOCS_PER_TERM, 
                    const unsigned int TERM_PER_DOC):
    cdef vector[unsigned int] term_indices
    term_indices.resize(height)
    cdef vector[unsigned int] doc_indices
    doc_indices.resize(width)
    cdef vector[double] scores
	
    cdef bool res = preprocess_tf(dereference(A.getThis()), term_indices, doc_indices, scores,
                             MAX_ITER, DOCS_PER_TERM, TERM_PER_DOC)
    if not res:
        return res
    else:
        return (term_indices, doc_indices, scores)

def PyWriteMtxFile(PyDoubleTermFrequencyMatrix M, const string& file_path, scores, const unsigned int precision):
    cdef vector[double] cscores
    for i in range(0, len(scores)):
        cscores.push_back(scores[i])
    cdef const double *pScore0 = &cscores[0]
    return M.getThis().WriteMtxFile(file_path, pScore0, precision)

#########################################################################################
def WriteStringsToFile(filepath, stringarr, valid_indices, N):
    f = open(filepath, "w")
    for i in range(0, N):
        index = valid_indices[i]
        f.write(stringarr[index]+"\n")
    f.close()
