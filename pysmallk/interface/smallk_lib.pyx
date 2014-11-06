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

# nmf_main.pyx: numpy arrays -> extern from "nmf_main.h"
# /Users/barrydrake/development/python/spyder_proj/repositories/xdata2/nmf/cython/interface/nmf_lib.pyx
# 3 steps:
# cython f.pyx  -> f.c
# link: python f-setup.py build_ext --inplace  -> f.so, a dynamic library
# py test-f.py: import f gets f.so, f.fpy below calls fc()

import numpy
cimport numpy as np
import cython

from libc.stdlib cimport malloc, free
from libc.string cimport memcpy, memset
from libc.string cimport strdup, strcpy
from libcpp.string cimport string
from libcpp cimport bool
from cython.operator import dereference

from libcpp.vector cimport vector

import argparse
import sys


cdef extern from "smallk.hpp" namespace "smallk":
    cdef enum Algorithm: 
        MU     # Lee & Seung, multiplicative updating
        HALS   # Cichocki & Pan, hierarchical alternating least squares
        RANK2  # Kuang and Park, rank2 specialization
        BPP    # Kim and Park, Block Principle Pivoting method for inner NNLS

    cdef enum OutputFormat:
        XML    # PG_i / PG_1 (ratio of the ith PG to that of iteration 1)
        JSON   # relative change in the Frobenius norm of W

    ctypedef char *stdStringR "const std::string &"

    void Initialize(int& argc, char **&argv)
    int IsInitialized()
    void Finalize()
    
    # version info
    unsigned int GetMajorVersion()
    unsigned int GetMinorVersion()
    unsigned int GetPatchLevel()
    string GetVersionString()
    unsigned int GetOutputPrecision()
    double GetNmfTolerance()
    unsigned int GetMaxIter()
    unsigned int GetMinIter()
    unsigned int GetMaxThreads()
    void Reset()
    string GetOutputDir()
    unsigned int GetMaxTerms()
    OutputFormat GetOutputFormat()
    double GetHierNmf2Tolerance()



    void LoadMatrix(stdStringR filepath)
    void Nmf(int k, Algorithm algorithm, stdStringR initfile_w, stdStringR initfile_h)


    void SetOutputPrecision(const unsigned int num_digits)
    void SetMinIter(const unsigned int min_iterations)
    void SetNmfTolerance(const double tol)
    void SetMaxIter(const unsigned int max_iterations)
    void SetMaxThreads(const unsigned int mt)
    void SetHierNmf2Tolerance(const double tol)
    void SetMaxTerms(const unsigned int max_terms)
    void SetOutputDir(const string& output_dir)
    bool IsMatrixLoaded()

    void LoadMatrix(double *buffer, 
                    unsigned int ldim, 
                    unsigned int height, 
                    unsigned int width)
    void LoadMatrix(unsigned int height, 
                    unsigned int width,
                    unsigned int nz, 
                    vector[double]& data,
                    vector[unsigned int]& row_indices,
                    vector[unsigned int]& col_offsets)
    const double* LockedBufferW(unsigned int& ldim, 
                                unsigned int& height,
                                unsigned int& width)
    const double* LockedBufferH(unsigned int& ldim,
                                unsigned int& height,
                                unsigned int& width)

    void HierNmf2(const unsigned int num_clusters)
    void SetOutputFormat(const OutputFormat form)
    void LoadDictionary(const string& filepath)
    void LoadDictionary(const vector[string]& terms)

               

cdef extern from "nmf.hpp":
    cdef enum Result:
        OK              = 0
        NOTINITIALIZED  = -1
        INITIALIZED     = -2
        BAD_PARAM       = -3
        FAILURE         = -4
        SIZE_TOO_LARGE  = -5

    cdef enum NmfAlgorithm:
        MU
        BPP
        HALS
        RANK2

    cdef enum NmfProgressAlgorithm:
        PG_RATIO
        DELTA_FNORM

    cdef struct NmfStats:
        NmfStats()
        unsigned long long elapsed_us
        int iteration_count

    cdef struct NmfOptions:
        double tol
        NmfAlgorithm algorithm
        NmfProgressAlgorithm prog_est_algorithm
        int height
        int width
        int k
        int min_iter
        int max_iter
        int tolcount
        int max_threads
        bool verbose
        bool normalize

    void NmfInitialize(int argc, char* argv[])
    Result NmfIsInitialized()
    void NmfFinalize()


cdef extern from "flat_clust.hpp":
    cdef struct FlatClustOptions:
        NmfOptions nmf_opts
        int maxterms
        int num_clusters
        bool verbose
    Result FlatClust(const NmfOptions& options,
                     double* buf_A, int ldim_A,
                     double* buf_W, int ldim_W,
                     double* buf_H, int ldim_H,
                     NmfStats& stats)
    Result FlatClustSparse(const NmfOptions& options,
                           const unsigned int height,
                           const unsigned int width,
                           const unsigned int nz, 
                           const unsigned int* col_offsets,
                           const unsigned int* row_indices,
                           const double* data,
                           double* buf_W, int ldim_W,
                           double* buf_H, int ldim_H,
                           NmfStats& stats)
    void FlatClustSparseTest(const NmfOptions& options,
                           const unsigned int height,
                           const unsigned int width,
                           const unsigned int nz, 
                           const unsigned int* col_offsets,
                           const unsigned int* row_indices,
                           const double* data,
                           double* buf_W, int ldim_W,
                           double* buf_H, int ldim_H,
                           NmfStats& stats)


cdef extern from "sparse_matrix_decl.hpp":
    ctypedef unsigned int* const_uint_ptr "const unsigned int*"
    ctypedef double* const_double_ptr "const double*"

    cdef cppclass SparseMatrix[T]:
        SparseMatrix()
        unsigned int Height()
        unsigned int Width()
        unsigned int Size()

        const unsigned int* LockedColBuffer() 
        const unsigned int* LockedRowBuffer() 
        const T* LockedDataBuffer() # BUG: const T* no good 
        void Reserve(const unsigned int height, const unsigned int width, const unsigned int nzmax)
        void BeginLoad()
        void EndLoad()
        void Load(const unsigned int row, const unsigned int col, const T& value)


cdef extern from "file_loader.hpp":
    bool IsDense(const string& file_path)
    bool IsSparse(const string& file_path)
    bool LoadSparseMatrix(const string& file_path, 
                          SparseMatrix[double]& A,
                          unsigned int& height,
                          unsigned int& width,
                          unsigned int& nz)
    bool LoadDenseMatrix(const string& file_path,
                         vector[double]& data,
                         unsigned int& height,
                         unsigned int& nz)

cdef extern from "random.hpp":
    cdef cppclass Random:
        Random()
        void SeedFromTime()
        double RandomDouble(const double& center, const double& radius)

cdef extern from "matrix_generator.hpp":
    bool RandomMatrix(vector[double]& buf, 
                      const unsigned int height,
                      const unsigned int width,
                      Random& rng,
                      const double rng_center,
                      const double rng_radius)

cdef extern from "delimited_file.hpp":
    bool LoadDelimitedFile(vector[double]& buf,
                           unsigned int& height,
                           unsigned int& width, 
                           const string& filename,
                           const char DELIM)
    bool LoadMatrixArray[T](vector[vector[T]]& buf, 
                            const unsigned int height,
                            const unsigned int width,
                            const string& filename,
                            const char DELIM)

cdef extern from "file_format.hpp":
    cdef enum FileFormat:
        CSV
        XML
        JSON

cdef extern from "assignments.hpp":
    bool WriteAssignmentsFile(const vector[int]& labels, 
                             const string& filepath)

    void ComputeAssignments[T](vector[int]& assignments, 
                       const T* buf_h,
                       const unsigned int ldim_h,
                       const unsigned int k,
                       const unsigned int n)

cdef extern from "terms.hpp":
    void TopTerms(const int maxterms,
             const double* buf_w, 
             const unsigned int ldim,
             const unsigned int height,
             const unsigned int width,
             vector[int]& term_indices)

cdef extern from "flat_clust_output.hpp":
    void FlatClustWriteResults(const string& assignfilepath,
                          const string& resultfilepath,
                          const vector[int]& assignments,
                          const vector[string]& dictionary,
                          const vector[int]& term_indices,
                          const FileFormat form,
                          const unsigned int maxterms,
                          const unsigned int num_docs,
                          const unsigned int num_clusters)


cdef extern from "clust.hpp":
    cdef struct ClustStats:
        ClustStats()
        int nmf_count
        int max_count

    cdef struct ClustOptions:
        NmfOptions nmf_opts
        int maxterms
        double unbalanced
        int trial_allowance
        int num_clusters
        bool verbose
        bool flat
    Result Clust(const ClustOptions& options,
                 double* buf_A, int ldim_A,
                 double* buf_w, double* buf_h,
                 vector[vector[double]]& w_initializers,
                 vector[vector[double]]& h_initializers,
                 vector[int]& assignments, 
                 Tree& tree,
                 ClustStats& stats)
    Result ClustSparse(const ClustOptions& options,
                       const SparseMatrix[double]& A,
                       double* buf_w, double* buf_h,
                       vector[vector[double]]& w_initializers,
                       vector[vector[double]]& h_initializers,
                       vector[int]& assignments,
                       Tree& tree,
                       ClustStats& stats)

cdef extern from "tree.hpp":
    cdef cppclass Tree:
        bool Write(const string& filepath,
                   const FileFormat& form,
                   const vector[string]& dictionary)
        void Print(const vector[string]& dictionary)



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
    cdef bool LoadMatrixMarketFile(const string& file_path,
                                   SparseMatrix[double]& A, 
                                   unsigned int& height,
                                   unsigned int& width,
                                   unsigned int& nnz)
    void LoadSparseMatrix(const unsigned int height, 
                    const unsigned int width,
                    const unsigned int nz, 
                    const vector[double]& data,
                    const vector[unsigned int]& row_indices,
                    const vector[unsigned int]& col_offsets,
                    SparseMatrix[double]& A)


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


cdef extern from "delimited_file.hpp":
    cdef bool WriteDelimitedFile(const double* buffer, 
                        const unsigned int ldim,
                        const unsigned int height, 
                        const unsigned int width,
                        const string& filename,
                        const unsigned int precision,
                        const char DELIM)

cdef extern from "term_frequency_matrix.hpp":
    cdef cppclass TermFrequencyMatrix:
        TermFrequencyMatrix(const SparseMatrix[double]& S, bool boolean_mode)
        unsigned int Width()
        unsigned int Height()
        unsigned int Size()
        bool WriteMtxFile(const string& file_path, const double* scores, 
                          const unsigned int precision)
        const unsigned int* LockedColBuffer()
        const TFData* LockedTFDataBuffer()
    cdef struct TFData:
        TFData()
        unsigned int row
        unsigned int count 


@cython.boundscheck(False)
@cython.wraparound(False)


#----------------------------------------------------------------------------------


#data structures

cdef class TermFreqMatrix:
    cdef TermFrequencyMatrix *thisptr
    def __cinit__(self, Sparse A, bool b_mode):
        self.thisptr = new TermFrequencyMatrix(dereference(A.get()), b_mode)
    def Width(self):
        return self.thisptr.Width()
    def Height(self):
        return self.thisptr.Height()
    cdef TermFrequencyMatrix* get(self):
        return self.thisptr


cdef class Rand:
    cdef Random* thisptr
    def __cinit__(self):
        self.thisptr = new Random()
    def seed_from_time(self):
        self.thisptr.SeedFromTime()
    cdef Random* get(self):
        return self.thisptr
    def double(self, center, radius):
        return self.thisptr.RandomDouble(center, radius)

def random(const unsigned int height, const unsigned int width,
                   Rand r, const double r_center, const double r_radius):
    cdef vector[double] buf
    cdef bool ok = RandomMatrix(buf, height, width, dereference(r.get()), r_center, r_radius)
    return buf

cdef class Sparse:
    cdef SparseMatrix[double] *thisptr
    def __cinit__(self):
        self.thisptr = new SparseMatrix[double]()
    def Height(self):
        return self.thisptr.Height()
    def Width(self):
        return self.thisptr.Width()
    def Size(self):
        return self.thisptr.Size()

    cdef const unsigned int* LockedColBuffer(self):
        return self.thisptr.LockedColBuffer()
    cdef const unsigned int* LockedRowBuffer(self):
        return self.thisptr.LockedRowBuffer()
    cdef const double* LockedDataBuffer(self):
        return self.thisptr.LockedDataBuffer()

    cdef SparseMatrix[double]* get(self):
        return self.thisptr
    def Reserve(self, height, width, nzmax):
        self.thisptr.Reserve(height, width, nzmax)
    def BeginLoad(self):
        self.thisptr.BeginLoad()
    def EndLoad(self):
        self.thisptr.EndLoad()
    def Load(self, row, col, value):
        self.thisptr.Load(row, col, value)


cdef class TermFreqData:
    cdef TFData* thisptr
    def __cinit__(self):
        self.thisptr = <TFData*>malloc(sizeof(TFData))
        self.thisptr.row = 0
        self.thisptr.count = 0
    cdef TFData* get(self):
        return self.thisptr
    property row:
        def __get__(self):
            return self.thisptr.row
    property count:
        def __get__(self):
            return self.thisptr.count


cdef class Stats:
    cdef NmfStats* thisptr
    def __cinit__(self):
        self.thisptr = <NmfStats*>malloc(sizeof(NmfStats))
        self.thisptr.elapsed_us = 0
        self.thisptr.iteration_count = 0
    cdef NmfStats* get(self):
        return self.thisptr
    property elapsed_us:
        def __get__(self):
            return self.thisptr.elapsed_us
    property iteration_count:
        def __get__(self):
            return self.thisptr.iteration_count

cdef class TreeResults:
    cdef Tree* thisptr
    def __cinit__(self):
        self.thisptr = new Tree()

    def write(self, const string& filepath, const FileFormat& form, const vector[string]& dictionary):
        return self.thisptr.Write(filepath, form, dictionary)

    cdef Tree* get(self):
        return self.thisptr

cdef class ClusterStats:
    cdef ClustStats* thisptr
    def __cinit__(self):
        self.thisptr = <ClustStats*>malloc(sizeof(ClustStats))
        self.thisptr.nmf_count = 0
        self.thisptr.max_count = 0
    cdef ClustStats* get(self):
        return self.thisptr
    property nmf_count:
        def __get__(self):
            return self.thisptr.nmf_count
    property max_count:
        def __get__(self):
            return self.thisptr.max_count

#----------------------------------------------------------------------------------


# Helper functions

def is_dense(string infile):
    return IsDense(infile)

def is_sparse(string infile):
    return IsSparse(infile)

def get_alg(alg_name):
    if (alg_name == 'MU'):
        return 0
    elif (alg_name == 'HALS'):
        return 2
    elif (alg_name == 'RANK2'):
        return 3
    elif (alg_name == 'BPP'):
        return 1

def get_outputformat(form):
    if form == "XML":
        return XML
    elif form == "JSON":
        return JSON

#----------------------------------------------------------------------------------


# Write output functions

def write_mtx(string& filepath, Sparse S, unsigned int precision):
    res = WriteMatrixMarketFile(filepath, dereference(S.get()), precision)
    return res

def write_delimited(vector[double] A, string& filename, unsigned int precision, unsigned int ldim, unsigned int height, unsigned int width):
    res = WriteDelimitedFile(&A[0], ldim, height, width, filename, precision, ',')
    return res

def write_termfreq(TermFreqMatrix M, const string& file_path, scores, const unsigned int precision):
    cdef vector[double] cscores
    for i in range(0, len(scores)):
        cscores.push_back(scores[i])
    cdef const double *pScore0 = &cscores[0]
    return M.get().WriteMtxFile(file_path, pScore0, precision)

def write_strings(filepath, stringarr, valid_indices, N):
    f = open(filepath, "w")
    for i in range(0, N):
        index = valid_indices[i]
        f.write(stringarr[index]+"\n")
    f.close()

#----------------------------------------------------------------------------------

# Loading input

def load_matrix_internal(filepath=None, height=None, width=None, delim=None, num_initializers=None,
                buffer=None, nz=None, row_indices=None, col_offsets=None, matrix=None, column_major=False):
    if filepath != None:
        sparse = is_sparse(filepath)
        dense = is_dense(filepath)

    if filepath != None and sparse:
        m, h, w = load_sparse_file(filepath, Sparse())
    elif num_initializers != None:
        if delim != None:
            m, h, w = load_matrixarray(height, width, filepath, num_initializers, DELIM=delim)
        else:
            m, h, w = load_matrixarray(height, width, filepath, num_initializers)
    elif row_indices != None and col_offsets != None:
        print 'about to call python load_sparse'
        m, h, w = load_sparse(height, width, nz, buffer, row_indices, col_offsets, Sparse())
    elif buffer != None and height != None and width != None:
        m, h, w = buffer, height, width
    elif matrix != None:
        m, h, w = load_numpy(matrix, column_major)
    elif height != None and width != None and filepath != None:
        if delim != None:
            m, h, w = load_delimited(height, width, filepath, DELIM=delim)
        else:
            m, h, w = load_delimited(height, width, filepath)
    elif filepath != None and dense:
        m, h, w = load_dense_file(filepath)
    return m, h, w


def load_delimited(unsigned int height, unsigned int width,
                        const string& filename, const char DELIM = ','):
    cdef vector[double] buf
    cdef bool ok = LoadDelimitedFile(buf, height, width, filename, DELIM)
    return buf, height, width

def load_sparse_file(const string& file_path, Sparse A):
    cdef unsigned int m = 0, n = 0, nnz = 0
    cdef bool ok = LoadSparseMatrix(file_path, dereference(A.get()), m, n, nnz)
    return A, m, n

def load_sparse(const unsigned int height, 
                    const unsigned int width,
                    const unsigned int nz, 
                    const vector[double]& data,
                    const vector[unsigned int]& row_indices,
                    const vector[unsigned int]& col_offsets, Sparse A):
  LoadSparseMatrix(height, width, nz, data, row_indices, col_offsets, dereference(A.get()))
  return A, height, width

def load_dense_file(const string& file_path):
    cdef unsigned int m = 0, n = 0
    cdef vector[double] buf_a
    cdef bool ok = LoadDenseMatrix(file_path, buf_a, m, n)
    return buf_a, m, n

def load_matrixarray(const unsigned int matrix_height, const unsigned int matrix_width, 
                      const string& filename, unsigned int size, const char DELIM = ','):
    cdef vector[vector[double]] buf
    buf.resize(size)
    cdef bool ok = LoadMatrixArray[double](buf, matrix_height, matrix_width, filename, DELIM)
    return buf, matrix_height, matrix_width

def validate_load(np.ndarray[double, ndim=2, mode="c"] matrix): # matrix is typed as numpy array of doubles
    return matrix.reshape(1, matrix.shape[0]*matrix.shape[1])[0], matrix.shape[0], matrix.shape[1]

def load_numpy(matrix, column_major):
    if column_major:
        matrix_column_major = matrix
    else:
        matrix_column_major = matrix.T
    while True:
        try:
            D, height, width = validate_load(matrix_column_major) # 2, buffer numpy array; see 1. above
            return D, height, width
        except ValueError:
            if not matrix_column_major.flags['C_CONTIGUOUS']:
                print 'Array not C-contiguous. Rebuilding with contiguous memory allocation...'
                matrix_column_major = numpy.ascontiguousarray(matrix_column_major)
            else:
                print 'Array not float or double type. Rebuilding with correct type...'
                matrix_column_major = matrix_column_major.astype(numpy.float)


#----------------------------------------------------------------------------------

# Initialize and finalize functions

def init():
    av = sys.argv
    ac = len(av)
    cdef char **c_arr = <char**>malloc((ac+1) * sizeof(char*))
    cdef char* c_string
    cdef char* char_str

    # argv[argc] must be null terminated 
    memset(c_arr, '\0', (ac+1)*sizeof(char*))
    for i in xrange(0, ac):
        py_string = av[i]
        # convert to byte
        py_byte_string = py_string.encode('ascii')
        # convert to char*
        char_str = py_byte_string
        # Now allocate memory
        data_len = len(char_str) + 1
        c_arr[i] = <char *>malloc(data_len)
        # Null-terminate the string
        memset(c_arr[i], '\0', data_len*sizeof(char))
        # Copy it
        memcpy(c_arr[i], char_str, data_len*sizeof(char))

    NmfInitialize(ac, c_arr)
    free(c_arr)

def isinit():
    return NmfIsInitialized() == INITIALIZED

def finalize():
    NmfFinalize()


#----------------------------------------------------------------------------------


cdef class SmallkAPI:

    def __init__(self):
        self.init()
        if not isinit():
            print 'ERROR'
 
    # Creates a parser that takes the same inputs and has the same defaults as the binary smallk tool.
    # Returns a dictionary containing the command line inputs or their default values. The fields in 
    # that dictionary are: matrixfile, kval, algorithm, tol_val, infile_W, infile_H, outdir, outprecision, 
    # maxiter, miniter, maxthreads.
    def parser(self):
        parser = argparse.ArgumentParser(description="Run NMF via python binding")
        parser.add_argument('--matrixfile', action='store', required=True, metavar='matrixfile')
        parser.add_argument('--k', action='store', required=True, type=int, metavar='k')
        parser.add_argument('--dictfile', action='store', required=False, metavar='dictfile')
        parser.add_argument('--hiernmf2', action='store', required=False, metavar='hiernmf2', default=0, choices=[0,1])
        parser.add_argument('--algorithm', action='store', required=False, default='BPP', metavar='algorithm', choices=['MU','HALS','RANK2','BPP'])
        parser.add_argument('--stopping', action='store', required=False, metavar='stopping', default='PG_RATIO', choices=['PG_RATIO','DELTA'])
        parser.add_argument('--tol', action='store', type=float, required=False, metavar='tol', default=0.005)
        parser.add_argument('--tolcount', action='store', type=int, required=False, metavar='tolcount', default=1)
        parser.add_argument('--infile_W', action='store', required=False, metavar='infile_W', default="")
        parser.add_argument('--infile_H', action='store', required=False, metavar='infile_H', default="")
        parser.add_argument('--outfile_W', action='store', required=False, metavar='outfile_W', default="w.csv")
        parser.add_argument('--outfile_H', action='store', required=False, metavar='outfile_H', default="h.csv")
        parser.add_argument('--outprecision', action='store', type=int, required=False, metavar="outprecision", default=6)
        parser.add_argument('--maxiter', action='store', type=int, required=False, metavar="maxiter", default=5000)
        parser.add_argument('--miniter', action='store', type=int, required=False, metavar="miniter", default=5)
        parser.add_argument('--maxthreads', action='store', type=int, required=False, metavar="maxthreads", default=8)
        parser.add_argument('--maxterms', action='store', type=int, required=False, metavar="maxterms", default=5)
        parser.add_argument('--normalize', action='store', type=int, required=False, metavar="normalize", default=1)
        parser.add_argument('--verbose', action='store', type=int, required=False, metavar="verbose", default=1)

        args = parser.parse_args()
        return args


    def get_outputformat(self, format):
        if format == "XML":
            return 0
        elif format == "JSON":
            return 1

    def get_major_version(self):
        return GetMajorVersion()

    def get_minor_version(self):
        return GetMinorVersion()

    def get_patch_level(self):
        return GetPatchLevel()

    def get_version_string(self):
        return GetVersionString()


    # Load the input matrix. Four combinations are available, which provide for 
    # loading via a file path, direct loading of a dense matrix, direct loading 
    # of a sparse matrix, or loading a numpy matrix. For example:

    # File path: load_matrix(filepath=’/path/to/file/’)
    # Sparse matrix: load_matrix(height=matrix_height, width=matrix_width, nz=non_zero_count, buffer=non_zero_elements, row_indices=rows, col_offsets=cols)
    # Dense matrix: load_matrix(buffer=matrix_buffer, height=matrix_height, width=matrix_width)
    # Numpy matrix: load_matrix(matrix=matrix)
    def load_matrix(self, filepath=None, height=None, width=None, delim=None, buffer=None,
                matrix=None, nz=None, row_indices=None, col_offsets=None, column_major=False):

        if row_indices != None and col_offsets != None:
            self.load_sparse_buffer(height, width, nz, buffer, row_indices, col_offsets)
        elif buffer != None and height != None and width != None:
            self.load_dense_buffer(buffer, height, width)
        elif matrix != None:
            self.load_numpy(matrix, column_major)
        if filepath != None:
            self.load_matrix_file(filepath)

    def load_matrix_file(self, char* filepath):
        LoadMatrix(filepath)

    def load_dense_buffer(self, vector[double]& buffer, unsigned int height, unsigned int width):
        LoadMatrix(&buffer[0], height, height, width)
      
    def load_sparse_buffer(self, unsigned int height, 
                      unsigned int width,
                      unsigned int nz, 
                      vector[double]& buffer,
                      vector[unsigned int]& row_indices,
                      vector[unsigned int]& col_offsets):
        LoadMatrix(height, width, nz, buffer, row_indices, col_offsets)

    def validate_load(self, np.ndarray[double, ndim=2, mode="c"] matrix): # matrix is typed as numpy array of doubles
        LoadMatrix(&matrix[0,0], matrix.shape[0], matrix.shape[0], matrix.shape[1])


    def load_numpy(self, matrix, column_major):
        if column_major:
            matrix_column_major = matrix
        else:
            matrix_column_major = matrix.T

        while True:
            try:
                self.validate_load(matrix_column_major) # 2, buffer numpy array; see 1. above
                break
            except ValueError:
                if not matrix_column_major.flags['C_CONTIGUOUS']:
                    print 'Array not C-contiguous. Rebuilding with contiguous memory allocation...'
                    matrix_column_major = numpy.ascontiguousarray(matrix_column_major)
                else:
                    print 'Array not float or double type. Rebuilding with correct type...'
                    matrix_column_major = matrix_column_major.astype(numpy.float)


    # Runs NMF on the loaded matrix using the supplied algorithm and implementation details.
    def nmf(self, unsigned int k, algorithm, char* infile_W="", char* infile_H="",
            precision=4, min_iter=5, max_iter=5000,
            tol=0.005, max_threads=8, outdir="."):      
        #set inputs
        SetOutputPrecision(precision)
        SetMinIter(min_iter)
        SetMaxIter(max_iter)
        SetNmfTolerance(tol)
        SetMaxThreads(max_threads)
        SetOutputDir(outdir)

        Nmf(k, get_alg(algorithm), infile_W, infile_H)

    # Returns a dictionary of the supplied inputs to the nmf function.
    def get_inputs(self):
        inputs = {}
        inputs['precision'] = GetOutputPrecision()
        inputs['min_iter'] = GetMinIter()
        inputs['max_iter'] = GetMaxIter()
        inputs['tol'] = GetNmfTolerance()
        inputs['max_threads'] = GetMaxThreads()
        inputs['outdir'] = GetOutputDir()
        inputs['format'] = "XML" if GetOutputFormat() == 0 else "JSON"
        return inputs

    def is_matrix_loaded(self):
        return IsMatrixLoaded()


    def init(self):
        av = sys.argv
        ac = len(av)

        cdef char **c_arr = <char**>malloc((ac+1) * sizeof(char*))
        #cdef char **c_arr = <char**>malloc((ac+1) * sizeof(char))
        cdef char* c_string
        cdef char* char_str

        # argv[argc] must be null terminated 
        memset(c_arr, '\0', (ac+1)*sizeof(char*))
        #memset(c_arr, '\0', (ac+1)*sizeof(char))
        for i in xrange(0, ac):
        #for i in xrange(ac):
            py_string = av[i]

            # convert to byte
            py_byte_string = py_string.encode('ascii')

            # convert to char*
            char_str = py_byte_string

            # Now allocate memory
            data_len = len(char_str) + 1
            c_arr[i] = <char *>malloc(data_len)

            # Null-terminate the string
            memset(c_arr[i], '\0', data_len*sizeof(char))
            #memset(c_arr[i], '\0', (data_len+1)*sizeof(char))

            # Copy it
            memcpy(c_arr[i], char_str, data_len*sizeof(char))

        # Initialize
        Initialize(ac, c_arr)

        # free up all allocated memory
        for i in xrange(ac):
            #del(c_arr[i]) # added by bld
            free(c_arr[i])
        #del(c_arr) # added by bld
        free(c_arr)

    def isinit(self):
        return IsInitialized();


    def finalize(self): 
        Finalize()

    # Returns the output H matrix as a numpy array.
    def get_H(self):
        cdef unsigned int ldim = 0
        cdef unsigned int height = 0
        cdef unsigned int width = 0
        h_pointer = LockedBufferH(ldim, height, width)
        h_all = []
        for i in range(height*width):
            h_all.append(h_pointer[i])
        h_all = numpy.array(h_all).reshape(width, height)
        return h_all.T

    # Returns the output W matrix as a numpy array.
    def get_W(self):
        cdef unsigned int ldim = 0
        cdef unsigned int height = 0
        cdef unsigned int width = 0
        w_pointer = LockedBufferW(ldim, height, width)
        w_all = []
        for i in range(height*width):
            w_all.append(w_pointer[i])
        w_all = numpy.array(w_all).reshape(width, height)
        return w_all.T

    # Runs HierNMF2 on the loaded matrix, using the provided running parameters. 
    # There are two options for loading the dictionary: passing a filepath or a 
    # list containing the dictionary.
        
    # File path: hiernmf2(5, dict_filepath=’/path/to/dictionary/’)
    # List: hiernmf2(5, dictionary=[list, containing, terms, from, dictionary])
    def hiernmf2(self, num_clusters, dict_filepath='', dictionary=None, 
        format="XML", maxterms=5, hiernmf2tolerance=0.0001):
        if dict_filepath != '':
            self.load_dictionary(dict_filepath)
        elif dictionary != None:
            self.load_dictionary_from_buffer(dictionary)

        SetHierNmf2Tolerance(hiernmf2tolerance)
        SetMaxTerms(maxterms)
        SetOutputFormat(self.get_outputformat(format))
        HierNmf2(num_clusters)
  

    def load_dictionary(self, const string filepath):
        LoadDictionary(filepath)

    def load_dictionary_from_buffer(self, vector[string]& dictionary):
        LoadDictionary(dictionary)

#----------------------------------------------------------------------------------


cdef class Clustering:

    cdef dict opts
    cdef Sparse matrix
    cdef vector[double] dense_matrix
    cdef TreeResults tree
    cdef ClusterStats stats
    cdef Stats nmfstats

    cdef vector[double] buf_w, buf_h
    cdef ClustOptions clust_opts
    cdef unsigned int height, width, k, required_size, num_initializers
    cdef Rand rng
    cdef vector[vector[double]] w_initializers, h_initializers
    cdef vector[double] w_init, h_init
    cdef vector[int] term_indices, assignments, assignments_flat
    cdef vector[double] w, h
    cdef np.ndarray dictionary
    cdef unsigned int maxterms
    cdef bool sparse

    def __init__(self):
        init()
        if not isinit():
            print 'ERROR'


    def get_alg(self, alg_name):
        if (alg_name == 'MU'):
            return 0
        elif (alg_name == 'HALS'):
            return 1
        elif (alg_name == 'RANK2'):
            return 2
        elif (alg_name == 'BPP'):
            return 3

    def finalize(self):
        finalize()

    # Load the input matrix. Four combinations are available, which provide for 
    # loading via a file path, direct loading of a dense matrix, direct loading 
    # of a sparse matrix, or loading a numpy matrix. For example:

    # File path: load_matrix(filepath=’/path/to/file/’)
    # Sparse matrix: load_matrix(height=matrix_height, width=matrix_width, nz=non_zero_count, buffer=non_zero_elements, row_indices=rows, col_offsets=cols)
    # Dense matrix: load_matrix(buffer=matrix_buffer, height=matrix_height, width=matrix_width)
    # Numpy matrix: load_matrix(matrix=matrix)
    def load_matrix(self, **kwargs):
        #assuming is sparse for now
        if 'filepath' in kwargs:
            if is_sparse(kwargs['filepath']):
                self.matrix, self.height, self.width = load_matrix_internal(filepath=kwargs['filepath'])
                self.sparse = True
            else:
                self.dense_matrix, self.height, self.width = load_matrix_internal(filepath=kwargs['filepath'])

        else:
            if 'row_indices' in kwargs:
                self.sparse = True
                self.matrix, self.height, self.width = load_matrix_internal(**kwargs)
            else:
                self.sparse = False
                self.dense_matrix, self.height, self.width = load_matrix_internal(**kwargs)

    # Load the dictionary by providing the path to the file, filepath.
    def load_dictionary(self, dictfile='', dictionary=None):
        if dictfile != '':
            self.dictionary = numpy.genfromtxt(dictfile, dtype="str")
        else:
            self.dictionary = numpy.array(dictionary)


    def compute_flat_assignments(self, vector[int]& flat):
        ComputeAssignments[double](flat, &(self.h[0]), self.k, self.k, self.width)
        return flat


    def topterms(self):
        cdef vector[int] temp_indices
        temp_indices.resize(self.maxterms*self.k)
        TopTerms(self.maxterms, &(self.w[0]), self.height, self.height, self.k, temp_indices)
        self.term_indices = temp_indices


    def write_flatclust(self, string& assignfile, string& resultfile, vector[int]& assignments, vector[string]& dictionary, 
                            vector[int]& term_indices, unsigned int n, const FileFormat form):

        FlatClustWriteResults(assignfile, resultfile, assignments, dictionary, term_indices, 
                             form, self.maxterms, n, self.k)

    # Return the top term indices for each cluster. The length of the returned array is maxterms*k, 
    # with the first maxterms elements belonging to the first cluster, the second maxterms elements 
    # belonging to the second cluster, etc.  
    def get_flat_top_terms(self):
        return self.term_indices

    # Return the list of cluster assignments for each document.
    def get_assignments(self):
        return self.assignments_flat

#----------------------------------------------------------------------------------
# TODO: NMF
#----------------------------------------------------------------------------------

cdef class Flatclust(Clustering):

    def __init__(self):
        super(Flatclust, self).__init__()
        self.nmfstats = Stats()
        self.tree = TreeResults()

    # Creates a parser that takes the same inputs and has the same defaults as the binary smallk tool.
    # Returns a dictionary containing the command line inputs or their default values. The fields in 
    # that dictionary are: matrixfile, dictfile, clusters, algorithm, tol, infile_W, infile_H, outdir, 
    # maxterms, verbose, format, assignfile, treefile, maxiter, miniter, maxthreads.
    def parser(self):

        parser = argparse.ArgumentParser()
        parser.add_argument("--matrixfile", action="store", required=True,    metavar="matrixfile")
        parser.add_argument("--dictfile",   action="store", required=True,    metavar="dictfile")
        parser.add_argument("--clusters",   action="store", required=True,    metavar="clusters",   type=int)
        parser.add_argument("--algorithm",  action="store", required=False,   metavar="algorithm",    default="BPP", choices=["HALS", "RANK2", "BPP"])
        parser.add_argument("--infile_W",   action="store", required=False,   metavar="infile_W",     default="")
        parser.add_argument("--infile_H",   action="store", required=False,   metavar="infile_H",     default="")
        parser.add_argument("--tol",        action="store", required=False,   metavar="tol",        type=float,  default=0.0001)
        parser.add_argument("--outdir",     action="store", required=False,   metavar="outdir",       default="")
        parser.add_argument("--miniter",    action="store", required=False,   metavar="miniter",    type=int,  default=5)
        parser.add_argument("--maxiter",    action="store", required=False,   metavar="maxiter",    type=int,  default=5000)
        parser.add_argument("--maxterms",   action="store", required=False,   metavar="maxterms",     default=5)
        parser.add_argument("--maxthreads", action="store", required=False,   metavar="maxthreads",   default=8)
        parser.add_argument("--verbose",    action="store", required=False,   metavar="verbose",      default=True)
        parser.add_argument("--format",     action="store", required=False,   metavar="format",       default="XML")
        parser.add_argument("--assignfile", action="store", required=False,   metavar="assignfile",   default="assignments")
        parser.add_argument("--treefile", action="store", required=False,   metavar="treefile",   default="tree")
        args = parser.parse_args()
        return args

    # Run flat clustering on the loaded matrix, using the provided input options
    def cluster(self, k,
        infile_W='',
        infile_H='',
        algorithm="BPP",
        maxterms=5,
        verbose=True,
        min_iter=5,
        max_iter=5000,
        max_threads=8,
        tol=0.0001):
 
        self.k = k
        self.maxterms = maxterms

        rng = Rand()
        rng.seed_from_time()
        RNG_CENTER = 0.5
        RNG_RADIUS = 0.5
        num_initializers = 1
        if not infile_W:
            w_init = random(self.height, k, rng, RNG_CENTER, RNG_RADIUS)
        else:
            w, height, width = load_matrix_internal(height=self.height, width=k, filepath=infile_W, num_initializers=num_initializers)
            w_init = w[0]

        if not infile_H:
            h_init = random(k, self.width, rng, RNG_CENTER, RNG_RADIUS)
        else:
            h, height, width = load_matrix_internal(height=k, width=self.width, filepath=infile_H, num_initializers=num_initializers)
            h_init = h[0]
        cdef NmfOptions nmf_opts
        nmf_opts.tol = tol
        nmf_opts.algorithm = self.get_alg(algorithm)
        nmf_opts.prog_est_algorithm = PG_RATIO
        nmf_opts.height = self.height
        nmf_opts.width = self.width
        if algorithm == "RANK2":
            nmf_opts.k = 2
        else:
            nmf_opts.k = k
        nmf_opts.min_iter = min_iter
        nmf_opts.max_iter = max_iter
        nmf_opts.tolcount = 1
        nmf_opts.max_threads = max_threads
        nmf_opts.verbose = verbose
        nmf_opts.normalize = True

        if self.sparse:
            self.w, self.h = self.flatclust_sparse(nmf_opts, self.matrix, w_init, h_init, self.height, self.k, self.nmfstats)
        else:
            self.w, self.h = self.flatclust(nmf_opts, self.dense_matrix, w_init, h_init, self.height, self.height, self.k, self.nmfstats)
        self.assignments_flat = self.compute_flat_assignments(self.assignments_flat)
        self.topterms()


    def flatclust_sparse(self, nmf_opts, Sparse A, 
                      vector[double]& buf_w, vector[double]& buf_h,
                      unsigned int ldim_w, unsigned int ldim_h, 
                      Stats stats): 
        FlatClustSparse(nmf_opts, A.Height(), A.Width(), A.Size(),
                                      A.LockedColBuffer(),
                                      A.LockedRowBuffer(),
                                      A.LockedDataBuffer(),
                                      &(buf_w[0]), ldim_w,
                                      &(buf_h[0]), ldim_h,
                                      dereference(stats.get()))
    

        return buf_w, buf_h

    def flatclust(self, nmf_opts, 
                vector[double]& buf_a, vector[double]& buf_w, vector[double]& buf_h, 
                unsigned int ldim_a, unsigned int ldim_w, unsigned int ldim_h, 
                Stats stats):

    
        FlatClust(nmf_opts, &(buf_a[0]), ldim_a,
                                &(buf_w[0]), ldim_w,
                                &(buf_h[0]), ldim_h,
                                dereference(stats.get()))

        return buf_w, buf_h


    # Writes the assignment file and treefile in the provided format and to the provided output directory.
    def write_output(self, assignfile, treefile, outdir='./', format='XML'):
        print 'Writing output files...'
        if format == "XML":
            tree = outdir + treefile + '_' + str(self.k) + '.xml'
        else:
            tree = outdir + treefile + '_' + str(self.k) + '.json'
        assign = outdir + assignfile + '_' + str(self.k) + '.csv'
        self.write_flatclust(assign, tree, self.assignments_flat, self.dictionary, 
            self.term_indices, self.width, get_outputformat(format))

#----------------------------------------------------------------------------------


cdef class Hierclust(Clustering):

    cdef FileFormat format
    cdef unsigned int flat

    def __init__(self):
        super(Hierclust, self).__init__()
        self.stats = ClusterStats()
        self.tree = TreeResults()

    # Creates a parser that takes the same inputs and has the same defaults as the binary smallk tool.
    # Returns a dictionary containing the command line inputs or their default values. The fields in 
    # that dictionary are: matrixfile, dictfile, clusters, tol, infile_W, infile_H, outdir, maxterms, 
    # verbose, format, assignfile, treefile, maxiter, miniter, maxthreads, unbalanced, trial_allowance, flat.
    def parser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--matrixfile", action="store", required=True,  metavar="matrixfile")
        parser.add_argument("--dictfile",   action="store", required=True,  metavar="dictfile")
        parser.add_argument("--clusters",   action="store", required=True,  metavar="clusters",   type=int)
        parser.add_argument("--infile_W",   action="store", required=False, metavar="infile_W",     default="")
        parser.add_argument("--infile_H",   action="store", required=False, metavar="infile_H",     default="")
        parser.add_argument("--tol",        action="store", required=False, metavar="tol",        type=float,  default=0.0001)
        parser.add_argument("--outdir",     action="store", required=False, metavar="outdir",       default="")
        parser.add_argument("--miniter",    action="store", required=False, metavar="miniter",    type=int,  default=5)
        parser.add_argument("--maxiter",    action="store", required=False, metavar="maxiter",    type=int,  default=5000)
        parser.add_argument("--maxterms",   action="store", required=False, metavar="maxterms",   type=int,  default=5)
        parser.add_argument("--maxthreads", action="store", required=False, metavar="maxthreads", type=int,  default=8)
        parser.add_argument("--unbalanced", action="store", required=False, metavar="unbalanced", type=float,  default=0.1)
        parser.add_argument("--trial_allowance", action="store", required=False, metavar="trial_allowance", type=int, default=3)
        parser.add_argument("--flat",       action="store", required=False, metavar="flat",  type=int,       default=0)
        parser.add_argument("--verbose",    action="store", required=False, metavar="verbose",      default=True)
        parser.add_argument("--format",     action="store", required=False, metavar="format",       default="XML", choices=["XML", "JSON"])
        parser.add_argument("--treefile",  action="store", required=False, metavar="treefile",    default="tree")
        parser.add_argument("--assignfile", action="store", required=False, metavar="assignfile",   default="assignments")
        args = parser.parse_args()
        return args

    # Run hierarchical clustering on the loaded matrix, using the provided input options.
    def cluster(self, k, 
        infile_W='',
        infile_H='',
        maxterms=5,
        unbalanced=0.1,
        trial_allowance=3,
        verbose=True,
        flat=0,
        min_iter=5,
        max_iter=5000,
        max_threads=8,
        tol=0.0001):

        self.k = k
        self.flat = flat
        self.maxterms = maxterms

        rng = Rand()
        rng.seed_from_time()
        RNG_CENTER = 0.5
        RNG_RADIUS = 0.5
        required_size = self.height*2
        num_initializers = 2*k
        if not infile_W:
            w_initializers = [[0]*required_size for i in xrange(num_initializers)] 
            for i in xrange(num_initializers):
                w_initializers[i] = random(self.height, 2, rng, RNG_CENTER, RNG_RADIUS)
        else:
            w_initializers, width, height = load_matrix_internal(height=self.height, width=2, filepath=infile_W, num_initializers=num_initializers)
        required_size = 2*self.width
        if not infile_H:
            h_initializers = [[0]*required_size for i in xrange(num_initializers)] 
            for i in xrange(num_initializers):
                h_initializers[i] = random(2, self.width, rng, RNG_CENTER, RNG_RADIUS)
        else:
            h_initializers, width, height = load_matrix_internal(height=2, width=self.width, filepath=infile_H, num_initializers=num_initializers)

        #set opts
        self.clust_opts.maxterms                      = int(maxterms)
        self.clust_opts.unbalanced                    = float(unbalanced)
        self.clust_opts.trial_allowance               = int(trial_allowance)
        self.clust_opts.num_clusters                  = int(k)
        self.clust_opts.verbose                       = verbose
        self.clust_opts.flat                          = int(flat)
        self.clust_opts.nmf_opts.tol                  = float(tol)
        self.clust_opts.nmf_opts.algorithm            = RANK2
        self.clust_opts.nmf_opts.prog_est_algorithm   = PG_RATIO
        self.clust_opts.nmf_opts.height               = self.height
        self.clust_opts.nmf_opts.width                = self.width
        self.clust_opts.nmf_opts.k                    = int(k)
        self.clust_opts.nmf_opts.min_iter             = min_iter
        self.clust_opts.nmf_opts.max_iter             = max_iter
        self.clust_opts.nmf_opts.tolcount             = 1
        self.clust_opts.nmf_opts.max_threads          = max_threads
        self.clust_opts.nmf_opts.verbose              = verbose
        self.clust_opts.nmf_opts.normalize            = False


        if self.sparse:
            self.w, self.h, self.assignments = self.clust_sparse(self.matrix, w_initializers, h_initializers, self.tree, self.stats, self.height, self.width, k)
        else:
            self.w, self.h, self.assignments = self.clust_dense(self.dense_matrix, self.height, w_initializers, h_initializers, self.tree, self.stats, self.height, self.width, k)
        if self.flat == 1:
            self.assignments_flat = self.compute_flat_assignments(self.assignments_flat)
        self.topterms()

    # Return the top term indices for each cluster. The length of the returned array is maxterms*k, 
    # with the first maxterms elements belonging to the first cluster, the second maxterms elements 
    # belonging to the second cluster, etc.  
    def get_flat_top_terms(self):
        if self.flat == 1:
            return self.term_indices
        else:
            print 'ERROR: To get top terms indices, rerun hierarchical clusting with flat=1'
            return None

    def clust_dense(self, vector[double]& buf_a, unsigned int ldim_a,
            vector[vector[double]]& w_initializers, vector[vector[double]]& h_initializers,
            TreeResults tree, ClusterStats stats,
            unsigned int m, unsigned int n, unsigned int num_clusters):
        self.buf_w.resize(m*num_clusters)
        self.buf_h.resize(n*num_clusters)
        # cdef vector[int] assignments

        # cdef Result res = Clust(clust_opts, &(buf_a[0]), ldim_a, &(buf_w[0]), &(buf_h[0]),
        Clust(self.clust_opts, &(buf_a[0]), ldim_a, &(self.buf_w[0]), &(self.buf_h[0]),
                                w_initializers, h_initializers, self.assignments, 
                                dereference(tree.get()), dereference(stats.get()))
        return self.buf_w, self.buf_h, self.assignments

    def clust_sparse(self, Sparse A, 
                  vector[vector[double]]& w_initializers, vector[vector[double]]& h_initializers,
                  TreeResults tree, ClusterStats stats,
                  unsigned int m, unsigned int n, unsigned int num_clusters):
        # cdef vector[double] buf_w, buf_h
        self.buf_w.resize(m*num_clusters)
        self.buf_h.resize(n*num_clusters)
        
        ClustSparse(self.clust_opts, dereference(A.get()), &(self.buf_w[0]), &(self.buf_h[0]), 
                                      w_initializers, h_initializers, self.assignments, 
                                      dereference(tree.get()), dereference(stats.get()))

        return self.buf_w, self.buf_h, self.assignments

    # Writes the assignment file and treefile in the provided format and to the provided output directory.
    def write_output(self, assignfile, treefile, outdir='./', format='XML'):

        print 'Writing output files...'
        if format == 'XML':
            tree = outdir + treefile + '_' + str(self.k) + '.xml'
        elif format == 'JSON':
            tree = outdir + treefile + '_' + str(self.k) + '.json'
        assign = outdir + assignfile + '_' + str(self.k) + '.csv'
        if self.flat == 1:
            self.write_assignments(self.assignments_flat, assign)
            self.write_flatclust(assign, tree, self.assignments_flat, self.dictionary, 
                self.term_indices, self.width, get_outputformat(format))
        else:
            self.write_assignments(self.assignments, assign)
            self.tree.write(tree, get_outputformat(format), self.dictionary)


    def write_assignments(self, vector[int]& labels, const string& filepath):
        return WriteAssignmentsFile(labels, filepath)

    def get_assignments(self):
        if self.flat:
            return self.assignments_flat
        else:
            return self.assignments



#----------------------------------------------------------------------------------

cdef class Matrixgen:

    cdef Rand rng
    cdef bool is_sparse
    cdef Sparse S
    cdef vector[double] A
    cdef unsigned int height
    cdef unsigned int width



    def __init__(self):
        self.rng = Rand()
        self.rng.seed_from_time()
        self.is_sparse = False

    # Creates a parser that takes the same inputs and has the same defaults as the binary matrixgen tool.
    # Returns a dictionary containing the command line inputs or their default values. The fields in 
    # that dictionary are: height, width, filename, type, rng_center, rng_radius, precision, and nz_per_col.
    def parser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--height",         action="store", required=True,  metavar="height", type=int)
        parser.add_argument("--width",          action="store", required=True, metavar="width", type=int)
        parser.add_argument("--filename",       action="store", required=True, metavar="filename", type=str)
        parser.add_argument("--type",           action="store", required=False, metavar="type", default='UNIFORM',
                        choices=['UNIFORM', 'DENSE_DIAG', 'SPARSE_DIAG','IDENTITY', 'ONES', 'ZEROS', 'SPARSE'])
        parser.add_argument("--rng_center",     action="store", required=False, metavar="rng_center", default=0.5, type=float)
        parser.add_argument("--rng_radius",     action="store", required=False, metavar="rng_radius", default=0.5, type=float)
        parser.add_argument("--precision",      action="store", required=False, metavar="precision", default=6, type=int)
        parser.add_argument("--nz_per_col",     action="store", required=False, metavar="nz_per_col", default=1, type=int)

        args = parser.parse_args()
        return args

    # Writes the generated matrix to the filename with the specified precision.
    def write_output(self, filename, precision=6):
        if self.is_sparse:
            if not write_mtx(filename, self.S, precision):
                print 'Matrixgen error - sparse matrix file write failed.'
            else:
                print 'Matrix succesfully written.'
        else:
            if not write_delimited(self.A, filename, precision, self.height, self.height, self.width):
                print 'Matrixgen error - file write failed.'
            else: 
                print 'Matrix succesfully written.'

    # Generates a uniform matrix of height m and width n with the RNG attributes of center and radius.
    def uniform(self, unsigned int m, unsigned int n, double center=0.5, double radius=0.5):
        cdef unsigned int c = 0
        cdef vector[double] A
        cdef unsigned int r
        self.height = m
        self.width = n
        while c != n:
            r = 0
            while r != m:
                A.push_back(self.rng.double(center, radius))
                r += 1
            c += 1
        self.A = A
        
        return A
    
    # Generates a dense diagonal matrix of height m and width n with the RNG attributes of center and radius.
    def densediag(self, unsigned int m, unsigned int n, double center=0.5, double radius=0.5):
        cdef unsigned int c = 0
        cdef vector[double] A
        cdef unsigned int r
        self.height = m
        self.width = n
        while c != n:
            r = 0
            while r != m:
                A.push_back(0.0)
                r += 1
            c += 1

        c = 0
        while c != n:
            A[c*m + c] = self.rng.double(center, radius)
            c += 1
        self.A = A
        return A

    # Generates an identify matrix of height m and width n.
    def identity(self, unsigned int m, unsigned int n):
        cdef unsigned int c = 0
        cdef vector[double] A
        cdef unsigned int r
        self.height = m
        self.width = n
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
        self.A = A
        return A
    
    # Generates a sparse diagonal matrix of width n with the RNG attributes of center and radius.  
    def sparsediag(self, unsigned int n, double center=0.5, double radius=0.5):
        S = Sparse()
        self.is_sparse = True
        S.Reserve(n, n, n)
        S.BeginLoad()
        cdef unsigned int c = 0
        self.width = n
        while c != n:
            S.Load(c, c, self.rng.double(center, radius))
            c += 1
        S.EndLoad()
        self.S = S
        return S

    # Generates an ones matrix of height m and width n.
    def ones(self, unsigned int m, unsigned int n):
        cdef unsigned int c = 0
        cdef unsigned int r
        cdef vector[double] A
        self.height = m
        self.width = n
        while c != n:
            r = 0
            while r != m:
                A.push_back(1.0)
                r += 1
            c += 1
        self.A = A
        return A

    # Generates a zeros matrix of height m and width n.
    def zeros(self, unsigned int m, unsigned int n):
        cdef unsigned int c = 0
        cdef unsigned int r
        cdef vector[double] A
        self.height = m
        self.width = n
        while c != n:
            r = 0
            while r != m:
                A.push_back(0.0)
                r += 1
            c += 1
        self.A = A
        return A

    # Generates an identify matrix of height m and width n and nz non-zero elements.
    def sparse(self, unsigned int m, unsigned int n, unsigned int nz):
        self.is_sparse = True
        self.height = m
        self.width = n
        S = Sparse()
        RandomSparseMatrix(dereference(self.rng.get()), dereference(S.get()), nz, m, m, n, n)
        self.S = S
        return S

#----------------------------------------------------------------------------------

cdef class Preprocessor:

    cdef unsigned int maxiter, docsperterm, termsperdoc, precision, boolean_mode
    cdef unsigned int height, width
    cdef TermFreqMatrix tfm
    cdef vector[unsigned int] term_indices, doc_indices
    cdef unsigned int col_offsets_pointer
    cdef vector[unsigned int] reduced_term_indices
    cdef vector[double] reduced_scores
    cdef vector[double] scores
    cdef Sparse A
    cdef list dictionary, documents
    # cdef np.ndarray dictionary, documents
    cdef TFData tfdataArray
    cdef list row_indices, counts, col_offsets, term_ind, doc_ind

    def __init__(self):
        pass

    # Creates a parser that takes the same inputs and has the same defaults as the binary 
    # preprocessor tool.Returns a dictionary containing the command line inputs or their 
    # default values. The fields in that dictionary are: indir, outdir, docs_per_term, 
    # terms_per_doc, maxiter, precision, boolean_mode.
    def parser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--indir",         action="store", required=True,  metavar="indir")
        parser.add_argument("--outdir",        action="store", required=False, metavar="outdir",        default="./")
        parser.add_argument("--docs_per_term", action="store", required=False, metavar="docs_per_term", default=3)
        parser.add_argument("--terms_per_doc", action="store", required=False, metavar="terms_per_doc", default=5)
        parser.add_argument("--maxiter",       action="store", required=False, metavar="maxiter",       default=1000)
        parser.add_argument("--precision",     action="store", required=False, metavar="precision",     default=4)
        parser.add_argument("--boolean_mode",  action="store", required=False, metavar="boolean_mode",     default=0)
        args = parser.parse_args()
        return args

    # Writes the results of the preprocessing to the filesystem. The matrix is written 
    # to matrix_filepath, the dictionary to dict_filepath, and the documents to docs_filepath. 
    # The outputs are written with the specified precision.    
    def write_output(self, matrix_filepath, dict_filepath, docs_filepath, precision=4):
        write_termfreq(self.tfm, matrix_filepath, self.scores, precision)
        print 'about to write dictionary'
        print len(self.dictionary)
        print len(self.term_ind)
        print self.height
        print max(self.term_ind)
        write_strings(dict_filepath, self.dictionary, self.term_ind, self.height)
        print 'about to write documents'
        write_strings(docs_filepath, self.documents, self.doc_ind, self.width)


    # Preprocesses the matrix, dictionary, and documents provided.
    def preprocess(self,
        maxiter=1000,
        docsperterm=3,
        termsperdoc=5,
        boolean_mode=0):
        self.tfm = TermFreqMatrix(self.A, boolean_mode)
        cdef vector[unsigned int] term_indices
        term_indices.resize(self.height)
        cdef vector[unsigned int] doc_indices
        doc_indices.resize(self.width)
        cdef vector[double] scores
        
        cdef bool res = preprocess_tf(dereference(self.tfm.get()), term_indices, doc_indices, scores,
                                 maxiter, docsperterm, termsperdoc)
        if not res:
            print 'ERROR: preprocess()'
            return None
        else:
            self.width = self.tfm.Width()
            self.height = self.tfm.Height()
            col_offsets_pointer = self.tfm.get().LockedColBuffer()
            tfdata_pointer = self.tfm.get().LockedTFDataBuffer()
            self.row_indices = []
            # self.counts = []
            self.col_offsets = []
            self.term_ind = []
            self.doc_ind = []

        for i in range(len(scores)):
            # self.counts.append(tfdata_pointer[i].count) #would be used for count matrix
            self.row_indices.append(tfdata_pointer[i].row)

        for i in range(self.width + 1):
            self.col_offsets.append(col_offsets_pointer[i])

        for i in range(self.width):
            self.doc_ind.append(doc_indices[i])

        for i in range(self.height):
            self.term_ind.append(term_indices[i])

        self.scores = scores


    # Load the input matrix. Two combinations are available, which provide for loading via a file path
    # or direct loading of a dense matrix. For example:

    # File path: load_inputmatrix(filepath=’/path/to/file/’)
    # Dense matrix: load_inputmatrix(height=matrix_height, width=matrix_width, nz=non_zero_count, 
    #     buffer=non_zero_elements, row_indices=rows, col_offsets=cols)

    def load_inputmatrix(self, filepath=None, height=None, width=None, nz=None, buffer=None, row_indices=None, col_offsets=None):
        if filepath != None:
            self.A, self.height, self.width = load_matrix_internal(filepath=filepath)
        elif row_indices != None:
            self.A, self.height, self.width = load_matrix_internal(height=height, width=width, nz=nz, buffer=buffer, row_indices=row_indices, col_offsets=col_offsets)


    # Load the dictionary, either by providing the path to the file, filepath, 
    # or by providing a list of the terms, dictionary.
    def load_dictionary(self, filepath=None, dictionary=None):
        if filepath != None:
            with open(filepath) as f:
                self.dictionary = f.read().split('\n')
        else:
            self.dictionary = dictionary

    # Load the documents, either by providing the path to the file, filepath, 
    # or by providing a list of the documents, documents.
    def load_documents(self, filepath=None, documents=None):
        if filepath != None:
            with open(filepath) as f:
                self.documents = f.read().split('\n')
        else:
            self.documents = documents

    # Returns the pruned list of document ids.
    def get_reduced_documents(self):

        return [self.documents[i] for i in self.doc_ind]

    # Returns the pruned dictionary.
    def get_reduced_dictionary(self):
        return [self.dictionary[i] for i in self.term_ind]

    # Returns the list of non-zero values in the pruned sparse matrix.
    def get_reduced_scores(self):
        return self.scores

    # Returns the list of row_indices for the pruned sparse matrix.
    def get_reduced_row_indices(self):
        return self.row_indices

    # Returns the list of row_indices for the pruned sparse matrix.
    def get_reduced_col_offsets(self):
        return self.col_offsets



#----------------------------------------------------------------------------------

#instantiate the objects
smallkapi = SmallkAPI()
flatclust = Flatclust()
matrixgen = Matrixgen()
preprocessor = Preprocessor()
hierclust = Hierclust()



