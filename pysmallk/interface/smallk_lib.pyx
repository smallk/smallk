#####################################################################################
# -*- coding: utf-8 -*-
# Copyright 2013,2014,2015,2016 Georgia Institute of Technology.
#
# Licensed under the Apache License, Version 2.0 (the â€œLicenseâ€); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by
# applicable law or agreed to in writing, software distributed under the License
# is distributed on as â€œAS ISâ€ BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the specific language
# governing permissions and limitations under the License.
#
#####################################################################################

# imports Python functions
import numpy as np

# imports the Numpy C API
cimport numpy as np

import cython
from cython.operator import dereference

from libc.stdlib cimport malloc, free
from libc.string cimport memcpy, memset
from libc.string cimport strdup, strcpy

from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector

import argparse
import sys


#####################################################################################
#
#                                   C++ functions
#
#####################################################################################

cdef extern from "smallk.hpp" namespace "smallk":
    cdef enum Algorithm: 
        MU
        HALS
        RANK2
        BPP
    cdef enum OutputFormat:
        XML
        JSON
    ctypedef char *stdStringR "const std::string &"
    void Initialize(int& argc, char **&argv) except +
    int IsInitialized()
    void Finalize()
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
    void LoadMatrix(stdStringR filepath) except +
    void Nmf(int k, Algorithm algorithm, stdStringR initfile_w, stdStringR initfile_h) except +
    void SetOutputPrecision(const unsigned int num_digits) except +
    void SetMinIter(const unsigned int min_iterations) except +
    void SetNmfTolerance(const double tol) except +
    void SetMaxIter(const unsigned int max_iterations) except +
    void SetMaxThreads(const unsigned int mt) except +
    void SetHierNmf2Tolerance(const double tol) except +
    void SetMaxTerms(const unsigned int max_terms) except +
    void SetOutputDir(const string& output_dir) except +
    bool IsMatrixLoaded()
    void LoadMatrix(double *buffer, unsigned int ldim, unsigned int height, unsigned int width) except +
    void LoadMatrix(unsigned int height, unsigned int width, unsigned int nz, vector[double]& data,
        vector[unsigned int]& row_indices, vector[unsigned int]& col_offsets) except +
    const double* LockedBufferW(unsigned int& ldim, unsigned int& height, unsigned int& width) except +
    const double* LockedBufferH(unsigned int& ldim, unsigned int& height, unsigned int& width) except +
    void HierNmf2(const unsigned int num_clusters) except +
    void SetOutputFormat(const OutputFormat form) except +
    void LoadDictionary(const string& filepath) except +
    void LoadDictionary(const vector[string]& terms) except +

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
    void NmfInitialize(int argc, char* argv[]) except +
    Result NmfIsInitialized()
    void NmfFinalize()

cdef extern from "flat_clust.hpp":
    cdef struct FlatClustOptions:
        NmfOptions nmf_opts
        int maxterms
        int num_clusters
        bool verbose
    Result FlatClust(const NmfOptions& options, double* buf_A, int ldim_A, double* buf_W, int ldim_W,
        double* buf_H, int ldim_H, NmfStats& stats) except +
    Result FlatClustSparse(const NmfOptions& options, const unsigned int height, const unsigned int width,
        const unsigned int nz, const unsigned int* col_offsets, const unsigned int* row_indices,
        const double* data, double* buf_W, int ldim_W, double* buf_H, int ldim_H, NmfStats& stats) except +
    void FlatClustSparseTest(const NmfOptions& options, const unsigned int height, const unsigned int width,
        const unsigned int nz, const unsigned int* col_offsets, const unsigned int* row_indices,
        const double* data, double* buf_W, int ldim_W, double* buf_H, int ldim_H, NmfStats& stats) except +

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
        void Reserve(const unsigned int height, const unsigned int width, const unsigned int nzmax) except +
        void BeginLoad()
        void EndLoad()
        void Load(const unsigned int row, const unsigned int col, const T& value) except +

cdef extern from "file_loader.hpp":
    bool IsDense(const string& file_path) except +
    bool IsSparse(const string& file_path) except +
    bool LoadSparseMatrix(const string& file_path, SparseMatrix[double]& A, unsigned int& height,
        unsigned int& width, unsigned int& nz) except +
    bool LoadDenseMatrix(const string& file_path, vector[double]& data, unsigned int& height,
        unsigned int& nz) except +

cdef extern from "random.hpp":
    cdef cppclass Random:
        Random()
        void SeedFromTime()
        double RandomDouble(const double& center, const double& radius) except +

cdef extern from "matrix_generator.hpp":
    bool RandomMatrix(vector[double]& buf, const unsigned int height, const unsigned int width,
        Random& rng, const double rng_center, const double rng_radius) except +

cdef extern from "delimited_file.hpp":
    bool LoadDelimitedFile(vector[double]& buf, unsigned int& height, unsigned int& width, 
        const string& filename, const char DELIM) except +
    bool LoadMatrixArray[T](vector[vector[T]]& buf, const unsigned int height,
        const unsigned int width, const string& filename, const char DELIM) except +

cdef extern from "file_format.hpp":
    cdef enum FileFormat:
        CSV
        XML
        JSON
        TXT

cdef extern from "assignments.hpp":
    bool WriteAssignmentsFile(const vector[unsigned int]& labels, const string& filepath) except +
    void ComputeAssignments[T](vector[unsigned int]& assignments, const T* buf_h,
        const unsigned int ldim_h, const unsigned int k, const unsigned int n) except +
    void ComputeFuzzyAssignments[T](vector[float]& probabilities, const T* buf_h,  
        const unsigned int ldim_h, const unsigned int k, const unsigned int n) except +

cdef extern from "terms.hpp":
    void TopTerms(const int maxterms, const double* buf_w, const unsigned int ldim,
        const unsigned int height, const unsigned int width, vector[int]& term_indices) except +

cdef extern from "flat_clust_output.hpp":
    void FlatClustWriteResults(const string& assignfilepath, const string& fuzzyfilepath,
        const string& resultfilepath, const vector[unsigned int]& assignments,
        const vector[float]& probabilities, const vector[string]& dictionary,
        const vector[int]& term_indices, const FileFormat form, const unsigned int maxterms,
        const unsigned int num_docs, const unsigned int num_clusters) except +

cdef extern from "clust.hpp":
    cdef struct ClustStats:
        ClustStats()
        int nmf_count
        int max_count
    cdef cppclass ClustOptions:
        NmfOptions nmf_opts
        int maxterms
        double unbalanced
        int trial_allowance
        int num_clusters
        bool verbose
        bool flat
        string initdir
    Result Clust(const ClustOptions& options, double* buf_A, int ldim_A, double* buf_w, double* buf_h,
        Tree& tree, ClustStats& stats, Random& rng) except +
    Result ClustSparse(const ClustOptions& options, const SparseMatrix[double]& A,
        double* buf_w, double* buf_h, Tree& tree, ClustStats& stats, Random& rng) except +
    Result ClustTest(const ClustOptions& options, const SparseMatrix[double]& A, Tree& tree, ClustStats& stats, 
        Random& rng, double* buf_w, double* buf_h,) except +

cdef extern from "hierclust_writer.hpp":
    cdef cppclass IHierclustWriter:
        IHierclustWriter()

cdef extern from "hierclust_writer_factory.hpp":
    IHierclustWriter* CreateHierclustWriter(const FileFormat& format) except +

cdef extern from "tree.hpp":
    cdef cppclass Tree[double]:
        Tree()
        bool WriteTree(IHierclustWriter* writer, const string& filepath, const vector[string]& dictionary) except +
        void Print()
        bool WriteAssignments(const string& filepath) except +
        vector[unsigned int]& Assignments()

cdef extern from "sparse_matrix_ops.hpp":
    cdef cppclass RandomSparseMatrix:
        void RandomSparseMatrix(Random& rng, SparseMatrix[double]& A, 
            const unsigned int nonzeros_per_column, const unsigned int min_height, 
            const unsigned int max_height, const unsigned int min_width, const unsigned int max_width)

cdef extern from "sparse_matrix_io.hpp":
    cdef bool WriteMatrixMarketFile(const string& file_path, const SparseMatrix[double]& S,
        const unsigned int precision) except +
    cdef bool LoadMatrixMarketFile(const string& file_path, SparseMatrix[double]& A, 
        unsigned int& height, unsigned int& width, unsigned int& nnz) except +
    void LoadSparseMatrix(const unsigned int height, const unsigned int width, const unsigned int nz, 
        const vector[double]& data, const vector[unsigned int]& row_indices,
        const vector[unsigned int]& col_offsets, SparseMatrix[double]& A) except +

cdef extern from "preprocess.hpp":
    cdef unsigned int TT = 0xFFFFFFFF
    cdef unsigned int FF = 0
    bool preprocess_tf(TermFrequencyMatrix& A, vector[unsigned int]& term_indices,
        vector[unsigned int]& doc_indices, vector[double]& scores, const unsigned int MAX_ITER,
        const unsigned int DOCS_PER_TERM, const unsigned int TERM_PER_DOC) except +

cdef extern from "delimited_file.hpp":
    cdef bool WriteDelimitedFile(const double* buffer, const unsigned int ldim, const unsigned int height, 
        const unsigned int width, const string& filename, const unsigned int precision, const char DELIM) except +

cdef extern from "term_frequency_matrix.hpp":
    cdef cppclass TermFrequencyMatrix:
        TermFrequencyMatrix(const SparseMatrix[double]& S, bool boolean_mode) except +
        unsigned int Width()
        unsigned int Height()
        unsigned int Size()
        bool WriteMtxFile(const string& file_path, const double* scores, const unsigned int precision) except +
        const unsigned int* LockedColBuffer()
        const TFData* LockedTFDataBuffer()
    cdef struct TFData:
        TFData()
        unsigned int row
        unsigned int count

#####################################################################################
#
#                              Cython data structures
#
#####################################################################################

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
    
    cdef Tree[double] *thisptr
    
    cdef IHierclustWriter *writer
    
    def __cinit__(self):
        self.thisptr = new Tree[double]()
    
    def write(self, const string& filepath, const vector[string]& dictionary, const FileFormat& form):
        return self.thisptr.WriteTree(CreateHierclustWriter(form), filepath, dictionary)

    def write_assignments(self, const string& filepath):
        return self.thisptr.WriteAssignments(filepath)

    def get_assignments(self, vector[unsigned int]& assignments):
        assignments = self.thisptr.Assignments()
        return assignments

    cdef Tree[double]* get(self):
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

#####################################################################################
#
#                              Python common functions
#
#####################################################################################

# Performance enhancement:
# Cython is free to assume that indexing operations ([]-operator) in the code 
# will not cause any IndexErrors to be raised. 
@cython.boundscheck(False)

# Performance enhancement:
# Cython will not ensure that python indexing is not used.
@cython.wraparound(False)

#-------------------------------- private functions ---------------------------------#

def _random(const unsigned int height, const unsigned int width, Rand r, const double r_center, const double r_radius):
    cdef vector[double] buf
    cdef bool ok = RandomMatrix(buf, height, width, dereference(r.get()), r_center, r_radius)
    return buf

def _is_dense(string infile):
    return IsDense(infile)

def _is_sparse(string infile):
    return IsSparse(infile)

def _get_alg(alg_name):
    if (alg_name.lower() == 'mu'):
        return 0
    elif (alg_name.lower() == 'hals'):
        return 2
    elif (alg_name.lower() == 'rank2'):
        return 3
    elif (alg_name.lower() == 'bpp'):
        return 1

def _get_outputformat(form):
    if form.lower() == "xml":
        return 1
    elif form.lower() == "json":
        return 2

def _write_mtx(string& filepath, Sparse S, unsigned int precision):
    res = WriteMatrixMarketFile(filepath, dereference(S.get()), precision)
    return res

def _write_delimited(vector[double] A, string& filename, unsigned int precision, unsigned int ldim, 
    unsigned int height, unsigned int width):
    res = WriteDelimitedFile(&A[0], ldim, height, width, filename, precision, ',')
    return res

def _write_termfreq(TermFreqMatrix M, const string& file_path, scores, const unsigned int precision):
    cdef vector[double] cscores
    for i in range(0, len(scores)):
        cscores.push_back(scores[i])
    cdef const double *pScore0 = &cscores[0]
    return M.get().WriteMtxFile(file_path, pScore0, precision)

def _write_strings(filepath, stringarr, valid_indices, N):
    f = open(filepath, "w")
    for i in range(0, N):
        index = valid_indices[i]
        f.write(stringarr[index]+"\n")
    f.close()

def _load_matrix_internal(filepath="", height=0, width=0, delim="", buffer=[], nz=0, row_indices=[], 
    col_offsets=[], matrix=[], column_major=False, sparse_matrix=None):
    if sparse_matrix != None:
        m, h, w = _load_sparse_obj(sparse_matrix, height, width)
    if filepath != "":
        sparse = _is_sparse(filepath)
        dense = _is_dense(filepath)

    if filepath != "" and sparse:
        m, h, w = _load_sparse_file(filepath, Sparse())
    elif height != 0 and width != 0 and filepath != "":
        if delim != "":
            m, h, w = _load_matrixarray(height, width, filepath, DELIM=delim)
        else:
            m, h, w = _load_matrixarray(height, width, filepath)
    elif len(row_indices) > 0 and len(col_offsets) > 0:
        m, h, w = _load_sparse(height, width, nz, buffer, row_indices, col_offsets, Sparse())
    elif len(buffer) > 0 and height != 0 and width != 0:
        m, h, w = buffer, height, width
    elif len(matrix) > 0:
        m, h, w = _load_numpy(matrix, column_major)
    elif height != 0 and width != 0 and filepath != "":
        if delim != "":
            m, h, w = _load_delimited(height, width, filepath, DELIM=delim)
        else:
            m, h, w = _load_delimited(height, width, filepath)
    elif filepath != "" and dense:
        m, h, w = _load_dense_file(filepath)
    return m, h, w

def _load_delimited(unsigned int height, unsigned int width, const string& filename, const char DELIM = ','):
    cdef vector[double] buf
    cdef bool ok = LoadDelimitedFile(buf, height, width, filename, DELIM)
    return buf, height, width

def _load_sparse_file(const string& file_path, Sparse A):
    cdef unsigned int m = 0, n = 0, nnz = 0
    cdef bool ok = LoadSparseMatrix(file_path, dereference(A.get()), m, n, nnz)
    return A, m, n

def _load_sparse_obj(Sparse A, unsigned int height, unsigned int width):
    return A, height, width

def _load_sparse(const unsigned int height, const unsigned int width, const unsigned int nz, 
    const vector[double]& data, const vector[unsigned int]& row_indices,
    const vector[unsigned int]& col_offsets, Sparse A):
    LoadSparseMatrix(height, width, nz, data, row_indices, col_offsets, dereference(A.get()))
    return A, height, width

def _load_dense_file(const string& file_path):
    cdef unsigned int m = 0, n = 0
    cdef vector[double] buf_a
    cdef bool ok = LoadDenseMatrix(file_path, buf_a, m, n)
    return buf_a, m, n

def _load_matrixarray(const unsigned int matrix_height, const unsigned int matrix_width, const string& filename,
    const char DELIM = ','):
    cdef vector[vector[double]] buf
    buf.resize(1)
    cdef bool ok = LoadMatrixArray[double](buf, matrix_height, matrix_width, filename, DELIM)
    return buf, matrix_height, matrix_width

# matrix is typed as numpy array of doubles
def _validate_load(np.ndarray[double, ndim=2, mode="c"] matrix): 
    return matrix.reshape(1, matrix.shape[0]*matrix.shape[1])[0], matrix.shape[0], matrix.shape[1]

def _load_numpy(matrix, column_major):
    if column_major:
        matrix_column_major = matrix
    else:
        matrix_column_major = matrix.T
    while True:
        try:
            D, height, width = _validate_load(matrix_column_major)
            return D, height, width
        except ValueError:
            if not matrix_column_major.flags['C_CONTIGUOUS']:
                print 'Array not C-contiguous. Rebuilding with contiguous memory allocation...'
                matrix_column_major = np.ascontiguousarray(matrix_column_major)
            else:
                print 'Array not float or double type. Rebuilding with correct type...'
                matrix_column_major = matrix_column_major.astype(np.float)

def _init():
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

def _isinit():
    return NmfIsInitialized() == INITIALIZED

def _finalize():
    NmfFinalize()

#####################################################################################
#
#                              Python SmallK functions
#
#####################################################################################

#-------------------------------- public functions ---------------------------------#


cdef class SmallkAPI:

    cdef bool _dictionary_loaded

    def __init__(self):
        self._init()
        if not _isinit():
            print 'ERROR'
        self._dictionary_loaded = False
 

    # \brief Returns the parsed arguments for the default command line application
    # \return The dictionary containing the parsed arguments
    def parser(self):
        parser = argparse.ArgumentParser(description="Run NMF via python binding")
        parser.add_argument('--matrixfile', action='store', 
            required=True, metavar='matrixfile')
        parser.add_argument('--k', action='store', 
            required=True, type=int, metavar='k')
        parser.add_argument('--dictfile', action='store', 
            required=False, metavar='dictfile', default="")
        parser.add_argument('--hiernmf2', action='store', 
            required=False, metavar='hiernmf2', default=0, choices=[0,1])
        parser.add_argument('--algorithm', action='store', 
            required=False, default='BPP', metavar='algorithm', choices=['MU','HALS','RANK2','BPP'])
        parser.add_argument('--stopping', action='store', 
            required=False, metavar='stopping', default='PG_RATIO', choices=['PG_RATIO','DELTA'])
        parser.add_argument('--tol', action='store', type=float, 
            required=False, metavar='tol', default=0.005)
        parser.add_argument('--tolcount', action='store', type=int, 
            required=False, metavar='tolcount', default=1)
        parser.add_argument('--infile_W', action='store', 
            required=False, metavar='infile_W', default="")
        parser.add_argument('--infile_H', action='store', 
            required=False, metavar='infile_H', default="")
        parser.add_argument('--outfile_W', action='store', 
            required=False, metavar='outfile_W', default="w.csv")
        parser.add_argument('--outfile_H', action='store', 
            required=False, metavar='outfile_H', default="h.csv")
        parser.add_argument('--outprecision', action='store', type=int, 
            required=False, metavar="outprecision", default=6)
        parser.add_argument('--maxiter', action='store', type=int, 
            required=False, metavar="maxiter", default=5000)
        parser.add_argument('--miniter', action='store', type=int, 
            required=False, metavar="miniter", default=5)
        parser.add_argument('--maxthreads', action='store', type=int, 
            required=False, metavar="maxthreads", default=8)
        parser.add_argument('--maxterms', action='store', type=int, 
            required=False, metavar="maxterms", default=5)
        parser.add_argument('--normalize', action='store', type=int, 
            required=False, metavar="normalize", default=1)
        parser.add_argument('--verbose', action='store', type=int, 
            required=False, metavar="verbose", default=1)
        args = parser.parse_args()
        return args


    # \brief Returns the major version of SmallK
    # \return unsigned int representing major version
    def get_major_version(self):
        return GetMajorVersion()

    # \brief Returns the minor version of SmallK
    # \return unsigned int representing minor version
    def get_minor_version(self):
        return GetMinorVersion()

    # \brief Returns the patch level of SmallK
    # \return unsigned int representing patch level
    def get_patch_level(self):
        return GetPatchLevel()

    # \brief Returns a string representation of the version of SmallK
    # \return string representing the major and minor version and patch level
    def get_version_string(self):
        return GetVersionString()

    # \brief Load an input matrix
    #
    # To load a matrix from a file:
    #   \param[kwarg] filepath      The path to the input matrix
    #
    # To load a sparse matrix from python:
    #   \param[kwarg] height        The height of the sparse matrix
    #   \param[kwarg] width         The width of the sparse matrix
    #   \param[kwarg] nz            The number of non-zeros in the sparse matrix
    #   \param[kwarg] buffer        List of doubles containing the non-zero elements of the sparse matrix
    #   \param[kwarg] row_indices   List of integers representing the row indices of the sparse matrix
    #   \param[kwarg] col_offsets   List of integers representing the column offsets of the sparse matrix
    #
    # To load a dense matrix from python:
    #   \param[kwarg] height        The height of the dense matrix
    #   \param[kwarg] width         The width of the dense matrix
    #   \param[kwarg] buffer        List of doubles containing the elements of the dense matrix
    #
    # To load a numpy matrix from python:
    #   \param[kwarg] matrix        The numpy matrix
    #   \param[kwarg] column_major  Boolean for whether or not the matrix is column major (optional)
    def load_matrix(self, filepath="", height=0, width=0, delim="", buffer=[], matrix=[], 
        nz=0, row_indices=[], col_offsets=[], column_major=False, sparse_matrix=None):
        if len(row_indices) > 0 and len(col_offsets) > 0:
            self._load_sparse_buffer(height, width, nz, buffer, row_indices, col_offsets)
        # elif sparse_matrix != None:
        #     self.load_sparse(sparse_matrix, height, width)
        elif len(buffer) > 0 and height != 0 and width != 0:
            self._load_dense_buffer(buffer, height, width)
        elif len(matrix) > 0:
            self._load_numpy(matrix, column_major)
        if filepath != "":
            self._load_matrix_file(filepath)

    # \brief Runs NMF on the loaded matrix using the supplied algorithm and implementation details
    # \param[in]    k           The desired number of clusters
    # \param[in]    algorithm   The desired NMF algorithm
    # \param[kwarg] infile_W    Initialization for W (optional)
    # \param[kwarg] infile_H    Initialization for H (optional)
    # \param[kwarg] precision   Precision for calcuations (optional)
    # \param[kwarg] min_iter    Minimum number of iterations (optional)
    # \param[kwarg] max_iter    Maximum number of iterations (optional)
    # \param[kwarg] tol         Tolerance for determing convergence (optional)
    # \param[kwarg] max_threads Maximum number of threads to use (optional)
    # \param[kwarg] outdir      Output directory for files (optional)
    def nmf(self, unsigned int k, algorithm, char* infile_W="", char* infile_H="", precision=4, min_iter=5, 
        max_iter=5000, tol=0.005, max_threads=8, outdir="."):
        if self.is_matrix_loaded():
            try:
                SetOutputPrecision(precision)
                SetMinIter(min_iter)
                SetMaxIter(max_iter)
                SetNmfTolerance(tol)
                SetMaxThreads(max_threads)
                SetOutputDir(outdir)
                Nmf(k, _get_alg(algorithm), infile_W, infile_H)
            except:
                raise
        else:
            print 'Error: No matrix loaded, not running NMF.'

    # \brief Returns a dictionary of the supplied inputs to the nmf function
    # \return The dictionary containing the supplied inputs
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

    # \brief Indicates whether or not a matrix has been loaded
    # \return boolean representing if a matrix is loaded
    def is_matrix_loaded(self):
        return IsMatrixLoaded()

    # \brief Cleans up the elemental and smallk environment
    def finalize(self): 
        Finalize()

    # \brief Returns the output H matrix
    # \return Numpy array of the output H matrix
    def get_H(self):
        cdef unsigned int ldim = 0
        cdef unsigned int height = 0
        cdef unsigned int width = 0
        h_pointer = LockedBufferH(ldim, height, width)
        h_all = []
        for i in range(height*width):
            h_all.append(h_pointer[i])
        h_all = np.array(h_all).reshape(width, height)
        return h_all.T

    # \brief Returns the output W matrix
    # \return Numpy array of the output W matrix
    def get_W(self):
        cdef unsigned int ldim = 0
        cdef unsigned int height = 0
        cdef unsigned int width = 0
        w_pointer = LockedBufferW(ldim, height, width)
        w_all = []
        for i in range(height*width):
            w_all.append(w_pointer[i])
        w_all = np.array(w_all).reshape(width, height)
        return w_all.T

    # \brief Runs HierNMF2 on the loaded matrix
    # \param[in]    k           The desired number of clusters
    # \param[kwarg] format      Output format, XML or JSON (optional)
    # \param[kwarg] maxterms    Maximum number of terms (optional)
    # \param[kwarg] tol         Tolerance to use for determining convergence (optional)
    def hiernmf2(self, k, format="XML", maxterms=5, tol=0.0001):
        if self._dictionary_loaded and self.is_matrix_loaded():
            SetHierNmf2Tolerance(tol)
            SetMaxTerms(maxterms)
            SetOutputFormat(self._get_outputformat(format))
            HierNmf2(k)
        else:
            print 'Error: No dictionary loaded.'

    # \brief Loads a dictionary to use for computing top terms
    # \param[kwarg] filepath    The filepath for the dictionary
    # OR
    # \param[kwarg] dictionary  List containing the dictionary strings
    def load_dictionary(self, filepath="", dictionary=[]):
        if filepath != "":
            self._load_dictionary_from_file(filepath)
        elif len(dictionary) > 0:
            self._load_dictionary_from_buffer(dictionary)
        else:
            print 'Error: Invalid dictionary arguments.'

#-------------------------------- private functions ---------------------------------#

    def _load_matrix_file(self, char* filepath):
        LoadMatrix(filepath)

    def _load_dense_buffer(self, vector[double]& buffer, unsigned int height, unsigned int width):
        LoadMatrix(&buffer[0], height, height, width)
      
    def _load_sparse_buffer(self, unsigned int height, unsigned int width, unsigned int nz, vector[double]& buffer,
        vector[unsigned int]& row_indices, vector[unsigned int]& col_offsets):
        LoadMatrix(height, width, nz, buffer, row_indices, col_offsets)

    def _get_outputformat(self, format):
        if format.lower() == "xml":
            return 0
        elif format.lower() == "json":
            return 1

    def _validate_load(self, np.ndarray[double, ndim=2, mode="c"] matrix):
        LoadMatrix(&matrix[0,0], matrix.shape[0], matrix.shape[0], matrix.shape[1])

    def _load_numpy(self, matrix, column_major):
        if column_major:
            matrix_column_major = matrix
        else:
            matrix_column_major = matrix.T

        while True:
            try:
                self._validate_load(matrix_column_major)
                break
            except ValueError:
                if not matrix_column_major.flags['C_CONTIGUOUS']:
                    print 'Array not C-contiguous. Rebuilding with contiguous memory allocation...'
                    matrix_column_major = np.ascontiguousarray(matrix_column_major)
                else:
                    print 'Array not float or double type. Rebuilding with correct type...'
                    matrix_column_major = matrix_column_major.astype(np.float)
    def _init(self):
        av = sys.argv
        ac = len(av)
        cdef char **c_arr = <char**>malloc((ac+1) * sizeof(char*))
        cdef char* c_string
        cdef char* char_str
        memset(c_arr, '\0', (ac+1)*sizeof(char*))
        for i in xrange(0, ac):
            py_string = av[i]
            py_byte_string = py_string.encode('ascii')
            char_str = py_byte_string
            data_len = len(char_str) + 1
            c_arr[i] = <char *>malloc(data_len)
            memset(c_arr[i], '\0', data_len*sizeof(char))
            memcpy(c_arr[i], char_str, data_len*sizeof(char))
        Initialize(ac, c_arr)
        for i in xrange(ac):
            free(c_arr[i])
        free(c_arr)

    def _isinit(self):
        return IsInitialized();

    def _load_dictionary_from_file(self, const string filepath):
        LoadDictionary(filepath)
        self._dictionary_loaded = True

    def _load_dictionary_from_buffer(self, vector[string]& dictionary):
        LoadDictionary(dictionary)
        self._dictionary_loaded = True


#####################################################################################
#
#                         Python common clustering functions
#
#####################################################################################

#-------------------------------- public functions ---------------------------------#

cdef class Clustering:

    cdef dict opts
    cdef Sparse matrix
    cdef vector[double] dense_matrix
    cdef TreeResults tree
    cdef ClusterStats stats
    cdef Stats nmfstats
    cdef vector[double] buf_w, buf_h
    cdef ClustOptions clust_opts
    cdef unsigned int height, width, k, required_size
    cdef Rand rng
    cdef vector[vector[double]] w_initializers, h_initializers
    cdef vector[double] w_init, h_init
    cdef vector[int] term_indices, assignments
    cdef vector[unsigned int] assignments_flat
    cdef vector[float] probabilities
    cdef vector[double] w, h
    cdef np.ndarray dictionary
    cdef unsigned int maxterms
    cdef bool _sparse
    cdef string initdir

    def __init__(self):
        _init()
        if not _isinit():
            print 'ERROR'

    # \brief Cleans up the elemental and smallk environment
    def finalize(self):
        _finalize()

    # \brief Load an input matrix
    #
    # To load a matrix from a file:
    #   \param[kwarg] filepath      The path to the input matrix
    #
    # To load a sparse matrix from python:
    #   \param[kwarg] height        The height of the sparse matrix
    #   \param[kwarg] width         The width of the sparse matrix
    #   \param[kwarg] nz            The number of non-zeros in the sparse matrix
    #   \param[kwarg] buffer        List of doubles containing the non-zero elements of the sparse matrix
    #   \param[kwarg] row_indices   List of integers representing the row indices of the sparse matrix
    #   \param[kwarg] col_offsets   List of integers representing the column offsets of the sparse matrix
    #
    # To load a sparse matrix from Matrixgen:
    #   \param[kwarg] height        The height of the sparse matrix
    #   \param[kwarg] width         The width of the sparse matrix
    #   \param[kwarg] sparse_matrix The sparse matrix returned from Matrixgen
    #
    # To load a dense matrix from python:
    #   \param[kwarg] height        The height of the dense matrix
    #   \param[kwarg] width         The width of the dense matrix
    #   \param[kwarg] buffer        List of doubles containing the elements of the dense matrix
    #
    # To load a numpy matrix from python:
    #   \param[kwarg] matrix        The numpy matrix
    #   \param[kwarg] column_major  Boolean for whether or not the matrix is column major (optional)
    def load_matrix(self, **kwargs):
        if 'filepath' in kwargs:
            if _is_sparse(kwargs['filepath']):
                self.matrix, self.height, self.width = _load_matrix_internal(filepath=kwargs['filepath'])
                self._sparse = True
            else:
                self.dense_matrix, self.height, self.width = _load_matrix_internal(filepath=kwargs['filepath'])
        else:
            if 'row_indices' in kwargs or 'sparse_matrix' in kwargs:
                self._sparse = True
                self.matrix, self.height, self.width = _load_matrix_internal(**kwargs)
            else:
                self._sparse = False
                self.dense_matrix, self.height, self.width = _load_matrix_internal(**kwargs)

    # \brief Loads a dictionary to use for computing top terms
    # \param[kwarg] filepath    The filepath for the dictionary
    # OR
    # \param[kwarg] dictionary  List containing the dictionary strings
    def load_dictionary(self, filepath='', dictionary=[]):
        if filepath != '':
            with open(filepath) as dictionary:
                dictionary = dictionary.read().split("\n")
                dictionary.pop()
            self.dictionary = np.array(dictionary)
        elif len(dictionary) > 0:
            self.dictionary = np.array(dictionary)
        else:
            print 'Error: Invalid dictionary.'

    # \brief Return the top term indices for each cluster
    #
    # The length of the returned array is maxterms*k, with the first maxterms elements belonging 
    # to the first cluster, the second maxterms elements belonging to the second cluster, etc.
    #
    # \return List of the term_indices
    def get_top_indices(self):
        return self.term_indices

    # \brief Return the list of cluster assignments for each document
    # \return List of the assignments
    def get_assignments(self):
        return self.assignments_flat

    # \brief Return the top terms for each cluster
    #
    # The length of the returned array is maxterms*k, with the first maxterms elements belonging 
    # to the first cluster, the second maxterms elements belonging to the second cluster, etc.
    #
    # \return List of the top terms
    def get_top_terms(self, filepath="", dictionary=[]):
        if filepath != "":
            with open(filepath) as dictionary:
                dictionary = dictionary.read().split("\n")
                dictionary.pop()
        terms_indices = self.get_flat_top_indices()
        if len(terms_indices) > 0:
            return [x for i, x in enumerate(dictionary) if i == idx for idx in terms_indices]

#-------------------------------- private functions ---------------------------------#

    def _get_alg(self, alg_name):
        if (alg_name.lower() == 'mu'):
            return 0
        elif (alg_name.lower() == 'hals'):
            return 1
        elif (alg_name.lower() == 'rank2'):
            return 2
        elif (alg_name.lower() == 'bpp'):
            return 3

    def _compute_assignments(self, vector[unsigned int]& flat):
        ComputeAssignments[double](flat, &(self.h[0]), self.k, self.k, self.width)
        return flat

    def _compute_fuzzy_assignments(self, vector[float]& probabilities):
        ComputeFuzzyAssignments[double](probabilities, &(self.h[0]), self.k, self.k, self.width)
        return probabilities

    def _topterms(self):
        cdef vector[int] temp_indices
        temp_indices.resize(self.maxterms*self.k)
        TopTerms(self.maxterms, &(self.w[0]), self.height, self.height, self.k, temp_indices)
        self.term_indices = temp_indices

    def _write_flatclust(self, string& assignfile, string& fuzzyfile, string& resultfile, 
        vector[unsigned int]& assignments,  vector[float]& probabilities, vector[string]& dictionary, 
        vector[int]& term_indices, unsigned int width, const FileFormat form):
        FlatClustWriteResults(assignfile, fuzzyfile, resultfile, assignments, probabilities, dictionary, 
            term_indices, form, self.maxterms, width, self.k)

#####################################################################################
#
#                            Python flatclust functions
#
#####################################################################################

#-------------------------------- public functions ---------------------------------#

cdef class Flatclust(Clustering):

    def __init__(self):
        super(Flatclust, self).__init__()

    # \brief Returns the parsed arguments for the default command line application
    #
    # The command line arguemnts are the same as those for the C++ binary application flatclust
    #
    # \return The dictionary containing the parsed arguments
    def parser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--matrixfile", action="store", 
            required=True,    metavar="matrixfile")
        parser.add_argument("--dictfile",   action="store", 
            required=True,    metavar="dictfile")
        parser.add_argument("--clusters",   action="store", 
            required=True,    metavar="clusters",   type=int)
        parser.add_argument("--algorithm",  action="store", 
            required=False,   metavar="algorithm",    default="BPP", choices=["HALS", "RANK2", "BPP"])
        parser.add_argument("--infile_W",   action="store", 
            required=False,   metavar="infile_W",     default="")
        parser.add_argument("--infile_H",   action="store", 
            required=False,   metavar="infile_H",     default="")
        parser.add_argument("--tol",        action="store", 
            required=False,   metavar="tol",        type=float,  default=0.0001)
        parser.add_argument("--outdir",     action="store", 
            required=False,   metavar="outdir",       default="")
        parser.add_argument("--miniter",    action="store", 
            required=False,   metavar="miniter",    type=int,  default=5)
        parser.add_argument("--maxiter",    action="store", 
            required=False,   metavar="maxiter",    type=int,  default=5000)
        parser.add_argument("--maxterms",   action="store", 
            required=False,   metavar="maxterms",     type = int, default=5)
        parser.add_argument("--maxthreads", action="store", 
            required=False,   metavar="maxthreads",   type = int, default=8)
        parser.add_argument("--verbose",    action="store", 
            required=False,   metavar="verbose",      default=True)
        parser.add_argument("--format",     action="store", 
            required=False,   metavar="format",       default="XML")
        parser.add_argument("--assignfile", action="store", 
            required=False,   metavar="assignfile",   default="assignments")
        parser.add_argument("--treefile", action="store", 
            required=False,   metavar="treefile",   default="tree")
        parser.add_argument("--fuzzyfile", action="store", 
            required=False,   metavar="fuzzyfile",   default="assignments_fuzzy")
        args = parser.parse_args()
        return args

    # \brief Runs NMF on the loaded matrix using the supplied algorithm and implementation details
    # \param[in]    k           The desired number of clusters
    # \param[kwarg] infile_W    Initialization for W (optional)
    # \param[kwarg] infile_H    Initialization for H (optional)
    # \param[kwarg] algorithm   The desired NMF algorithm (optional)
    # \param[kwarg] maxterms    Maximum number of terms per cluster (optional)
    # \param[kwarg] verbose     Boolean for whether or not to be verbose (optional)
    # \param[kwarg] min_iter    Minimum number of iterations (optional)
    # \param[kwarg] max_iter    Maximum number of iterations (optional)
    # \param[kwarg] max_threads Maximum number of threads to use (optional)
    # \param[kwarg] tol         Tolerance for determing convergence (optional)
    def cluster(self, k, infile_W='', infile_H='', algorithm="BPP", maxterms=5, verbose=True, min_iter=5,
        max_iter=5000, max_threads=8, tol=0.0001):
        self.nmfstats = Stats()
        self.tree = TreeResults()
        self.k = k
        self.maxterms = maxterms
        rng = Rand()
        rng.seed_from_time()
        RNG_CENTER = 0.5
        RNG_RADIUS = 0.5
        if not infile_W:
            w_init = _random(self.height, k, rng, RNG_CENTER, RNG_RADIUS)
        else:
            w, height, width = _load_matrix_internal(height=self.height, width=k, filepath=infile_W)
            w_init = w[0]

        if not infile_H:
            h_init = _random(k, self.width, rng, RNG_CENTER, RNG_RADIUS)
        else:
            h, height, width = _load_matrix_internal(height=k, width=self.width, filepath=infile_H) 
            h_init = h[0]
        cdef NmfOptions nmf_opts
        nmf_opts.tol = tol
        nmf_opts.algorithm = self._get_alg(algorithm)
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

        if self._sparse:
            self.w, self.h = self._flatclust_sparse(nmf_opts, self.matrix, w_init, h_init, self.height, self.k, 
                self.nmfstats)
        else:
            self.w, self.h = self._flatclust(nmf_opts, self.dense_matrix, w_init, h_init, self.height, self.height, 
                self.k, self.nmfstats)
        self.assignments_flat = self._compute_assignments(self.assignments_flat)
        self.probabilities = self._compute_fuzzy_assignments(self.probabilities)
        self._topterms()

    # \brief Writes the flatclust results to files
    # \param[in]        assignfile     The filepath for writing assignments
    # \param[in]        fuzzyfile      The filepath for writing fuzzy assignments
    # \param[in]        treefile       The filepath for the tree results
    # \param[kwargs]    outdir         The output directory for the output files (optional)
    # \param[kwargs]    format         The output format JSON or XML (optional)
    def write_output(self, assignfile, fuzzyfile, treefile, outdir='./', format='XML'):
        print 'Writing output files...'
        if format == "XML":
            if 'xml' not in treefile:
                tree = outdir + treefile + "_" + str(self.k) + '.xml'
            else:
                tree = outdir + treefile
        else:
            if 'json' not in treefile:
                tree = outdir + treefile + "_" + str(self.k) + '.json'
            else:
                tree = outdir + treefile
        if 'csv' in assignfile:
            assign = outdir + assignfile 
        else:
            assign = outdir + assignfile + "_" + str(self.k) + '.csv'

        if 'csv' in fuzzyfile:
            fuzzy = outdir + fuzzyfile 
        else:
            fuzzy = outdir + fuzzyfile + "_" + str(self.k) + '.csv'
        self._write_flatclust(assign, fuzzy, tree, self.assignments_flat, self.probabilities, self.dictionary, 
            self.term_indices, self.width, _get_outputformat(format))

#-------------------------------- private functions ---------------------------------#

    def _flatclust_sparse(self, nmf_opts, Sparse A, vector[double]& buf_w, vector[double]& buf_h, unsigned int ldim_w, 
        unsigned int ldim_h, Stats stats): 
        FlatClustSparse(nmf_opts, A.Height(), A.Width(), A.Size(), A.LockedColBuffer(), A.LockedRowBuffer(),
            A.LockedDataBuffer(), &(buf_w[0]), ldim_w, &(buf_h[0]), ldim_h, dereference(stats.get()))
        return buf_w, buf_h

    def _flatclust(self, nmf_opts, vector[double]& buf_a, vector[double]& buf_w, vector[double]& buf_h, 
                unsigned int ldim_a, unsigned int ldim_w, unsigned int ldim_h, Stats stats):
        FlatClust(nmf_opts, &(buf_a[0]), ldim_a, &(buf_w[0]), ldim_w, &(buf_h[0]), ldim_h, dereference(stats.get()))
        return buf_w, buf_h

#####################################################################################
#
#                            Python hierclust functions
#
#####################################################################################

#-------------------------------- public functions ---------------------------------#

cdef class Hierclust(Clustering):
    cdef FileFormat format
    cdef unsigned int flat

    def __init__(self):
        super(Hierclust, self).__init__()

    # \brief Returns the parsed arguments for the default command line application
    #
    # The command line arguemnts are the same as those for the C++ binary application hierclust
    #
    # \return The dictionary containing the parsed arguments
    def parser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--matrixfile", action="store", 
            required=True,  metavar="matrixfile")
        parser.add_argument("--dictfile",   action="store", 
            required=True,  metavar="dictfile")
        parser.add_argument("--clusters",   action="store", 
            required=True,  metavar="clusters",   type=int)
        parser.add_argument("--tol",        action="store", 
            required=False, metavar="tol",        type=float,  default=0.0001)
        parser.add_argument("--outdir",     action="store", 
            required=False, metavar="outdir",       default="")
        parser.add_argument("--miniter",    action="store", 
            required=False, metavar="miniter",    type=int,  default=5)
        parser.add_argument("--maxiter",    action="store", 
            required=False, metavar="maxiter",    type=int,  default=5000)
        parser.add_argument("--maxterms",   action="store", 
            required=False, metavar="maxterms",   type=int,  default=5)
        parser.add_argument("--maxthreads", action="store", 
            required=False, metavar="maxthreads", type=int,  default=8)
        parser.add_argument("--unbalanced", action="store", 
            required=False, metavar="unbalanced", type=float,  default=0.1)
        parser.add_argument("--trial_allowance", action="store", 
            required=False, metavar="trial_allowance", type=int, default=3)
        parser.add_argument("--flat",       action="store", 
            required=False, metavar="flat",  type=int,       default=0)
        parser.add_argument("--verbose",    action="store", 
            required=False, metavar="verbose",      default=True)
        parser.add_argument("--format",     action="store", 
            required=False, metavar="format",       default="XML", choices=["XML", "JSON"])
        parser.add_argument("--treefile",  action="store", 
            required=False, metavar="treefile",    default="tree")
        parser.add_argument("--initdir",  action="store", 
            required=False, metavar="initdir",    default="")
        parser.add_argument("--assignfile", action="store", 
            required=False, metavar="assignfile",   default="assignments")
        parser.add_argument("--fuzzyfile", action="store", 
            required=False,   metavar="fuzzyfile",   default="assignments_fuzzy")
        args = parser.parse_args()
        return args

    # \brief Runs HierNMF2 on the loaded matrix using the supplied algorithm and implementation details
    # \param[in]    k               The desired number of clusters
    # \param[kwarg] initdir         Initialization matrices (optional)
    # \param[kwarg] maxterms        Maximum number of terms per cluster (optional)
    # \param[kwarg] unbalanced      Unbalanced parameter (optional)
    # \param[kwarg] trial_allowance Number of trials to use (optional)
    # \param[kwarg] verbose         Boolean for whether or not to be verbose (optional)
    # \param[kwarg] flat            Whether or not to flatten the results (optional)
    # \param[kwarg] min_iter        Minimum number of iterations (optional)
    # \param[kwarg] max_iter        Maximum number of iterations (optional)
    # \param[kwarg] max_threads     Maximum number of threads to use (optional)
    # \param[kwarg] tol             Tolerance for determing convergence (optional)
    def cluster(self, k, initdir='', maxterms=5, unbalanced=0.1, trial_allowance=3, verbose=True, flat=0,
        min_iter=5, max_iter=5000, max_threads=8, tol=0.0001):
        self.stats = ClusterStats()
        self.tree = TreeResults()
        self.k = k
        self.flat = flat
        self.maxterms = maxterms
        rng = Rand()
        rng.seed_from_time()
        RNG_CENTER = 0.5
        RNG_RADIUS = 0.5
        required_size = self.height*2

        self.clust_opts.maxterms                      = int(maxterms)
        self.clust_opts.initdir                       = initdir
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

        if self._sparse:
            self.w, self.h = self._clust_sparse(self.matrix, self.tree, self.stats, self.height, self.width, 
                k, rng, initdir)
        else:
            self.w, self.h, self.assignments = self._clust_dense(self.dense_matrix, self.height, self.tree, 
                self.stats, self.height, self.width, k, rng)
        if self.flat == 1:
            self.assignments_flat = self._compute_assignments(self.assignments_flat)
            self.probabilities = self._compute_fuzzy_assignments(self.probabilities)
            self._topterms()

    # \brief Return the top term indices for each cluster
    #
    # The length of the returned array is maxterms*k, with the first maxterms elements belonging 
    # to the first cluster, the second maxterms elements belonging to the second cluster, etc.
    #
    # \return List of the term_indices
    def get_top_indices(self):
        if self.flat == 1:
            return self.term_indices
        else:
            print 'ERROR: To get top terms indices, rerun hierarchical clusting with flat=1'
            return []

    # \brief Writes the flatclust results to files
    # \param[in]        assignfile     The filepath for writing assignments
    # \param[in]        fuzzyfile      The filepath for writing fuzzy assignments
    # \param[in]        treefile       The filepath for the tree results
    # \param[kwargs]    outdir         The output directory for the output files (optional)
    # \param[kwargs]    format         The output format JSON or XML (optional)
    def write_output(self, assignfile, treefile, fuzzyfile, outdir='./', format='XML'):
        print 'Writing output files...'
        if format == "XML":
            if 'xml' not in treefile:
                tree = outdir + treefile + "_" + str(self.k) + '.xml'
            else:
                tree = outdir + treefile
        else:
            if 'json' not in treefile:
                tree = outdir + treefile + "_" + str(self.k) + '.json'
            else:
                tree = outdir + treefile

        if 'csv' in assignfile:
            assign = outdir + assignfile
        else:
            assign = outdir + assignfile + "_" + str(self.k) + '.csv'

        if 'csv' in fuzzyfile:
            fuzzy = outdir + fuzzyfile
        else:
            fuzzy = outdir + fuzzyfile + "_" + str(self.k) + '.csv'
        self.tree.write_assignments(assign)
        self.tree.write(tree, self.dictionary, _get_outputformat(format))
        if self.flat == 1:
            self._write_flatclust(assign, fuzzy, tree, self.assignments_flat, self.probabilities, self.dictionary, 
                self.term_indices, self.width, _get_outputformat(format))

    # \brief Return the list of cluster assignments for each document
    # \return List of the assignments
    def get_assignments(self):
        cdef vector[unsigned int] assignments
        cdef list assignments_python
        assignments_python = []
        if self.flat:
            return self.assignments_flat
        else:
            for each in self.tree.get_assignments(assignments):
                if each == 4294967295:
                    each = -1
                assignments_python.append(each)
            return assignments_python

#-------------------------------- private functions ---------------------------------#

    def _write_assignments(self, vector[unsigned int]& labels, const string& filepath):
        return WriteAssignmentsFile(labels, filepath)

    def _clust_dense(self, vector[double]& buf_a, unsigned int ldim_a, TreeResults tree, ClusterStats stats,
            unsigned int m, unsigned int n, unsigned int num_clusters, Rand rng):
        self.buf_w.resize(m*num_clusters)
        self.buf_h.resize(n*num_clusters)
        Clust(self.clust_opts, &(buf_a[0]), ldim_a, &(self.buf_w[0]), &(self.buf_h[0]), 
                                dereference(tree.get()), dereference(stats.get()), dereference(rng.get()))
        return self.buf_w, self.buf_h, self.assignments

    def _clust_sparse(self, Sparse A, TreeResults tree, ClusterStats stats, unsigned int m, unsigned int n, 
        unsigned int num_clusters, Rand rng, string initdir):
        self.buf_w.resize(m*num_clusters)
        self.buf_h.resize(n*num_clusters)
        res = ClustSparse(self.clust_opts, dereference(A.get()), &(self.buf_w[0]), &(self.buf_h[0]), 
            dereference(tree.get()), dereference(stats.get()), dereference(rng.get()))
        return self.buf_w, self.buf_h

#####################################################################################
#
#                            Python matrixgen functions
#
#####################################################################################

#-------------------------------- public functions ---------------------------------#

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

    # \brief Returns the parsed arguments for the default command line application
    #
    # The command line arguments are the same as those for the C++ binary application matrixgen
    #
    # \return The dictionary containing the parsed arguments
    def parser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--height",         action="store", 
            required=True,  metavar="height", type=int)
        parser.add_argument("--width",          action="store", 
            required=True, metavar="width", type=int)
        parser.add_argument("--filename",       action="store", 
            required=True, metavar="filename", type=str)
        parser.add_argument("--type",           action="store", 
            required=False, metavar="type", default='UNIFORM', 
            choices=['UNIFORM', 'DENSE_DIAG', 'SPARSE_DIAG','IDENTITY', 'ONES', 'ZEROS', 'SPARSE'])
        parser.add_argument("--rng_center",     action="store", 
            required=False, metavar="rng_center", default=0.5, type=float)
        parser.add_argument("--rng_radius",     action="store", 
            required=False, metavar="rng_radius", default=0.5, type=float)
        parser.add_argument("--precision",      action="store", 
            required=False, metavar="precision", default=6, type=int)
        parser.add_argument("--nz_per_col",     action="store", 
            required=False, metavar="nz_per_col", default=1, type=int)
        args = parser.parse_args()
        return args

    # \brief Writes the generated matrix to file
    # \param[in]        filename     The filepath for writing the matrix
    # \param[kwarg]     precision    The precision with which to write the matrix
    def write_output(self, filename, precision=6):
        if self.is_sparse:
            if not _write_mtx(filename, self.S, precision):
                print 'Matrixgen error - sparse matrix file write failed.'
            else:
                print 'Matrix succesfully written.'
        else:
            if not _write_delimited(self.A, filename, precision, self.height, self.height, self.width):
                print 'Matrixgen error - file write failed.'
            else: 
                print 'Matrix succesfully written.'

    # \brief Generates a uniform matrix
    # \param[in]        m       The desired height
    # \param[in]        n       The desired width
    # \param[kwarg]     center  Center with which to initialize the RNG 
    # \param[kwarg]     radius  Radius with which to initialize the RNG 
    # \return A tuple of the list of values, the height, and the width
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
        return A, m, n
    
    # \brief Generates a dense diagonal matrix
    # \param[in]        m       The desired height
    # \param[in]        n       The desired width
    # \param[kwarg]     center  Center with which to initialize the RNG 
    # \param[kwarg]     radius  Radius with which to initialize the RNG 
    # \return A tuple of the list of values, the height, and the width
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
        return A, m, n

    # \brief Generates an identify matrix of height m and width n.
    # \param[in]        m       The desired height
    # \param[in]        n       The desired width
    # \return A tuple of the list of values, the height, and the width
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
        return A, m, n
    
    # \brief Generates a sparse diagonal matrix of width n with the RNG attributes of center and radius
    # \param[in]        n       The desired width
    # \param[kwarg]     center  Center with which to initialize the RNG 
    # \param[kwarg]     radius  Radius with which to initialize the RNG 
    # \return A tuple of the list of values, the height, and the width
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
        return S, n, n

    # \brief Generates a matrix of ones of height m and width n
    # \param[in]        m       The desired height
    # \param[in]        n       The desired width
    # \return A tuple of the list of values, the height, and the width
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
        return A, m, n

    # \brief Generates a matrix of zeros of height m and width n
    # \param[in]        m       The desired height
    # \param[in]        n       The desired width
    # \return A tuple of the list of values, the height, and the width
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
        return A, m, n

    # \brief Generates a random sparse matrix
    # \param[in]        m       The desired height
    # \param[in]        n       The desired width
    # \param[in]        nz      The desired non-zeros
    # \return A tuple of the list of values, the height, and the width
    def sparse(self, unsigned int m, unsigned int n, unsigned int nz):
        self.is_sparse = True
        self.height = m
        self.width = n
        S = Sparse()
        RandomSparseMatrix(dereference(self.rng.get()), dereference(S.get()), nz, m, m, n, n)
        self.S = S
        return S, m, n

#####################################################################################
#
#                            Python preprocessor functions
#
#####################################################################################

#-------------------------------- public functions ---------------------------------#

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
    cdef TFData tfdataArray
    cdef list row_indices, counts, col_offsets, term_ind, doc_ind

    def __init__(self):
        pass

    # \brief Returns the parsed arguments for the default command line application
    #
    # The command line arguemnts are the same as those for the C++ binary application preprocessor
    #
    # \return The dictionary containing the parsed arguments
    def parser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--indir",         action="store", 
            required=True,  metavar="indir")
        parser.add_argument("--outdir",        action="store", 
            required=False, metavar="outdir",        default="./")
        parser.add_argument("--docs_per_term", action="store", 
            required=False, metavar="docs_per_term", default=3)
        parser.add_argument("--terms_per_doc", action="store", 
            required=False, metavar="terms_per_doc", default=5)
        parser.add_argument("--maxiter",       action="store", 
            required=False, metavar="maxiter",       default=1000)
        parser.add_argument("--precision",     action="store", 
            required=False, metavar="precision",     default=4)
        parser.add_argument("--boolean_mode",  action="store", 
            required=False, metavar="boolean_mode",     default=0)
        args = parser.parse_args()
        return args

    # \brief Writes the preprocessor results to files
    # \param[in]        matrix_filepath     The filepath for writing the matrix
    # \param[in]        dict_filepath       The filepath for writing the dictionary
    # \param[in]        docs_filepath       The filepath for the documents
    # \param[kwargs]    precision           The precision with which to write the outputs (optional)
    def write_output(self, matrix_filepath, dict_filepath, docs_filepath, precision=4):
        _write_termfreq(self.tfm, matrix_filepath, self.scores, precision)
        _write_strings(dict_filepath, self.dictionary, self.term_ind, self.height)
        _write_strings(docs_filepath, self.documents, self.doc_ind, self.width)

    # \brief Preprocesses the matrix
    # \param[kwarg] maxiter      The maximum number of iterations (optional)
    # \param[kwarg] docsperterm  The number of documents required per term (optional)
    # \param[kwarg] termsperdoc  The number of terms requried per document (optional)
    # \param[kwarg] boolean_mode All nonzero matrix elements will be treated as if they had 
    #                            the value 1.0  (optional)
    def preprocess(self, maxiter=1000, docsperterm=3, termsperdoc=5, boolean_mode=0):
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

    # \brief Load an input matrix
    #
    # To load a matrix from a file:
    #   \param[kwarg] filepath      The path to the input matrix
    #
    # To load a sparse matrix from Matrixgen:
    #   \param[kwarg] height        The height of the sparse matrix
    #   \param[kwarg] width         The width of the sparse matrix
    #   \param[kwarg] sparse_matrix The sparse matrix returned from Matrixgen
    #
    # To load a sparse matrix from python:
    #   \param[kwarg] height        The height of the sparse matrix
    #   \param[kwarg] width         The width of the sparse matrix
    #   \param[kwarg] nz            The number of non-zeros in the sparse matrix
    #   \param[kwarg] buffer        List of doubles containing the non-zero elements of the sparse matrix
    #   \param[kwarg] row_indices   List of integers representing the row indices of the sparse matrix
    #   \param[kwarg] col_offsets   List of integers representing the column offsets of the sparse matrix
    def load_matrix(self, filepath="", height=0, width=0, nz=0, buffer=[], row_indices=[], 
        col_offsets=[], sparse_matrix=None):
        if filepath != "":
            self.A, self.height, self.width = _load_matrix_internal(filepath=filepath)
        elif len(row_indices) > 0:
            self.A, self.height, self.width = _load_matrix_internal(height=height, width=width, 
                nz=nz, buffer=buffer, row_indices=row_indices, col_offsets=col_offsets)

    # \brief Loads a dictionary
    # \param[kwarg] filepath    The filepath for the dictionary
    # OR
    # \param[kwarg] dictionary  List containing the dictionary strings
    def load_dictionary(self, filepath="", dictionary=[]):
        if filepath != "":
            with open(filepath) as f:
                self.dictionary = f.read().split('\n')
        else:
            self.dictionary = dictionary

    # \brief Loads the documents file
    # \param[kwarg] filepath    The filepath for the documents
    # OR
    # \param[kwarg] documents   List containing the docuemnts strings
    def load_documents(self, filepath="", documents=[]):
        if filepath != "":
            with open(filepath) as f:
                self.documents = f.read().split('\n')
        else:
            self.documents = documents
 
    # \brief Returns the reduced documents
    # \return The documents in a list
    def get_reduced_documents(self):
        return [self.documents[i] for i in self.doc_ind]

    # \brief Returns the reduced dictionary
    # \return The dictionary in a list
    def get_reduced_dictionary(self):
        return [self.dictionary[i] for i in self.term_ind]

    # \brief Returns the non-zero scores from the reduced matrix
    # \return The scores in a list
    def get_reduced_scores(self):
        return self.scores

    # \brief Returns the row indices for the reduced matrix
    # \return The row indices in a list 
    def get_reduced_row_indices(self):
        return self.row_indices

    # \brief Returns the column offsets for the reduced matrix
    # \return The column offsets in a list
    def get_reduced_col_offsets(self):
        return self.col_offsets

    # \brief Loads the additional field file
    # \param[kwarg] filepath    The filepath for the field
    # OR
    # \param[kwarg] values      List containing the field strings
    def get_reduced_field(self, filepath="", values=[]):
        if filepath != "":
            with open(filepath) as f:
                values = f.read().split('\n')
        return [values[i] for i in self.doc_ind]

