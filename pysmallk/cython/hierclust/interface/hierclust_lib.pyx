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
import numpy 
cimport numpy
from libc.stdlib cimport malloc, free
from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector
from cython.operator import dereference
import time
from libc.string cimport memcpy, memset


cdef extern from "random.hpp":
    cdef cppclass Random:
        Random()
        void SeedFromTime()

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
        HALS
        RANK2
        BPP

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

cdef extern from "file_format.hpp":
    cdef enum FileFormat:
        CSV
        XML
        JSON

cdef extern from "sparse_matrix_decl.hpp":
    cdef cppclass SparseMatrix[T]:
        SparseMatrix()
        unsigned int Height()
        unsigned int Width()
        unsigned int Size()

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
        
cdef extern from "command_line.hpp":
    cdef struct CommandLineOptions:
        ClustOptions clust_opts
        string infile_A
        string infile_W
        string infile_H
        string dictfile
        string outdir
        string treefile
        string assignfile
        FileFormat format
        bool show_help

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

cdef extern from "matrix_generator.hpp":
    bool RandomMatrix(vector[double]& buf, 
                      const unsigned int height,
                      const unsigned int width,
                      Random& rng,
                      const double rng_center,
                      const double rng_radius)

cdef extern from "delimited_file.hpp":
    bool LoadMatrixArray[T](vector[vector[T]]& buf, 
                            const unsigned int height,
                            const unsigned int width,
                            const string& filename,
                            const char DELIM)

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
    void FlatClustWriteResults(const string& outdir,
                               const vector[int]& assignments,
                               const vector[string]& dictionary,
                               const vector[int]& term_indices,
                               const FileFormat format,
                               const unsigned int maxterms,
                               const unsigned int num_docs,
                               const unsigned int num_clusters)

####################################################################################
@cython.boundscheck(False)
@cython.wraparound(False)
####################################################################################

cdef class PyRandom:
    cdef Random* thisptr
    def __cinit__(self):
        self.thisptr = new Random()
    def SeedFromTime(self):
        self.thisptr.SeedFromTime()
    cdef Random* getThis(self):
        return self.thisptr

cdef class PyDoubleSparseMatrix:
    cdef SparseMatrix[double] *thisptr
    def __cinit__(self):
        self.thisptr = new SparseMatrix[double]()
    def Height(self):
        return self.thisptr.Height()
    def Width(self):
        return self.thisptr.Width()
    def Size(self):
        return self.thisptr.Size()
    cdef SparseMatrix[double]* getThis(self):
        return self.thisptr

cdef class PyTree:
    cdef Tree* thisptr
    def __cinit__(self):
        self.thisptr = new Tree()
    def Write(self, const string& filepath, const FileFormat& form, const vector[string]& dictionary):
        return self.thisptr.Write(filepath, form, dictionary)
    cdef Tree* getThis(self):
        return self.thisptr

cdef class PyClustStats:
    cdef ClustStats* thisptr
    def __cinit__(self):
        self.thisptr = <ClustStats*>malloc(sizeof(ClustStats))
        self.thisptr.nmf_count = 0
        self.thisptr.max_count = 0
    cdef ClustStats* getThis(self):
        return self.thisptr
    property nmf_count:
        def __get__(self):
            return self.thisptr.nmf_count
    property max_count:
        def __get__(self):
            return self.thisptr.max_count


def PyNmfInitialize(ac, av):
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

def PyNmfIsInitialized():
    return NmfIsInitialized() == INITIALIZED

def PyNmfFinalize():
    NmfFinalize()

def PyIsDense(string infile):
    return IsDense(infile)
def PyIsSparse(string infile):
    return IsSparse(infile)
def PyLoadSparseMatrix(const string& file_path, PyDoubleSparseMatrix A):
    cdef unsigned int m = 0, n = 0, nnz = 0
    cdef bool ok = LoadSparseMatrix(file_path, dereference(A.getThis()), m, n, nnz)
    return (ok, m, n, nnz)
def PyLoadDenseMatrix(const string& file_path):
    cdef unsigned int m = 0, n = 0
    cdef vector[double] buf_a
    cdef bool ok = LoadDenseMatrix(file_path, buf_a, m, n)
    return (ok, m, n, buf_a)

def PyRandomMatrix(const unsigned int height, const unsigned int width,
                   PyRandom r,
                   const double r_center, const double r_radius, unsigned int required_size):
    cdef vector[double] buf
    buf.resize(required_size)
    cdef bool ok = RandomMatrix(buf, height, width, dereference(r.getThis()), r_center, r_radius)
    return (ok, buf)

def PyLoadMatrixArray(const unsigned int matrix_height, const unsigned int matrix_width, 
                      const string& filename, unsigned int num_initializers, const char DELIM = ','):
    cdef vector[vector[double]] buf
    buf.resize(num_initializers)
    cdef bool ok = LoadMatrixArray[double](buf, matrix_height, matrix_width, filename, DELIM)
    return (ok, buf)

def PyLoadMatrixArray(const unsigned int matrix_height, const unsigned int matrix_width, vector[double]& buf):
    pass

def PyClust(opts, vector[double]& buf_a, unsigned int ldim_a,
            vector[vector[double]]& w_initializers, vector[vector[double]]& h_initializers,
            PyTree tree, PyClustStats stats,
            unsigned int m, unsigned int n, unsigned int num_clusters):
    cdef vector[double] buf_w, buf_h
    buf_w.resize(m*num_clusters)
    buf_h.resize(n*num_clusters)
    cdef vector[int] assignments

    cdef ClustOptions clust_opts
    temp = opts["clust_opts"]
    clust_opts.maxterms                      = temp["maxterms"]
    clust_opts.unbalanced                    = temp["unbalanced"]
    clust_opts.trial_allowance               = temp["trial_allowance"]
    clust_opts.num_clusters                  = temp["num_clusters"]
    clust_opts.verbose                       = temp["verbose"]
    clust_opts.flat                          = temp["flat"]
    clust_opts.nmf_ops.tol                   = temp["nmf_opts"]["tol"]
    clust_opts.nmf_opts.algorithm            = temp["nmf_ops"]["algorithm"]
    clust_opts.nmf_opts.prog_est_algorithm   = temp["nmf_opts"]["prog_est_algorithm"]
    clust_opts.nmf_opts.height               = temp["nmf_opts"]["height"]
    clust_opts.nmf_opts.width                = temp["nmf_opts"]["width"]
    clust_opts.nmf_opts.k                    = temp["nmf_opts"]["k"]
    clust_opts.nmf_opts.min_iter             = temp["nmf_opts"]["min_iter"]
    clust_opts.nmf_opts.max_iter             = temp["nmf_opts"]["max_iter"]
    clust_opts.nmf_opts.tolcount             = temp["nmf_opts"]["tolcount"]
    clust_opts.nmf_opts.max_threads          = temp["nmf_opts"]["max_threads"]
    clust_opts.nmf_opts.verbose              = temp["nmf_opts"]["verbose"]
    clust_opts.nmf_opts.normalize            = temp["nmf_opts"]["normalize"]

    cdef Result res = Clust(clust_opts, &(buf_a[0]), ldim_a, &(buf_w[0]), &(buf_h[0]),
                            w_initializers, h_initializers, assignments, 
                            dereference(tree.getThis()), dereference(stats.getThis()))
    return (res == OK, buf_w, buf_h, assignments)

def PyClustSparse(opts, PyDoubleSparseMatrix A, 
                  vector[vector[double]]& w_initializers, vector[vector[double]]& h_initializers,
                  #w_initializers, h_initializers,
                  PyTree tree, PyClustStats stats,
                  unsigned int m, unsigned int n, unsigned int num_clusters):
    cdef vector[double] buf_w, buf_h
    buf_w.resize(m*num_clusters)
    buf_h.resize(n*num_clusters)
    cdef vector[int] assignments
        
    cdef ClustOptions clust_opts
    temp = opts["clust_opts"]
    clust_opts.maxterms                       = temp["maxterms"]
    clust_opts.unbalanced                     = temp["unbalanced"]
    clust_opts.trial_allowance                = temp["trial_allowance"]
    clust_opts.num_clusters                   = temp["num_clusters"]
    clust_opts.verbose                        = temp["verbose"]
    clust_opts.flat                           = temp["flat"]
    clust_opts.nmf_opts.tol                   = temp["nmf_opts"]["tol"]
    clust_opts.nmf_opts.algorithm             = temp["nmf_opts"]["algorithm"]
    clust_opts.nmf_opts.prog_est_algorithm    = temp["nmf_opts"]["prog_est_algorithm"]
    clust_opts.nmf_opts.height                = temp["nmf_opts"]["height"]
    clust_opts.nmf_opts.width                 = temp["nmf_opts"]["width"]
    clust_opts.nmf_opts.k                     = temp["nmf_opts"]["k"]
    clust_opts.nmf_opts.min_iter              = temp["nmf_opts"]["min_iter"]
    clust_opts.nmf_opts.max_iter              = temp["nmf_opts"]["max_iter"]
    clust_opts.nmf_opts.tolcount              = temp["nmf_opts"]["tolcount"]
    clust_opts.nmf_opts.max_threads           = temp["nmf_opts"]["max_threads"]
    clust_opts.nmf_opts.verbose               = temp["nmf_opts"]["verbose"]
    clust_opts.nmf_opts.normalize             = temp["nmf_opts"]["normalize"]
    
    cdef Result res = ClustSparse(clust_opts, dereference(A.getThis()), &(buf_w[0]), &(buf_h[0]), 
                                  w_initializers, h_initializers, assignments, 
                                  #w_inits, h_inits, assignments,
                                  dereference(tree.getThis()), dereference(stats.getThis()))
    return (res == OK, buf_w, buf_h, assignments)

def PyComputeAssignments(vector[double]& buf_h, unsigned int ldim_h, unsigned int k, unsigned int n):
    cdef vector[int] assignments_flat
    ComputeAssignments[double](assignments_flat, &(buf_h[0]), ldim_h, k, n)
    return assignments_flat

def PyTopTerms(const int maxterms, vector[double]& buf_w, const unsigned int ldim, 
               const unsigned int height, const unsigned int width, unsigned int num_clusters):
    cdef vector[int] term_indices
    term_indices.resize(maxterms, num_clusters)
    TopTerms(maxterms, &(buf_w[0]), ldim, height, width, term_indices)
    return term_indices

def PyWriteAssignmentsFile(vector[int]& labels, const string& filepath):
    return WriteAssignmentsFile(labels, filepath)

def PyFlatClustWriteResults(const string& outdir, const vector[int]& assignments, const vector[string]& dictionary,
                            const vector[int]& indices, const FileFormat form, const unsigned int maxterms,
                            const unsigned int num_docs, const unsigned int num_clusters):
    FlatClustWriteResults(outdir, assignments, dictionary, indices, form, maxterms, num_docs, num_clusters)

def PyParseCommandLine(args):
    opts = {}
    opts["clust_opts"] = {}
    opts["clust_opts"]["nmf_opts"] = {}
    opts["clust_opts"]["nmf_opts"]["height"]    = 0
    opts["clust_opts"]["nmf_opts"]["width"]     = 0
    opts["clust_opts"]["nmf_opts"]["k"]         = 0
    opts["clust_opts"]["nmf_opts"]["min_iter"]  = args.miniter
    opts["clust_opts"]["nmf_opts"]["max_iter"]  = args.maxiter
    opts["clust_opts"]["nmf_opts"]["tol"]       = args.tol
    opts["clust_opts"]["nmf_opts"]["tolcount"]  = 1
    opts["clust_opts"]["nmf_opts"]["verbose"]   = (args.verbose != 0)  
    opts["clust_opts"]["nmf_opts"]["normalize"] = False
    opts["clust_opts"]["nmf_opts"]["algorithm"] = RANK2
    
    
    opts["clust_opts"]["nmf_opts"]["prog_est_algorithm"] = PG_RATIO
    
    opts["clust_opts"]["maxterms"] = args.maxterms
    opts["clust_opts"]["trial_allowance"] = args.trial_allowance
    opts["clust_opts"]["unbalanced"] = args.unbalanced
    opts["clust_opts"]["num_clusters"] = args.clusters
    opts["clust_opts"]["verbose"] = True
    opts["clust_opts"]["flat"] = (args.flat != 0)
                                                                                                        
    opts["infile_A"] = args.matrixfile
    opts["infile_W"] = args.infile_W
    opts["infile_H"] = args.infile_H
    opts["dictfile"] = args.dictfile
    opts["outdir"] = args.outdir
    if not args.treefile:
        opts["treefile"] = args.outdir+"tree_%d."%args.clusters+args.format.lower()
    else:
        opts["treefile"] = args.treefile

    if not args.assignfile:
        opts["assignfile"] = args.outdir+"assignments_%d.csv"%args.clusters
    else:
        opts["assignfile"] = args.assignfile
    
    opts["show_help"] = False
    if args.format == "XML":
        opts["format"] = XML
    elif args.format == "JSON":
        opts["format"] = JSON

    # see command_line.cpp for correct source; going with this for now
    opts["clust_opts"]["nmf_opts"]["max_threads"] = args.maxthreads

    return opts
