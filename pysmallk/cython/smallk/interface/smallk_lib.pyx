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

# nmf_main.pyx: numpy arrays -> extern from "nmf_main.h"
# /Users/barrydrake/development/python/spyder_proj/repositories/xdata2/nmf/cython/interface/nmf_lib.pyx
# 3 steps:
# cython f.pyx  -> f.c
# link: python f-setup.py build_ext --inplace  -> f.so, a dynamic library
# py test-f.py: import f gets f.so, f.fpy below calls fc()

cimport numpy as np
import cython

from libc.stdlib cimport malloc, free
from libc.string cimport memcpy, memset
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "smallk.hpp" namespace "smallk":
    cdef enum Algorithm: 
        MU     # Lee & Seung, multiplicative updating
        HALS   # Cichocki & Pan, hierarchical alternating least squares
        RANK2  # Kuang and Park, rank2 specialization
        BPP   

    cdef enum OutputFormat:
        XML    # PG_i / PG_1 (ratio of the ith PG to that of iteration 1)
        JSON   # relative change in the Frobenius norm of W

    ctypedef char *stdStringR "const std::string &"

    void Initialize(int& argc, char **&argv)
    int IsInitialized()
    void Finalize()
    
    # version info
    unsigned int GetMajorVersion()
    unsigned int getMinorVersion()
    unsigned int GetPatchLevel()
    string GetVersionString()

    void LoadMatrix(stdStringR filepath)
    void Nmf(int k, Algorithm algorithm, stdStringR initfile_w, stdStringR initfile_h)

    void SetOutputPrecision(const unsigned int num_digits)
    void SetMinIter(const unsigned int min_iterations)
    void SetNmfTolerance(const double tol)
    void SetMaxIter(const unsigned int max_iterations)
    void SetMaxThreads(const unsigned int mt)
    void SetOutputDir(const string& output_dir)

#    void SetOutputFormat(const OutputFormat form)
#    void LoadDictionary(const string& filepath)
#    bool HierNmf2(const unsigned int num_clusters)
    
# end cdef extern

@cython.boundscheck(False)
@cython.wraparound(False)

def get_algorithm(alg_name):
    if (alg_name == 'MU'):
        return MU
    elif (alg_name == 'HALS'):
        return HALS
    elif (alg_name == 'RANK2'):
        return RANK2
    elif (alg_name == 'BPP'):
        return BPP

def get_outputformat(form):
    if form == "XML":
        return XML
    elif form == "JSON":
        return JSON


def py_initialize(ac, av):
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

    # Initialize
    Initialize(ac, c_arr)

    # free up all allocated memory
    for i in xrange(ac):
        free(c_arr[i])
    free(c_arr)

def py_isInitialized():
    return IsInitialized();

def py_loadMatrix(char* filepath):
    LoadMatrix(filepath)

def py_nmf(k, algorithm, char* initfile_w, char* initfile_h):
    '''
    Calls non-distributed NMF
      k           Rank k
      algorithm   Algorithm used for NMF
                    'MU' | 'HALS' | 'RANK2' | 'BPP'
      initfile_w  Path to W Matrix file
      initfile_h  Path to H Matrix file
    '''
    cdef int c_k = k
    Nmf(c_k, algorithm, initfile_w, initfile_h)

def py_finalize(): 
    Finalize()

def PySetOutputPrecision(const unsigned int num_digits):
    SetOutputPrecision(num_digits)

def PySetMinIter(const unsigned int min_iterations):
    SetMinIter(min_iterations)

def PySetNmfTolerance(const double tol):
    SetNmfTolerance(tol)

def PySetMaxIter(const unsigned int max_iterations):
    SetMaxIter(max_iterations)

def PySetMaxThreads(const unsigned int mt):
    SetMaxThreads(mt)

def PySetOutputDir(const string& output_dir):
    SetOutputDir(output_dir)

#unused for running regular nmf
#def PySetOutputFormat(const OutputFormat form):
#    SetOutputFormat(form)
#
#def PyLoadDictionary(const string& filepath):
#    LoadDictionary(filepath)
#
#def PyHierNmf2(const unsigned int num_clusters):
#    HierNmf2(num_clusters)
# end smallklibpy

