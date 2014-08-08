import cython
import numpy 
cimport numpy
from libc.stdlib cimport malloc, free
from libcpp.string cimport string


from libcpp cimport bool
from libcpp.vector cimport vector
from cython.operator import dereference
from libc.string cimport strdup, strcpy
from libc.string cimport memcpy, memset


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
                           
cdef extern from "command_line.hpp":
    cdef struct CommandLineOptions:
        FlatClustOptions clust_opts
        string infile_A
        string infile_W
        string infile_H
        string dictfile
        string outdir
        string clustfile
        string assignfile
        FileFormat format
        bool show_help
    #void ParseCommandLine(int argc, char* argv[], CommandLineOptions& opts)

cdef extern from "nmf.hpp":
    void NmfInitialize(int argc, char* argv[])
    Result NmfIsInitialized()
    void NmfFinalize()

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

cdef extern from "file_format.hpp":
    cdef enum FileFormat:
        CSV
        XML
        JSON

cdef extern from "assignments.hpp":
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

################################################################################
@cython.boundscheck(False)
@cython.wraparound(False)
################################################################################

cdef class PyRandom:
    cdef Random* thisptr
    def __cinit__(self):
        self.thisptr = new Random()
    def SeedFromTime(self):
        self.thisptr.SeedFromTime()
    cdef Random* getThis(self):
        return self.thisptr

def PyRandomMatrix(const unsigned int height, const unsigned int width,
                   PyRandom r,
                   const double r_center, const double r_radius):
    cdef vector[double] buf
    cdef bool ok = RandomMatrix(buf, height, width, dereference(r.getThis()), r_center, r_radius)
    return (ok, buf)

def PyLoadDelimitedFile(unsigned int height, unsigned int width,
                        const string& filename, const char DELIM = ','):
    cdef vector[double] buf
    cdef bool ok = LoadDelimitedFile(buf, height, width, filename, DELIM)
    return (ok, buf)



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
    cdef const unsigned int* LockedColBuffer(self):
        return self.thisptr.LockedColBuffer()
    cdef const unsigned int* LockedRowBuffer(self):
        return self.thisptr.LockedRowBuffer()
    cdef const double* LockedDataBuffer(self):
        return self.thisptr.LockedDataBuffer()
    cdef SparseMatrix[double]* getThis(self):
        return self.thisptr

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

cdef class PyNmfStats:
    cdef NmfStats* thisptr
    def __cinit__(self):
        self.thisptr = <NmfStats*>malloc(sizeof(NmfStats))
        self.thisptr.elapsed_us = 0
        self.thisptr.iteration_count = 0
    cdef NmfStats* getThis(self):
        return self.thisptr
    property elapsed_us:
        def __get__(self):
            return self.thisptr.elapsed_us
    property iteration_count:
        def __get__(self):
            return self.thisptr.iteration_count

cdef class PyCommandLineOptions:
    cdef CommandLineOptions* thisptr
    def __cinit__(self):
        self.thisptr = <CommandLineOptions*>malloc(sizeof(CommandLineOptions))
    cdef CommandLineOptions* getThis(self):
        return self.thisptr

# 
# BUG: SEE BUG LINE
# 
def PyParseCommandLine(ac, av, PyCommandLineOptions opts, args):
    #cdef char **c_arr = to_cstring_array(av)
    #ParseCommandLine(ac, c_arr, dereference(opts.getThis()))
    #free(c_arr)
    #cdef int user_max_threads = -1
    #(opts.getThis()).clust_opts.nmf_opts.height = 0
    #(opts.getThis()).clust_opts.nmf_opts.width = 0
    #(opts.getThis()).clust_opts.nmf_opts.k = 0
    (opts.getThis()).clust_opts.nmf_opts.min_iter = 5
    print (opts.getThis()).clust_opts.nmf_opts.min_iter
    #(opts.getThis()).clust_opts.nmf_opts.max_iter = 5000
    #(opts.getThis()).clust_opts.nmf_opts.tol = 0.0001
    #(opts.getThis()).clust_opts.nmf_opts.tolcount = 1
    #(opts.getThis()).clust_opts.nmf_opts.verbose = True
    #(opts.getThis()).clust_opts.nmf_opts.normalize = True
    #(opts.getThis()).clust_opts.nmf_opts.algorithm = BPP

    #(opts.getThis()).clust_opts.nmf_opts.prog_est_algorithm = PG_RATIO

    #(opts.getThis()).clust_opts.maxterms = 5
    #(opts.getThis()).clust_opts.num_clusters = 0
    #(opts.getThis()).clust_opts.verbose = True
    
    # BUG: none of the below works, Cython no like strings
    #ParseCommandLine(ac, av, opts, args)
    #cdef string temp = ""
    #(opts.getThis()).infile_A = temp 
    #(opts.getThis()).infile_W = ""
    #(opts.getThis()).infile_H = ""
    #(opts.getThis()).dictfile = ""
    #(opts.getThis()).outdir = ""
    #(opts.getThis()).clustfile = ""
    #(opts.getThis()).assignfile = ""
    #(opts.getThis()).show_help = False
    #(opts.getThis()).format = XML

    #opts.getThis().infile_A = string(PyString_AsString(args.matrixfile))

def PyFlatClustSparse(opts, PyDoubleSparseMatrix A, 
                      vector[double]& buf_w, vector[double]& buf_h,
                      unsigned int ldim_w, unsigned int ldim_h, 
                      PyNmfStats stats):
    cdef NmfOptions nmf_opts
    temp = opts["clust_opts"]["nmf_opts"]
    nmf_opts.tol = temp["tol"]
    nmf_opts.algorithm = temp["algorithm"]
    nmf_opts.prog_est_algorithm = temp["prog_est_algorithm"]
    nmf_opts.height = temp["height"]
    nmf_opts.width = temp["width"]
    nmf_opts.k = temp["k"]
    nmf_opts.min_iter = temp["min_iter"]
    nmf_opts.max_iter = temp["max_iter"]
    nmf_opts.tolcount = temp["tolcount"]
    nmf_opts.max_threads = temp["max_threads"]
    nmf_opts.verbose = temp["verbose"]
    nmf_opts.normalize = temp["normalize"]

    #cdef vector[double] buf_w_vec
    #for w in buf_w:
    #    buf_w_vec.push_back(w)
    #cdef vector[double] buf_h_vec
    #for h in buf_h:
    #    buf_h_vec.push_back(h)
    cdef Result res = FlatClustSparse(nmf_opts, A.Height(), A.Width(), A.Size(),
                                      A.LockedColBuffer(),
                                      A.LockedRowBuffer(),
                                      A.LockedDataBuffer(),
                                      &(buf_w[0]), ldim_w,
                                      &(buf_h[0]), ldim_h,
                                      dereference(stats.getThis()))
    #return (res == OK, buf_w_vec, buf_h_vec)
    return res == OK
def PyFlatClust(opts, 
                vector[double]& buf_a, vector[double]& buf_w, vector[double]& buf_h, 
                unsigned int ldim_a, unsigned int ldim_w, unsigned int ldim_h, 
                PyNmfStats stats):
    cdef NmfOptions nmf_opts
    temp = opts["clust_opts"]["nmf_opts"]
    nmf_opts.tol = temp["tol"]
    nmf_opts.algorithm = temp["algorithm"]
    nmf_opts.prog_est_algorithm = temp["prog_est_algorithm"]
    nmf_opts.height = temp["height"]
    nmf_opts.width = temp["width"]
    nmf_opts.k = temp["k"]
    nmf_opts.min_iter = temp["min_iter"]
    nmf_opts.max_iter = temp["max_iter"]
    nmf_opts.tolcount = temp["tolcount"]
    nmf_opts.max_threads = temp["max_threads"]
    nmf_opts.verbose = temp["verbose"]
    nmf_opts.normalize = temp["normalize"]
    
    #cdef vector[double] buf_a_vec = buf_a
    #for a in buf_a:
    #    buf_a_vec.push_back(a)
    #cdef vector[double] buf_w_vec
    #for w in buf_w:
    #    buf_w_vec.push_back(w)
    #cdef vector[double] buf_h_vec
    #for h in buf_h:
    #    buf_h_vec.push_back(h)

    cdef Result res = FlatClust(nmf_opts, 
                                &(buf_a[0]), ldim_a,
                                &(buf_w[0]), ldim_w,
                                &(buf_h[0]), ldim_h,
                                dereference(stats.getThis()))
    return (res == OK, buf_w, buf_h)
    #return res == OK

def PyFlatClustAndWriteResults(opts, string& afile, string& cfile, vector[string]& dictionary, 
                            vector[double]& buf_h, vector[double]& buf_w, 
                            ldim_h, ldim_w, m, n, k, PyNmfStats stats):
    cdef vector[int] assignments
    assignments.resize(n)
    cdef vector[int] term_indices
    term_indices.resize(opts["clust_opts"]["maxterms"]*k)

    #cdef vector[double] buf_h_vec
    #for h in buf_h:
    #    buf_h_vec.push_back(h)
    #cdef vector[double] buf_w_vec
    #for w in buf_w:
    #    buf_w_vec.push_back(w)
    #cdef vector[string] dictionary_vec
    #for d in dictionary:
    #    dictionary_vec.push_back(d)


    import time
    start = time.time()
    ComputeAssignments[double](assignments, &(buf_h[0]), ldim_h, k, n)
    TopTerms(opts["clust_opts"]["maxterms"], &(buf_w[0]), ldim_w, m, k, term_indices)

    FlatClustWriteResults(afile, cfile, 
                         assignments, dictionary, term_indices, 
                         opts["format"], opts["clust_opts"]["maxterms"],
                         n, opts["clust_opts"]["num_clusters"])
    
    end = time.time()
    elapsed = (end-start)*1000000
    elapsed += stats.elapsed_us
    if opts["clust_opts"]["verbose"]:
        print "Elapsed wall clock time: ", elapsed/1000000.0


# workaround: pure python
def PyParseCommandLine2(args):
    opts = {}
    opts["clust_opts"] = {}
    opts["clust_opts"]["nmf_opts"] = {}
    opts["clust_opts"]["nmf_opts"]["height"]       = 0
    opts["clust_opts"]["nmf_opts"]["width"]        = 0
    opts["clust_opts"]["nmf_opts"]["k"]            = args.clusters
    opts["clust_opts"]["nmf_opts"]["min_iter"]     = args.miniter
    opts["clust_opts"]["nmf_opts"]["max_iter"]     = args.maxiter
    opts["clust_opts"]["nmf_opts"]["tol"]          = args.tol
    opts["clust_opts"]["nmf_opts"]["tolcount"]     = 1
    opts["clust_opts"]["nmf_opts"]["verbose"]      = args.verbose  
    opts["clust_opts"]["nmf_opts"]["normalize"]    = True
    if args.algorithm == "HALS":
        opts["clust_opts"]["nmf_opts"]["algorithm"] = HALS 
    elif args.algorithm == "RANK2":
        opts["clust_opts"]["nmf_opts"]["algorithm"] = RANK2
    elif args.algorithm == "BPP":
        opts["clust_opts"]["nmf_opts"]["algorithm"] = BPP

    opts["clust_opts"]["nmf_opts"]["prog_est_algorithm"] = PG_RATIO

    opts["clust_opts"]["maxterms"] = 5
    opts["clust_opts"]["num_clusters"] = args.clusters
    opts["clust_opts"]["verbose"] = True

    opts["infile_A"] = args.matrixfile
    opts["infile_W"] = args.infile_W
    opts["infile_H"] = args.infile_H
    opts["dictfile"] = args.dictfile
    opts["outdir"] = args.outdir
    if not args.clustfile:
        opts["clustfile"] = args.outdir+"clusters_%d."%args.clusters+args.format.lower()
    else:
        opts["clustfile"] = args.clustfile
    
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
    
    if args.algorithm == "RANK2":
        if args.clusters != 2:
            args.clusters = 2    
            opts["clust_opts"]["nmf_opts"]["num_clusters"] = 2
            opts["clust_opts"]["nmf_opts"]["k"]            = 2
    return opts
