// Copyright 2014 Georgia Institute of Technology
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <string>
#include <vector>
#include <thread>
#include <cstdint>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include "smallk.hpp"
#include "timer.hpp"
#include "utils.hpp"
#include "random.hpp"
#include "elemental.hpp"
#include "size_check.hpp"
#include "file_format.hpp"
#include "file_loader.hpp"
#include "run_hier_nmf2.hpp"
#include "delimited_file.hpp"
#include "matrix_generator.hpp"
#include "flat_clust_output.hpp"

namespace smallk
{
using namespace smallk;

using std::cout;
using std::cerr;
using std::endl;

typedef double R;

static Random rng;
static bool matrix_loaded = false, dict_loaded = false;
static bool is_sparse = false;
static SparseMatrix<R> A;
static std::vector<R> buf_a, buf_w, buf_h;
static unsigned int m=0u, n=0u, k=0u, nnz=0u;
static unsigned int ldim_a=0u, ldim_w=0u, ldim_h=0u;

static double nmf_tolerance;
static double hier_nmf2_tolerance;
static unsigned int max_iter;
static unsigned int min_iter;
static unsigned int max_threads;
static unsigned int maxterms;
static unsigned int outprecision;
static OutputFormat clustfile_format;
static std::string outdir;
static std::string matrix_filepath, dict_filepath;

std::vector<std::vector<R> > w_initializers;
std::vector<std::vector<R> > h_initializers;
std::vector<std::string> dictionary;

static NmfOptions   nmf_opts;
static ClustOptions hierclust_opts;  // change to HierclustOptions - TBD

static const std::string DEFAULT_FILENAME_W("w.csv");
static const std::string DEFAULT_FILENAME_H("h.csv");

void PrintNmfOpts(NmfOptions opts);
void PrintHierclustOpts(ClustOptions opts);
bool HierNmf2Internal(bool generate_flat,
                      const unsigned int num_clusters,
                      const std::string& outdir,
                      const OutputFormat clustfile_format);

//-----------------------------------------------------------------------------
void Reset()
{
    min_iter = 5;
    max_iter = 5000;
    nmf_tolerance = 0.005;
    hier_nmf2_tolerance = 0.0001;
    
    // use two threads if the HW cannot provide concurrency info
    unsigned int hw_threads = std::thread::hardware_concurrency();
    if (0 == hw_threads)
        hw_threads = 2;
    max_threads = hw_threads;

    maxterms = 5;
    outprecision = 6;
    clustfile_format = OutputFormat::JSON;

    // default outdir is the current directory
    outdir = std::string("");

    dict_loaded   = false;
    matrix_loaded = false;
    dict_filepath.clear();
    matrix_filepath.clear();

    A.Clear();
    buf_a.clear();
    buf_w.clear();
    buf_h.clear();
    m = n = k = nnz = ldim_a = ldim_w = ldim_h = 0u;
}

//-----------------------------------------------------------------------------
void Initialize(int& argc, char**& argv)
{
    Reset();
    rng.SeedFromTime();
    elem::Initialize(argc, argv);
}

//-----------------------------------------------------------------------------
bool IsInitialized()
{
    return elem::Initialized();
}

//-----------------------------------------------------------------------------
void Finalize()
{
    // any cleanup - TBD

    elem::Finalize();
}

//-----------------------------------------------------------------------------
unsigned int GetMajorVersion()
{
    return SMALLK_MAJOR_VERSION;
}

//-----------------------------------------------------------------------------
unsigned int GetMinorVersion()
{
    return SMALLK_MINOR_VERSION;
}

//-----------------------------------------------------------------------------
unsigned int GetPatchLevel()
{
    return SMALLK_PATCH_LEVEL;
}

//-----------------------------------------------------------------------------
std::string GetVersionString()
{
    std::ostringstream version;
    version << GetMajorVersion() << "." 
            << GetMinorVersion() << "." 
            << GetPatchLevel();
    return version.str();
}

//-----------------------------------------------------------------------------
void SeedRNG(const int seed)
{
    rng.SeedFromInt(seed);
}

//-----------------------------------------------------------------------------
void LoadMatrix(const std::string& filepath)
{
    if (filepath.empty())
        throw std::runtime_error("smallk error (LoadMatrix): matrix filename is invalid.");

    cout << "Loading matrix..." << endl;

    matrix_loaded = false;
    if (IsSparse(filepath))
    {
        if (!LoadSparseMatrix(filepath, A, m, n, nnz))
        {
            matrix_filepath.clear();
            std::ostringstream msg;
            msg << "smallk error (LoadMatrix): load failed for file " << filepath;
            throw std::runtime_error(msg.str());
        }
        
        is_sparse = true;
    }
    else // (IsDense(filepath))
    {
        bool ok = LoadDenseMatrix(filepath, buf_a, m, n);
        if (!ok || (buf_a.size() < m*n))
        {
            matrix_filepath.clear();
            std::ostringstream msg;
            msg << "smallk error (LoadMatrix): load failed for file " << filepath;
            throw std::runtime_error(msg.str());
        }
        
        is_sparse = false;
        ldim_a = m;
    }
    
    matrix_loaded = true;
    matrix_filepath = filepath;
}

//-----------------------------------------------------------------------------
bool IsMatrixLoaded()
{
    return matrix_loaded;
}

//-----------------------------------------------------------------------------
std::string GetOutputDir()
{
    return outdir;
}

//-----------------------------------------------------------------------------
void SetOutputDir(const std::string& output_dir)
{
    std::ostringstream msg;
    msg << "smallk error (SetOutputDir): ";

    // Check to see if this is an absolute or relative path.  Construct the
    // fully-qualified path, then check to see if the specified directory
    // actually exists.
    
    std::string full_path;
    if (!output_dir.empty() && ('/' != output_dir[0]))
    {
        // a relative path was given, so get current dir
        if (!GetCurrentDirectory(full_path))
        {
            msg << "could not determine current directory.";
            throw std::runtime_error(msg.str());
        }

        full_path = EnsureTrailingPathSep(full_path);
    }

    full_path += output_dir;
    if (!DirectoryExists(full_path))
    {
        msg << "the directory ";
        msg << "\"" << full_path << "\"";
        msg << " does not exist.";
        throw std::logic_error(msg.str());
    }

    outdir = EnsureTrailingPathSep(full_path);
}

//-----------------------------------------------------------------------------
double GetNmfTolerance()         {return nmf_tolerance;}
double GetHierNmf2Tolerance()    {return hier_nmf2_tolerance;}
unsigned int GetMaxIter()        {return max_iter;}
unsigned int GetMinIter()        {return min_iter;}
unsigned int GetMaxThreads()     {return max_threads;}
unsigned int GetMaxTerms()       {return maxterms;}

//-----------------------------------------------------------------------------
void SetNmfTolerance(const double tol)
{
    if ( (tol <= 0.0) || (tol >= 1.0))
    {
        std::ostringstream msg;
        msg << "smallk error (SetNmfTolerance): ";
        msg << "require 0.0 < tol < 1.0";
        throw std::logic_error(msg.str());
    }

    nmf_tolerance = tol;
}

//-----------------------------------------------------------------------------
void SetHierNmf2Tolerance(const double tol)
{
    if ( (tol <= 0.0) || (tol >= 1.0))
    {
        std::ostringstream msg;
        msg << "smallk error (SetHierNmf2Tolerance): ";
        msg << "require 0.0 < tol < 1.0";
        throw std::logic_error(msg.str());
    }

    hier_nmf2_tolerance = tol;
}

//-----------------------------------------------------------------------------
void SetMaxIter(const unsigned int max_iterations)
{
    max_iter = max_iterations;
}

//-----------------------------------------------------------------------------
void SetMinIter(const unsigned int min_iterations)
{
    min_iter = min_iterations;
}

//-----------------------------------------------------------------------------
void SetMaxThreads(const unsigned int mt)
{
    unsigned int hw_threads = std::thread::hardware_concurrency();
    max_threads = std::min(mt, hw_threads);
}

//-----------------------------------------------------------------------------
void SetMaxTerms(const unsigned int max_terms)
{
    maxterms = max_terms;
}

//-----------------------------------------------------------------------------
unsigned int GetOutputPrecision()
{
    return outprecision;
}

//-----------------------------------------------------------------------------
void SetOutputPrecision(const unsigned int num_digits)
{
    outprecision = num_digits;
}

//-----------------------------------------------------------------------------
bool Nmf(const unsigned int kval, 
         const Algorithm algorithm,
         const std::string& csv_file_w,
         const std::string& csv_file_h)
{
    if (!matrix_loaded)
        throw std::logic_error("smallk error (NMF): no matrix has been loaded.");

    // Check the sizes of matrix W(m, k) and matrix H(k, n) and make sure 
    // they don't overflow Elemental's default signed int index type.

    if (!SizeCheck<int>(m, kval))
        throw std::logic_error("smallk error (Nmf): mxk matrix W is too large.");
    
    if (!SizeCheck<int>(kval, n))
        throw std::logic_error("smallk error (Nmf): kxn matrix H is too large.");

    k = kval;

    // convert to the 'NmfAlgorithm' type in nmf.hpp
    switch (algorithm)
    {
    case Algorithm::MU:
        nmf_opts.algorithm = NmfAlgorithm::MU;
        break;
    case Algorithm::HALS:
        nmf_opts.algorithm = NmfAlgorithm::HALS;
        break;
    case Algorithm::RANK2:
        nmf_opts.algorithm = NmfAlgorithm::RANK2;
        break;
    case Algorithm::BPP:
        nmf_opts.algorithm = NmfAlgorithm::BPP;
        break;
    default:
        throw std::logic_error("smallk error (NMF): unknown NMF algorithm.");
    }

    // set k == 2 for Rank2 algorithm
    if (NmfAlgorithm::RANK2 == nmf_opts.algorithm)
        k = 2;

    ldim_w = m;
    ldim_h = k;

    if (buf_w.size() < m*k)
        buf_w.resize(m*k);
    if (buf_h.size() < k*n)
        buf_h.resize(k*n);
    
    // initialize matrices W and H
    bool ok;
    unsigned int height_w = m, width_w = k, height_h = k, width_h = n;

    cout << "Initializing matrix W..." << endl;
    if (csv_file_w.empty())
        ok = RandomMatrix(&buf_w[0], ldim_w, m, k, rng);
    else
        ok = LoadDelimitedFile(buf_w, height_w, width_w, csv_file_w);
    if (!ok)
    {
        std::ostringstream msg;
        msg << "smallk error (Nmf): load failed for file ";
        msg << "\"" << csv_file_w << "\"";
        throw std::runtime_error(msg.str());
    }

    if ( (height_w != m) || (width_w != k))
    {
        cerr << "\tdimensions of matrix W are " << height_w
             << " x " << width_w << endl;
        cerr << "\texpected " << m << " x " << k << endl;
        throw std::logic_error("smallk error (Nmf): non-conformant matrix W.");
    }

    cout << "Initializing matrix H..." << endl;
    if (csv_file_h.empty())
        ok = RandomMatrix(&buf_h[0], ldim_h, k, n, rng);
    else
        ok = LoadDelimitedFile(buf_h, height_h, width_h, csv_file_h);

    if (!ok)
    {
        std::ostringstream msg;
        msg << "smallk error (Nmf): load failed for file ";
        msg << "\"" << csv_file_h << "\"";
        throw std::runtime_error(msg.str());
    }
    
    if ( (height_h != k) || (width_h != n))
    {
        cerr << "\tdimensions of matrix H are " << height_h
             << " x " << width_h << endl;
        cerr << "\texpected " << k << " x " << n << endl;
        throw std::logic_error("smallk error (Nmf): non-conformant matrix H.");
    }    

    // The ratio of projected gradient norms doesn't seem to work very well
    // with MU.  We frequently observe a 'leveling off' behavior and the 
    // convergence is even slower than usual.  So for MU use the relative
    // change in the Frobenius norm of W as the stopping criterion, which
    // always seems to behave well, even though it is on shaky theoretical
    // ground.

    if (NmfAlgorithm::MU == nmf_opts.algorithm)
        nmf_opts.prog_est_algorithm = NmfProgressAlgorithm::DELTA_FNORM;
    else
        nmf_opts.prog_est_algorithm = NmfProgressAlgorithm::PG_RATIO;

    nmf_opts.tol         = nmf_tolerance;
    nmf_opts.height      = m;
    nmf_opts.width       = n;
    nmf_opts.k           = k;
    nmf_opts.min_iter    = min_iter;
    nmf_opts.max_iter    = max_iter;
    nmf_opts.tolcount    = 1;
    nmf_opts.max_threads = max_threads;
    nmf_opts.verbose     = true;
    nmf_opts.normalize   = true;

    // display all params to user
    PrintNmfOpts(nmf_opts);

    NmfStats stats;
    NmfResult result;
    if (is_sparse)
    {
        result = NmfSparse(nmf_opts, 
                           A.Height(), A.Width(), A.Size(),
                           A.LockedColBuffer(),
                           A.LockedRowBuffer(),
                           A.LockedDataBuffer(),
                           &buf_w[0], ldim_w,
                           &buf_h[0], ldim_h,
                           stats);
    }
    else
    {
        result = Nmf(nmf_opts,
                     &buf_a[0], ldim_a,
                     &buf_w[0], ldim_w,
                     &buf_h[0], ldim_h,
                     stats);
    }

    cout << "Elapsed wall clock time: ";
    cout << ElapsedTime(stats.elapsed_us) << endl;
    cout << endl;

    if (NmfResult::OK != result)
        throw std::runtime_error("smallk error (Nmf): NMF solver failure.");

    // write the computed W and H factors to disk

    std::string outfile_w, outfile_h;
    if (outdir.empty())
    {
        outfile_w = DEFAULT_FILENAME_W;
        outfile_h = DEFAULT_FILENAME_H;
    }
    else
    {
        outfile_w = outdir + DEFAULT_FILENAME_W;
        outfile_h = outdir + DEFAULT_FILENAME_H;
    }

    cout << "Writing output files..." << endl;
    
    if (!WriteDelimitedFile(&buf_w[0], ldim_w, m, k, outfile_w, outprecision))
        throw std::runtime_error("smallk error (Nmf): could not write W result.");
    
    if (!WriteDelimitedFile(&buf_h[0], ldim_h, k, n, outfile_h, outprecision))
        throw std::runtime_error("smallk error (Nmf): could not write H result.");

    return true;
}

//-----------------------------------------------------------------------------
const double* LockedBufferW(unsigned int& ldim, 
                            unsigned int& height,
                            unsigned int& width)
{
    ldim = m;
    height = m;
    width = k;
    return &buf_w[0];
}

//-----------------------------------------------------------------------------
const double* LockedBufferH(unsigned int& ldim,
                            unsigned int& height,
                            unsigned int& width)
{
    ldim = k;
    height = k;
    width = n;
    return &buf_h[0];
}

//-----------------------------------------------------------------------------
void LoadDictionary(const std::string& filepath)
{
    cout << "loading dictionary..." << endl;
  
    dictionary.clear();
    dict_loaded = false;
    if (!LoadStringsFromFile(filepath, dictionary))
    {
        dict_filepath.clear();
        std::ostringstream msg;
        msg << "smallk error (LoadDictionary): load failed for file " << filepath;
        throw std::runtime_error(msg.str());
    }

    dict_filepath = filepath;
    dict_loaded = true;
}

//-----------------------------------------------------------------------------
OutputFormat GetOutputFormat()
{
    return clustfile_format;
}

//-----------------------------------------------------------------------------
void SetOutputFormat(const OutputFormat format)
{
    clustfile_format = format;
}

//-----------------------------------------------------------------------------
void HierNmf2Init(const unsigned int num_clusters)
{
    // The hierclust code requires 2*num_clusters separate W and H init 
    // matrices, each of which has size mx2 (W) or 2xn (H).  Make sure that
    // these sizes do not overflow Elemental's default index type (int).

    if (!SizeCheck<int>(m, 2))
        throw std::logic_error("smallk error (HierNmf2): matrix height too large.");
    if (!SizeCheck<int>(2, n))
        throw std::logic_error("smallk error (HierNmf2): matrix width too large.");

    // setup random initializers
    unsigned int num_initializers = 2*num_clusters;

    if (w_initializers.size() < num_initializers)
        w_initializers.resize(num_initializers);
    if (h_initializers.size() < num_initializers)
        h_initializers.resize(num_initializers);

    cout << "creating random W initializers..." << endl;

    unsigned int required_size = m*2;
    for (unsigned int i=0; i<num_initializers; ++i)
    {
        w_initializers[i].resize(required_size);
        RandomMatrix(w_initializers[i], m, 2, rng);
    }

    cout << "creating random H initializers..." << endl;
    
    // no initializer file, so use random init
    required_size = 2*n;
    for (unsigned int i=0; i<num_initializers; ++i)
    {
        h_initializers[i].resize(required_size);
        RandomMatrix(h_initializers[i], 2, n, rng);
    }
}

//-----------------------------------------------------------------------------
void HierNmf2Internal(bool generate_flat,
                      const unsigned int num_clusters,
                      const OutputFormat clustfile_format)
{
    if (!matrix_loaded)
        throw std::logic_error("smallk error (HierNmf2): no matrix has been loaded.");

    if (!dict_loaded)
        throw std::logic_error("smallk error (HierNmf2): no dictionary has been loaded.");

    HierNmf2Init(num_clusters);

    std::string output_dir = EnsureTrailingPathSep(outdir);

    // set options
    hierclust_opts.nmf_opts.tol                = hier_nmf2_tolerance;
    hierclust_opts.nmf_opts.algorithm          = NmfAlgorithm::RANK2;
    hierclust_opts.nmf_opts.prog_est_algorithm = NmfProgressAlgorithm::PG_RATIO;
    hierclust_opts.nmf_opts.height             = m;
    hierclust_opts.nmf_opts.width              = n;
    hierclust_opts.nmf_opts.k                  = 2;
    hierclust_opts.nmf_opts.min_iter           = min_iter;
    hierclust_opts.nmf_opts.max_iter           = max_iter;
    hierclust_opts.nmf_opts.tolcount           = 1;
    hierclust_opts.nmf_opts.max_threads        = max_threads;
    hierclust_opts.nmf_opts.verbose            = false;
    hierclust_opts.nmf_opts.normalize          = true;
    hierclust_opts.maxterms                    = maxterms;
    hierclust_opts.unbalanced                  = 0.1;
    hierclust_opts.trial_allowance             = 3;
    hierclust_opts.num_clusters                = num_clusters;
    hierclust_opts.verbose                     = true;
    hierclust_opts.flat                        = generate_flat;

    // convert to internal FileFormat enum
    FileFormat format = FileFormat::JSON;
    if (OutputFormat::XML == clustfile_format)
        format = FileFormat::XML;

    // the assignments file is always in CSV format
    std::ostringstream assign_fname;
    assign_fname << "assignments_" << num_clusters;
    std::string hier_assignfile = 
        output_dir + AppendExtension(assign_fname.str(), FileFormat::CSV);

    // the tree file (clustfile) is in the user-specified format
    std::ostringstream tree_fname;
    tree_fname << "tree_" << num_clusters;
    std::string hier_treefile = 
        output_dir + AppendExtension(tree_fname.str(), format);

    // init complete, so print all options
    PrintHierclustOpts(hierclust_opts);
    
    // run the hierarchical clustering code

    Tree tree;
    ClustStats stats;
    std::vector<int> assignments, assignments_flat;
    std::vector<int> term_indices(maxterms * num_clusters);

    Timer timer;
    timer.Start();

    bool ok = RunHierNmf2(m, n, A, buf_a, 
                          w_initializers, 
                          h_initializers,
                          assignments, 
                          assignments_flat, 
                          term_indices, tree,
                          stats, hierclust_opts);

    timer.Stop();
    double elapsed = timer.ReportMilliseconds();

    cout << "\nElapsed wall clock time: ";
    if (elapsed < 1000.0)
        cout << elapsed << " ms." << endl;
    else
        cout << elapsed*0.001 << " s." << endl;

    if (!ok)
    {
        throw std::runtime_error("smallk error (HierNMF2): HierNMF2 fatal error.");
    }
    else
    {
        int num_converged = stats.nmf_count - stats.max_count;
        cout << num_converged << "/" << stats.nmf_count << " factorizations"
             << " converged." << endl << endl;

        // write results

        if (hierclust_opts.verbose)
            cout << "Writing output files..." << endl;

        if (!WriteAssignmentsFile(assignments, hier_assignfile))
            cerr << "\terror writing assignments file" << endl;

        if (!tree.Write(hier_treefile, format, dictionary))
            cerr << "\terror writing hierarchical results file" << endl;

        if (hierclust_opts.flat)
        {
            FlatClustWriteResults(output_dir,
                                  assignments_flat,
                                  dictionary, term_indices,
                                  format,
                                  hierclust_opts.maxterms, n,
                                  hierclust_opts.num_clusters);
        }
    }
}

//-----------------------------------------------------------------------------
void HierNmf2(const unsigned int num_clusters)
{
    HierNmf2Internal(false, num_clusters, clustfile_format);
}

//-----------------------------------------------------------------------------
void HierNmf2WithFlat(const unsigned int num_clusters)
{
    HierNmf2Internal(true, num_clusters, clustfile_format);
}

//-----------------------------------------------------------------------------
void PrintNmfOpts(NmfOptions opts)
{
    cout << "\n                parameters: \n" << endl;
    cout << "\t         algorithm: ";
    switch (opts.algorithm)
    {
    case NmfAlgorithm::MU:
         cout << "Multiplicative Updating";
         break;
    case NmfAlgorithm::HALS:
        cout << "HALS";
        break;
    case NmfAlgorithm::RANK2:
        cout << "Rank 2";
        break;
    case NmfAlgorithm::BPP:
        cout << "Nonnegative Least Squares with Block Principal Pivoting";
        break;
    default:
        cout << "*** UNKNOWN NMF ALGORITHM ***";
        break;
    }
    cout << endl;

    cout << "\tstopping criterion: ";
    switch (opts.prog_est_algorithm)
    {
    case NmfProgressAlgorithm::PG_RATIO:
        cout << "Ratio of Projected Gradients";
        break;
    case NmfProgressAlgorithm::DELTA_FNORM:
        cout << "Relative Change in the F-norm of W";
        break;
    default:
        cout << "*** UNKNOWN PROGRESS ESTIMATION ALGORITHM ***";
        break;
    }
    cout << endl;

    cout << "\t            height: " << opts.height << endl;
    cout << "\t             width: " << opts.width << endl;
    cout << "\t                 k: " << opts.k << endl;
    cout << "\t           miniter: " << opts.min_iter << endl;
    cout << "\t           maxiter: " << opts.max_iter << endl;
    cout << "\t               tol: " << opts.tol << endl;
    cout << "\t        matrixfile: " << matrix_filepath << endl;
    cout << "\t        maxthreads: " << opts.max_threads << endl;
    cout << endl;
}

//-----------------------------------------------------------------------------
void PrintHierclustOpts(ClustOptions opts)
{
    cout << "\n\t        parameters: \n" << endl;
    cout << "\t            height: " << opts.nmf_opts.height << endl;
    cout << "\t             width: " << opts.nmf_opts.width << endl;
    cout << "\t        matrixfile: " << matrix_filepath << endl;
    cout << "\t          dictfile: " << dict_filepath << endl;
//    cout << "\t          treefile: " 
    cout << "\t               tol: " << opts.nmf_opts.tol << endl;
    cout << "\t           miniter: " << opts.nmf_opts.min_iter << endl;
    cout << "\t           maxiter: " << opts.nmf_opts.max_iter << endl;
    cout << "\t          maxterms: " << opts.maxterms << endl;
    cout << "\t        maxthreads: " << opts.nmf_opts.max_threads << endl;
}

} // namespace smallk
