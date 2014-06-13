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

#include <thread>
#include <string>
#include <limits>
#include <cassert>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <getopt.h>
#include "utils.hpp"
#include "constants.hpp"
#include "file_format.hpp"
#include "command_line.hpp"
#include "thread_utils.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;

static option longopts[] = 
{
    { "matrixfile",      required_argument,   NULL,    'a' },
    { "dictfile",        required_argument,   NULL,    'b' },
    { "clusters",        required_argument,   NULL,    'c' },
    { "tol",             required_argument,   NULL,    'd' },
    { "outdir",          required_argument,   NULL,    'e' },
    { "miniter",         required_argument,   NULL,    'f' },
    { "maxiter",         required_argument,   NULL,    'g' },
    { "help",            no_argument,         NULL,    'h' },
    { "algorithm",       required_argument,   NULL,    'i' },
    { "verbose",         required_argument,   NULL,    'k' },
    { "maxthreads",      required_argument,   NULL,    'l' },
    { "maxterms",        required_argument,   NULL,    'm' },
    { "infile_W",        required_argument,   NULL,    'n' },
    { "infile_H",        required_argument,   NULL,    'o' },
//    { "stopping",        required_argument,   NULL,    'p' },
    { "clustfile",       required_argument,   NULL,    'q' },
    { "assignfile",      required_argument,   NULL,    'r' },
    { "format",          required_argument,   NULL,    's' },
    { 0, 0, 0, 0}
};

//-----------------------------------------------------------------------------
void PrintOpts(const CommandLineOptions& opts)
{
    cout << "\n     Command line options: \n" << endl;

    cout << "\t            height: " << opts.clust_opts.nmf_opts.height << endl;
    cout << "\t             width: " << opts.clust_opts.nmf_opts.width << endl;
    cout << "\t        matrixfile: " << opts.infile_A << endl;
    cout << "\t          infile_W: " << opts.infile_W << endl;
    cout << "\t          infile_H: " << opts.infile_H << endl;
    cout << "\t          dictfile: " << opts.dictfile << endl;
    cout << "\t        assignfile: " << opts.assignfile << endl;

    std::string format = STRING_XML;
    if (FileFormat::JSON == opts.format)
        format = STRING_JSON;
    cout << "\t            format: " << format << endl;

    cout << "\t         clustfile: " << opts.clustfile << endl;
    cout << "\t         algorithm: ";
    switch (opts.clust_opts.nmf_opts.algorithm)
    {
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

    cout << "\t          clusters: " << opts.clust_opts.num_clusters << endl;
    cout << "\t               tol: " << opts.clust_opts.nmf_opts.tol << endl;
    cout << "\t            outdir: " << opts.outdir << endl;
    cout << "\t           miniter: " << opts.clust_opts.nmf_opts.min_iter << endl;
    cout << "\t           maxiter: " << opts.clust_opts.nmf_opts.max_iter << endl;
    cout << "\t          maxterms: " << opts.clust_opts.maxterms << endl;
    cout << "\t        maxthreads: " << opts.clust_opts.nmf_opts.max_threads << endl;
    cout << "\t           verbose: " << opts.clust_opts.verbose << endl;
    cout << endl;
}

//-----------------------------------------------------------------------------
void ShowHelp(const std::string& program_name)
{
    cout << endl;
    cout << "Usage: " << program_name << endl;
    cout << "        --matrixfile <filename>      Filename of the matrix to be factored." << endl;
    cout << "                                     Either CSV format for dense or MatrixMarket format ";
    cout << "for sparse." << endl;
    cout << "        --dictfile <filename>        The name of the dictionary file." << endl;
    cout << "        --clusters <integer>         The number of clusters to generate." << endl;
    cout << "        [--algorithm  BPP]           The NMF algorithm to use: " << endl;
    cout << "                                         HALS:  hierarchical alternating least squares" << endl;
    cout << "                                         RANK2: rank2 with optimal active set selection" << endl;
    cout << "                                                (for two clusters only)" << endl;
    cout << "                                         BPP:   block principal pivoting" << endl;
    cout << "        [--infile_W  (empty)]        Dense matrix to initialize W, CSV file." << endl;
    cout << "                                     The matrix has m rows and 'clusters' columns." << endl;
    cout << "                                     If unspecified, W will be randomly initialized." << endl;
    cout << "        [--infile_H  (empty)]        Dense matrix to initialize H, CSV file. " << endl;
    cout << "                                     The matrix has 'clusters' rows and n columns." << endl;
    cout << "                                     If unspecified, H will be randomly initialized. " << endl;    
    cout << "        [--tol  0.0001]              Tolerance value for the progress metric. " << endl;
    cout << "        [--outdir  (empty)]          Output directory.  If unspecified, results will be " << endl;
    cout << "                                     written to the current directory." << endl;
    cout << "        [--miniter  5]               Minimum number of iterations to perform." << endl;
    cout << "        [--maxiter  5000]            Maximum number of  iterations to perform. " << endl;
    cout << "        [--maxterms  5]              Number of terms per node. " << endl;
 
    unsigned int hw_threads = GetMaxThreadCount();
    cout << "        [--maxthreads " << setw(3) << hw_threads << "]           Upper limit to thread count. " 
         << endl;

    cout << "        [--verbose  1]               Whether to print updates to the screen." << endl;
    cout << "                                         1 == yes, 0 == no" << endl;
    cout << "        [--format  XML]              Format of the output file containing the tree." << endl;
    cout << "                                         XML: XML format" << endl;
    cout << "                                         JSON: JavaScript Object Notation" << endl;
    cout << "        [--clustfile clusters_N.ext] Name of the output XML file containing the tree." << endl;
    cout << "                                     N is the number of clusters for this run." << endl;
    cout << "                                     The string 'ext' depends on the desired format." << endl;
    cout << "                                     This filename is relative to the outdir." << endl;
    cout << "        [--assignfile assignments_N.csv]  Name of the file containing final assignments." << endl;
    cout << "                                          N is the number of clusters for this run." << endl;
    cout << "                                          This filename is relative to the outdir." << endl;
    cout << endl;
}

//-----------------------------------------------------------------------------
bool ParseCommandLine(int argc, char* argv[], CommandLineOptions& opts)
{
    std::string tmp;
    int user_max_threads = -1;

    // set nmf_opts defaults
    opts.clust_opts.nmf_opts.height       = 0;
    opts.clust_opts.nmf_opts.width        = 0;
    opts.clust_opts.nmf_opts.k            = 0;
    opts.clust_opts.nmf_opts.min_iter     = 5;
    opts.clust_opts.nmf_opts.max_iter     = 5000;
    opts.clust_opts.nmf_opts.tol          = 0.0001;
    opts.clust_opts.nmf_opts.tolcount     = 1;
    opts.clust_opts.nmf_opts.verbose      = true;  // print NMF updates
    opts.clust_opts.nmf_opts.normalize    = true;
    opts.clust_opts.nmf_opts.algorithm    = NmfAlgorithm::BPP;

    // this is mandatory
    opts.clust_opts.nmf_opts.prog_est_algorithm = 
        NmfProgressAlgorithm::PG_RATIO;

    // max_threads will be set below

    // set clust_opts defaults
    opts.clust_opts.maxterms              = 5;
    opts.clust_opts.num_clusters          = 0;
    opts.clust_opts.verbose               = true;
    
    // set command line opts defaults
    opts.infile_A    = std::string("");
    opts.infile_W    = std::string("");
    opts.infile_H    = std::string("");
    opts.dictfile    = std::string("");
    opts.outdir      = std::string("");
    opts.clustfile    = std::string("");
    opts.assignfile  = std::string("");
    opts.show_help   = false;
    opts.format      = FileFormat::XML;

    char c;
    int index;
    while (-1 != (c = getopt_long(argc, argv, 
                                  ":a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s",
                                  longopts, &index)))
    {
        switch (c)
        {
        case 'a':  // matrixfile
            opts.infile_A = std::string(optarg);
            break;
        case 'b':  // dictfile
            opts.dictfile = std::string(optarg);
            break;
        case 'c':  // clusters
            opts.clust_opts.num_clusters = atoi(optarg);
            opts.clust_opts.nmf_opts.k   = atoi(optarg);
            break;
        case 'd':  // tol
            opts.clust_opts.nmf_opts.tol = atof(optarg);
            break;
        case 'e':  // outdir
            opts.outdir = std::string(optarg);
            break;
        case 'f':  // miniter
            opts.clust_opts.nmf_opts.min_iter = atoi(optarg);
            break;
        case 'g':  // maxiter
            opts.clust_opts.nmf_opts.max_iter = atoi(optarg);
            break;
        case 'h':  // help
            opts.show_help = true;
            break;
        case 'i':  // algorithm
            tmp = std::string(optarg);
            ToUpper(tmp);
            if (STRING_HALS == tmp)
                opts.clust_opts.nmf_opts.algorithm = NmfAlgorithm::HALS;
            else if (STRING_RANK2 == tmp)
                opts.clust_opts.nmf_opts.algorithm = NmfAlgorithm::RANK2;
            else if (STRING_BPP == tmp)
                opts.clust_opts.nmf_opts.algorithm = NmfAlgorithm::BPP;
            else
                InvalidValue(tmp);
            break;
        case 'k':  // verbose
            opts.clust_opts.verbose = (0 != atoi(optarg));
            break;
        case 'l':  // maxthreads
            user_max_threads = atoi(optarg);
            break;
        case 'm':  // maxterms
            opts.clust_opts.maxterms = atoi(optarg);
            break;
        case 'n':  // infile_W
            opts.infile_W = std::string(optarg);
            break;
        case 'o':  // infile_H
            opts.infile_H = std::string(optarg);
            break;
        case 'q':  // clustfile
            opts.clustfile = std::string(optarg);
            break;
        case 'r':  // assignfile
            opts.assignfile = std::string(optarg);
            break;
        case 's':  // format
            tmp = std::string(optarg);
            ToUpper(tmp);
            if (STRING_XML == tmp)
                opts.format = FileFormat::XML;
            else if (STRING_JSON == tmp)
                opts.format = FileFormat::JSON;
            else
                InvalidValue(tmp);
            break;
        case ':':  // missing option argument
            assert(optind >= 1);
            assert(optind <= argc);
            cerr << "missing argument for option " << argv[optind-1] << endl;
            return false;
            break;
        case '?':  // invalid option
        default:
            assert(optind >= 1);
            assert(optind <= argc);
            cerr << "invalid option: " << argv[optind-1] << endl;
            return false;
            break;
        }
    }

    // if no command line arg, user wants help
    if (1 == argc)
        opts.show_help = true;

    // found --help on the command line
    if (opts.show_help)
        return false;

    // adjust the thread count to not exceed HW capabilities
    int hw_max_threads = GetMaxThreadCount();

    if (user_max_threads <= 0)
        user_max_threads = hw_max_threads;
    opts.clust_opts.nmf_opts.max_threads = std::min(user_max_threads, hw_max_threads);

    // set verbosity
    if (!opts.clust_opts.verbose)
        opts.clust_opts.nmf_opts.verbose = false;

    // check for required options
    if (opts.infile_A.empty())
    {
        cerr << "required command line argument --matrixfile not found" << endl;
        return false;
    }

    if (opts.dictfile.empty())
    {
        cerr << "required command line argument --dictfile not found" << endl;
        return false;
    }

    if (0 == opts.clust_opts.num_clusters)
    {
        cerr << "required command line argument --clusters not found" << endl;
        return false;
    }

    unsigned int num_clusters = opts.clust_opts.num_clusters;
    std::string output_dir = EnsureTrailingPathSep(opts.outdir);
    
    // construct the name of the assignment file if empty
    std::string filename = opts.assignfile;
    if (filename.empty())
    {
        std::ostringstream fname;
        fname << "assignments_" << num_clusters;
        filename = AppendExtension(fname.str(), FileFormat::CSV);
    }
    
    opts.assignfile = output_dir + filename;
    
    // construct the name of the result file if empty
    filename = opts.clustfile;
    if (opts.clustfile.empty())
    {
        std::ostringstream fname;
        fname << "clusters_" << num_clusters;
        filename = AppendExtension(fname.str(), opts.format);
    }
    
    opts.clustfile = output_dir + filename;

    // force two clusters if RANK2 algorithm
    if (NmfAlgorithm::RANK2 == opts.clust_opts.nmf_opts.algorithm)
    {
        if (2 != opts.clust_opts.num_clusters)
        {
            if (opts.clust_opts.verbose)
                cout << "warning: forcing num_clusters=2 for RANK2 algorithm" << endl;

            opts.clust_opts.num_clusters = 2;
            opts.clust_opts.nmf_opts.k   = 2;
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
bool IsValid(const CommandLineOptions& opts)
{
    // a nonempty output directory must actually exist
    if (!opts.outdir.empty() && !DirectoryExists(opts.outdir))
    {
        cerr << "the specified output directory \"";
        cerr << opts.outdir << "\"";
        cerr << " does not exist" << endl;
        return false;
    }

    // number of clusters must be > 0
    if (opts.clust_opts.num_clusters <= 0)
    {
        cerr << "value for --clusters must be a positive integer" << endl;
        return false;
    }

    // 0 < tolerance < 1.0
    if ( (opts.clust_opts.nmf_opts.tol <= 0.0) || (opts.clust_opts.nmf_opts.tol >= 1.0))
    {
        cerr << "tolerance must be in the interval (0.0, 1.0)" << endl;
        return false;
    }

    // min iterations > 0
    if (opts.clust_opts.nmf_opts.min_iter <= 0)
    {
        cerr << "miniter must be a positive integer" << endl;
        return false;
    }

    // max iterations > 0
    if (opts.clust_opts.nmf_opts.max_iter <= 0)
    {
        cerr << "maxiter must be a positive integer" << endl;
        return false;
    }

    // maxterms > 0
    if (opts.clust_opts.maxterms <= 0)
    {
        cerr << "maxterms must be a positive integer" << endl;
        return false;
    }

    if ( (NmfProgressAlgorithm::PG_RATIO != opts.clust_opts.nmf_opts.prog_est_algorithm) &&
         (NmfProgressAlgorithm::DELTA_FNORM != opts.clust_opts.nmf_opts.prog_est_algorithm))
    {
        cerr << "clustlib error: unknown stopping criterion " << endl;
        return false;
    }

    assert(opts.clust_opts.num_clusters == opts.clust_opts.nmf_opts.k);

    return true;
}
