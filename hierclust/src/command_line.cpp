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
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <getopt.h>
#include "nmf.hpp"
#include "clust.hpp"
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
    { "trial_allowance", required_argument,   NULL,    'i' },
    { "unbalanced",      required_argument,   NULL,    'j' },
    { "verbose",         required_argument,   NULL,    'k' },
    { "maxthreads",      required_argument,   NULL,    'l' },
    { "maxterms",        required_argument,   NULL,    'm' },
    { "initdir",         required_argument,   NULL,    'n' },
    { "treefile",        required_argument,   NULL,    'q' },
    { "assignfile",      required_argument,   NULL,    'r' },
    { "flat",            required_argument,   NULL,    's' },
    { "format",          required_argument,   NULL,    't' },
    { 0, 0, 0, 0}
};

//-----------------------------------------------------------------------------
void PrintOpts(const CommandLineOptions& opts)
{
    cout << "\n     Command line options: \n" << endl;

    cout << "\t            height: " << opts.clust_opts.nmf_opts.height << endl;
    cout << "\t             width: " << opts.clust_opts.nmf_opts.width << endl;
    cout << "\t        matrixfile: " << opts.infile_A << endl;
    cout << "\t           initdir: " << opts.clust_opts.initdir << endl;
    cout << "\t          dictfile: " << opts.dictfile << endl;
    cout << "\t        assignfile: " << opts.assignfile << endl;

    std::string format = STRING_XML;
    if (FileFormat::JSON == opts.format)
        format = STRING_JSON;
    cout << "\t            format: " << format << endl;
    cout << "\t          treefile: " << opts.treefile << endl;
    cout << "\t          clusters: " << opts.clust_opts.num_clusters << endl;
    cout << "\t               tol: " << opts.clust_opts.nmf_opts.tol << endl;
    cout << "\t            outdir: " << opts.outdir << endl;
    cout << "\t           miniter: " << opts.clust_opts.nmf_opts.min_iter << endl;
    cout << "\t           maxiter: " << opts.clust_opts.nmf_opts.max_iter << endl;
    cout << "\t          maxterms: " << opts.clust_opts.maxterms << endl;
    cout << "\t        maxthreads: " << opts.clust_opts.nmf_opts.max_threads << endl;
    cout << "\t        unbalanced: " << opts.clust_opts.unbalanced << endl;
    cout << "\t   trial_allowance: " << opts.clust_opts.trial_allowance << endl;
    cout << "\t              flat: " << opts.clust_opts.flat << endl;
    cout << "\t           verbose: " << opts.clust_opts.verbose << endl;
    cout << endl;
}

//-----------------------------------------------------------------------------
void ShowHelp(const std::string& program_name)
{
    cout << endl;
    cout << "Usage: " << program_name << endl;
    cout << "        --matrixfile <filename>     Filename of the matrix to be factored." << endl;
    cout << "                                    Either CSV format for dense or MatrixMarket format ";
    cout << "for sparse." << endl;
    cout << "        --dictfile <filename>       The name of the dictionary file." << endl;
    cout << "        --clusters <integer>        The number of clusters to generate." << endl;
    cout << "        [--initdir  (empty)]        Directory of initializers for all Rank2 factorizations." << endl;
    cout << "                                    If unspecified, random init will be used. " << endl;    
    cout << "        [--tol  0.0001]             Tolerance value for each factorization. " << endl;
    cout << "        [--outdir  (empty)]         Output directory.  If unspecified, results will be " << endl;
    cout << "                                    written to the current directory." << endl;
    cout << "        [--miniter  5]              Minimum number of iterations to perform." << endl;
    cout << "        [--maxiter  5000]           Maximum number of  iterations to perform. " << endl;
    cout << "        [--maxterms  5]             Number of terms per node. " << endl;

    unsigned int hw_threads = GetMaxThreadCount();
    cout << "        [--maxthreads  " << setw(3) << hw_threads << "]         Upper limit to thread count. " 
         << endl;

    cout << "        [--unbalanced  0.1]         Threshold for determining leaf node imbalance. " << endl;
    cout << "        [--trial_allowance  3]      Number of split attempts. " << endl;
    cout << "        [--flat  0]                 Whether to generate a flat clustering result. " << endl;
    cout << "                                        1 == yes, 0 == no" << endl;
    cout << "        [--verbose  1]              Whether to print updates to the screen." << endl;
    cout << "                                        1 == yes, 0 == no" << endl;
    cout << "        [--format  XML]             Format of the output file containing the tree." << endl;
    cout << "                                        XML: XML format" << endl;
    cout << "                                        JSON: JavaScript Object Notation" << endl;
    cout << "        [--treefile  tree_N.ext]    Name of the output file containing the tree." << endl;
    cout << "                                    N is the number of clusters for this run." << endl;
    cout << "                                    The string 'ext' depends on the desired format." << endl;
    cout << "                                    This filename is relative to the outdir." << endl;
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
    opts.clust_opts.nmf_opts.verbose      = false;  // nmf is silent
    opts.clust_opts.nmf_opts.normalize    = false;  // rank2 normalizes on each iter
    opts.clust_opts.nmf_opts.algorithm    = NmfAlgorithm::RANK2;
    opts.clust_opts.initdir               = std::string("");
    
    // stopping criterion is always the ratio of projected gradients
    opts.clust_opts.nmf_opts.prog_est_algorithm = 
        NmfProgressAlgorithm::PG_RATIO;

    // max_threads will be set below

    // set clust_opts defaults
    opts.clust_opts.maxterms              = 5;
    opts.clust_opts.trial_allowance       = 3;
    opts.clust_opts.unbalanced            = 0.1;
    opts.clust_opts.num_clusters          = 0;
    opts.clust_opts.verbose               = true;
    opts.clust_opts.flat                  = false;
    
    // set command line opts defaults
    opts.infile_A    = std::string("");
    opts.dictfile    = std::string("");
    opts.outdir      = std::string("");
    opts.treefile    = std::string("");
    opts.assignfile  = std::string("");
    opts.show_help   = false;
    opts.format      = FileFormat::XML;

    char c;
    int index;
    while (-1 != (c = getopt_long(argc, argv, 
                                  ":a:b:c:d:e:f:g:h:i:j:k:l:m:n:p:q:r:s:t",
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
        case 'i':  // trial_allowance
            opts.clust_opts.trial_allowance = atoi(optarg);
            break;
        case 'j':  // unbalanced
            opts.clust_opts.unbalanced = atof(optarg);
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
        case 'n':  // initdir
            opts.clust_opts.initdir = std::string(optarg);
            break;
        case 'q':  // treefile
            opts.treefile = std::string(optarg);
            break;
        case 'r':  // assignfile
            opts.assignfile = std::string(optarg);
            break;
        case 's':  // flat
            opts.clust_opts.flat = (0 != atoi(optarg));
            break;
        case 't':  // format
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

    // adjust thread count for user specifications
    if (user_max_threads <= 0)
        user_max_threads = hw_max_threads;
    opts.clust_opts.nmf_opts.max_threads = std::min(user_max_threads, hw_max_threads);

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

    // add a trailing path separator to the initdir, if any
    if (!opts.clust_opts.initdir.empty())
    {
        std::string initdir = EnsureTrailingPathSep(opts.clust_opts.initdir);
        opts.clust_opts.initdir = initdir;
    }
    
    std::string output_dir = EnsureTrailingPathSep(opts.outdir);
    
    // construct the name of the assignment file if unspecified by user
    std::string filename = opts.assignfile;
    if (filename.empty())
    {
        std::ostringstream fname;
        fname << "assignments_" << num_clusters;
        filename = AppendExtension(fname.str(), FileFormat::CSV);
    }
    
    opts.assignfile = output_dir + filename;
    
    // construct the name of the tree file if unspecified by user
    filename = opts.treefile;
    if (opts.treefile.empty())
    {
        std::ostringstream fname;
        fname << "tree_" << num_clusters;
        filename = AppendExtension(fname.str(), opts.format);
    }
    
    opts.treefile = output_dir + filename;
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

    // specified initdir must exist
    if (!opts.clust_opts.initdir.empty() && !DirectoryExists(opts.clust_opts.initdir))
    {
        cerr << "the specified init directory \"";
        cerr << opts.clust_opts.initdir << "\"";
        cerr << " does not exist" << endl;
        return false;
    }
    
    // validate opts.clust_opts, but ignore the matrix checks since the 
    // matrices have yet to be loaded
    if (!IsValid(opts.clust_opts, false))
        return false;

    return true;
}
