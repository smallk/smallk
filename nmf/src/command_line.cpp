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

#include <limits>
#include <thread>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <getopt.h>
#include "nmf.hpp"
#include "utils.hpp"
#include "constants.hpp"
#include "command_line.hpp"
#include "thread_utils.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;

static option longopts[] = 
{
    { "matrixfile",     required_argument,   NULL,    'a' },
    { "k",              required_argument,   NULL,    'b' },
    { "algorithm",      required_argument,   NULL,    'c' },
    { "stopping",       required_argument,   NULL,    'd' },
    { "tol",            required_argument,   NULL,    'e' },
    { "tolcount",       required_argument,   NULL,    'f' },
    { "infile_W",       required_argument,   NULL,    'g' },
    { "infile_H",       required_argument,   NULL,    'h' },
    { "outfile_W",      required_argument,   NULL,    'i' },
    { "outfile_H",      required_argument,   NULL,    'j' },
    { "miniter",        required_argument,   NULL,    'k' },
    { "maxiter",        required_argument,   NULL,    'l' },
    { "outprecision",   required_argument,   NULL,    'm' },
    { "maxthreads",     required_argument,   NULL,    'n' },
    { "normalize",      required_argument,   NULL,    'o' },
    { "verbose",        required_argument,   NULL,    'p' },
    { "help",           no_argument,         NULL,    'q' },
    { 0, 0, 0, 0}
};


//-----------------------------------------------------------------------------
void PrintOpts(const CommandLineOptions& opts)
{
    cout << "\n      Command line options: \n" << endl;
    cout << "\t         algorithm: ";
    switch (opts.nmf_opts.algorithm)
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
    switch (opts.nmf_opts.prog_est_algorithm)
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

    cout << "\t            height: " << opts.nmf_opts.height << endl;
    cout << "\t             width: " << opts.nmf_opts.width << endl;
    cout << "\t                 k: " << opts.nmf_opts.k << endl;
    cout << "\t           miniter: " << opts.nmf_opts.min_iter << endl;
    cout << "\t           maxiter: " << opts.nmf_opts.max_iter << endl;
    cout << "\t               tol: " << opts.nmf_opts.tol << endl;
    cout << "\t          tolcount: " << opts.nmf_opts.tolcount << endl;
    cout << "\t           verbose: " << opts.nmf_opts.verbose << endl;
    cout << "\t         normalize: " << opts.nmf_opts.normalize << endl;
    cout << "\t      outprecision: " << opts.output_precision << endl;
    cout << "\t        matrixfile: " << opts.infile_A << endl;

    std::string inW = opts.infile_W;
    if (inW.empty())
        inW = STRING_RANDOM;
    cout << "\t          infile_W: " << opts.infile_W << endl;

    std::string inH = opts.infile_H;
    if (inH.empty())
        inH = STRING_RANDOM;
    cout << "\t          infile_H: " << opts.infile_H << endl;

    cout << "\t         outfile_W: " << opts.outfile_W << endl;
    cout << "\t         outfile_H: " << opts.outfile_H << endl;
    cout << "\t        maxthreads: " << opts.nmf_opts.max_threads << endl;
    cout << endl;
}

//-----------------------------------------------------------------------------
void ShowHelp(const std::string& program_name)
{
    cout << endl;
    cout << "Usage: " << program_name << endl;
    cout << "        --matrixfile <filename>  Filename of the matrix to be factored." << endl;
    cout << "                                 Either CSV format for dense or MatrixMarket format ";
    cout << "for sparse." << endl;
    cout << "        --k <integer value>      The common dimension for factors W and H." << endl;
    cout << "        [--algorithm  BPP]       NMF algorithms: " << endl;
    cout << "                                     MU:    multiplicative updating " << endl;
    cout << "                                     HALS:  hierarchical alternating least squares" << endl;
    cout << "                                     RANK2: rank2 with optimal active set selection" << endl;
    cout << "                                     BPP:   block principal pivoting" << endl;
    cout << "        [--stopping  PG_RATIO]   Stopping criterion: " << endl;
    cout << "                                     PG_RATIO: Ratio of projected gradients" << endl;
    cout << "                                     DELTA:    Change in relative F-norm of W" << endl;
    cout << "        [--tol  0.005]           Tolerance for the selected stopping criterion." << endl;
    cout << "        [--tolcount  1]          Tolerance count; declare convergence after this many " << endl;
    cout << "                                 iterations with metric < tolerance; default is to " << endl;
    cout << "                                 declare convergence on the first such iteration." << endl;
    cout << "        [--infile_W  (empty)]    Dense mxk matrix to initialize W; CSV file." << endl;
    cout << "                                 If unspecified, W will be randomly initialized." << endl;
    cout << "        [--infile_H  (empty)]    Dense kxn matrix to initialize H; CSV file. " << endl;
    cout << "                                 If unspecified, H will be randomly initialized. " << endl;
    cout << "        [--outfile_W  w.csv]     Filename for the W matrix result." << endl;
    cout << "        [--outfile_H  h.csv]     Filename for the H matrix result." << endl;
    cout << "        [--miniter  5]           Minimum number of iterations to perform. " << endl;
    cout << "        [--maxiter  5000]        Maximum number of iterations to perform." << endl;
    cout << "        [--outprecision  6]      Write results with this many digits of precision." << endl;

    unsigned int hw_threads = GetMaxThreadCount();
    cout << "        [--maxthreads  " << setw(3) << hw_threads << "]      Upper limit to thread count. " 
         << endl;
    cout << "        [--normalize  1]         Whether to normalize W and scale H." << endl;
    cout << "                                     1 == yes, 0 == no " << endl;
    cout << "        [--verbose  1]           Whether to print updates to the screen. " << endl;
    cout << "                                     1 == print updates, 0 == silent " << endl;
    cout << endl;
}

//-----------------------------------------------------------------------------
bool ParseCommandLine(int argc, char* argv[], CommandLineOptions& opts)
{
    std::string tmp;
    int precision, user_max_threads;

    // set defaults
    opts.nmf_opts.algorithm          = NmfAlgorithm::BPP;
    opts.nmf_opts.height             = 0;
    opts.nmf_opts.width              = 0;
    opts.nmf_opts.k                  = 0;
    opts.nmf_opts.min_iter           = 5;
    opts.nmf_opts.max_iter           = 5000;
    opts.nmf_opts.verbose            = true;
    opts.nmf_opts.normalize          = true;
    opts.nmf_opts.tol                = 0.005;
    opts.nmf_opts.tolcount           = 1;
    opts.nmf_opts.prog_est_algorithm = NmfProgressAlgorithm::PG_RATIO;
    opts.show_help                   = false;
    opts.infile_A                    = std::string("");
    opts.infile_W                    = std::string("");
    opts.infile_H                    = std::string("");
    opts.outfile_W                   = DEFAULT_OUTFILE_W;
    opts.outfile_H                   = DEFAULT_OUTFILE_H;
    user_max_threads                 = -1;
    
    // The W and H matrices will be written to disk with this many
    // digits after the decimal point.
    opts.output_precision = 6;

    char c;
    int index;
    while (-1 != (c = getopt_long(argc, argv, 
                                  ":a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q",
                                  longopts, &index)))
    {
        switch (c)
        {
        case 'a':  // --matrixfile
            opts.infile_A = std::string(optarg);
            break;

        case 'b':  // --k
            opts.nmf_opts.k = atoi(optarg);
            break;

        case 'c':  // --algorithm
            tmp = std::string(optarg);
            ToUpper(tmp);
            if (STRING_MU == tmp)
                opts.nmf_opts.algorithm = NmfAlgorithm::MU;
            else if (STRING_HALS == tmp)
                opts.nmf_opts.algorithm = NmfAlgorithm::HALS;
            else if (STRING_RANK2 == tmp)
                opts.nmf_opts.algorithm = NmfAlgorithm::RANK2;
            else if (STRING_BPP == tmp)
                opts.nmf_opts.algorithm = NmfAlgorithm::BPP;
            else
                InvalidValue(tmp);
            break;

        case 'd':  // --stopping
            tmp = std::string(optarg);
            ToUpper(tmp);
            if (STRING_PG_RATIO == tmp)
                opts.nmf_opts.prog_est_algorithm = NmfProgressAlgorithm::PG_RATIO;
            else if (STRING_DELTA == tmp)
                opts.nmf_opts.prog_est_algorithm = NmfProgressAlgorithm::DELTA_FNORM;
            else
                InvalidValue(tmp);
            break;

        case 'e':  // --tol
            opts.nmf_opts.tol = atof(optarg);
            break;

        case 'f':  // --tolcount
            opts.nmf_opts.tolcount = atoi(optarg);
            break;

        case 'g':  // infile_W
            opts.infile_W = std::string(optarg);
            break;

        case 'h':  // infile_H
            opts.infile_H = std::string(optarg);
            break;

        case 'i':  // outfile_W
            opts.outfile_W = std::string(optarg);
            break;

        case 'j':  // outfile_H
            opts.outfile_H = std::string(optarg);
            break;

        case 'k':  // miniter
            opts.nmf_opts.min_iter = atoi(optarg);
            break;

        case 'l':  // maxiter
            opts.nmf_opts.max_iter = atoi(optarg);
            break;

        case 'm':  // outprecision
            precision = atoi(optarg);
            if (precision <= 0)
                precision = std::numeric_limits<float>::max_digits10;
            else if (precision >= std::numeric_limits<double>::max_digits10)
                precision = std::numeric_limits<double>::max_digits10;
            opts.output_precision = precision;
            break;

        case 'n':  // maxthreads
            user_max_threads = atoi(optarg);
            break;

        case 'o':  // normalize
            opts.nmf_opts.normalize = (0 != atoi(optarg));
            break;

        case 'p':  // verbose
            opts.nmf_opts.verbose = (0 != atoi(optarg));
            break;

        case 'q':  // help
            opts.show_help = true;
            break;
        case ':':
            // missing option argument
            cerr << "missing argument for option " << argv[optind-1] << endl;
            return false;
            break;
        case '?':
        default:
            // invalid option
            assert(optind >= 1);
            assert(optind <= argc);
            cerr << "invalid option: " << argv[optind-1] << endl;
            return false;
            break;
        }
    }

    // if no command line args, user wants help
    if (1 == argc)
        opts.show_help = true;

    // found --help on the command line
    if (opts.show_help)
        return false;

    // adjust the thread count to not exceed HW capabilities
    int hw_max_threads = GetMaxThreadCount();

    // if user did not specify the thread count, set it to HW max
    if (user_max_threads <= 0)
        user_max_threads = hw_max_threads;
    opts.nmf_opts.max_threads = std::min(user_max_threads, hw_max_threads);

    // check for required options
    if (opts.infile_A.empty())
    {
        cerr << "required command line argument --matrixfile not found" << endl;
        return false;
    }

    if ( (0 == opts.nmf_opts.k) && 
         (NmfAlgorithm::RANK2 != opts.nmf_opts.algorithm))
    {
        cerr << "required command line argument --k not found" << endl;
        return false;
    }

    // force k==2 if using RANK2
    if (NmfAlgorithm::RANK2 == opts.nmf_opts.algorithm)
    {
        if (2 != opts.nmf_opts.k)
        {
            cerr << "warning: forcing k=2 for RANK2 algorithm" << endl;
            opts.nmf_opts.k = 2;
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
bool IsValid(const CommandLineOptions& opts)
{
    // validate but ignore matrix, which has yet to be loaded
    return IsValid(opts.nmf_opts, false);
}
