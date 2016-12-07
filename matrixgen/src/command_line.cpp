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

#include <iostream>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include "command_line.hpp"
#include "matrix_types.hpp"
#include "utils.hpp"

using std::cout;
using std::cerr;
using std::endl;

//-----------------------------------------------------------------------------
void print_opts(const CommandLineOptions& opts)
{
    cout << "\n      Command line options: \n" << endl;
    cout << "\t            height: " << opts.height << endl;
    cout << "\t             width: " << opts.width << endl;
    cout << "\t          filename: " << opts.filename << endl;
    cout << "\t              type: " << opts.type << endl;
    cout << "\t        rng_center: " << opts.rng_center << endl;
    cout << "\t        rng_radius: " << opts.rng_radius << endl;
    cout << "\t         precision: " << opts.precision << endl;
    cout << "\t        nz_per_col: " << opts.nz_per_col << endl;
    cout << endl;
}

//-----------------------------------------------------------------------------
void print_usage(const std::string& program_name)
{
    cout << endl;
    cout << "Usage: " << program_name << endl;
    cout << "         --height <number of rows> " << endl;
    cout << "         --width  <number of cols> " << endl;
    cout << "         --filename <path> " << endl;
    cout << "        [--type  UNIFORM]  UNIFORM:     matrix with uniformly-distributed random entries" << endl;
    cout << "                           DENSE_DIAG:  dense diagonal matrix with uniform";
    cout << " random entries" << endl;
    cout << "                           SPARSE_DIAG: sparse diagonal matrix with uniform";
    cout << " random entries" << endl;
    cout << "                           IDENTITY:    identity matrix" << endl;
    cout << "                           ONES:        matrix of all ones" << endl;
    cout << "                           ZEROS:       matrix of all zeros" << endl;
    cout << "                           SPARSE:      sparse matrix with uniform random entries" 
         << endl;
    cout << "                                        specify 'nz_per_col' to control occupancy" << endl;
    cout << endl;
    cout << "        [--rng_center  0.5]   center of random numbers" << endl;
    cout << "        [--rng_radius  0.5]   radius of random numbers" << endl;
    cout << "        [--precision   6]     digits of precision" << endl;
    cout << "        [--nz_per_col  1]     (SPARSE only) nonzeros per column" << endl;
    cout << endl;
}

//-----------------------------------------------------------------------------
void invalid_arg(const std::string& arg)
{
    std::string msg("Invalid command-line argument: ");
    msg += arg;
    throw std::runtime_error(msg);
}

//-----------------------------------------------------------------------------
void invalid_value(const std::string& arg)
{
    std::string msg("Invalid value specified for command-line argument ");
    msg += arg;
    throw std::runtime_error(msg);
}

//-----------------------------------------------------------------------------
void parse_command_line(int argc, char* argv[], CommandLineOptions& opts)
{
    std::string tmp;

    // set defaults
    opts.height        = 0;
    opts.width         = 0;
    opts.filename      = std::string("");
    opts.type          = std::string("UNIFORM");
    opts.rng_center    = 0.5;
    opts.rng_radius    = 0.5;
    opts.precision     = 6;
    opts.nz_per_col    = 1;
    
    for (int k=1; k<argc; k += 2)
    {
        if ( ('-' == argv[k][0]) && ('-' == argv[k][1]))
        {
            char c = argv[k][2];
            if ('h' == c)
            {
                // --height
                int height = atoi(argv[k+1]);
                if (height <= 0)
                    invalid_value(std::string(argv[k]));
                
                opts.height = height;
            }
            else if ('w' == c)
            {
                // --width
                int width = atoi(argv[k+1]);
                if (width <= 0)
                    invalid_value(std::string(argv[k]));
                
                opts.width = width;
            }
            else if ('f' == c)
            {
                // --filename
                opts.filename = std::string(argv[k+1]);
            }
            else if ('t' == c)
            {
                // --type
                opts.type = std::string(argv[k+1]);
            }
            else if ('p' == c)
            {
                // --precision
                int precision = atoi(argv[k+1]);
                if (precision <= 0)
                    invalid_value(std::string(argv[k]));

                // set an upper bound on the precision
                opts.precision = precision;
                if (opts.precision > std::numeric_limits<double>::max_digits10)
                    opts.precision = std::numeric_limits<double>::max_digits10;
            }
            else if ('n' == c)
            {
                // --nz_per_col
                int nz_per_col = atoi(argv[k+1]);
                if (nz_per_col <= 0)
                    invalid_value(std::string(argv[k]));

                opts.nz_per_col = nz_per_col;
            }
            else if ('r' == c)
            {
                char c2 = argv[k][6];
                if ('c' == c2)
                {
                    // --rng_center
                    opts.rng_center = atof(argv[k+1]);
                }
                else if ('r' == c2)
                {
                    // --rng_radius
                    double rng_radius = atoi(argv[k+1]);
                    if (rng_radius <= 0.0)
                        invalid_value(std::string(argv[k]));

                    opts.rng_radius = rng_radius;
                }
                else
                {
                    invalid_arg(std::string(argv[k]));
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
bool validate_opts(const CommandLineOptions& opts)
{
    // filename must be specified
    if (opts.filename.empty())
    {
        cerr << "matrixgen error - required argument /'--filename/' not found" << endl;
        return false;
    }
    
    // height must be positive
    if (opts.height <= 0)
    {
        cerr << "matrixgen error - height must be a positive integer" << endl;
        return false;
    }

    // width must be positive
    if (opts.width <= 0)
    {
        cerr << "matrixgen error - width must be a positive integer" << endl;
        return false;
    }

    // check for supported type
    std::string type = opts.type;
    if ( !CaseInsensitiveEquals(TYPE_UNIFORM, type)    &&
         !CaseInsensitiveEquals(TYPE_DENSEDIAG, type)  &&
         !CaseInsensitiveEquals(TYPE_SPARSEDIAG, type) &&
         !CaseInsensitiveEquals(TYPE_IDENTITY, type)   &&
         !CaseInsensitiveEquals(TYPE_ONES, type)       &&
         !CaseInsensitiveEquals(TYPE_ZEROS, type)      &&
         !CaseInsensitiveEquals(TYPE_SPARSE, type))
    {
        cerr << "matrixgen error - unknown matrix type " << type << endl;
        return false;
    }

    // diagonal and identity matrices must be square
    if ( CaseInsensitiveEquals(TYPE_DENSEDIAG, type)  ||
         CaseInsensitiveEquals(TYPE_SPARSEDIAG, type) ||
         CaseInsensitiveEquals(TYPE_IDENTITY, type))
    {
        if (opts.height != opts.width)
        {
            cerr << "matrixgen error - height and width must be equal for the specified matrix type" 
                 << endl;
            return false;
        }
    }

    // sparse matrices must have fewer nonzeros per column than the number of rows
    if (CaseInsensitiveEquals(TYPE_SPARSE, type))
    {
        if (opts.nz_per_col > opts.height)
        {
            cerr << "matrixgen error - nz_per_col cannot exceed the height" << endl;
            return false;
        }
    }

    return true;
}
