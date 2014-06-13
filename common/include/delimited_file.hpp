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

#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <stdexcept>

// specialization for integer data
bool WriteDelimitedFile(const int* buffer,
                        const unsigned int ldim,
                        const unsigned int height,
                        const unsigned int width,
                        const std::string& filename,
                        const char DELIM = ',');

// returns true if the file has an extension supported by this reader
bool IsDelimitedFile(const std::string& filename);

// skip any initial blank lines and/or commented lines in the file
std::streampos
SkipBlankLinesAndComments(std::ifstream& infile,
                          std::string& line);

void GetDimensions(std::ifstream& infile,
                   std::string& line,
                   unsigned int& height,
                   unsigned int& width,
                   const char DELIM = ',');

//-----------------------------------------------------------------------------
template <typename T>
bool WriteDelimitedFile(const T* buffer, 
                        const unsigned int ldim,
                        const unsigned int height, 
                        const unsigned int width,
                        const std::string& filename,
                        const unsigned int precision,
                        const char DELIM = ',')
{
    std::ofstream outfile(filename);
    if (!outfile)
        return false;

    outfile << std::scientific;
    outfile.precision(precision);

    // write in row-major order, to match Matlab
    for (unsigned int r=0; r != height; ++r)
    {
        for (unsigned int c=0; c != (width-1); ++c)
            outfile << buffer[c*ldim + r] << DELIM;

        outfile << buffer[(width-1)*ldim + r] << std::endl;
    }

    outfile.close();
    return true;
}

//-----------------------------------------------------------------------------
template <typename T>
bool LoadDelimitedFile(std::vector<T>& buffer,
                       unsigned int& height,
                       unsigned int& width,
                       const std::string& filename,
                       const char DELIM = ',')
{
    // This file reader skips blank lines and comment lines at the beginning
    // of the file.  Comment lines are identified by their first character;
    // the comment chars are found in delimited_file.cpp.
    
    std::ifstream infile(filename);
    if (!infile)
        return false;

    std::string line;
    std::streampos start_pos = SkipBlankLinesAndComments(infile, line);

    // check for empty file
    if (line.empty())
        return false;

    GetDimensions(infile, line, height, width, DELIM);
    
    // set the buffer to the correct size
    buffer.resize(height * width);

    // clear the flags and rewind to where the data begins
    infile.clear();
    infile.seekg(start_pos);

    char dummy;
    unsigned int r = 0;
    while (std::getline(infile, line))
    {
        if (infile.eof())
            break;

        // extract the next row of data and store in column-major order
        std::istringstream data(line);
        for (unsigned int c=0; c != width; ++c)
        {
            // write row r into the buffer; also extract delimiter
            data >> buffer[c*height + r];
            data >> dummy;
        }

        ++r;
    }

    assert(height == r);
    infile.close();
    return true;
}

//-----------------------------------------------------------------------------
template <typename T>
bool LoadMatrixArray(std::vector<std::vector<T> >& matrices,
                     const unsigned int matrix_height,
                     const unsigned int matrix_width,
                     const std::string& filename,
                     const char DELIM = ',')
{
    // Load an array of matrices, each of dimension 
    // matrix_height x matrix_width, from a single delimited file.  The 
    // 'matrices' array is an array of data buffers into which the file 
    // data is written.  The leading dimension is equal to matrix_height.
    
    // The data file can contain more data than is necessary to load the
    // matrices.  The first 'matrix_width' columns will be read, and the
    // first 'num_matrices * matrix_height' lines will be read.

    const unsigned int num_matrices = matrices.size();
    const unsigned int matrix_size = matrix_height * matrix_width;

    std::ifstream infile(filename);
    if (!infile)
        return false;

    std::string line;
    std::streampos start_pos = SkipBlankLinesAndComments(infile, line);

    // check for empty file
    if (line.empty())
        return false;

    unsigned int file_height, file_width;
    GetDimensions(infile, line, file_height, file_width, DELIM);
    
    // clear the flags and rewind to where the data begins
    infile.clear();
    infile.seekg(start_pos);

    bool horizontally_packed = false;
    bool vertically_packed   = false;

    // ensure that the file contains enough data
    if (matrix_width < matrix_height)
    {
        unsigned int expected_file_width = matrix_width * num_matrices;
        if (file_width < expected_file_width)
        {
            std::cerr << "\tLoadMatrixArray: width mismatch." << std::endl;
            std::cerr << "\tExpected " << expected_file_width << " elements "
                      << "per line, found " << file_width << "." << std::endl;
            return false;
        }
        if (file_height < matrix_height)
        {
            std::cerr << "\tLoadMatrixArray: height mismatch." << std::endl;
            std::cerr << "\tExpected " << matrix_height << " lines "
                      << "of data, found " << file_height << "." << std::endl;
            return false;
        }
        horizontally_packed = true;
    }
    else if (matrix_height < matrix_width)
    {
        unsigned int expected_file_height = matrix_height * num_matrices;
        if (file_height < expected_file_height)
        {
            std::cerr << "\tLoadMatrixArray: height mismatch."  << std::endl;
            std::cerr << "\tExpected " << expected_file_height << " lines "
                      << "of data, found " << file_height << "." << std::endl;
            return false;
        }
        if (file_width < matrix_width)
        {
            std::cerr << "\tLoadMatrixArray: width mismatch." << std::endl;
            std::cerr << "\tExpected " << matrix_width << " elements "
                      << "per line, found " << file_width << "." << std::endl;
            return false;
        }
        vertically_packed = true;
    }
    else
    {
        // matrix_width == matrix_height
        if ( (file_width == matrix_width * num_matrices) &&
             (file_height == matrix_height))
        {
            horizontally_packed = true;
        }
        else if ( (file_height == matrix_height * num_matrices) &&
                  (file_width == matrix_width))
        {
            vertically_packed = true;
        }
        else
        {
            std::cerr << "\tLoadMatrixArray: invalid init file " << filename 
                      << std::endl;
            return false;
        }
    }

    assert(horizontally_packed || vertically_packed);

    // resize the buffers if necessary
    for (unsigned int i=0; i<num_matrices; ++i)
    {
        unsigned int size = matrices[i].size();
        if (size < matrix_size)
            matrices[i].resize(matrix_size);
    }

    char dummy;
    unsigned int matrices_loaded;

    if (horizontally_packed)
    {
        unsigned int expected_width = num_matrices * matrix_width;
        for (unsigned int r=0; r<matrix_height; ++r)
        {
            if (!std::getline(infile, line))
                break;
            if (infile.eof())
                break;

            assert(line.size() >= expected_width);

            matrices_loaded = 0;
            std::istringstream data(line);
            for (unsigned int q=0; q<num_matrices; ++q)
            {
                // write row r into each buffer
                for (unsigned int c=0; c<matrix_width; ++c)
                {
                    data >> (matrices[q])[c*matrix_height + r];
                    data >> dummy; // delimiter
                }
                ++matrices_loaded;
            }

            assert(num_matrices == matrices_loaded);
        }
    }
    else
    {
        matrices_loaded = 0;
        for (unsigned int q=0; q<num_matrices; ++q)
        {
            for (unsigned int r=0; r<matrix_height; ++r)
            {
                if (!std::getline(infile, line))
                    break;
                if (infile.eof())
                    break;

                assert(line.size() >= matrix_width);

                std::istringstream data(line);

                // write row r in to the buffer for matrix q
                for (unsigned int c=0; c<matrix_width; ++c)
                {
                    data >> (matrices[q])[c*matrix_height + r];
                    data >> dummy; // delimiter
                }
            }
            ++matrices_loaded;
        }

        assert(num_matrices == matrices_loaded);
    }
        
    infile.close();

    if (num_matrices != matrices_loaded)
    {
        std::cerr << "\tLoadMatrixArray: loaded " << matrices_loaded 
                  << " matrices; expected " << num_matrices << std::endl;
        return false;
    }

    return true;
}


