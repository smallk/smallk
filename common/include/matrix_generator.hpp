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

#include <vector>
#include <thread>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "random.hpp"
#include "thread_utils.hpp"

namespace MatrixGenData
{
    // use sequential code to initialize matrices with fewer than this many
    // elements; use parallel code for matrices larger than this
    static unsigned int CUTOFF = 32768;
}

//-----------------------------------------------------------------------------
template <typename T>
bool RandomMatrix(std::vector<T>& buf,
                  const unsigned int height, // number of data rows, <= ldim
                  const unsigned int width,  // number of data cols
                  Random& rng,
                  const T rng_center = T(0.5), 
                  const T rng_radius = T(0.5))
{
    // Fill a matrix of dimension height x width with random numbers.
    // The matrix is stored in the buffer 'buf', which will be resized if
    // needed.

    typename std::vector<T>::size_type required_size = height*width;
    if (buf.size() < required_size)
        buf.resize(required_size);

    unsigned int max_threads = GetMaxThreadCount();

    if ( (required_size < MatrixGenData::CUTOFF) || (1 == max_threads))
        RandomMatrixSequential(&buf[0], height, height, width, rng, rng_center, rng_radius);
    else
        RandomMatrixParallel(&buf[0], height, height, width, rng, rng_center, rng_radius);
    
    return true;
}

//-----------------------------------------------------------------------------
template <typename T>
bool RandomMatrix(T* buf,
                  const unsigned int ldim,   // leading dimension of buffer
                  const unsigned int height, // number of data rows, <= ldim
                  const unsigned int width,  // number of data cols
                  Random& rng,
                  const T rng_center = T(0.5), 
                  const T rng_radius = T(0.5))
{
    // Fill a column-major matrix with random numbers.
    // Assumes the buffer is large enough.

    uint64_t size = height * width;
    unsigned int max_threads = GetMaxThreadCount();

    // if fewer than 16k elements use sequential code to avoid thread overhead
    if ( (size < MatrixGenData::CUTOFF) || (1 == max_threads))
        RandomMatrixSequential(buf, ldim, height, width, rng, rng_center, rng_radius);
    else
        RandomMatrixParallel(buf, ldim, height, width, rng, rng_center, rng_radius);

    return true;
}

//-----------------------------------------------------------------------------
//
//                   H E L P E R    F U N C T I O N S
//
//  Not to be called from user code.  Users should instead call one of the
//  'RandomMatrix' functions above.
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
template <typename T>
void ThreadFuncRank2W(const int rng_seed,
                      const unsigned int num_threads,
                      const unsigned int index,
                      T* buf,
                      const unsigned int ldim,                      
                      const unsigned int height,
                      const T rng_center,
                      const T rng_radius)                     
{
    // This thread function randomly initializes rank2 'W' matrices of
    // dimension height x 2.  Each thread initializes a maximum of
    // 'height/num_threads' rows in each column.

    Random rng;
    rng.SeedFromInt(rng_seed);
    
    unsigned int num_elts = height / num_threads;
    unsigned int start = index * num_elts;
    unsigned int end = std::min(height, (index+1) * num_elts);

    // column 0
    for (unsigned int r=start; r<end; ++r)
        buf[r] = (T) rng.RandomDouble(rng_center, rng_radius);
    
    // column 1
    unsigned int offset = ldim;
    for (unsigned int r=start; r<end; ++r)
        buf[offset + r] = (T) rng.RandomDouble(rng_center, rng_radius);
}

//-----------------------------------------------------------------------------
template <typename T>
void ThreadFuncInitTall(const int rng_seed,
                        const unsigned int num_threads,
                        const unsigned int index,
                        T* buf,
                        const unsigned int ldim,                         
                        const unsigned int height,
                        const unsigned int width,
                        const T rng_center,
                        const T rng_radius)
{
    // This thread function randomly initializes matrices having height > width.
    // Each thread initializes 'height/num_threads' rows in each column.

    assert(height > width);
    
    Random rng;
    rng.SeedFromInt(rng_seed);
    
    unsigned int num_rows = height / num_threads;
    unsigned int start = index * num_rows;
    unsigned int end = std::min(height, (index+1) * num_rows);

    for (unsigned int c=0; c<width; ++c)
    {
        unsigned int offset = c*ldim;
        for (unsigned int r=start; r<end; ++r)
        {
            buf[offset + r] = (T) rng.RandomDouble(rng_center, rng_radius);
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void ThreadFuncRank2H(const int rng_seed,
                      const unsigned int num_threads,
                      const unsigned int index,
                      T* buf,
                      const unsigned int ldim,                      
                      const unsigned int width,
                      const T rng_center,
                      const T rng_radius)
{
    // This thread function randomly initializes rank2 'H' matrices of
    // dimension 2 x width.  Each thread initializes a maximum of
    // 'width/num_threads' columns.

    Random rng;
    rng.SeedFromInt(rng_seed);
    
    unsigned int num_cols = width / num_threads;
    unsigned int start = index * num_cols;
    unsigned int end = std::min(width, (index+1) * num_cols);

    for (unsigned int c=start; c<end; ++c)
    {
        unsigned int offset = c*ldim;

        // row 0
        buf[offset + 0] = (T) rng.RandomDouble(rng_center, rng_radius);

        // row 1
        buf[offset + 1] = (T) rng.RandomDouble(rng_center, rng_radius);
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void ThreadFuncInitWide(const int rng_seed,
                        const unsigned int num_threads,
                        const unsigned int index,
                        T* buf,
                        const unsigned int ldim,                         
                        const unsigned int height,
                        const unsigned int width,
                        const T rng_center,
                        const T rng_radius)
{
    // This thread function randomly initializes matrices having height <= width.
    // Each thread initializes a maximum of 'width/num_threads' columns.

    assert(height <= width);

    Random rng;
    rng.SeedFromInt(rng_seed);
    
    unsigned int num_cols = width / num_threads;
    unsigned int start = index * num_cols;
    unsigned int end = std::min(width, (index+1) * num_cols);

    for (unsigned int c=start; c<end; ++c)
    {
        unsigned int offset = c*ldim;
        for (unsigned int r=0; r<height; ++r)
        {
            buf[offset + r] = (T) rng.RandomDouble(rng_center, rng_radius);
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void RandomMatrixSequential(T* buf,
                            const unsigned int ldim,
                            const unsigned int height,
                            const unsigned int width,
                            Random& rng,
                            const T rng_center,
                            const T rng_radius)
{
    // Fills a column-major matrix with random numbers.
    // IMPORTANT: this function assumes the buffer is large enough.

    for (unsigned int c=0; c<width; ++c)
    {
        unsigned int col_offset = c*ldim;
        for (unsigned int r=0; r<height; ++r)
        {
            buf[r + col_offset] = (T) rng.RandomDouble(rng_center, rng_radius);
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void RandomMatrixParallel(T* buf,
                          const unsigned int ldim,
                          const unsigned int height,
                          const unsigned int width,
                          Random& rng,
                          const T rng_center,
                          const T rng_radius)
{
    std::vector<std::thread> threads;
    unsigned int max_threads = GetMaxThreadCount();    

    // the rng supplied as argument is used to generate a random seed
    // for each thread's private rng
    
    if (height <= width)
    {
        if (2 == height)
        {
            // special case for rank2 H matrices
            for (unsigned int k=0; k<max_threads; ++k)
            {
                threads.push_back(std::thread(ThreadFuncRank2H<T>,
                                              rng.RandomInt(),
                                              max_threads, k,
                                              buf, ldim, width,
                                              rng_center, rng_radius));
            }
        }
        else
        {
            for (unsigned int k=0; k<max_threads; ++k)
            {
                threads.push_back(std::thread(ThreadFuncInitWide<T>,
                                              rng.RandomInt(),
                                              max_threads, k,
                                              buf, ldim, height, width,
                                              rng_center, rng_radius));
            }
        }
    }
    else
    {
        if (2 == width)
        {
            // special case for rank2 W matrices
            for (unsigned int k=0; k<max_threads; ++k)
            {
                threads.push_back(std::thread(ThreadFuncRank2W<T>,
                                              rng.RandomInt(),
                                              max_threads, k,
                                              buf, ldim, height,
                                              rng_center, rng_radius));
            }
        }
        else
        {
            for (unsigned int k=0; k<max_threads; ++k)
            {
                threads.push_back(std::thread(ThreadFuncInitTall<T>,
                                              rng.RandomInt(),
                                              max_threads, k,
                                              buf, ldim, height, width,
                                              rng_center, rng_radius));
            }
        }
    }

    // wait for threads to finish
    std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
}
