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

enum class NmfResult
{
    OK             =  0,
    NOTINITIALIZED = -1,
    INITIALIZED    = -2,
    BAD_PARAM      = -3,
    FAILURE        = -4,
    SIZE_TOO_LARGE = -5
};

enum class NmfAlgorithm
{
    MU,       // Lee & Seung, multiplicative updating
    HALS,     // Cichocki & Pan, hierarchical alternating least squares
    RANK2,    // Kuang and Park, rank2 specialization
    BPP       // Kim and Park, block principal pivoting
};

// progress estimation algorithms (stopping criterion)
enum class NmfProgressAlgorithm
{
    PG_RATIO,    // PG_i / PG_1 (ratio of the ith PG to that of iteration 1)
    DELTA_FNORM  // relative change in the Frobenius norm of W
};

struct NmfStats
{
    NmfStats()
    {
        elapsed_us = 0u;
        iteration_count = 0;
    }

    unsigned long long elapsed_us;
    int iteration_count;
};

struct NmfOptions
{
    double tol;
    NmfAlgorithm algorithm;
    NmfProgressAlgorithm prog_est_algorithm;
    int height;
    int width;
    int k;
    int min_iter;
    int max_iter;
    int tolcount;
    int max_threads;
    bool verbose;
    bool normalize;
};

void NmfInitialize(int argc, char* argv []);
NmfResult NmfIsInitialized();
void NmfFinalize();

bool IsValid(const NmfOptions& opts, bool validate_matrix = true);

NmfResult Nmf(const NmfOptions& options,
              double* buf_A, int ldim_A,
              double* buf_W, int ldim_W,
              double* buf_H, int ldim_H,
              NmfStats& stats);

NmfResult NmfSparse(const NmfOptions& options,
                    const unsigned int height,       // height of sparse matrix
                    const unsigned int width,        // width of sparse matrix
                    const unsigned int nz,           // number of nonzeros
                    const unsigned int* col_offsets,
                    const unsigned int* row_indices,
                    const double* data,
                    double* buf_W, int ldim_W,
                    double* buf_H, int ldim_H,
                    NmfStats& stats);

