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

#include <string>
#include "nmf.hpp"

struct FlatClustOptions
{
    NmfOptions nmf_opts;

    int maxterms;
    int num_clusters;
    bool verbose;
};

bool IsValid(const FlatClustOptions& opts, bool validate_matrix = true);

Result FlatClust(const NmfOptions& options,
                 double* buf_A, int ldim_A,
                 double* buf_W, int ldim_W,
                 double* buf_H, int ldim_H,
                 NmfStats& stats);

Result FlatClustSparse(const NmfOptions& options,
                       const unsigned int height,       // height of sparse matrix
                       const unsigned int width,        // width of sparse matrix
                       const unsigned int nz,           // number of nonzeros
                       const unsigned int* col_offsets,
                       const unsigned int* row_indices,
                       const double* data,
                       double* buf_W, int ldim_W,
                       double* buf_H, int ldim_H,
                       NmfStats& stats);
