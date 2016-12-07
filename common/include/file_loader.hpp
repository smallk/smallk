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
#include "sparse_matrix_decl.hpp"
#include "sparse_matrix_impl.hpp"
#include "sparse_matrix_io.hpp"
#include "delimited_file.hpp"

// Given a file name, determine if it contains a 
// sparse or dense matrix in one of the supported formats.
bool IsDense (const std::string& file_path);
bool IsSparse(const std::string& file_path);

//-----------------------------------------------------------------------------
template <typename T>
bool LoadSparseMatrix(const std::string& file_path,
                      SparseMatrix<T>& A,
                      unsigned int& height,
                      unsigned int& width,
                      unsigned int& nz)
{
    if (IsMatrixMarketFile(file_path))
        return LoadMatrixMarketFile(file_path, A, height, width, nz);

    // other sparse formats TBD

    return false;
}

//-----------------------------------------------------------------------------
template <typename T>
bool LoadDenseMatrix(const std::string& file_path,
                     std::vector<T>& data,
                     unsigned int& height,
                     unsigned int& width)
{
    if (IsDelimitedFile(file_path))
        return LoadDelimitedFile(data, height, width, file_path);

    // other dense formats TBD

    return false;
}
