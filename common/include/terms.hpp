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
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include "dense_matrix.hpp"

//-----------------------------------------------------------------------------
template <typename T>
void TopTerms(const int maxterms,
              const DenseMatrix<T>& V, // column vector (single col of W)
              std::vector<int>& sort_indices,
              std::vector<int>& term_indices)
{
    // Sort the row indices for topic vector V into decreasing order.
    // Compare data elements to rearrange the indices.

    int height = V.Height();

    if (sort_indices.size() < static_cast<unsigned int>(height))
        throw std::runtime_error("TopTerms: index array too small");

    if (term_indices.size() < static_cast<unsigned int>(maxterms))
        throw std::runtime_error("TopTerms: term array too small");

    const T* data = V.LockedBuffer();

    // initialize the indices for the sort
    for (int q=0; q<height; ++q)
        sort_indices[q] = q;

    std::sort(&sort_indices[0], &sort_indices[0] + height,
              [&data](int i1, int i2) {return data[i1] > data[i2];});

    size_t max_terms = std::min(maxterms, height);
    for (size_t q=0; q<max_terms; ++q)
    {
        int index = sort_indices[q];
        assert(index >= 0);
        assert(index < height);
        term_indices[q] = index;
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void TopTerms(const int maxterms,
              const T* buf_w,
              const unsigned int ldim,
              const unsigned int height, 
              const unsigned int width,
              std::vector<int>& term_indices)
{
    if (height < width)
        throw std::logic_error("TopTerms: height of W buffer must be >= width");

    unsigned int max_terms = maxterms;
    if (max_terms > height)
        max_terms = height;

    std::vector<unsigned int> sort_indices(height);

    // The term indices for each column will be packed into the term_indices
    // array.  Indices for column 0 will occupy elements 0..maxterms-1.
    // Indices for the next column will occupy the next maxterm elements, etc.

    if (term_indices.size() < maxterms*width)
        term_indices.resize(maxterms*width);

    for (unsigned int c=0; c<width; ++c)
    {
        unsigned int term_offset = c*maxterms;

        // initialize indices for this column's sort
        for (unsigned int q=0; q<height; ++q)
            sort_indices[q] = q;

        // offset to column c in buf_w
        unsigned int col_offset = c*height;
        const T* data = &(buf_w[col_offset]);

        std::sort(&sort_indices[0], &sort_indices[0] + height,
                  [&data](unsigned int i1, unsigned int i2) {return data[i1] > data[i2];});

        for (unsigned int q=0; q<max_terms; ++q)
        {
            unsigned int index = sort_indices[q];
            assert(index >= 0u);
            assert(index < height);
            term_indices[term_offset + q] = index;
        }
    }
}
