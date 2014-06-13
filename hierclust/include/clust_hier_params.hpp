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

// This file defines a set of parameters needed
// inside hierarchical clustering based on rank-2 NMF.
// Should be included in:
//   clust_hier.hpp
//   clust_hier_sparse.hpp

#pragma once

#include <list>
#include <vector>
#include "nmf.hpp"
#include "dense_matrix.hpp"

template <typename T,
          template <typename> class MatrixType,
          template <typename, template <typename> class Matrix> class Solver>
class ClustHierParams
{
public:
    
    static T CLUST_OUTLIER_LABEL;
    static T CLUST_LEAF_PRIORITY;

    NmfOptions nmf_opts;
	int trial_allowance;
	T unbalanced;

	std::vector<T>& labels;
	std::vector<T>& priority;
	std::list< DenseMatrix<T> >& W_buffer;
    Solver<T, MatrixType>& solver;

	ClustHierParams(const NmfOptions& _nmf_opts,
                    int _trial_allowance,
                    T _unbalanced,
                    std::vector<T>& _labels,
                    std::vector<T>& _priority,
                    std::list<DenseMatrix<T> >& _W_buffer,
                    Solver<T, MatrixType>& _solver)
        : nmf_opts(_nmf_opts),
          trial_allowance(_trial_allowance),
          unbalanced(_unbalanced),
          labels(_labels),
          priority(_priority),
          W_buffer(_W_buffer),
          solver(_solver)
      {}
};

template <typename T, 
          template <typename> class MatrixType,
          template <typename, template <typename> class Matrix> class Solver>
T ClustHierParams<T, MatrixType, Solver>::CLUST_OUTLIER_LABEL = T(-1);

template <typename T,
          template <typename> class MatrixType,
          template <typename, template <typename> class Matrix> class Solver>
T ClustHierParams<T, MatrixType, Solver>::CLUST_LEAF_PRIORITY = T(-2);

