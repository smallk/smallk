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
#include <string>

bool Postprocess(const std::vector<std::vector<int> >& clusters,
                 const std::vector<std::vector<double> >& result_Ws,
                 const std::vector<std::string>& dictionary,
                 const int height,
                 const int max_leaf_nodes,
                 const int num_top_terms, 
                 const int term_count, // term-frequency matrix height  
                 const int doc_count); // term-frequency matrix width

