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
#include <vector>
#include <fstream>

class IHierclustWriter
{
    // Interface for writing hierclust results.
public:

    virtual ~IHierclustWriter() {}

    virtual void WriteHeader     (std::ofstream& outfile, const int doc_count)    = 0;
    virtual void WriteNodeBegin  (std::ofstream& outfile, const int node_id)      = 0;
    virtual void WriteLevel      (std::ofstream& outfile, const int level)        = 0;
    virtual void WriteLabel      (std::ofstream& outfile, const int label)        = 0;
    virtual void WriteParentId   (std::ofstream& outfile, const int parent_id)    = 0;
    virtual void WriteParentLabel(std::ofstream& outfile, const int parent_label) = 0;
    virtual void WriteLeftChild  (std::ofstream& outfile, 
                                  const bool is_left_child, 
                                  const int lc_label)                             = 0;
    virtual void WriteRightChild (std::ofstream& outfile, const int rc_label)     = 0;
    virtual void WriteDocCount   (std::ofstream& outfile, const int count)        = 0;
    virtual void WriteTopTerms   (std::ofstream& outfile, 
                                  const std::vector<int>& term_indices,
                                  const std::vector<std::string>& dictionary)     = 0;
    virtual void WriteNodeEnd    (std::ofstream& outfile)                         = 0;
    virtual void WriteFooter     (std::ofstream& outfile)                         = 0;
};
