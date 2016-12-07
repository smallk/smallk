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

#include "hierclust_writer.hpp"

class HierclustJsonWriter : public IHierclustWriter
{
public:

    HierclustJsonWriter() : nodes_written_(0u) {}

    void WriteHeader     (std::ofstream& outfile, const int doc_count);
    void WriteNodeBegin  (std::ofstream& outfile, const int node_id);
//    void WriteLevel      (std::ofstream& outfile, const int level);
//    void WriteLabel      (std::ofstream& outfile, const int label);
    void WriteParentId   (std::ofstream& outfile, const int parent_id);
//    void WriteParentLabel(std::ofstream& outfile, const int parent_label);
    void WriteLeftChild  (std::ofstream& outfile, 
                          const bool is_left_child, 
                          const int lc_label);
    void WriteRightChild (std::ofstream& outfile, const int rc_label);
    void WriteDocCount   (std::ofstream& outfile, const int count);
    void WriteTopTerms   (std::ofstream& outfile, 
                          const std::vector<int>& term_indices,
                          const std::vector<std::string>& dictionary);
    void WriteNodeEnd    (std::ofstream& outfile);
    void WriteFooter     (std::ofstream& outfile);

private:

    unsigned int nodes_written_;
};

