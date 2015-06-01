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
#include "file_format.hpp"

void FlatClustWriteResults(const std::string& assignfilepath,
                           const std::string& fuzzyfilepath,
                           const std::string& resultfilepath,
                           const std::vector<unsigned int>& assignments,
                           const std::vector<float>& probabilities,
                           const std::vector<std::string>& dictionary,
                           const std::vector<int>& term_indices,
                           const FileFormat format,
                           const unsigned int maxterms,
                           const unsigned int num_docs,
                           const unsigned int num_clusters);

void FlatClustWriteResults(const std::string& outdir, 
                           const std::vector<unsigned int>& assignments,
                           const std::vector<float>& probabilities,
                           const std::vector<std::string>& dictionary,
                           const std::vector<int>& term_indices,
                           const FileFormat format, 
                           const unsigned int maxterms,
                           const unsigned int num_docs,
                           const unsigned int num_clusters);
