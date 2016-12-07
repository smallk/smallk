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

class IFlatclustWriter
{
    // Interface for writing flatclust results.
public:

    virtual ~IFlatclustWriter() {}

    virtual void WriteHeader   (std::ofstream& outfile, const int doc_count) = 0;
    virtual void WriteNodeBegin(std::ofstream& outfile, const int node_id)   = 0;
    virtual void WriteDocCount (std::ofstream& outfile, const int count)     = 0;
    virtual void WriteTopTerms (std::ofstream& outfile,
                                const unsigned int offset,
                                const unsigned int maxterms,
                                const std::vector<int>& term_indices,
                                const std::vector<std::string>& dictionary)  = 0;
    virtual void WriteNodeEnd  (std::ofstream& outfile)                      = 0;
    virtual void WriteFooter   (std::ofstream& outfile)                      = 0;
};

//-----------------------------------------------------------------------------
class FlatclustXmlWriter : public IFlatclustWriter
{
public:

    void WriteHeader     (std::ofstream& outfile, const int doc_count);
    void WriteNodeBegin  (std::ofstream& outfile, const int node_id);
    void WriteDocCount   (std::ofstream& outfile, const int count);
    void WriteTopTerms   (std::ofstream& outfile, 
                          const unsigned int offset,
                          const unsigned int maxterms,
                          const std::vector<int>& term_indices,
                          const std::vector<std::string>& dictionary);
    void WriteNodeEnd    (std::ofstream& outfile);
    void WriteFooter     (std::ofstream& outfile);
};

//-----------------------------------------------------------------------------
class FlatclustJsonWriter : public IFlatclustWriter
{
public:

    FlatclustJsonWriter() : nodes_written_(0u) {}

    void WriteHeader     (std::ofstream& outfile, const int doc_count);
    void WriteNodeBegin  (std::ofstream& outfile, const int node_id);
    void WriteDocCount   (std::ofstream& outfile, const int count);
    void WriteTopTerms   (std::ofstream& outfile,
                          const unsigned int offset,
                          const unsigned int maxterms,
                          const std::vector<int>& term_indices,
                          const std::vector<std::string>& dictionary);
    void WriteNodeEnd    (std::ofstream& outfile);
    void WriteFooter     (std::ofstream& outfile);

private:

    unsigned int nodes_written_;
};

