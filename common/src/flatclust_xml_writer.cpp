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

#include <sstream>
#include <iostream>
#include "flatclust_writer.hpp"

using std::cout;
using std::endl;

// indentation
const std::string S4("    ");
const std::string S8  = S4 + S4;
const std::string S12 = S4 + S8;
const std::string S16 = S4 + S12;

//-----------------------------------------------------------------------------
void FlatclustXmlWriter::WriteHeader(std::ofstream& outfile, 
                                     const int doc_count)
{
    std::ostringstream count_str;
    count_str << doc_count;

    outfile << "<?xml version=\"1.0\"?>" << endl;
    outfile << "<DataSet id=\"" << count_str.str() << "\">" << endl;
}

//-----------------------------------------------------------------------------
void FlatclustXmlWriter::WriteNodeBegin(std::ofstream& outfile, 
                                        const int node_id)
{
    outfile << S4 << "<node id=\"" << node_id << "\">" << endl;
}

//-----------------------------------------------------------------------------
void FlatclustXmlWriter::WriteDocCount(std::ofstream& outfile, 
                                       const int count)
{
    outfile << S8 << "<doc_count>" << count << "</doc_count>" << endl;
}

//-----------------------------------------------------------------------------
void FlatclustXmlWriter::WriteTopTerms(std::ofstream& outfile,
                                       const unsigned int offset,
                                       const unsigned int maxterms,
                                       const std::vector<int>& term_indices,
                                       const std::vector<std::string>& dictionary)
{
    outfile << S8 << "<top_terms>" << endl;
    for (size_t q=0; q != maxterms; ++q)
    {
        int index = term_indices[offset + q];
        outfile << S12 << "<term name=\"" << dictionary[index] << "\"/>" << endl;
    }
    outfile << S8 << "</top_terms>" << endl;
}

//-----------------------------------------------------------------------------
void FlatclustXmlWriter::WriteNodeEnd(std::ofstream& outfile)
{
    outfile << S4 << "</node>" << endl;
}

//-----------------------------------------------------------------------------
void FlatclustXmlWriter::WriteFooter(std::ofstream& outfile)
{
    outfile << "</DataSet>" << endl;
}
