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
#include "hierclust_xml_writer.hpp"

using std::cout;
using std::endl;

// indentation
const std::string S4("    ");
const std::string S8  = S4 + S4;
const std::string S12 = S4 + S8;
const std::string S16 = S4 + S12;

//-----------------------------------------------------------------------------
void HierclustXmlWriter::WriteHeader(std::ofstream& outfile, 
                                     const int doc_count)
{
    outfile << "<?xml version=\"1.0\"?>" << endl;
    outfile << "<DataSet id=\"" << doc_count << "\">" << endl;
}

//-----------------------------------------------------------------------------
void HierclustXmlWriter::WriteNodeBegin(std::ofstream& outfile, 
                                        const int node_id)
{
    outfile << S4 << "<node id=\"" << node_id << "\">" << endl;
}
/*
//-----------------------------------------------------------------------------
void HierclustXmlWriter::WriteLevel(std::ofstream& outfile, 
                                    const int level)
{
    outfile << S8 << "<level>" << level << "</level>" << endl;
}

//-----------------------------------------------------------------------------
void HierclustXmlWriter::WriteLabel(std::ofstream& outfile, 
                                    const int label)
{
    outfile << S8 << "<label>" << label << "</label>" << endl;
}
*/
//-----------------------------------------------------------------------------
void HierclustXmlWriter::WriteParentId(std::ofstream& outfile, 
                                       const int parent_id)
{
    outfile << S8 << "<parent_id>" << parent_id << "</parent_id>" << endl;
}
/*
//-----------------------------------------------------------------------------
void HierclustXmlWriter::WriteParentLabel(std::ofstream& outfile, 
                                          const int parent_label)
{
    outfile << S8 << "<parent_label>" << parent_label 
            << "</parent_label>" << endl;
}
*/
//-----------------------------------------------------------------------------
void HierclustXmlWriter::WriteLeftChild(std::ofstream& outfile, 
                                        const bool is_left_child,
                                        const int lc_id)
{
    outfile << S8 << "<left_child>" << (is_left_child ? "true" : "false") 
            << "</left_child>" << endl;

    outfile << S8 << "<left_child_id>" << lc_id
            << "</left_child_id>" << endl;
}

//-----------------------------------------------------------------------------
void HierclustXmlWriter::WriteRightChild(std::ofstream& outfile, 
                                         const int rc_id)
{
    outfile << S8 << "<right_child_id>" << rc_id
            << "</right_child_id>" << endl;
}

//-----------------------------------------------------------------------------
void HierclustXmlWriter::WriteDocCount(std::ofstream& outfile, 
                                       const int count)
{
    outfile << S8 << "<doc_count>" << count << "</doc_count>" << endl;
}

//-----------------------------------------------------------------------------
void HierclustXmlWriter::WriteTopTerms(std::ofstream& outfile, 
                                       const std::vector<int>& term_indices,
                                       const std::vector<std::string>& dictionary)
{
    outfile << S8 << "<top_terms>" << endl;
    for (size_t q=0; q != term_indices.size(); ++q)
    {
        int index = term_indices[q];
        outfile << S12 << "<term name=\"" << dictionary[index] << "\"/>" << endl;
    }
    outfile << S8 << "</top_terms>" << endl;
}

//-----------------------------------------------------------------------------
void HierclustXmlWriter::WriteNodeEnd(std::ofstream& outfile)
{
    outfile << S4 << "</node>" << endl;
}

//-----------------------------------------------------------------------------
void HierclustXmlWriter::WriteFooter(std::ofstream& outfile)
{
    outfile << "</DataSet>" << endl;
}
