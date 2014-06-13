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

#include <map>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "utils.hpp"
#include "assignments.hpp"
#include "flatclust_writer.hpp"
#include "flat_clust_output.hpp"

using std::cout;
using std::cerr;
using std::endl;

//-----------------------------------------------------------------------------
IFlatclustWriter* CreateFlatclustWriter(const FileFormat& format)
{
    IFlatclustWriter* writer = nullptr;

    switch (format)
    {
    case FileFormat::XML:
        writer = new FlatclustXmlWriter();
        break;
    case FileFormat::JSON:
        writer = new FlatclustJsonWriter();
        break;
    default:
        throw std::logic_error("CreateFlatclustWriter: invalid file format");
        break;
    }

    return writer;
}

//-----------------------------------------------------------------------------
void FlatClustWriteResults(const std::string& assignfilepath,
                           const std::string& resultfilepath,
                           const std::vector<int>& assignments,
                           const std::vector<std::string>& dictionary,
                           const std::vector<int>& term_indices,
                           const FileFormat format,
                           const unsigned int maxterms,
                           const unsigned int num_docs,
                           const unsigned int num_clusters)
{
    if (term_indices.size() < num_clusters*maxterms)
        throw std::logic_error("FlatClustWriteResults: term count is invalid");

    // count docs for each node
    std::map<int, int> doc_counts;
    for (unsigned int i=0; i != assignments.size(); ++i)
    {
        int cluster_id = assignments[i];
        auto it = doc_counts.find(cluster_id);
        if (doc_counts.end() != it)
        {
            // already in map
            it->second += 1;
        }
        else
        {
            // add to map
            doc_counts.insert(std::make_pair(cluster_id, 1));
        }
    }

    if (doc_counts.size() != num_clusters)
    {
        cout << "Warning: only " << doc_counts.size() 
             << " clusters received an assignment." << endl << endl;
    }

    if (!WriteAssignmentsFile(assignments, assignfilepath))
        cerr << "\terror writing flat assignments file" << endl;
    
    IFlatclustWriter* writer = CreateFlatclustWriter(format);
    if (nullptr == writer)
        throw std::logic_error("FlatClustWriteResults: invalid output format");

    std::ofstream outfile(resultfilepath);
    if (!outfile)
    {
        delete writer;
        cerr << "FlatClustWriteResults: could not open output file " 
             << resultfilepath << endl;
        return;
    }

    writer->WriteHeader(outfile, num_docs);

    for (unsigned int i=0; i<num_clusters; ++i)
    {
        writer->WriteNodeBegin(outfile, i);

        auto it = doc_counts.find(i);
        if (doc_counts.end() == it)
        {
            writer->WriteDocCount(outfile, 0u);
        }
        else
        {
            writer->WriteDocCount(outfile, doc_counts[i]);
            writer->WriteTopTerms(outfile, i*maxterms, maxterms, 
                                  term_indices, dictionary);
        }

        writer->WriteNodeEnd(outfile);
    }

    writer->WriteFooter(outfile);

    delete writer;
    outfile.close();
}

//-----------------------------------------------------------------------------
void FlatClustWriteResults(const std::string& outdir,
                           const std::vector<int>& assignments,
                           const std::vector<std::string>& dictionary,
                           const std::vector<int>& term_indices,
                           const FileFormat format,
                           const unsigned int maxterms,
                           const unsigned int num_docs,
                           const unsigned int num_clusters)
{
    std::string filename;
    std::string output_dir = EnsureTrailingPathSep(outdir);
    
    // construct flat assignments file name
    std::ostringstream assign_name;
    assign_name << "assignments_flat_" << num_clusters;
    filename = AppendExtension(assign_name.str(), FileFormat::CSV);
    std::string flat_assignfile = output_dir + filename;
    
    // construct flat result file name
    std::ostringstream result_name;
    result_name << "clusters_" << num_clusters;
    filename = AppendExtension(result_name.str(), format);
    std::string flat_resultfile = output_dir + filename;

    FlatClustWriteResults(flat_assignfile, flat_resultfile, assignments, 
                          dictionary, term_indices, format, maxterms, 
                          num_docs, num_clusters);
}
