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

#include <iostream>
#include <fstream>
#include <cassert>
#include <set>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <cstdlib>
#include "postprocess.hpp"
#include "vector_utils.hpp"
#include "node.hpp"

using std::cout;
using std::cerr;
using std::endl;

//-----------------------------------------------------------------------------
void Setdiff1dMinIndex(const std::vector<int>& set1, 
                       const std::vector<int>& set2,
                       int& val)
{
    // Find the smallest element of set1 that is not in set2 and return it.
    // The python code does a full 'setdiff1d' operation, but that is not
    // necessary, since all we need is a single element.

    if (0 == set1.size())
        throw std::runtime_error("postprocess::Setdiff1dMinIndex: error - set1 is empty");

    // if set2 is empty, return the smallest entry in set 1
    if (0 == set2.size())
    {
        val = set1[0];
        return;
    }

    // search set2 for each element of set 1
    for (size_t i=0; i != set1.size(); ++i)
    {
        int elt = set1[i];        
        bool in_set_2 = std::binary_search(set2.begin(), set2.end(), elt);
        if (!in_set_2)
        {
            val = elt;
            return;
        }
    }

    throw std::runtime_error("postprocess::Setdiff1dMinIndex: error - setdiff is the empty set");
}

//-----------------------------------------------------------------------------
void UniqueNonnegIndices(const std::vector<int>& cluster,
                         std::vector<int>& result)
{
    size_t i;
    result.clear();

    // insert the indices into a set; this also sorts the unique indices
    std::set<int> unique;
    for (i=0; i != cluster.size(); ++i)
    {
        int val = cluster[i];
        if (val >= 0)
            unique.insert(val);
    }

    i=0;
    result.resize(unique.size());
    for (auto it=unique.begin(); it != unique.end(); ++it, ++i)
        result[i] = *it;
}

//-----------------------------------------------------------------------------
void BuildIndexMap(std::map<int, int>& M, std::vector<int>& unique_indices)
{
    M.clear();
    for (size_t i=0, offset=0; i != unique_indices.size(); ++i, ++offset)
        M[unique_indices[i]] = offset;
}

//-----------------------------------------------------------------------------
int GetOldIndex(const int new_index)
{
    int mask = 1;
    while (mask <= new_index)
        mask <<= 1;
    return new_index - (mask >> 1);
}

//-----------------------------------------------------------------------------
bool SetClusterOffsets(const int cluster_index,
                       const std::map<int, int>& index_map,
                       const std::vector<int>& cluster,
                       std::vector<int>& cluster_mapped)
{
    cluster_mapped.clear();
    for (size_t q=0; q != cluster.size(); ++q)
    {
        if (-1 == cluster[q])
            continue;

        auto it = index_map.find(cluster[q]);
        if (it == index_map.end())
        {
            cerr << "\nPostprocess: unexpected document index " 
                 << cluster[q] << " in cluster " << cluster_index << endl;
            return false;
        }

        cluster_mapped.push_back(it->second);
    }
    
    return true;
}

//-----------------------------------------------------------------------------
void Write(std::ofstream& outfile,
           Node* tree, const Node& node,
           std::map<std::pair<int, int>, std::vector<std::string> >& word_map)
{
    char item_buf[64], name_buf[64];

    if (!node.isroot_)
    {
        sprintf(item_buf, "%d", node.item_count_);
        if (-1 != node.left_child_)
        {
            outfile << "  <left_child doc_count=\"" << item_buf << "\">" << endl;
        }
        else
        {
            sprintf(name_buf, "Label: %d", node.index_);
            outfile << "  <right_child doc_count=\"" << item_buf 
                    << "\" name=\"" << name_buf << "\">" << endl;
        }
        
        // print words here
        outfile << "    <topics>" << endl;
        auto p = std::make_pair(node.k_-2, node.voc_index_);
        std::vector<std::string>& words = word_map[p];
        for (auto it=words.begin(); it != words.end(); ++it)
            outfile << "      <topic name=\"" << (*it) << "\"/>" << endl;
        outfile << "    </topics>" << endl;

        if (-1 != node.left_child_)
            outfile << "  </left_child>" << endl;
        else
            outfile << "  </right_child>" << endl;
    }

    if (-1 != node.left_child_)
    {
        assert(-1 != node.right_child_);
        Write(outfile, tree, tree[node.left_child_], word_map);
        Write(outfile, tree, tree[node.right_child_], word_map);
    }
}

//-----------------------------------------------------------------------------
bool WriteXML(const std::string& filename, 
              const int doc_count,
              Node* tree,
              std::map<std::pair<int, int>, std::vector<std::string> >& word_map)
{
    std::ofstream outfile(filename);
    if (!outfile)
        return false;

    char buf[64];
    sprintf(buf, "%d", doc_count);
    outfile << "<?xml version=\"1.0\"?>" << endl;
    outfile << "<DataSet id=\"" << buf << "\">" << endl;

    Write(outfile, tree, tree[0], word_map);

    outfile << "</DataSet>" << endl;
    outfile.close();
    return true;
}

//-----------------------------------------------------------------------------
void TopRankedWords(Node* tree, 
                    const unsigned int node_count,
                    const int height, 
                    const int num_top_terms,
                    const std::vector<std::string>& dictionary,
                    const std::vector< std::vector<double> >& result_Ws,
                    std::map<std::pair<int, int>, std::vector<std::string> >& word_map)
{
    // Walk the tree and extract the dictionary info.  Build a map of the 
    // strings requried by each unique combination of k-value and column. 
    // The k-value and column index will be stored as a pair, which will serve
    // as the key to the dictionary.

    std::vector<int> indices(height);
    std::vector<std::string> words(num_top_terms);
    
    word_map.clear();
    for (unsigned int i=0; i<node_count; ++i)
    {
        // don't print the
        if (tree[i].isroot_)
            continue;

        int w_index = tree[i].k_ - 2;
        assert(w_index >= 0);

        // compute the width of this W matrix
        size_t num_elts = result_Ws[w_index].size();
        int w_width = num_elts / height;
        assert(w_width >= 2);

        int col = tree[i].voc_index_;
        assert(col >= 0);

        if (0 != (num_elts % height))
        {
            cerr << "TopRankedWords: error - matrix W[" << w_index
                 << "] has an invalid size." << endl;
            cerr << "num_elts: " << num_elts << ", height: " 
                 << height << ", w_width: " << w_width << endl;
            return;
        }
        if (col >= w_width)
        {
            cerr << "TopRankedWords: error - the column index for W[" << w_index
                 << "] is too large. " << endl;
            cerr << "W matrix width: " << w_width << ", col index: " << col << endl;
            return;
        }

        auto p = std::make_pair(w_index, col);
        auto it = word_map.find(p);
        if (it == word_map.end())
        {
            // add the dictionary words for this node

            // data buffer for the W matrix at index 'w_index'
            const double* wdata = &(result_Ws[w_index])[0];

            // data for this column begins at this offset
            size_t col_offset = col*height;

            // column data starts here
            const double* data = &wdata[col_offset];

            // setup the indices for the sort
            for (int q=0; q<height; ++q)
                indices[q] = q;

            // sort the indices in this column by comparing data elements;
            // sort in DECREASING order, so highest-ranked terms come first
            std::sort(&indices[0], &indices[0] + height,
                      [&data](int i1, int i2) {return data[i1] > data[i2];});

            // store the top words in the vector for this key
            word_map[p].clear();
            size_t max_terms = std::min(num_top_terms, height);
            for (size_t q=0; q<max_terms; ++q)
            {
                assert(indices[q] >= 0);
                assert(indices[q] < height);
                word_map[p].push_back(dictionary[indices[q]]);
            }
        }
    }
}

//-----------------------------------------------------------------------------
bool Postprocess(const std::vector<std::vector<int> >& clusters,
                 const std::vector<std::vector<double> >& result_Ws,
                 const std::vector<std::string>& dictionary,
                 const int height,
                 const int max_leaf_nodes,
                 const int num_top_terms,
                 const int term_count, const int doc_count)
{
    // this code performs (max_leaf_nodes-1) 'splits'
    
    std::map<int, int> index_map;
    std::vector<int> unique_indices, cluster_mapped;
    std::vector< std::vector<int> > counts(max_leaf_nodes-1);

    //for (size_t i=0; i != clusters.size(); ++i)
    for (int k=2; k<=max_leaf_nodes; ++k)
    {
        int i = k-2;
        const std::vector<int>& cluster = clusters[i];

        // find the sorted, unique, nonnegative indices in the current cluster
        UniqueNonnegIndices(cluster, unique_indices);

        // build a map of offsets for each unique index
        BuildIndexMap(index_map, unique_indices);

        // store the offset in unique_indices for all members of the cluster
        if (!SetClusterOffsets(i, index_map, cluster, cluster_mapped))
            return false;

        // Compute the number of times each unique index occurs in each cluster; 
        // these are the document counts.  Store a different vector for each
        // cluster.  Invariant: sum of the values in each list == total number
        // of documents in the data set (width of term-frequency matrix).
        counts[i].resize(unique_indices.size());
        Histogram(cluster_mapped.begin(), cluster_mapped.end(), counts[i].begin(), counts[i].end());

        // for the W-to-index conversion, only need to process the top terms
    }

    unique_indices.clear();

    Node node;
    int node_alloc_index = 0;
    int root_node_offset=0, old_node_offset, new_node_offset; 
    std::vector<Node> tree;
    std::vector<int> cur_indices = {0};
    int max_level = 0, new_index;
    int num_splits = static_cast<int>(counts.size());
    for (int one_split=0; one_split < num_splits; ++one_split)
    {
        const std::vector<int>& cluster = clusters[one_split];
        std::vector<int>& count = counts[one_split];

        UniqueNonnegIndices(cluster, unique_indices);

        // find smallest element of unique_indices not in cur_indices
        Setdiff1dMinIndex(unique_indices, cur_indices, new_index);
        BuildIndexMap(index_map, unique_indices);

        int old_index = GetOldIndex(new_index);
        int old_idx = index_map[old_index];
        int new_idx = index_map[new_index];

        //cout << "old_index [" << one_split << "]: " << old_index << endl;
        //cout << "old_idx   [" << one_split << "]: " << old_idx << endl;
        //cout << "new_idx   [" << one_split << "]: " << new_idx << endl;
        
        if (0 == one_split)
        {
            int test_count = count[old_idx] + count[new_idx];
            assert(doc_count == test_count);
            root_node_offset = node_alloc_index++;
            tree.push_back(node);
            tree[root_node_offset].Init(true, root_node_offset, count[old_idx] + count[new_idx]);
        }

        old_node_offset = node_alloc_index++;
        new_node_offset = node_alloc_index++;
        tree.push_back(node);
        tree.push_back(node);
        tree[old_node_offset].Init(false, old_index, count[old_idx], one_split+2, old_idx);
        tree[new_node_offset].Init(false, new_index, count[new_idx], one_split+2, new_idx);
        int level = Insert(&tree[0], root_node_offset, old_node_offset, new_node_offset, new_index);
        if (level > max_level)
            max_level = level;

        cur_indices = unique_indices;
        std::sort(cur_indices.begin(), cur_indices.end());
    }

    std::map<std::pair<int, int>, std::vector<std::string> > word_map;
    TopRankedWords(&tree[0], tree.size(), height, num_top_terms, dictionary, result_Ws, word_map);

    char buf[8];
    sprintf(buf, "%d", max_leaf_nodes);
    std::string filename = std::string("DataSet_") + 
                           std::string(buf) + std::string(".xml");;
    WriteXML(filename, doc_count, &tree[0], word_map);

    return true;
}
