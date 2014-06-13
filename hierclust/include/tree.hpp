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

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include "file_format.hpp"
#include "hierclust_writer.hpp"

class TreeNode
{
public:

    int label;
    int index;
    bool is_left_child;
    int doc_count;
    int parent_label;
    int parent_index;
    int left_child_label;  // in the next level
    int right_child_label; // in the next level

    // indices of top-ranked terms
    std::vector<int> term_indices;
};

class Tree
{
public:

    void Clear();

    // return the index of the node on the given level with the given label
    bool NodeIndex(const int level, const int label, int& index);

    void InsertNode(const int level, const TreeNode node)
    {
        tree_.insert(std::make_pair(level, node));
    }

    bool UpdateChildLabels(const int level, const int label,
                           const int new_lc_label, const int new_rc_label);

    bool Update(const int level, const int label, 
                const int new_lc_label, const int new_rc_label,
                const std::vector<int>& top_term_indices);

    // incrementally update the term count histogram for leaf nodes
    void TermDistribution(std::map<int, int>& histogram);

    void Print(const std::vector<std::string>& dictionary);
//    bool WriteXml(const std::string& filepath,
//                  const std::vector<std::string>& dictionary);

    bool Write(const std::string& filepath,
               const FileFormat& format, 
               const std::vector<std::string>& dictionary);

private:

    // multimapmap<level, Treenode>
    std::multimap<int, TreeNode> tree_;

    std::multimap<int, TreeNode>::iterator Find(const int level, const int label);

    // void Write(std::ofstream& outfile, 
    //            const std::vector<std::string>& dictionary,
    //            std::multimap<int, TreeNode>::iterator& it);

    void WriteNodes(IHierclustWriter* writer,
                    std::ofstream& outfile,
                    const std::vector<std::string>& dictionary,
                    std::multimap<int, TreeNode>::iterator& it);

    bool IsLeafNode(const TreeNode& node)
    {
        return ( (-1 == node.left_child_label) && 
                 (-1 == node.right_child_label));
    }
};

