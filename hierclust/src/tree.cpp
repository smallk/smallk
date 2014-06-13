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

#include <cassert>
#include <iostream>
#include <fstream>
#include "tree.hpp"
#include "hierclust_writer_factory.hpp"

using std::cout;
using std::cerr;
using std::endl;

//-----------------------------------------------------------------------------
void Tree::Clear()
{
    tree_.clear();
}

//-----------------------------------------------------------------------------
bool Tree::NodeIndex(const int level, const int label, int& index)
{
    for (auto it = tree_.lower_bound(level);
         it != tree_.upper_bound(level);
         ++it)
    {
        if (label == it->second.label)
        {
            index = it->second.index;
            return true;
        }
    }

    return false;
}

//-----------------------------------------------------------------------------
bool Tree::UpdateChildLabels(const int level, const int label,
                             const int new_lc_label, const int new_rc_label)
{
    for (auto it = tree_.lower_bound(level); 
         it != tree_.upper_bound(level); 
         ++it)
    {
        if (label == it->second.label)
        {
            it->second.left_child_label = new_lc_label;
            it->second.right_child_label = new_rc_label;
            return true;
        }
    }

    return false;
}

//-----------------------------------------------------------------------------
bool Tree::Update(const int level, const int label, 
                  const int new_lc_label, const int new_rc_label,
                  const std::vector<int>& top_term_indices)
{
    for (auto it = tree_.lower_bound(level); 
         it != tree_.upper_bound(level); 
         ++it)
    {
        if (label == it->second.label)
        {
            it->second.left_child_label = new_lc_label;
            it->second.right_child_label = new_rc_label;
            it->second.term_indices = top_term_indices;
            return true;
        }
    }

    return false;
}
    
//-----------------------------------------------------------------------------
void Tree::Print(const std::vector<std::string>& dictionary)
{
    size_t dict_size = dictionary.size();

    for (auto it=tree_.begin(); it != tree_.end(); ++it)
    {
        int level = it->first;
        std::cout << "\t[" << level << "] label: " 
                  << it->second.label << std::endl;
        std::cout << "\t[" << level << "] index: " 
                  << it->second.index << std::endl;
        std::cout << "\t[" << level << "] is_left_child: " 
                  << (it->second.is_left_child ? "true" : "false") << std::endl;
        std::cout << "\t[" << level << "] doc_count: " 
                  << it->second.doc_count << std::endl;
        std::cout << "\t[" << level << "] parent_label: " 
                  << it->second.parent_label << std::endl;
        std::cout << "\t[" << level << "] parent_index: " 
                  << it->second.parent_index << std::endl;
        std::cout << "\t[" << level << "] left_child_label: " 
                  << it->second.left_child_label << std::endl;
        std::cout << "\t[" << level << "] right_child_label: " 
                  << it->second.right_child_label << std::endl;
        for (size_t q=0; q != it->second.term_indices.size(); ++q)
        {
            size_t term_index = it->second.term_indices[q];
            assert(term_index >= 0);
            assert(term_index < dict_size);
            std::string term = dictionary[term_index];
            std::cout << "\t\t" << term << std::endl;
        }
    }
}

//-----------------------------------------------------------------------------
std::multimap<int, TreeNode>::iterator 
Tree::Find(const int level, const int label)
{
    // return an iterator pointing to the node at the given level and label
    std::multimap<int, TreeNode>::iterator it = tree_.end();

    
    for (it = tree_.lower_bound(level);
         it != tree_.upper_bound(level);
         ++it)
    {
        if (label == it->second.label)
        {
            break;
        }
    }

    return it;
}

//-----------------------------------------------------------------------------
void Tree::TermDistribution(std::map<int, int>& histogram)
{
    for (auto it = tree_.begin(); it != tree_.end(); ++it)
    {
        if ( (-1 != it->second.left_child_label) &&
             (-1 != it->second.right_child_label))
        {
            // not a leaf node
            continue;
        }

        for (auto itv = it->second.term_indices.begin(); 
             itv != it->second.term_indices.end(); ++itv)
        {
            int term_index = *itv;
            auto itm = histogram.find(term_index);
            if (histogram.end() == itm)
            {
                // new term; add to histogram with count of 1
                histogram.insert(std::make_pair(term_index, 1));
            }
            else
            {
                // increment the count for this term
                itm->second += 1;
            }
        }
    }
}

//-----------------------------------------------------------------------------
bool Tree::Write(const std::string& filepath,
                 const FileFormat& format,
                 const std::vector<std::string>& dictionary)
{
    IHierclustWriter* writer = CreateHierclustWriter(format);
    if (nullptr == writer)
        throw std::logic_error("Tree::Write: invalid output format.");

    std::ofstream outfile(filepath);
    if (!outfile)
    {
        delete writer;
        cerr << "Tree::Write: could not open output file " << filepath << endl;
        return false;
    }

    // count the total number of documents in the leaf nodes
    int doc_count = 0;
    for (auto it=tree_.begin(); it != tree_.end(); ++it)
    {
        if (IsLeafNode(it->second))
            doc_count += it->second.doc_count;
    }

    std::multimap<int, TreeNode>::iterator it = tree_.begin();

    writer->WriteHeader(outfile, doc_count);
    WriteNodes(writer, outfile, dictionary, it);
    writer->WriteFooter(outfile);
    outfile.close();
    delete writer;
    return true;
}

//-----------------------------------------------------------------------------
void Tree::WriteNodes(IHierclustWriter* writer,
                      std::ofstream& outfile,
                      const std::vector<std::string>& dictionary,
                      std::multimap<int, TreeNode>::iterator& it)
{
    if (it == tree_.end())
        return;

    // skip the root node
    if (-1 != it->second.parent_index)
    {
        writer->WriteNodeBegin(outfile, it->second.index);
        writer->WriteLevel(outfile, it->first);
        writer->WriteLabel(outfile, it->second.label);
        writer->WriteParentId(outfile, it->second.parent_index);
        writer->WriteParentLabel(outfile, it->second.parent_label);
        writer->WriteLeftChild(outfile, it->second.is_left_child, it->second.left_child_label);
        writer->WriteRightChild(outfile, it->second.right_child_label);
        writer->WriteDocCount(outfile, it->second.doc_count);
        writer->WriteTopTerms(outfile, it->second.term_indices, dictionary);
        writer->WriteNodeEnd(outfile);
    }

    int child_level       = it->first + 1;
    int left_child_label  = it->second.left_child_label;
    int right_child_label = it->second.right_child_label;

    // write XML for the left child
    if (-1 != left_child_label)
    {
        it = Find(child_level, left_child_label);
        assert(it != tree_.end());
        WriteNodes(writer, outfile, dictionary, it);
    }

    if (-1 != right_child_label)
    {
        it = Find(child_level, right_child_label);
        assert(it != tree_.end());
        WriteNodes(writer, outfile, dictionary, it);
    }
}

