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
#include <set>
#include <limits>
#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "file_format.hpp"
#include "dense_matrix.hpp"
#include "vector_utils.hpp"
#include "hierclust_writer.hpp"

template <typename T>
class TreeNode
{
public:

    TreeNode()  {is_valid = false;}
    
    T            priority;
    unsigned int parent_index;
    unsigned int left_child_index;
    unsigned int right_child_index;
    bool         is_valid;
    bool         is_left_child;
    
    // an mx1 vector, where m is the number of terms in the dictionary
    DenseMatrix<T> topic_vector;

    // indices of top-ranked terms
    std::vector<int> term_indices;

    // the set of document indices assigned to this node
    std::vector<unsigned int> docs;
};

//-----------------------------------------------------------------------------
template <typename T>
class Tree
{
public:
    Tree() {}

    // indices of the most recently created left and right child nodes
    unsigned int LeftChildIndex()  {return index0_;}
    unsigned int RightChildIndex() {return index1_;}        

    // topic vectors of the most recently created left and right child nodes
    DenseMatrix<T>& LeftChildTopicVector()  {return nodes_[index0_].topic_vector;}
    DenseMatrix<T>& RightChildTopicVector() {return nodes_[index1_].topic_vector;}

    // doc set of the most recently created left and right child nodes
    std::vector<unsigned int>& LeftChildDocs()  {return nodes_[index0_].docs;}
    std::vector<unsigned int>& RightChildDocs() {return nodes_[index1_].docs;}

    std::vector<unsigned int>& Outliers()    {return outliers_;}    
    std::vector<unsigned int>& Assignments() {return assignments_;}
    
    void Init(const unsigned int num_clusters,
              const unsigned int node_count,
              const unsigned int term_count,
              const unsigned int doc_count);

    // return the min and max leaf node priorities, along with the index
    // of the max priority leaf node
    void MinMaxLeafPriorities(T& min_priority,
                              T& max_priority,
                              unsigned int& max_priority_index);

    // split the root node
    void SplitRoot(const DenseMatrix<T>* W, const DenseMatrix<T>* H);
    
    // split the node at the specified index
    void Split(const unsigned int split_index,
               const DenseMatrix<T>* W,
               const DenseMatrix<T>* H);

    void SetNodePriority(const unsigned int node_index, const T priority);

    void ComputeTopTerms(const unsigned int max_terms);
    void ComputeAssignments();

    // Build an mxk W initializer matrix for flat clustering;
    // use the leaf node topic vectors generated from a hierclust run
    bool FlatclustInitW(DenseMatrix<T>& Winit,
                        const unsigned int m,  // height (document count)
                        const unsigned int k); // width  (num clusters)

    bool WriteAssignments(const std::string& filepath);
    
    bool WriteTree(IHierclustWriter* writer,
                   const std::string& filepath,
                   const std::vector<std::string>& dictionary);
    
    void Print();
    
private:

    static const unsigned int NONE;
    static const T            MAX_PRIORITY;
    
    // root node is not stored in the tree structure
    TreeNode<T> root_;

    // The tree structure is encoded in this vector.  The children of the root
    // node have indices 0 and 1, so the 0th entry of this vector contains
    // info for a CHILD of the root node.  There is no entry for the root node.
    std::vector< TreeNode<T> > nodes_;

    // leaf node indicators
    std::vector<bool> is_leaf_;

    // number of nodes currently in the tree (not necessarily contiguous
    // in the nodes_ array); the root node is not included in this count
    unsigned int active_nodes_;

    // indices of the most recent left and right child nodes
    unsigned int index0_, index1_;

    // original number of docs in the root node
    unsigned int total_docs_;

    // sum of all docs in the leaf nodes
    // (leaf_doc_count_ + outliers.size() == total_docs_)
    unsigned int leaf_doc_count_;
    
    // document indices not assigned to any cluster
    std::vector<unsigned int> outliers_;

    // assignments of documents to clusters
    // document j is assigned to the cluster at index 'assignments_[j]'
    std::vector<unsigned int> assignments_;
    
    void PartitionDocs(const DenseMatrix<T>* H);
    void UpdateTopicVectors(const DenseMatrix<T>* W);
    
    void WriteNodes(IHierclustWriter* writer,
                    std::ofstream& outfile,
                    const std::vector<std::string>& dictionary);
};

//-----------------------------------------------------------------------------
template <typename T>
void Tree<T>::Init(const unsigned int num_clusters,
                   const unsigned int node_count,
                   const unsigned int term_count,
                   const unsigned int doc_count)
{
    total_docs_ = doc_count;
    
    // initialize the root node
    root_.parent_index      = Tree<T>::NONE;
    root_.left_child_index  = Tree<T>::NONE;
    root_.right_child_index = Tree<T>::NONE;
    root_.priority          = Tree<T>::MAX_PRIORITY;
    root_.is_valid          = true;
    
    // allocate all remaining nodes
    nodes_.resize(node_count);

    // allocate the topic vector matrix (mx1) in each node
    for (unsigned int q=0; q<node_count; ++q)
        nodes_[q].topic_vector.Resize(term_count, 1);

    // setup the leaf node indicators
    is_leaf_.resize(node_count);
    for (unsigned int q=0; q<is_leaf_.size(); ++q)
        is_leaf_[q] = false;

    active_nodes_ = 0;
}

//-----------------------------------------------------------------------------
template <typename T>
void Tree<T>::MinMaxLeafPriorities(T& min_priority,
                                   T& max_priority,
                                   unsigned int& max_priority_index)
{
    // Find the min and max priority leaf nodes; return these priorities
    // along with the index of the max priority leaf node.
    
    min_priority = std::numeric_limits<T>::max();
    max_priority = std::numeric_limits<T>::lowest();

    for (unsigned int q=0; q<is_leaf_.size(); ++q)
    {
        // skip if this is not a leaf node
        if (!is_leaf_[q])
            continue;
        
        T node_priority = nodes_[q].priority;
        if ( (node_priority > T(0)) && (node_priority < min_priority))
            min_priority = node_priority;
        
        if (node_priority > max_priority)
        {
            max_priority = node_priority;
            max_priority_index = q;
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void Tree<T>::SplitRoot(const DenseMatrix<T>* W,
                        const DenseMatrix<T>* H)
{
    // children of the root node are at indices 0 and 1
    index0_ = 0;
    index1_ = 1;
    
    root_.left_child_index  = index0_;
    root_.right_child_index = index1_;

    nodes_[0].parent_index      = Tree<T>::NONE;
    nodes_[0].left_child_index  = Tree<T>::NONE;
    nodes_[0].right_child_index = Tree<T>::NONE;
    nodes_[0].is_valid          = true;
    nodes_[0].is_left_child     = true;
    is_leaf_[0] = true;

    nodes_[1].parent_index      = Tree<T>::NONE;
    nodes_[1].left_child_index  = Tree<T>::NONE;
    nodes_[1].right_child_index = Tree<T>::NONE;
    nodes_[1].is_valid          = true;
    nodes_[1].is_left_child     = false;
    is_leaf_[1] = true;

    active_nodes_ += 2;

    // Partition the documents assigned to the root node.  The root node
    // contains all documents, so the 'root_.docs' set does not need to be
    // explicitly created.
    
    unsigned int h_width = H->Width();
    for (unsigned int c=0; c<h_width; ++c)
    {
        if (H->Get(0, c) > H->Get(1, c))
            nodes_[0].docs.push_back(c);
        else
            nodes_[1].docs.push_back(c);
    }

    UpdateTopicVectors(W);
}

//-----------------------------------------------------------------------------
template <typename T>
void Tree<T>::Split(const unsigned int node_index,
                    const DenseMatrix<T>* W,
                    const DenseMatrix<T>* H)
{    
    assert(node_index < nodes_.size());

    // use two more nodes
    index0_ = active_nodes_ + 0;
    index1_ = active_nodes_ + 1;
    active_nodes_ += 2;

    assert(index0_ < nodes_.size());
    assert(index1_ < nodes_.size());
    
    // assign child indices to the split node
    nodes_[node_index].left_child_index  = index0_;
    nodes_[node_index].right_child_index = index1_;

    // the split node is no longer a leaf
    is_leaf_[node_index] = false;

    // set the left child's parent and make it a leaf
    nodes_[index0_].parent_index      = node_index;
    nodes_[index0_].left_child_index  = Tree<T>::NONE;
    nodes_[index0_].right_child_index = Tree<T>::NONE;
    nodes_[index0_].is_valid          = true;
    nodes_[index0_].is_left_child     = true;
    is_leaf_[index0_] = true;
    
    // set the right child's parent and make it a leaf
    nodes_[index1_].parent_index      = node_index;
    nodes_[index1_].left_child_index  = Tree<T>::NONE;
    nodes_[index1_].right_child_index = Tree<T>::NONE;
    nodes_[index1_].is_valid          = true;
    nodes_[index1_].is_left_child     = false;
    is_leaf_[index1_] = true;

    // Partition the documents assigned to the node being split.
    unsigned int h_width = H->Width();
    std::vector<unsigned int>& source_docs = nodes_[node_index].docs;
    
    for (unsigned int c=0; c<h_width; ++c)
    {
        if (H->Get(0, c) > H->Get(1, c))
            nodes_[index0_].docs.push_back(source_docs[c]);
        else
            nodes_[index1_].docs.push_back(source_docs[c]);
    }

    UpdateTopicVectors(W);
}

//-----------------------------------------------------------------------------
template <typename T>
void Tree<T>::SetNodePriority(const unsigned int node_index,
                              const T priority)
{
    assert(node_index < nodes_.size());
    assert(is_leaf_[node_index]);

    nodes_[node_index].priority = priority;
}

//-----------------------------------------------------------------------------
template <typename T>
void Tree<T>::UpdateTopicVectors(const DenseMatrix<T>* W)
{
    // The left child of a newly-split node gets the left column of W, and
    // the right child gets the right column of W.

    int m = W->Height();
    DenseMatrix<T> source_col, dest_col;

    // W(:, 0) -> topic_vectors[index0_]
    LockedView(source_col, *W, 0, 0, m, 1);
    View(dest_col, nodes_[index0_].topic_vector, 0, 0, m, 1);
    Copy(source_col, dest_col);

    // W(:, 1) -> topic_vectors[index1_]
    LockedView(source_col, *W, 0, 1, m, 1);
    View(dest_col, nodes_[index1_].topic_vector, 0, 0, m, 1);
    Copy(source_col, dest_col);
}

//-----------------------------------------------------------------------------
template <typename T>
void Tree<T>::ComputeTopTerms(const unsigned int max_terms)
{
    // nodes[0] always exists
    std::vector<int> sort_indices(nodes_[0].topic_vector.Height());
    
    // compute 'max_terms' top-ranked terms for each node
    for (unsigned int q=0; q<nodes_.size(); ++q)
    {
        if (!nodes_[q].is_valid)
            continue;
        
        nodes_[q].term_indices.resize(max_terms);
        
        TopTerms(max_terms,
                 nodes_[q].topic_vector,
                 sort_indices,
                 nodes_[q].term_indices);
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void Tree<T>::ComputeAssignments()
{
    outliers_.clear();
    assignments_.resize(total_docs_);
    std::fill(assignments_.begin(), assignments_.end(), Tree<T>::NONE);

    leaf_doc_count_ = 0;
    unsigned int doc_count;
    for (unsigned int q=0; q<nodes_.size(); ++q)
    {
        if (!is_leaf_[q])
            continue;

        doc_count = nodes_[q].docs.size();
        leaf_doc_count_ += doc_count;

        // for all doc indices in this cluster...
        for (unsigned int k=0; k<doc_count; ++k)
        {
            unsigned int doc_index = nodes_[q].docs[k];

            // the document at index 'doc_index' belongs to cluster q
            assignments_[doc_index] = q;
        }
    }

    // any unassigned docs are outliers; save doc indices in outliers_
    for (unsigned int q=0; q<assignments_.size(); ++q)
    {
        if (Tree<T>::NONE == assignments_[q])
            outliers_.push_back(q);
    }

    // sum of leaf_docs + outliers == original doc count
    assert( (leaf_doc_count_ + outliers_.size()) == total_docs_);
}

//-----------------------------------------------------------------------------
template <typename T>
bool Tree<T>::FlatclustInitW(DenseMatrix<T>& Winit,
                             const unsigned int m, const unsigned int k)
{
    // Collect the topic vectors for all leaf nodes into the 'Winit' matrix.
    
    // count the leaf nodes
    unsigned int leaf_count = 0;
    for (unsigned int q=0; q<nodes_.size(); ++q)
    {
        if (is_leaf_[q])
            ++leaf_count;
    }

    if (k != leaf_count)
    {
        std::cerr << "Insufficient number of leaf nodes for flat clustering."
                  << std::endl;
        return false;
    }
    
    if (m != static_cast<unsigned int>(nodes_[0].topic_vector.Height()) )
    {
        std::cerr << "Invalid W matrix height for flat clustering."
                  << std::endl;
        return false;
    }

    unsigned int c = 0;
    DenseMatrix<T> src_view, dest_view;
    for (unsigned int q=0; q<nodes_.size(); ++q)
    {
        if (!is_leaf_[q])
            continue;

        // view of the topic vector
        LockedView(src_view, nodes_[q].topic_vector, 0, 0, m, 1);

        // view of the dest column
        View(dest_view, Winit, 0, c, m, 1);

        Copy(src_view, dest_view);
        ++c;
    }

    // must have copied k columns
    return (k == c);
}

//-----------------------------------------------------------------------------
template <typename T>
bool Tree<T>::WriteAssignments(const std::string& filepath)
{
    // Writes a two-section CSV file.  The first set of CSV values records the
    // assignment of documents to clusters.  For a set of n documents, n values
    // will be written.  Each value is either the node index of the leaf node
    // that contains the document, or a '-1' indicating an outlier.
    
    std::ofstream outfile(filepath);
    if (!outfile)
    {
        std::cerr << "Tree::WriteAssignments: could not open output file "
                  << filepath << std::endl;
        return false;
    }

    // write cluster assignments or -1 for outliers
    outfile << assignments_[0];
    for (unsigned int q=1; q<assignments_.size(); ++q)
    {
        unsigned int label = assignments_[q];
        outfile << ",";
        if (Tree<T>::NONE == label)
            outfile << static_cast<int>(-1); // force -1 to appear
        else
            outfile << label;
    }
    outfile << std::endl;

    // write a blank line to separate clusters and outliers
    outfile << std::endl;

    // write outliers
    if (outliers_.size() > 0)
    {
        outfile << outliers_[0];
        for (unsigned int q=1; q<outliers_.size(); ++q)
            outfile << ',' << outliers_[q];
        outfile << std::endl;
    }
    
    outfile.close();
    return true;
}

//-----------------------------------------------------------------------------
template <typename T>
bool Tree<T>::WriteTree(IHierclustWriter* writer,
                        const std::string& filepath,
                        const std::vector<std::string>& dictionary)
{
    std::ofstream outfile(filepath);
    if (!outfile)
    {
        std::cerr << "Tree::Write: could not open output file "
                  << filepath << std::endl;
        return false;
    }

    writer->WriteHeader(outfile, leaf_doc_count_);
    this->WriteNodes(writer, outfile, dictionary);
    writer->WriteFooter(outfile);
    outfile.close();
    return true;
}

//-----------------------------------------------------------------------------
template <typename T>
void Tree<T>::WriteNodes(IHierclustWriter* writer,
                         std::ofstream& outfile,
                         const std::vector<std::string>& dictionary)
{
    for (unsigned int q=0; q<nodes_.size(); ++q)
    {
        // node id == node label == index in the nodes_ array
        writer->WriteNodeBegin(outfile, q);
        writer->WriteParentId(outfile, nodes_[q].parent_index);
        writer->WriteLeftChild(outfile, nodes_[q].is_left_child, nodes_[q].left_child_index);
        writer->WriteRightChild(outfile, nodes_[q].right_child_index);
        writer->WriteDocCount(outfile, nodes_[q].docs.size());
        writer->WriteTopTerms(outfile, nodes_[q].term_indices, dictionary);
        writer->WriteNodeEnd(outfile);
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void Tree<T>::Print()
{
    std::cout << "\n\ncluster sizes: " << std::endl;
    std::cout << "\t";
    for (unsigned int q=0; q<nodes_.size(); ++q)
        std::cout << nodes_[q].docs.size() << "  ";
    std::cout << std::endl;

    std::cout << "cluster sums: " << std::endl;
    std::cout << "\t";
    for (unsigned int q=0; q<nodes_.size(); ++q)
        std::cout << Sum(nodes_[q].docs) << "  ";
    std::cout << std::endl;
    
    std::cout << "Node tree: " << std::endl;
    for (unsigned int q=0; q<nodes_.size(); ++q)
        std::cout << (Tree<T>::NONE == nodes_[q].left_child_index ? 0 : nodes_[q].left_child_index) << " ";
    std::cout << std::endl;

    for (unsigned int q=0; q<nodes_.size(); ++q)
        std::cout << (Tree<T>::NONE == nodes_[q].right_child_index ? 0 : nodes_[q].right_child_index)  << " ";
    std::cout << std::endl;

    std::cout << "cluster priorities: " << std::endl;
    for (unsigned int q=0; q<nodes_.size(); ++q)
        std::cout << nodes_[q].priority << "  ";
    std::cout << std::endl;

    int leaf_count = 0;
    std::cout << "leaf nodes: " << std::endl;
    for (unsigned int q=0; q<is_leaf_.size(); ++q)
    {
        if (is_leaf_[q])
        {
            ++leaf_count;
            std::cout << q << ", ";
        }
    }
    std::cout << std::endl;
    std::cout << "leaf node count: " << leaf_count << std::endl;

    // std::cout << "topic vectors: " << std::endl;
    // for (unsigned int q=0; q<nodes_.size(); ++q)
    // {
    //     std::cout << "Node " << q << ": " << std::endl;
    //     for (unsigned int j=0; j<5; ++j)
    //     {
    //         std::cout << "\t" << nodes_[q].topic_vector.Get(j, 0) << std::endl;
    //     }
    // }

    std::cout << "top term indices: " << std::endl;
    for (unsigned int q=0; q<nodes_.size(); ++q)
    {
        if (!is_leaf_[q])
            continue;

        std::cout << "Leaf node " << q << ": " << std::endl;
        for (unsigned int j=0; j<5; ++j)
            std::cout << "\t" << nodes_[q].term_indices[j] << std::endl;
    }

    std::cout << "Found " << outliers_.size() << " outliers." << std::endl;
}
