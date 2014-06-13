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

#include "node.hpp"

//-----------------------------------------------------------------------------
void Node::Init(bool isroot, 
                const int index, 
                const int item_count)
{
    Init(isroot, index, item_count, -1, -1);
}

//-----------------------------------------------------------------------------
void Node::Init(bool isroot, 
                const int index, 
                const int item_count, 
                const int k, 
                const int voc_index)
{
    isroot_      = isroot;
    index_       = index;
    item_count_  = item_count;
    k_           = k;
    voc_index_   = voc_index;
    left_child_  = -1;
    right_child_ = -1;
}

//-----------------------------------------------------------------------------
int Insert(Node* tree,
           const int root_index,
           const int old_node_index,
           const int new_node_index,
           const int new_index)
{
    int split_node_index = root_index;
    int mask = 1;
    int level = 0;

    while (-1 != tree[split_node_index].left_child_)
    {
        if (mask & new_index)
            split_node_index = tree[split_node_index].right_child_;
        else
            split_node_index = tree[split_node_index].left_child_;
        mask <<= 1;
        level += 1;
    }

    tree[split_node_index].left_child_ = old_node_index;
    tree[split_node_index].right_child_ = new_node_index;
    return level;
}

