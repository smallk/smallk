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

struct Node
{
    Node() 
    {
        left_child_ = -1;
        right_child_ = -1;
    };

    void Init(bool isroot, 
              const int index, 
              const int item_count);

    void Init(bool isroot, 
              const int index, 
              const int item_count, 
              const int k, 
              const int voc_index);

    bool isroot_;
    int index_;
    int item_count_;
    int k_;
    int voc_index_;
    int left_child_;  // node index of left child, -1 if no left child
    int right_child_; // node index of right child, -1 if no right child
};

int Insert(Node* tree,
           const int root_index,
           const int old_node_index,
           const int new_node_index,
           const int new_index);

