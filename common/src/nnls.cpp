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

#include "nnls.hpp"

//-----------------------------------------------------------------------------
void UpdatePassiveSet(BitMatrix& passive_set,
                      const int PBAR,
                      const unsigned int n,
                      const BitMatrix& not_opt_cols,
                      const BitMatrix& nonopt_set,
                      const BitMatrix& infeas_set,
                      const std::vector<int>& not_good,
                      std::vector<int>& P, 
                      std::vector<int>& Ninf)
{
    std::vector<int> linear_indices;

    BitMatrix cols1 = not_opt_cols & (not_good < Ninf);
    BitMatrix cols2 = not_opt_cols & (not_good >= Ninf) & (P >= 1);
    BitMatrix tmp = not_opt_cols & (~cols1) & (~cols2);
    std::vector<int> cols3_ix;
    tmp.Find(cols3_ix);

    cols1.Find(linear_indices);
    if (!linear_indices.empty())
    {
        // P(cols1) = PBAR;
        IndexedAssign(P, cols1, PBAR);
        IndexedAssign(Ninf, cols1, not_good);

        BitMatrix bitmask = ColumnwiseAND(nonopt_set, cols1);
        passive_set.SetBits(bitmask);

        BitMatrix bitmask2 = ColumnwiseAND(infeas_set, cols1);
        passive_set.ClearBits(bitmask2);
    }
        
    cols2.Find(linear_indices);
    if (!linear_indices.empty())
    {
        // P(cols2) = P(cols2) - 1
        IndexedUpdate(P, cols2, -1);

        BitMatrix bitmask = ColumnwiseAND(nonopt_set, cols2);
        passive_set.SetBits(bitmask);

        BitMatrix bitmask2 = ColumnwiseAND(infeas_set, cols2);
        passive_set.ClearBits(bitmask2);
    }

    if (!cols3_ix.empty())
    {
        for (unsigned int i=0; i < cols3_ix.size(); ++i)
        {
            int ix = cols3_ix[i];
            unsigned int r1 = nonopt_set.MaxRowIndex(ix);
            unsigned int r2 = infeas_set.MaxRowIndex(ix);
            unsigned int row_to_change = std::max(r1, r2);
            passive_set.ToggleBit(row_to_change, ix);
        }
    }
}
