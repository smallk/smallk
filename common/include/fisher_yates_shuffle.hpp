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

// Randomly rearrange the elements of an array.

#include <algorithm>
#include "random.hpp"

template <typename RndIt>
void FisherYatesShuffle(RndIt start, RndIt end, Random& rng)
{
    int n = end - start;
    for (int i=n-1; i>=1; --i)
    {
        int j = rng.RandomRangeInt(0, i);
        std::swap( *(start+i), *(start+j) );
    }
}

