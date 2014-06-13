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

//-----------------------------------------------------------------------------
template <typename T, typename RndIt>
void MeanStdDev(RndIt begin, RndIt end, T& mean, T& stddev)
{
    // online algorithm from 
    // http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance

    int n  = 0;
    mean = T(0);
    T var  = T(0);
    for (auto it=begin; it != end; ++it)
    {
        n += 1;
        T x = *it;
        T delta = x - mean;
        mean += delta/n;
        var += delta*(x-mean);
    }

    var /= (n-1);
    stddev = sqrt(var);
}

