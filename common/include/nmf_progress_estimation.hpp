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

#include "progress_estimator_generic.hpp"

//-----------------------------------------------------------------------------
template <typename T>
inline 
void ReportProgress(const int iter_count, const T progress_metric)
{
    const int START_COUNT = 10;
    const int SKIP_COUNT  = 10;

    if ( (iter_count < START_COUNT) ||      // show the first few iterations
         (0 == (iter_count % SKIP_COUNT)))  // then periodically thereafter
    {
        std::cout << iter_count << ":\tprogress metric:\t" 
                  << progress_metric << std::endl;
    }
}


