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

#include <thread>
#include <stdexcept>
#include <algorithm>
#include "openmp_pragma.hpp"

namespace ThreadUtils
{
    // The call to hardware_concurrency() could return 0, indicating that
    // the maximum number of supported threads is unknown.  In such cases
    // set the upper limit for the thread count to 2.  We really do not
    // want to run single-threaded.

    static unsigned int max_threads = 
        std::max(2u, std::thread::hardware_concurrency());
}

//-----------------------------------------------------------------------------
inline void SetMaxThreadCount(const unsigned int max_threads)
{    
    if (0 == max_threads)
        throw std::logic_error("SetMaxThreads: thread count must be greater than zero.");

    ThreadUtils::max_threads = max_threads;
    OPENMP_API_CALL(omp_set_num_threads(max_threads));
}

//-----------------------------------------------------------------------------
inline unsigned int GetMaxThreadCount()
{
    return ThreadUtils::max_threads;
}
