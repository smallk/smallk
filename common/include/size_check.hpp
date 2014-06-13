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

#include <limits>
#include <cstdint>

//-----------------------------------------------------------------------------
template <typename T>
bool FitsWithin(const uint64_t testval)
{
    // determine whether type T is large enough to contain 'testval'

    const uint64_t MAX_VAL = 
        static_cast<uint64_t>(std::numeric_limits<T>::max());

    return (testval <= MAX_VAL);
}

//-----------------------------------------------------------------------------
template <typename T>
bool SizeCheck(const unsigned int height, const unsigned int width)
{
    uint64_t required_size = static_cast<uint64_t>(height);
    required_size *= static_cast<uint64_t>(width);
    return FitsWithin<T>(required_size);
}
