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

#include <cstdint>

// Population count - count the number of '1' bits in an unsigned integer.

//-----------------------------------------------------------------------------
inline uint32_t PopulationCount(const uint32_t x)
{
    // This code is from Hacker's Delight, 2nd ed. by Henry Warren, p. 82.

    uint32_t y = x - ((x >> 1) & 0x55555555);
    y = (y & 0x33333333) + ((y >> 2) & 0x33333333);
    y = (y + (y >> 4)) & 0x0F0F0F0F;
    y =  y + (y >> 8);
    y =  y + (y >> 16);
    return y & 0x0000003F;
}

//-----------------------------------------------------------------------------
inline uint32_t PopulationCount(const uint64_t x)
{
    uint32_t lo = static_cast<uint32_t>(x);
    uint32_t hi = static_cast<uint32_t>(x >> 32);
    return PopulationCount(lo) + PopulationCount(hi);
}

