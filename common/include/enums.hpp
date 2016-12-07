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

// These are used internally and are analogs of those in Elemental.

enum Orientation
{
    NORMAL,
    TRANSPOSE
};

enum NormType
{
    MAX_NORM,
    ONE_NORM,
    INFINITY_NORM,
    FROBENIUS_NORM
};

enum LeftOrRight
{
    LEFT,
    RIGHT
};

enum UpperOrLower
{
    UPPER,
    LOWER
};

enum UnitOrNonUnit
{
    UNIT,
    NON_UNIT
};
