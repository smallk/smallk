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

#include <string>
#include "nmf.hpp"
#include "flat_clust.hpp"
#include "file_format.hpp"

struct CommandLineOptions
{
    FlatClustOptions clust_opts;

    std::string infile_A;
    std::string infile_W;
    std::string infile_H;
    std::string dictfile;
    std::string outdir;
    std::string clustfile;
    std::string assignfile;
    FileFormat format;

    bool show_help;
};

void ShowHelp(const std::string& program_name);
void PrintOpts(const CommandLineOptions& opts);
bool ParseCommandLine(int argc, char* argv[], CommandLineOptions& opts);
bool IsValid(const CommandLineOptions& opts);

