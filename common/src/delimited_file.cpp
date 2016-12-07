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

#include "delimited_file.hpp"
#include "utils.hpp"

// file types supported by this reader
const std::string FILE_EXTENSION_CSV("CSV");

// comment chars; any line beginning with one of these chars is skipped
static const std::vector<char> COMMENT_CHARS = {'#', '%'};

//-----------------------------------------------------------------------------
bool IsDelimitedFile(const std::string& filename)
{
    std::string ext = FileExtension(filename);
    if (FILE_EXTENSION_CSV == ext)
        return true;
    // others TBD

    return false;
}

//-----------------------------------------------------------------------------
std::streampos 
SkipBlankLinesAndComments(std::ifstream& infile,
                          std::string& line)
{
    std::streampos start_pos = infile.tellg();

    // skip any blank lines or lines with comments
    while (std::getline(infile, line))
    {
        if (line.empty())
        {
            start_pos = infile.tellg();
            continue;
        }

        bool is_comment = false;
        for (unsigned int k=0; k != COMMENT_CHARS.size(); ++k)
        {
            if (COMMENT_CHARS[k] == line[0])
            {
                is_comment = true;
                break;
            }
        }

        if (is_comment)
        {
            start_pos = infile.tellg();
            continue;
        }

        break;
    }

    return start_pos;
}

//-----------------------------------------------------------------------------
void GetDimensions(std::ifstream& infile,
                   std::string& line,
                   unsigned int& height,
                   unsigned int& width,
                   const char DELIM)
{
    // count the delimiters in the first line
    width = 0;
    for (size_t c=0; c != line.size(); ++c)
    {
        if (DELIM == line[c])
            ++width;
    }

    // width == (num_delimiters + 1)
    ++width;

    // scan the file and count lines
    height = 1;
    while (std::getline(infile, line))
    {
        if (infile.eof())
            break;

        ++height;
    }
}

//-----------------------------------------------------------------------------
bool WriteDelimitedFile(const int* buffer, 
                        const unsigned int ldim,
                        const unsigned int height, 
                        const unsigned int width,
                        const std::string& filename,
                        const char DELIM)
{
    std::ofstream outfile(filename);
    if (!outfile)
        return false;

    // write in row-major order, to match Matlab
    for (unsigned int r=0; r != height; ++r)
    {
        for (unsigned int c=0; c != (width-1); ++c)
            outfile << buffer[c*ldim + r] << DELIM;

        outfile << buffer[(width-1)*ldim + r] << std::endl;
    }

    outfile.close();
    return true;
}

