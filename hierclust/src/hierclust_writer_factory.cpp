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

#include <stdexcept>
#include "hierclust_writer_factory.hpp"
#include "hierclust_xml_writer.hpp"
#include "hierclust_json_writer.hpp"

//-----------------------------------------------------------------------------
IHierclustWriter* CreateHierclustWriter(const FileFormat& format)
{
    IHierclustWriter* writer = nullptr;

    switch (format)
    {
    case FileFormat::XML:
        writer = new HierclustXmlWriter();
        break;
    case FileFormat::JSON:
        writer = new HierclustJsonWriter();
        break;
    default:
        throw std::logic_error("CreateHierclustWriter: invalid file format");
        break;
    }

    return writer;
}
