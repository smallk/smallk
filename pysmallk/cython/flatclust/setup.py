# -*- coding: utf-8 -*-
# Copyright 2013,2014 Georgia Institute of Technology.
#
# Licensed under the Apache License, Version 2.0 (the “License”); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by
# applicable law or agreed to in writing, software distributed under the License
# is distributed on as “AS IS” BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the specific language
# governing permissions and limitations under the License.

import numpy 
from distutils.core import setup
import distutils.util as du 
from distutils.extension import Extension
from Cython.Distutils import build_ext

import os

# do a platform check
thisPlatform = du.get_platform()
print "\nPlatform is ", thisPlatform, "\n"

#edit these paths to reflect your system configuration
if 'macosx' in thisPlatform:
    os.environ['CC'] = '/usr/local/bin/g++-4.9'
    os.environ['CXX'] = '/usr/local/bin/g++-4.9'
    os.environ['CPP'] = '/usr/local/bin/g++-4.9'
    os.environ['CMAKE_CXX_COMPILER'] = '/usr/local/bin/g++-4.9'
elif 'linux' in thisPlatform:
    os.environ['CC'] = '/usr/bin/g++-4.8'
    os.environ['CXX'] = '/usr/bin/g++-4.8'
    os.environ['CPP'] = '/usr/bin/g++-4.8'
    os.environ['CMAKE_CXX_COMPILER'] = '/usr/bin/g++-4.8'
else:
    print "\nPlatform is not supported.", "\n"

smallk_path = "../../../"

elem_base = "/usr/local/elemental"
elem_ver = "0.84"
elem_build_mode = "HybridRelease"

elem_path = "%s/%s/%s" % (elem_base, elem_ver, elem_build_mode)
elem_include = "%s/include" % elem_path
elem_lib = "%s/lib" % elem_path

more_sources = [smallk_path + "flatclust/src/flat_clust.cpp",
                smallk_path + "flatclust/src/flat_clust_options.cpp",
                smallk_path + "flatclust/src/command_line.cpp",
                smallk_path + "common/src/utils.cpp",
                smallk_path + "common/src/nmf.cpp",
                smallk_path + "common/src/nmf_options.cpp",
                smallk_path + "common/src/bit_matrix.cpp",
                smallk_path + "common/src/bit_matrix_ops.cpp",
                smallk_path + "common/src/nnls.cpp",
                smallk_path + "common/src/spooky_v2.cpp",
                smallk_path + "common/src/xxhash.cpp",
                smallk_path + "common/src/file_loader.cpp",
                smallk_path + "common/src/matrix_market_file.cpp",
                smallk_path + "common/src/constants.cpp",
                smallk_path + "common/src/flat_clust_output.cpp",
                smallk_path + "common/src/flatclust_xml_writer.cpp",
                smallk_path + "common/src/flatclust_json_writer.cpp",
                smallk_path + "common/src/delimited_file.cpp",
                smallk_path + "common/src/assignments.cpp",
                smallk_path + "common/src/system_mac.cpp"]



#get the appropriate include paths for the platform
if 'macosx' in thisPlatform:
    links = ["-lgomp", "-lpthread", "-lelemental"]
elif 'linux' in thisPlatform:
    links = ["-lgomp", "-lpthread", "-lelemental", "-lopenblas"]
else:
    print "\nPlatform is not supported.", "\n"


ext_modules = [
    Extension(
        name = "libflatclust",
        sources = ["interface/flatclust_lib.pyx"]+more_sources,
        include_dirs = ["/usr/local/include", "/usr/include", smallk_path + "common/include", smallk_path + "flatclust/include", elem_include, numpy.get_include()],
        library_dirs = [elem_lib],
        libraries = ["mpi_cxx", "mpi", "elemental"],
        language = "c++",
        extra_compile_args = ["-fopenmp", "-std=c++11", "-D ELEM_VER=84"],
        extra_link_args = links
    )]
    
setup(
    name = "libflatclust",
    cmdclass = {"build_ext" : build_ext},
    ext_modules = ext_modules
)
