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

smallk_path = "../../../../libsmallk-1.4.0/"


more_sources = [smallk_path + "matrixgen/src/main.cpp",
                smallk_path + "matrixgen/src/command_line.cpp",
                smallk_path + "common/src/matrix_market_file.cpp",
                smallk_path + "common/src/utils.cpp",
                smallk_path + "common/src/system_mac.cpp",
                smallk_path + "common/src/constants.cpp"]




#get the appropriate include paths for the platform
if 'macosx' in thisPlatform:
    include = [smallk_path + "common/include/", smallk_path + "matrixgen/include", numpy.get_include()]
    library = []
elif 'linux' in thisPlatform:
    include = ["/usr/local/include", smallk_path + "common/include/", smallk_path + "matrixgen/include", numpy.get_include()]
    library = ["/usr/local/lib/", "/usr/lib/"]
else:
    print "\nPlatform is not supported.", "\n"

ext_modules = [
    Extension(
        name = "libmatrixgen",
        sources = ["interface/matrixgen_lib.pyx"] + more_sources,
        include_dirs = include,
        libraries = ["m"],
        library_dirs = library,
        language="c++",
        extra_compile_args=["-std=c++11", "-D ELEM_VER=84", "-fpermissive"]
    )
    ]

setup(
    name = 'libmatrixgen',
    cmdclass = {'build_ext' : build_ext},
    ext_modules = ext_modules
    )
