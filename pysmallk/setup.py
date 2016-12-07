# -*- coding: utf-8 -*-
# Copyright 2013,2014,2015,2016 Georgia Institute of Technology.
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
print "thisPlatform: ", thisPlatform

smallk_path = "../"

# Path to the .hpp interface file
hpp_path = smallk_path + "smallk/include"

# Path to archive file (.a file)
archive_path = smallk_path + "build/bin/"
archive_filename = "libsmallk.a"
archive_file_path = "%s/%s" %(archive_path, archive_filename)

elem_include = os.environ['ELEM_INCLUDE']
#elem_include = os.environ["EL_INC"]
elem_lib = os.environ['ELEM_LIB']
#elem_lib = os.environ["EL_LIB"]

#hierclust sources
more_sources = [smallk_path + "hierclust/src/clust.cpp",
                smallk_path + "hierclust/src/clust_options.cpp",
                smallk_path + "hierclust/src/hierclust_json_writer.cpp",
                smallk_path + "hierclust/src/hierclust_writer_factory.cpp",
                smallk_path + "hierclust/src/hierclust_xml_writer.cpp",
                smallk_path + "hierclust/src/node.cpp",
                smallk_path + "hierclust/src/postprocess.cpp",
                smallk_path + "hierclust/src/tree.cpp"]

#preprocesseor and matrixgen sources
more_sources += [smallk_path + "preprocessor/src/preprocess.cpp",
                smallk_path + "common/src/term_frequency_matrix.cpp"]


#smallk and flatclust sources
more_sources += [smallk_path + "flatclust/src/flat_clust.cpp",
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
                smallk_path + "common/src/error.cpp",
                smallk_path + "common/src/flat_clust_output.cpp",
                smallk_path + "common/src/flatclust_xml_writer.cpp",
                smallk_path + "common/src/flatclust_json_writer.cpp",
                smallk_path + "common/src/delimited_file.cpp",
                smallk_path + "common/src/assignments.cpp",
                smallk_path + "common/src/system_mac.cpp"]

#get the appropriate include paths for the platform
if 'macosx' in thisPlatform:
#    links = ["-lgomp","-lpthread","-lelemental"]
    links = ["-lgomp","-lpthread","-lEl"]
elif 'linux' in thisPlatform:
#    links = ["-lgomp","-lpthread","-lelemental", "-lopenblas"]
    links = ["-lgomp","-lpthread","-lEl", "-lopenblas"]
else:
    print "\\nPlatform is not supported.", "\\n"

# extension module section
ext_modules = [
    Extension(
        name="pysmallk",
        sources=["interface/smallk_lib.pyx"] + more_sources,
        include_dirs = ["/usr/local/include", smallk_path + "common/include", smallk_path + "flatclust/include", smallk_path + "hierclust/include", smallk_path + "matrixgen/include", smallk_path + "preprocessor/include", hpp_path,elem_include, numpy.get_include()],
#        libraries=["mpi_cxx","mpi","m", "elemental"],
        libraries=["mpichcxx","mpi","m", "El"],
        library_dirs=[archive_path,"/usr/local/lib/", elem_lib],
        extra_objects=[archive_file_path],
        language="c++",
        extra_compile_args=["-fopenmp","-std=c++11", "-D ELEM_VER=85", "-fpermissive"],
        extra_link_args=links)
    ]

# call the main setup code from disutils
setup(
    name = "pysmallk",
    cmdclass = {"build_ext": build_ext},
    ext_modules = ext_modules,
    version=1.0,
    description="smallk.github.io",
    author="Barry L. Drake, Daniel C. Lee, Jason A. Poovey",
    author_email="barry.drake@gtri.gatech.edu, daniel.lee@gtri.gatech.edu, jason.poovey@gtri.gatech.edu",
    packages="libsmallk_cython"
    )

