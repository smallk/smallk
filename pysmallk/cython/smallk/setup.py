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

smallk_path = "../../../../libsmallk-1.3.0/"

# Path to the .hpp interface file
hpp_path = '/usr/local/smallk/include'

# Path to archive file (.a file)
archive_path = '/usr/local/smallk/lib'
archive_filename = 'libsmallk.a'
archive_file_path = '%s/%s' %(archive_path, archive_filename)

# Set the elem base, version and build mode
elem_base = '/usr/local/elemental'
elem_ver = '0.83'
elem_build_mode = 'HybridRelease'

# Elemental path to include and lib
elem_path = '%s/%s/%s' %(elem_base, elem_ver, elem_build_mode)
elem_include = '%s/include' %elem_path
elem_lib = '%s/lib' %elem_path

#get the appropriate include paths for the platform
if 'macosx' in thisPlatform:
    library = [archive_path,"/usr/local/lib/", elem_lib]
    links = ["-lgomp","-lpthread","-lelemental"]
elif 'linux' in thisPlatform:
    library = [archive_path,"/usr/local/lib/", elem_lib]
    links = ["-lgomp","-lpthread","-lelemental", "-lopenblas"]
else:
    print "\nPlatform is not supported.", "\n"

# extension module section
ext_modules = [
    Extension(
        name="libsmallkpy",
        sources=["interface/smallk_lib.pyx"],
        include_dirs = ["/usr/local/include", hpp_path,elem_include, numpy.get_include()],
        libraries=["mpi_cxx","mpi","m"],
        library_dirs=library,
        extra_objects=[archive_file_path],
        language="c++",
        extra_compile_args=["-fopenmp","-std=c++11", "-D ELEM_VER=83"],
        extra_link_args=links)
    ]

# call the main setup code from disutils
setup(
    name = 'libsmallkpy',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
    version=0.01,
    description="smallk.github.io",
    author="Barry L. Drake, Daniel C. Lee, Jason A. Poovey",
    author_email="barry.drake@gtri.gatech.edu, daniel.lee@gtri.gatech.edu, jason.poovey@gtri.gatech.edu",
    packages="libsmallk_cython"
    )

