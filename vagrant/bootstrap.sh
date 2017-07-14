#!/usr/bin/env bash
echo "================= in bootstrap.sh ========================================================"
date
# Make sure the package information is up-to-date
apt-get update

mkdir -p /home/ubuntu/Downloads

# Compilers
apt-get install -y gcc-5
apt-get install -y g++-5
apt-get install -y gfortran-5
apt-get install -y gfortran
apt-get install -y unzip

apt-get install -y python-dev
apt-get install -y python-pip
pip install --upgrade pip
pip install numpy
apt-get install -y python-scipy
pip install cython

#mpich needed by Elemental 0.85
apt-get install -y mpich

# Source control
apt-get install -y git

# Configuration
apt-get install -y cmake

# libflame
git clone https://github.com/flame/libflame.git /home/ubuntu/libflame
cd /home/ubuntu/libflame
./configure --prefix=/usr/local/flame --with-cc=/usr/bin/gcc-5 --with-ranlib=/usr/bin/gcc-ranlib-5
make -j4
make install
chown -R ubuntu /home/ubuntu/libflame

# OpenBLAS and LAPACK
apt-get install -y libopenblas-dev
apt-get install -y libatlas-dev liblapack-dev

# Elemental
mkdir /usr/local/elemental
export ELEMENTAL_INSTALL_DIR=/usr/local/elemental
wget -O /home/ubuntu/Downloads/Elemental-0.85.tgz http://libelemental.org/pub/releases/Elemental-0.85.tgz
tar -zxvf /home/ubuntu/Downloads/Elemental-0.85.tgz -C /home/ubuntu
cd /home/ubuntu/Elemental-0.85

#Version 0.85 of Elemental has an error in one of its cmake files. The file is: 
#	Elemental-0.85/cmake/tests/CXX.cmake
sed -i'' '1s/include(FindCXXFeatures)/include_directories(FindCXXFeatures)/' ./cmake/tests/CXX.cmake

mkdir build_hybrid
mkdir build_pure

cd build_hybrid
cmake -D CMAKE_INSTALL_PREFIX=/usr/local/elemental/0.85/HybridRelease -D CMAKE_BUILD_TYPE=HybridRelease -D CMAKE_CXX_COMPILER=/usr/bin/g++-5 -D CMAKE_CXX_FLAGS="-std=c++11 -fPIC" -D CMAKE_C_COMPILER=/usr/bin/gcc-5 -D CMAKE_C_FLAG=-fPIC -D CXX_FLAGS="-std=c++11 -fPIC" -D CMAKE_Fortran_COMPILER=/usr/bin/gfortran-5 -D MATH_LIBS="/usr/local/flame/lib/libflame.a;-L/usr/local/lib/ -lopenblas;/lib/x86_64-linux-gnu/libm.so.6" –D ELEM_EXAMPLES=ON –D ELEM_TESTS=ON  ..
make -j4
make install

cd ../build_pure
cmake -D CMAKE_INSTALL_PREFIX=/usr/local/elemental/0.85/PureRelease -D CMAKE_BUILD_TYPE=PureRelease -D CMAKE_CXX_COMPILER=/usr/bin/g++-5 -D CMAKE_CXX_FLAGS="-std=c++11 -fPIC" -D CMAKE_C_COMPILER=/usr/bin/gcc-5 -D CMAKE_C_FLAG=-fPIC -D CXX_FLAGS="-std=c++11 -fPIC" -D CMAKE_Fortran_COMPILER=/usr/bin/gfortran-5 -D MATH_LIBS="/usr/local/flame/lib/libflame.a;-L/usr/local/lib/ -lopenblas;/lib/x86_64-linux-gnu/libm.so.6" –D ELEM_EXAMPLES=ON –D ELEM_TESTS=ON  ..
make -j4
make install

#chown -R vagrant.vagrant /home/ubuntu/Elemental-0.84-p1


