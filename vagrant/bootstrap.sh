#!/usr/bin/env bash
echo "================= in bootstrap.sh ========================================================"
date
# Make sure the package information is up-to-date
apt-get update

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
pip install cython==0.22

#mpich needed by Elemental 0.85
apt-get install -y mpich

# Source control
apt-get install -y git
git config --global core.autocrlf input

# Configuration
apt-get install -y cmake

echo "----------------------------- libflame -----------------------------"

# libflame
git clone https://github.com/flame/libflame.git /home/vagrant/libflame
cd /home/vagrant/libflame
./configure --prefix=/usr/local/flame --with-cc=/usr/bin/gcc-5 --with-ranlib=/usr/bin/gcc-ranlib-5 CFLAGS=-fPIC CXXFLAGS=-fPIC --enable-shared
make -j`nproc`
make install
chown -R vagrant.vagrant /home/vagrant/libflame

echo "----------------------------- OpenBLas and LAPACK -----------------------------"

# OpenBLAS and LAPACK
apt-get install -y libopenblas-dev
apt-get install -y libatlas-dev liblapack-dev

echo "----------------------------- Elemental -----------------------------"

# Elemental
mkdir /usr/local/elemental
export ELEMENTAL_INSTALL_DIR=/usr/local/elemental
mkdir -p /home/vagrant/Downloads
wget -O /home/vagrant/Downloads/Elemental-0.85.tgz http://libelemental.org/pub/releases/Elemental-0.85.tgz
tar -zxvf /home/vagrant/Downloads/Elemental-0.85.tgz -C /home/vagrant
cd /home/vagrant/Elemental-0.85

#Version 0.85 of Elemental has an error in one of its cmake files. The file is: 
#	Elemental-0.85/cmake/tests/CXX.cmake
sed -i'' '1s/include(FindCXXFeatures)/include_directories(FindCXXFeatures)/' ./cmake/tests/CXX.cmake

mkdir build_hybrid
mkdir build_pure

cd build_hybrid
cmake -D CMAKE_INSTALL_PREFIX=/usr/local/elemental/0.85/HybridRelease -D CMAKE_BUILD_TYPE=HybridRelease -D CMAKE_CXX_COMPILER=/usr/bin/g++-5 -D CMAKE_CXX_FLAGS="-std=c++11 -fPIC" -D CMAKE_C_COMPILER=/usr/bin/gcc-5 -D CMAKE_C_FLAGS=-fPIC -D CXX_FLAGS="-std=c++11 -fPIC" -D CMAKE_Fortran_COMPILER=/usr/bin/gfortran-5 -D MATH_LIBS="/usr/local/flame/lib/libflame.a;-L/usr/local/lib/ -lopenblas;/lib/x86_64-linux-gnu/libm.so.6" –D ELEM_EXAMPLES=ON –D ELEM_TESTS=ON  ..
make -j`nproc`
make install

cd ../build_pure
cmake -D CMAKE_INSTALL_PREFIX=/usr/local/elemental/0.85/PureRelease -D CMAKE_BUILD_TYPE=PureRelease -D CMAKE_CXX_COMPILER=/usr/bin/g++-5 -D CMAKE_CXX_FLAGS="-std=c++11 -fPIC" -D CMAKE_C_COMPILER=/usr/bin/gcc-5 -D CMAKE_C_FLAGS=-fPIC -D CXX_FLAGS="-std=c++11 -fPIC" -D CMAKE_Fortran_COMPILER=/usr/bin/gfortran-5 -D MATH_LIBS="/usr/local/flame/lib/libflame.a;-L/usr/local/lib/ -lopenblas;/lib/x86_64-linux-gnu/libm.so.6" –D ELEM_EXAMPLES=ON –D ELEM_TESTS=ON  ..
make -j`nproc`
make install

chown -R vagrant.vagrant /home/vagrant/Elemental-0.85

ln -s /usr/local/elemental/0.85/HybridRelease/lib/*.so /usr/lib/

echo "----------------------------- SmallK -----------------------------"

# SmallK
apt-get install -y libmetis-dev
cd /home/vagrant
git clone https://github.com/smallk/smallk_data.git
cp /vagrant/libsmallk-1.6.2.tar.gz /home/vagrant/
#cp /vagrant/smallk_data.zip /home/vagrant/ ##### sink this directory with
tar -zxvf /home/vagrant/libsmallk-1.6.2.tar.gz -C /home/vagrant
#unzip /home/vagrant/smallk_data.zip -d /home/vagrant/smallk_data
cd /home/vagrant/libsmallk-1.6.2

# make SITE_PACKAGES_DIR=/usr/local/lib/python2.7/dist-packages/ install

#in pysmallk/setup.py change include_dirs to have "usr/include/mpich" (escaping $.*/[\]^)
sed -i'' 's/include_dirs = \["\/usr\/local\/include",/include_dirs = \["\/usr\/local\/include", "\/usr\/include\/mpich",/' ./pysmallk/setup.py

make all PYSMALLK=1 ELEMVER=0.85 DATA_DIR=../smallk_data
#install it in: /usr/local/smallk/bin
make install PYSMALLK=1 ELEMVER=0.85
#add pysmallk.so to the path
#this replaces something like: 
#ln -s /usr/local/lib/python2.7/site-packages/pysmallk.so  /usr/lib/python2.7/lib-dynload/pysmallk.so
echo "export PYTHONPATH=/usr/local/lib/python2.7/site-packages" >> /home/vagrant/.bashrc
echo "export PATH=/usr/local/smallk/bin:$PATH" >> home/vagrant/.bashrc

#run tests
make check PYSMALLK=1 ELEMVER=0.85 DATA_DIR=../smallk_data

chown -R vagrant.vagrant /home/vagrant/libsmallk-1.6.2

#verify paths
#python -c "import sys; print '\n'.join(sys.path)"
#test that pysmallk can be imported from anywhere
#python -c "import numpy; import pysmallk"
