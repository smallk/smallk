#!/usr/bin/env bash

# Make sure the package information is up-to-date
apt-get update

# Compilers
apt-get install -y gcc-4.8
apt-get install -y g++-4.8
apt-get install -y gfortran-4.8
apt-get install -y gfortran
apt-get install -y unzip

apt-get install -y python-dev
apt-get install -y python-pip
pip install numpy
apt-get install -y python-scipy
pip install cython

# OpenMPI
mkdir /home/vagrant/Downloads
wget -O /home/vagrant/Downloads/openmpi-1.8.1.tar.gz http://www.open-mpi.org/software/ompi/v1.8/downloads/openmpi-1.8.1.tar.gz
tar -zxvf /home/vagrant/Downloads/openmpi-1.8.1.tar.gz -C /home/vagrant
cd /home/vagrant/openmpi-1.8.1
./configure --enable-mpi-thread-multiple --disable-dlopen --prefix=/usr/local  CC=/usr/bin/gcc-4.8 CXX=/usr/bin/g++-4.8 F77=usr/bin/gfortran-4.8 FC=/usr/bin/gfortran-4.8
make -j4
make install
chown -R vagrant.vagrant /home/vagrant/Downloads
chown -R vagrant.vagrant /home/vagrant/openmpi-1.8.1

# Source control
apt-get install -y git

# Configuration
apt-get install -y cmake

# libflame
git clone https://github.com/flame/libflame.git /home/vagrant/libflame
cd /home/vagrant/libflame
./configure --prefix=/usr/local/flame --with-cc=/usr/bin/gcc-4.8 --with-ranlib=/usr/bin/gcc-ranlib-4.8
make -j4
make install
chown -R vagrant.vagrant /home/vagrant/libflame

# OpenBLAS
wget -O /home/vagrant/Downloads/openblas-0.2.9.tar http://github.com/xianyi/OpenBLAS/tarball/v0.2.9
tar -zxvf /home/vagrant/Downloads/openblas-0.2.9.tar -C /home/vagrant
cd /home/vagrant/xianyi-OpenBLAS-347dded
make BINARY=64 USE_OPENMP=1
make PREFIX=/usr/local/openblas/0.2.9/ install
ln -s /usr/local/openblas/0.2.9/lib/libopenblas.so /usr/local/lib/libopenblas.so
ln -s /usr/local/openblas/0.2.9/lib/libopenblas.so.0 /usr/local/lib/libopenblas.so.0
chown -R vagrant.vagrant /home/vagrant/xianyi-OpenBLAS-347dded


# Elemental
mkdir /usr/local/elemental
export ELEMENTAL_INSTALL_DIR=/usr/local/elemental
wget -O /home/vagrant/Downloads/Elemental-0.84-p1.tgz http://libelemental.org/pub/releases/Elemental-0.84-p1.tgz
tar -zxvf /home/vagrant/Downloads/Elemental-0.84-p1.tgz -C /home/vagrant
cd /home/vagrant/Elemental-0.84-p1
mkdir build_hybrid
mkdir build_pure
cd build_hybrid
cmake -D CMAKE_INSTALL_PREFIX=/usr/local/elemental/0.84-p1/HybridRelease -D CMAKE_BUILD_TYPE=HybridRelease -D CMAKE_CXX_COMPILER=/usr/bin/g++-4.8 -D CMAKE_CXX_FLAGS=-fPIC -D CMAKE_C_COMPILER=/usr/bin/gcc-4.8 -D CMAKE_C_FLAG=-fPIC -D CXX_FLAGS="-std=c++11 -fPIC" -D CMAKE_Fortran_COMPILER=/usr/bin/gfortran-4.8 -D MATH_LIBS="/usr/local/flame/lib/libflame.a;-L/usr/local/lib/ -lopenblas;/lib/x86_64-linux-gnu/libm.so.6" –D ELEM_EXAMPLES=ON –D ELEM_TESTS=ON  ..
make -j4
make install

cd ../build_pure
cmake -D CMAKE_INSTALL_PREFIX=/usr/local/elemental/0.84-p1/PureRelease -D CMAKE_BUILD_TYPE=PureRelease -D CMAKE_CXX_COMPILER=/usr/bin/g++-4.8 -D CMAKE_CXX_FLAGS=-fPIC -D CMAKE_C_COMPILER=/usr/bin/gcc-4.8 -D CMAKE_C_FLAG=-fPIC -D CXX_FLAGS="-std=c++11 -fPIC" -D CMAKE_Fortran_COMPILER=/usr/bin/gfortran-4.8 -D MATH_LIBS="/usr/local/flame/lib/libflame.a;-L/usr/local/lib/ -lopenblas;/lib/x86_64-linux-gnu/libm.so.6" –D ELEM_EXAMPLES=ON –D ELEM_TESTS=ON  ..
make -j4
make install
chown -R vagrant.vagrant /home/vagrant/Elemental-0.84-p1


# Java
cp /vagrant/jdk-8u25-linux-x64.gz /home/vagrant/Downloads
tar -zxvf /home/vagrant/Downloads/jdk-8u25-linux-x64.gz -C /home/vagrant
chown -R vagrant.vagrant /home/vagrant/jdk1.8.0_25

#somehow add to .bashrc
echo "export JAVA_HOME=/home/vagrant/jdk1.8.0_25/jre" >> .bashrc
apt-get install -y ant
ln -s /home/vagrant/jdk1.8.0_25/jre/bin/java /usr/local/bin/java



# SmallK
#cp /vagrant/libsmallk-1.4.0.tar.gz /home/vagrant/
#cp /vagrant/xdata_data.zip /home/vagrant/
#tar -zxvf /home/vagrant/libsmallk-1.4.0.tar.gz -C /home/vagrant
#unzip /home/vagrant/xdata_data.zip -d /home/vagrant
#cd /home/vagrant/libsmallk-1.4.0

#make PYSMALLK=1 JAVA=/usr/local/bin/java ANT=/usr/bin/ant
# make SITE_PACKAGES_DIR=/usr/local/lib/python2.7/dist-packages/ install
#make check
#chown -R vagrant.vagrant /home/vagrant/libsmallk-1.4.0
