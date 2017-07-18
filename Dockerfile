FROM ubuntu:16.04

ENV SMALLK_SRC=/usr/local/src/smallk

# Switch to Leaseweb mirrors
RUN sed -i 's|archive.ubuntu.com/ubuntu|mirror.us.leaseweb.net/ubuntu|g' /etc/apt/sources.list

# Touch this file for a complete APT re-update
ADD docker/apt-reupdate /tmp

# Apt updates and core utils
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y wget nano vim supervisor sudo openssh-client git xterm build-essential

# smallk dependencies
RUN apt-get install -y gcc-5 g++-5 gfortran-5 gfortran unzip python-dev python-pip python-numpy \
                       python-scipy mpich cmake libopenblas-dev libatlas-dev liblapack-dev libmetis-dev
RUN pip install --upgrade pip
RUN pip install cython==0.22

# Bring in the entire source tree
ADD . $SMALLK_SRC

# Build and install libflame
WORKDIR $SMALLK_SRC/modules/libflame
RUN ./configure --prefix=/usr/local/flame --with-cc=/usr/bin/gcc-5 --with-ranlib=/usr/bin/gcc-ranlib-5 \
    CFLAGS=-fPIC CXXFLAGS=-fPIC
RUN make -j`nproc`
RUN make install

# Build and install libelemental
# Version 0.85 of Elemental has an error in one of its cmake files. The file is:
#	Elemental-0.85/cmake/tests/CXX.cmake
WORKDIR $SMALLK_SRC/modules/libelemental
RUN sed -i '1s/include(FindCXXFeatures)/include_directories(FindCXXFeatures)/' ./cmake/tests/CXX.cmake

RUN mkdir $SMALLK_SRC/modules/libelemental/build_hybrid
WORKDIR $SMALLK_SRC/modules/libelemental/build_hybrid
RUN cmake -D CMAKE_INSTALL_PREFIX=/usr/local/elemental/0.85/HybridRelease \
    -D CMAKE_BUILD_TYPE=HybridRelease -D CMAKE_CXX_COMPILER=/usr/bin/g++-5 \
    -D CMAKE_CXX_FLAGS="-std=c++11 -fPIC" -D CMAKE_C_COMPILER=/usr/bin/gcc-5 -D CMAKE_C_FLAGS=-fPIC \
    -D CXX_FLAGS="-std=c++11 -fPIC" -D CMAKE_Fortran_COMPILER=/usr/bin/gfortran-5 \
    -D MATH_LIBS="/usr/local/flame/lib/libflame.a;-L/usr/local/lib/ -lopenblas;/lib/x86_64-linux-gnu/libm.so.6" \
    -D ELEM_EXAMPLES=ON -D ELEM_TESTS=ON ..
RUN make -j`nproc`
RUN make install

RUN mkdir $SMALLK_SRC/modules/libelemental/build_pure
WORKDIR $SMALLK_SRC/modules/libelemental/build_pure
RUN cmake -D CMAKE_INSTALL_PREFIX=/usr/local/elemental/0.85/PureRelease \
    -D CMAKE_BUILD_TYPE=PureRelease -D CMAKE_CXX_COMPILER=/usr/bin/g++-5 \
    -D CMAKE_CXX_FLAGS="-std=c++11 -fPIC" -D CMAKE_C_COMPILER=/usr/bin/gcc-5 -D CMAKE_C_FLAGS=-fPIC \
    -D CXX_FLAGS="-std=c++11 -fPIC" -D CMAKE_Fortran_COMPILER=/usr/bin/gfortran-5 \
    -D MATH_LIBS="/usr/local/flame/lib/libflame.a;-L/usr/local/lib/ -lopenblas;/lib/x86_64-linux-gnu/libm.so.6" \
    -D ELEM_EXAMPLES=ON -D ELEM_TESTS=ON ..
RUN make -j`nproc`
RUN make install

RUN ln -s /usr/local/elemental/0.85/HybridRelease/lib/*.so /usr/lib/

ENV ELEMENTAL_INSTALL_DIR=/usr/local/elemental

# Build and install smallk

WORKDIR $SMALLK_SRC
RUN make all PYSMALLK=1 ELEMVER=0.85
RUN make install PYSMALLK=1 ELEMVER=0.85

RUN chown -R 1000:1000 $SMALLK_SRC

# Set up a standard user with the same UID/GID as the standard user on the host
RUN mkdir -p /home/docker
ADD docker/.bashrc.additional /home/docker/.bashrc.additional
RUN cat /root/.bashrc /home/docker/.bashrc.additional >> /home/docker/.bashrc
RUN chown -R 1000:1000 /home/docker

RUN echo 'docker:x:1000:1000:Docker,,,:/home/docker:/bin/bash' >> /etc/passwd
RUN echo 'docker:x:1000:' >> /etc/group

# Give the docker user passwordless sudo
RUN echo 'docker ALL=(ALL) NOPASSWD: ALL' > /etc/sudoers.d/docker

USER 1000
WORKDIR /home/docker

CMD ["/bin/bash"]
