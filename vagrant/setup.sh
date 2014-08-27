# Anaconda
wget -O /home/vagrant/Downloads/Anaconda-2.0.1-Linux-x86_64.sh http://09c8d0b2229f813c1b93-c95ac804525aac4b6dba79b00b39d1d3.r79.cf1.rackcdn.com/Anaconda-2.0.1-Linux-x86_64.sh
sudo chmod +x /home/vagrant/Downloads/Anaconda-2.0.1-Linux-x86_64.sh
/home/vagrant/Downloads/Anaconda-2.0.1-Linux-x86_64.sh -b

# prepare environment
export PATH=~/anaconda/bin:$PATH

sudo sed -i '/include/ i\include /usr/local/lib' /etc/ld.so.conf
sudo ldconfig

# matrixgen cython
cd /home/vagrant/smallk/pysmallk/cython/matrixgen
python setup.py build_ext --inplace
python test_matrixgen.py --height 10 --width 10 --filename test.txt --type UNIFORM
python test_matrixgen.py --height 10 --width 10 --filename test.txt --type DENSE_DIAG
python test_matrixgen.py --height 10 --width 10 --filename test.txt --type SPARSE_DIAG
python test_matrixgen.py --height 10 --width 10 --filename test.txt --type ONES
python test_matrixgen.py --height 10 --width 10 --filename test.txt --type ZEROS
python test_matrixgen.py --height 10 --width 10 --filename test.txt --type SPARSE


# preprocessor cython
cd ../preprocessor
python setup.py build_ext --inplace
python test_preprocessor.py --indir data/

# smallk cython
cd ../smallk
python setup.py build_ext --inplace
python test_smallk.py --k 16 --algorithm MU --matrixfile data/reuters.mtx

# hierclust cython
cd ../hierclust
python setup.py build_ext --inplace
python test_hierclust.py --matrixfile data/reuters.mtx --dictfile data/reuters_dictionary.txt --clusters 16

# flatclust cython
cd ../flatclust
python setup.py build_ext --inplace
python test_flatclust.py --matrixfile data/reuters.mtx --dictfile data/reuters_dictionary.txt --clusters 16


