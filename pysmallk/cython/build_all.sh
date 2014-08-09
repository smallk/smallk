# build_all.sh
# just run the command line with no args to build the .so shared library files
# These can be copied to /usr/local/smallk/lib for ease of use

# build flatclust
cd flatclust
./clean.sh
python setup.py build_ext --inplace

# build hierclust
cd ../hierclust
./clean.sh
python setup.py build_ext --inplace

# build matrixgen
cd ../matrixgen
./clean.sh
python setup.py build_ext --inplace

# build preprocessor
cd ../preprocessor
./clean.sh
python setup.py build_ext --inplace

# build smallk
cd ../smallk
./clean.sh
python setup.py build_ext --inplace

cd ..
