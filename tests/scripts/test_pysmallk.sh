# This script runs as the 'make check' target from the master makefile.  The 
# positional argument $1 is the top-level folder for the xdata_data project.

echo "*****************************************************"
echo "*             pysmallk - command line               *"
echo "*                                                   *"
echo "*           Testing the smallk interface.           *"
echo "*                                                   *"
echo "*****************************************************"

w_file=w.csv
h_file=h.csv

if [ -f "$w_file" ]
then
    rm -f "$w_file"
fi

if [ -f "$h_file" ]
then
    rm -f "$h_file"
fi

python pysmallk/tests/smallkapi.py --matrixfile $1/reuters.mtx --infile_W $1/nmf_init_w.csv --infile_H $1/nmf_init_h.csv --k 8 --algorithm 'BPP' --outprecision 6 --miniter 1 --tol 0.005 --maxthreads 8

if cmp -s w.csv $1/test/nmf_result_w.csv; then
    echo "W matrix test passed"
else
    echo "W matrix test failed"
fi

if cmp -s h.csv $1/test/nmf_result_h.csv; then
    echo "H matrix test passed"
else
    echo "H matrix test failed"
fi

python pysmallk/tests/smallkapi.py --matrixfile $1/reuters.mtx --dictfile $1/reuters_dictionary.txt --k 5 

rm -f "$w_file"
rm -f "$h_file"
rm -f assignments_5.csv
rm -f tree_5.xml


echo "*****************************************************"
echo "*              pysmallk - in memory                 *"
echo "*                                                   *"
echo "*           Testing the smallk interface.           *"
echo "*                                                   *"
echo "*****************************************************"

w_file=w.csv
h_file=h.csv

if [ -f "$w_file" ]
then
    rm -f "$w_file"
fi

if [ -f "$h_file" ]
then
    rm -f "$h_file"
fi

python pysmallk/tests/smallkapi_inmem.py --indir $1

rm -f "$w_file"
rm -f "$h_file"

echo "*****************************************************"
echo "*             pysmallk - command line               *"
echo "*                                                   *"
echo "*            Testing the preprocessor.              *"
echo "*                                                   *"
echo "*****************************************************"

python pysmallk/tests/preprocessor.py --indir $1

if cmp -s reduced_matrix.mtx $1/test/reduced_matrix_20news.mtx; then
    echo "preprocessor matrix test passed"
else
    echo "preprocessor matrix test failed"
fi

if cmp -s reduced_dictionary.txt $1/test/reduced_dictionary_20news.txt; then
    echo "preprocessor dictionary test passed"
else
    echo "preprocessor dictionary test failed"
fi

if cmp -s reduced_documents.txt $1/test/reduced_documents_20news.txt; then
    echo "preprocessor documents test passed"
else
    echo "preprocessor documents test failed"
fi

rm -f reduced_matrix.mtx
rm -f reduced_dictionary.txt
rm -f reduced_documents.txt

echo "*****************************************************"
echo "*              pysmallk - in memory                 *"
echo "*                                                   *"
echo "*            Testing the preprocessor.              *"
echo "*                                                   *"
echo "*****************************************************"

python pysmallk/tests/preprocessor_inmem.py --indir $1

if cmp -s reduced_dictionary.txt $1/test/reduced_dictionary_20news.txt; then
    echo "preprocessor dictionary test passed"
else
    echo "preprocessor dictionary test failed"
fi

if cmp -s reduced_documents.txt $1/test/reduced_documents_20news.txt; then
    echo "preprocessor documents test passed"
else
    echo "preprocessor documents test failed"
fi

rm -f reduced_matrix.mtx
rm -f reduced_dictionary.txt
rm -f reduced_documents.txt

echo "*****************************************************"
echo "*             pysmallk - command line               *"
echo "*                                                   *"
echo "*                Testing hierclust.                 *"
echo "*                                                   *"
echo "*****************************************************"

echo "-------------------------------"
echo "                               " 
echo "  Reuters matrix, 12 clusters  "
echo "                               "
echo "-------------------------------"

xml_file=tree_12.xml
assign_file=assignments_12.csv

if [ -f "$xml_file" ]
then
    rm -f "$xml_file"
fi

if [ -f "$assign_file" ]
then
    rm -f "$assign_file"
fi

python pysmallk/tests/hierclust.py --matrixfile $1/reuters.mtx --dictfile $1/reuters_dictionary.txt --clusters 12 --initdir $1/test/matrices.reuters --format XML

if cmp -s "$xml_file" $1/test/reuters_tree_12.xml; then
    echo "XML file test passed"
else
    echo "XML file test failed"
fi

if cmp -s "$assign_file" $1/test/reuters_assignments_12.csv; then
    echo "assignment file test passed"
else
    echo "assignment file test failed"
fi

rm -f "$xml_file"
rm -f "$assign_file"

echo "-------------------------------"
echo "                               " 
echo "  20News matrix, 15 clusters  "
echo "                               "
echo "-------------------------------"

xml_file=tree_15.xml
assign_file=assignments_15.csv

if [ -f "$xml_file" ]
then
    rm -f "$xml_file"
fi

if [ -f "$assign_file" ]
then
    rm -f "$assign_file"
fi

python pysmallk/tests/hierclust.py --matrixfile $1/news20.mtx --dictfile $1/news20_dictionary.txt --clusters 15 --initdir $1/test/matrices.20news --miniter 1

if cmp -s "$xml_file" $1/test/news20_tree_15.xml; then
    echo "XML file test passed"
else
    echo "XML file test failed"
fi

if cmp -s "$assign_file" $1/test/news20_assignments_15.csv; then
    echo "assignment file test passed"
else
    echo "assignment file test failed"
fi

rm -f "$xml_file"
rm -f "$assign_file"

echo "*****************************************************"
echo "*              pysmallk - in memory                 *"
echo "*                                                   *"
echo "*              Testing hierclust.                   *"
echo "*                                                   *"
echo "*****************************************************"


python pysmallk/tests/hierclust_inmem.py --indir $1


echo "*****************************************************"
echo "*             pysmallk - command line               *"
echo "*                                                   *"
echo "*              Testing flatclust.                   *"
echo "*                                                   *"
echo "*****************************************************"

xml_file=tree_16.xml
assign_file=assignments_16.csv
fuzzy_file=assignments_fuzzy_16.csv

if [ -f "$xml_file" ]
then
    rm -f "$xml_file"
fi

if [ -f "$assign_file" ]
then
    rm -f "$assign_file"
fi

if [ -f "$fuzzy_file" ]
then
    rm -f "$fuzzy_file"
fi

python pysmallk/tests/flatclust.py --matrixfile $1/rnd_256_256.csv --dictfile $1/reuters_dictionary.txt --clusters 16 --infile_W $1/flatclust_init_w.csv --infile_H $1/flatclust_init_h.csv --miniter 1 --algorithm HALS --maxiter 5000

if cmp -s "$xml_file" $1/test/flatclust_rnd_clusters_16.xml; then
    echo "XML file test passed"
else
    echo "XML file test failed"
fi

if cmp -s "$assign_file" $1/test/flatclust_rnd_assignments_16.csv; then
    echo "assignment file test passed"
else
    echo "assignment file test failed"
fi

if cmp -s "$fuzzy_file" $1/test/flatclust_rnd_assignments_fuzzy_16.csv; then
    echo "fuzzy assignment file test passed"
else
    echo "fuzzy assignment file test failed"
fi

rm -f "$xml_file"
rm -f "$assign_file"
rm -f "$fuzzy_file"

echo "*****************************************************"
echo "*              pysmallk - in memory                 *"
echo "*                                                   *"
echo "*              Testing flatclust.                   *"
echo "*                                                   *"
echo "*****************************************************"


python pysmallk/tests/flatclust_inmem.py --indir $1



