# This script runs as the 'make check' target from the master makefile.

echo "*****************************************************"
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

smallk/bin/smallk_tester data

if cmp -s w.csv data/test/nmf_result_w.csv; then
    echo "W matrix test passed"
else
    echo "W matrix test failed"
fi

if cmp -s h.csv data/test/nmf_result_h.csv; then
    echo "H matrix test passed"
else
    echo "H matrix test failed"
fi

rm -f "$w_file"
rm -f "$h_file"
rm -f assignments_5.csv
rm -f tree_5.xml

echo "*****************************************************"
echo "*                                                   *"
echo "*            Testing the preprocessor.              *"
echo "*                                                   *"
echo "*****************************************************"

preprocessor/bin/preprocess_tf --indir data

if cmp -s reduced_matrix.mtx data/test/reduced_matrix_20news.mtx; then
    echo "preprocessor matrix test passed"
else
    echo "preprocessor matrix test failed"
fi

if cmp -s reduced_dictionary.txt data/test/reduced_dictionary_20news.txt; then
    echo "preprocessor dictionary test passed"
else
    echo "preprocessor dictionary test failed"
fi

if cmp -s reduced_documents.txt data/test/reduced_documents_20news.txt; then
    echo "preprocessor documents test passed"
else
    echo "preprocessor documents test failed"
fi

rm -f reduced_matrix.mtx
rm -f reduced_dictionary.txt
rm -f reduced_documents.txt

echo "*****************************************************"
echo "*                                                   *"
echo "*            Testing the NMF routines.              *"
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

nmf/bin/nmf --matrixfile data/reuters.mtx --algorithm BPP --k 8 --infile_W data/nmf_init_w.csv --infile_H data/nmf_init_h.csv --outprecision 6 --miniter 1

if cmp -s w.csv data/test/nmf_result_w.csv; then
    echo "NMF W matrix test passed"
else
    echo "NMF W matrix test failed"
fi

if cmp -s h.csv data/test/nmf_result_h.csv; then
    echo "NMF H matrix test passed"
else
    echo "NMF H matrix test failed"
fi

rm -f "$w_file"
rm -f "$h_file"

echo "*****************************************************"
echo "*                                                   *"
echo "*                Testing hierclust.                 *"
echo "*                                                   *"
echo "*****************************************************"

xml_file=tree_5.xml
assign_file=assignments_5.csv

if [ -f "$xml_file" ]
then
    rm -f "$xml_file"
fi

if [ -f "$assign_file" ]
then
    rm -f "$assign_file"
fi

hierclust/bin/hierclust --matrixfile data/reuters.mtx --dictfile data/reuters_dictionary.txt --clusters 5 --infile_W data/hierclust_init_w.csv --infile_H data/hierclust_init_h.csv --miniter 1

if cmp -s "$xml_file" data/test/reuters_tree_5.xml; then
    echo "hierclust cluster file test passed"
else
    echo "hierclust cluster file test failed"
fi

if cmp -s "$assign_file" data/test/reuters_assignments_5.csv; then
    echo "hierclust assignment file test passed"
else
    echo "hierclust assignment file test failed"
fi

rm -f "$xml_file"
rm -f "$assign_file"

echo "*****************************************************"
echo "*                                                   *"
echo "*              Testing flatclust.                   *"
echo "*                                                   *"
echo "*****************************************************"

xml_file=clusters_16.xml
assign_file=assignments_16.csv

if [ -f "$xml_file" ]
then
    rm -f "$xml_file"
fi

if [ -f "$assign_file" ]
then
    rm -f "$assign_file"
fi

flatclust/bin/flatclust --matrixfile data/rnd_256_256.csv --dictfile data/reuters_dictionary.txt --clusters 16 --infile_W data/flatclust_init_w.csv --infile_H data/flatclust_init_h.csv --miniter 1 --algorithm HALS --maxiter 5000

if cmp -s "$xml_file" data/test/flatclust_rnd_clusters_16.xml; then
    echo "XML file test passed"
else
    echo "XML file test failed"
fi

if cmp -s "$assign_file" data/test/flatclust_rnd_assignments_16.csv; then
    echo "assignment file test passed"
else
    echo "assignment file test failed"
fi

rm -f "$xml_file"
rm -f "$assign_file"
