# Test flat clustering
#
# From the main xdata3 folder run the following command:
#
#       sh tests/scripts/test_flatclust.sh <path_to_xdata_data>

if [ $# != 1 ]; then
    echo "Usage: from the main xdata3 folder, run the following command:"
    echo " "
    echo "     sh tests/scripts/test_flatclust.sh <path_to_xdata_data>"
    echo " "
    exit 1
fi

xml_file=clusters_16.xml
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

flatclust/bin/flatclust --matrixfile $1/rnd_256_256.csv --dictfile $1/reuters_dictionary.txt --clusters 16 --infile_W $1/flatclust_init_w.csv --infile_H $1/flatclust_init_h.csv --miniter 1 --algorithm HALS --maxiter 5000

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
