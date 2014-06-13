# Test 1: hierarchical clustering

xml_file=tree_80.xml
assign_file=assignments_80.csv

if [ -f "$xml_file" ]
then
    rm -f "$xml_file"
fi

if [ -f "$assign_file" ]
then
    rm -f "$assign_file"
fi

./bin/hierclust --matrixfile ../data/test/reduced_matrix_20news.mtx --dictfile ../data/test/reduced_dictionary_20news.txt --clusters 80 --infile_W ../data/test/nmf_20news_w_init.csv --infile_H ../data/test/nmf_20news_h_init.csv --miniter 1

if cmp -s "$xml_file" ../data/test/nmf_20news_tree_80.xml; then
    echo "XML file test passed"
else
    echo "XML file test failed"
fi

if cmp -s "$assign_file" ../data/test/nmf_20news_assignments_80.csv; then
    echo "assignment file test passed"
else
    echo "assignment file test failed"
fi

rm -f "$xml_file"
rm -f "$assign_file"
