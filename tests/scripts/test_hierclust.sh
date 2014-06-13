# Test 1: hierarchical clustering
 
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

./bin/hierclust --matrixfile ../data/reuters.mtx --dictfile ../data/reuters_dictionary.txt --clusters 5 --infile_W ../data/hierclust_init_w.csv --infile_H ../data/hierclust_init_h.csv --miniter 1

if cmp -s "$xml_file" ../data/test/reuters_tree_5.xml; then
    echo "XML file test passed"
else
    echo "XML file test failed"
fi

if cmp -s "$assign_file" ../data/test/reuters_assignments_5.csv; then
    echo "assignment file test passed"
else
    echo "assignment file test failed"
fi

rm -f "$xml_file"
rm -f "$assign_file"

