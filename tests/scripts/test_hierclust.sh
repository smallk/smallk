# Test hierarchical clustering
#
# From the main xdata3 folder run the following command:
#
#       sh tests/scripts/test_hierclust.sh <path_to_xdata_data>

if [ $# != 1 ]; then
    echo "Usage: from the main xdata3 folder, run the following command:"
    echo " "
    echo "     sh tests/scripts/test_hierclust.sh <path_to_xdata_data>"
    echo " "
    exit 1
fi

# Reuters matrix, 12 clusters ##

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

hierclust/bin/hierclust --matrixfile $1/reuters.mtx --dictfile $1/reuters_dictionary.txt --clusters 12 --initdir $1/test/matrices.reuters --miniter 1

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

## 20News matrix, 15 clusters ##

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

hierclust/bin/hierclust --matrixfile $1/news20.mtx --dictfile $1/news20_dictionary.txt --clusters 15 --initdir $1/test/matrices.20news --miniter 1

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
