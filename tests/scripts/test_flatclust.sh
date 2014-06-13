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

./bin/flatclust --matrixfile ../data/rnd_256_256.csv --dictfile ../data/reuters_dictionary.txt --clusters 16 --infile_W ../data/flatclust_init_w.csv --infile_H ../data/flatclust_init_h.csv --miniter 1 --algorithm HALS --maxiter 5000

if cmp -s "$xml_file" ../data/test/flatclust_rnd_clusters_16.xml; then
    echo "XML file test passed"
else
    echo "XML file test failed"
fi

if cmp -s "$assign_file" ../data/test/flatclust_rnd_assignments_16.csv; then
    echo "assignment file test passed"
else
    echo "assignment file test failed"
fi

rm -f "$xml_file"
rm -f "$assign_file"
