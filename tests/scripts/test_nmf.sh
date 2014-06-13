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

./bin/nmf --matrixfile ../data/reuters.mtx --algorithm BPP --k 8 --infile_W ../data/nmf_init_w.csv --infile_H ../data/nmf_init_h.csv --outprecision 6 --miniter 1

if cmp -s w.csv ../data/test/nmf_result_w.csv; then
    echo "W matrix test passed"
else
    echo "W matrix test failed"
fi

if cmp -s h.csv ../data/test/nmf_result_h.csv; then
    echo "H matrix test passed"
else
    echo "H matrix test failed"
fi

rm -f "$w_file"
rm -f "$h_file"
