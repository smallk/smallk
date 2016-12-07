# Test NMF
#
# From the main xdata3 folder run the following command:
#
#       sh tests/scripts/test_nmf.sh <path_to_xdata_data>

if [ $# != 1 ]; then
    echo "Usage: from the main xdata3 folder, run the following command:"
    echo " "
    echo "     sh tests/scripts/test_nmf.sh <path_to_xdata_data>"
    echo " "
    exit 1
fi

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

nmf/bin/nmf --matrixfile $1/reuters.mtx --algorithm BPP --k 8 --infile_W $1/nmf_init_w.csv --infile_H $1/nmf_init_h.csv --outprecision 6 --miniter 1

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

rm -f "$w_file"
rm -f "$h_file"
