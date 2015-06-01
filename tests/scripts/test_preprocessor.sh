# Test the preprocessor
#
# From the main xdata3 folder run the following command:
#
#       sh tests/scripts/test_preprocessor.sh <path_to_xdata_data>

if [ $# != 1 ]; then
    echo "Usage: from the main xdata3 folder, run the following command:"
    echo " "
    echo "     sh tests/scripts/test_preprocessor.sh <path_to_xdata_data>"
    echo " "
    exit 1
fi

preprocessor/bin/preprocess_tf --indir $1

if cmp -s reduced_matrix.mtx $1/test/reduced_matrix_20news.mtx; then
    echo "matrix test passed"
else
    echo "matrix test failed"
fi

if cmp -s reduced_dictionary.txt $1/test/reduced_dictionary_20news.txt; then
    echo "dictionary file test passed"
else
    echo "dictionary file test failed"
fi

if cmp -s reduced_documents.txt $1/test/reduced_documents_20news.txt; then
    echo "documents file test passed"
else
    echo "documents file test failed"
fi

rm -f reduced_matrix.mtx
rm -f reduced_dictionary.txt
rm -f reduced_documents.txt
