./bin/preprocess_tf --indir ../data

if cmp -s reduced_matrix.mtx ../data/test/reduced_matrix_20news.mtx; then
    echo "matrix test passed"
else
    echo "matrix test failed"
fi

if cmp -s reduced_dictionary.txt ../data/test/reduced_dictionary_20news.txt; then
    echo "dictionary file test passed"
else
    echo "dictionary file test failed"
fi

if cmp -s reduced_documents.txt ../data/test/reduced_documents_20news.txt; then
    echo "documents file test passed"
else
    echo "documents file test failed"
fi

rm -f reduced_matrix.mtx
rm -f reduced_dictionary.txt
rm -f reduced_documents.txt
