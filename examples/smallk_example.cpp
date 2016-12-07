#include <string>
#include <limits>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include "smallk.hpp"

using std::cout;
using std::cerr;
using std::endl;

const std::string FILENAME_W("nmf_rank2_init_w.csv");
const std::string FILENAME_H("nmf_rank2_init_h.csv");
const std::string FILENAME_MATRIX("reuters.mtx");
const std::string FILENAME_DICT("reuters_dictionary.txt");

static const char PATH_SEP = '/';

void MsgBox(const std::string& msg);

//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    // This example uses several matrices from the data directory in the 
    // smallk distribution.  The path to the data directory must be provided 
    // on the command line.

    if (1 == argc)
    {
        cerr << "usage: " << argv[0] << " <path_to_data_dir> " << endl;
        return -1;
    }

    // add a trailing '/' to make path construction easier below
    std::string data_dir(argv[1]);
    if (data_dir.size() >= 1)
    {
        auto sz = data_dir.size();
        if (PATH_SEP != data_dir[sz-1])
            data_dir += PATH_SEP;
    }

    // construct fully-qualified paths to matrices and dictionary
    std::string filepath_w      = data_dir + FILENAME_W;
    std::string filepath_h      = data_dir + FILENAME_H;
    std::string filepath_dict   = data_dir + FILENAME_DICT;
    std::string filepath_matrix = data_dir + FILENAME_MATRIX;

    // always call Initialize first
    smallk::Initialize(argc, argv);
    assert(smallk::IsInitialized());

    // here's how to get the smallk library version
    cout << "Smallk major version: " << smallk::GetMajorVersion() << endl;
    cout << "Smallk minor version: " << smallk::GetMinorVersion() << endl;
    cout << "Smallk patch level:   " << smallk::GetPatchLevel() << endl;
    cout << "Smallk version string: " << smallk::GetVersionString() << endl;

    try
    {
        ///////////////////////////////////////////////////////////////////////
        //
        // Load the 12411 x 7984 Reuters matrix, which will be used for all 
        // examples below.
        //
        ///////////////////////////////////////////////////////////////////////

        smallk::LoadMatrix(filepath_matrix);
        assert(smallk::IsMatrixLoaded());
        
        ///////////////////////////////////////////////////////////////////////
        //
        // Perform a nonnegative factorization of the Reuters matrix using the
        // default NMF-BPP algorithm.  Denoting the Reuters matrix by A, the 
        // Nmf function computes matrices W and H such that
        //
        //                           A ~ WH
        //
        // If A is mxn, then W is mxk and H is kxn.  The value of k is used
        // as an argument to the factorization routine.  For this example we
        // let k == 32.  Since the BPP algorithm is the default algorithm it
        // does not need to be specified explicitly as an argument.  We will 
        // also not specify initializer matrices for W and H, so that the
        // smallk code will randomly initialize them.
        //
        // The computed matrices W and H will be written to files w.csv and
        // h.csv in the output folder, which defaults to the current directory.
        //
        ///////////////////////////////////////////////////////////////////////

        MsgBox("Running NMF-BPP using k=32");

        smallk::Nmf(32);

        ///////////////////////////////////////////////////////////////////////
        //
        // Let's try that again, this time using the NMF-HALS algorithm with a 
        // k-value of 16.  The algorithm must be explicitly specified as an
        // argument since it is not the default.
        //
        ///////////////////////////////////////////////////////////////////////

        MsgBox("Running NMF-HALS using k=16");

        smallk::Nmf(16, smallk::HALS);

        ///////////////////////////////////////////////////////////////////////
        //
        // How about another nonnegative factorization, this time with k==2?
        // The RANK2 algorithm is specialized for k==2, so we'll use it.
        // Also, we will use initializer matrices for W and H, instead of
        // the random initialization used previously.
        //
        ///////////////////////////////////////////////////////////////////////
        
        MsgBox("Running NMF-RANK2 with W and H initializers");

        smallk::Nmf(2, smallk::RANK2, filepath_w, filepath_h);

        ///////////////////////////////////////////////////////////////////////
        //
        // Let's repeat the Rank2 factorization with different parameters.  
        // We will use a tolerance of 1.0e-5 instead of the default 0.005.  
        // The solver should require more iterations with the tighter 
        // tolerance, which you can verify from the output.
        //
        ///////////////////////////////////////////////////////////////////////

        MsgBox("Repeating the previous run with tol = 1.0e-5");

        smallk::SetNmfTolerance(1.0e-5);
        smallk::Nmf(2, smallk::RANK2, filepath_w, filepath_h);
        
        ///////////////////////////////////////////////////////////////////////
        //
        // Now suppose we want to examine the elements of matrices W and H.
        // The smallk library provides a READONLY pointer to the contents of
        // these matrices, which can be obtained by calling either 
        // LockedBufferW or LockedBufferH.
        // 
        // Each matrix is stored in column-major order, in a single contiguous
        // buffer.  The height of the buffer is called the 'leading dimension',
        // and it can be different from the height of the actual matrix.
        // The values of the buffer's leading dimension, the height of the 
        // matrix, and the width of the matrix are returned by the call to
        // LockedBuffer{W,H}.
        //
        // To illustrate how to use this information, we will examine the
        // contents of the W matrix and find the maximum and minimum values.
        //
        ///////////////////////////////////////////////////////////////////////
        
        unsigned int ldim, height, width;
        double min_w = std::numeric_limits<double>::max(), max_w = 0.0;

        // get a READONLY pointer to the buffer for matrix W
        const double* buf_w = smallk::LockedBufferW(ldim, height, width);

        // Iterate over the elements of W, looking for the max and min.
        // Since the matrix is stored in column-major order, it is most
        // efficient to scan the elements column-by-column.  The variable
        // 'c' is the column index and 'r' is the row index.
        for (unsigned int c=0; c != width; ++c)
        {
            // offset to the 0th element in column c
            unsigned int offset = c*ldim;

            // scan all elements in column c
            for (unsigned int r=0; r != height; ++r)
            {
                // the value of W(r, c)
                double elt = buf_w[offset + r];

                // check min and max
                if (elt < min_w)
                    min_w = elt;
                if (elt > max_w)
                    max_w = elt;
            }
        }

        cout << "Minimum value in W matrix: " << min_w << "." << endl;
        cout << "Maximum value in W matrix: " << max_w << "." << endl;

        ///////////////////////////////////////////////////////////////////////
        //
        // Now let's run a hierarchical clustering problem on the Reuters 
        // matrix and generate a total of 5 clusters.  To do this, we need
        // to load the Reuters dictionary.  The Reuters matrix is still loaded
        // so no need to load that again.
        //
        // We will use the default JSON output format for the clustering
        // results. 
        //
        ///////////////////////////////////////////////////////////////////////

        MsgBox("Running HierNMF2 with 5 clusters, JSON format");
        
        smallk::LoadDictionary(filepath_dict);
        smallk::HierNmf2(5);
        
        // The call to HierNmf2 will generate two output files in the current 
        // directory: 'assignments_5.csv' and 'tree_5.json'.  The assignment 
        // file contains an integer label for each of the 7984 documents in the
        // Reuters data set.  Documents that share a common label are in
        // the same cluster.  Any documents assigned the label '-1' are 
        // outliers and were not assigned to any cluster.

        // The JSON file contains the complete hierarchical factorization tree.
        // The leaf nodes have -1 for their left and right children.  For each
        // node you can see information on where the node is located in the
        // tree, the top terms for that node, etc.  This file contains enough
        // information to unambiguously reconstruct the factorization tree.

        ///////////////////////////////////////////////////////////////////////
        //
        // Now let's repeat the hierarchical clustering run, this time 
        // generating 10 clusters with 12 top terms per node.  We will also
        // use XML format for the cluster result file instead of JSON.
        //
        // We do not need to load either the dictionary or the Reuters matrix, 
        // since they are already loaded.
        //
        ///////////////////////////////////////////////////////////////////////
        
        MsgBox("Running HierNMF2 with 10 clusters, 12 terms, XML format");

        smallk::SetMaxTerms(12);
        smallk::SetOutputFormat(smallk::XML);
        smallk::HierNmf2(10);

        // This time two more files will be generated: 'assignments_10.csv' and
        // 'tree_10.xml'.  The tree is likely to be different (since the 
        // hierNMF2 algorithm uses random initialization internally).  Also 
        // note that each node now has 12 terms per node.

        ///////////////////////////////////////////////////////////////////////
        //
        // For the final example, we will repeat the clustering run, but this
        // time we generate 18 clusters with 8 terms per node.  We will also 
        // generate a flat clustering result.
        //
        ///////////////////////////////////////////////////////////////////////

        MsgBox("Running HierNmf2 with 18 clusters, 8 terms, with flat");

        smallk::SetMaxTerms(8);
        smallk::HierNmf2WithFlat(18);

        // Four new output files will be generated.  They are:
        //
        // 'assignments_18.csv' : assignments from hierarchical clustering
        // 'assignments_flat_18.csv' : assignments from flat clustering
        // 'tree_18.xml' : the hierarchical factorization tree
        // 'clusters_18.xml' : flat clustering results
    }
    catch (std::exception& e)
    {
        cerr << e.what() << endl;
    }
    
    smallk::Finalize();
    return 0;
}

//-----------------------------------------------------------------------------
void MsgBox(const std::string& msg)
{
    // This function prints a message enclosed in a starred box.

    const std::string STARS("************************************************************");
    const size_t WIDTH = STARS.size();

    size_t len = msg.size();
    size_t pad = (WIDTH - 2 - len)/2;

    cout << endl;
    cout << STARS << endl;
    cout << "*";
    for (size_t i=0; i<(WIDTH-2); ++i)
        cout << " ";
    cout << "*" << endl;

    cout << "*";
    for (size_t i=0; i<pad; ++i)
        cout << " ";
    cout << msg;
    for (size_t i=(1+pad+len); i<WIDTH-1; ++i)
        cout << " ";
    cout << "*" << endl;

    cout << "*";
    for (size_t i=0; i<(WIDTH-2); ++i)
        cout << " ";
    cout << "*" << endl;
    cout << STARS << endl;
}
