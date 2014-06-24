///
/// \file smallk.hpp
///
/// \brief The smallk API header file.
///
/// Copyright 2014 Georgia Institute of Technology
///
/// Licensed under the Apache License, Version 2.0 (the "License");
/// you may not use this file except in compliance with the License.
/// You may obtain a copy of the License at
///
///     http://www.apache.org/licenses/LICENSE-2.0
///
/// Unless required by applicable law or agreed to in writing, software
/// distributed under the License is distributed on an "AS IS" BASIS,
/// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
/// See the License for the specific language governing permissions and
/// limitations under the License.

#pragma once

#include <string>
#include <vector>

#define SMALLK_MAJOR_VERSION 1
#define SMALLK_MINOR_VERSION 0
#define SMALLK_PATCH_LEVEL   0

namespace smallk
{
    ///
    /// The NMF algorithms supported by this API.
    ///
    enum Algorithm
    {
        MU,      ///< Multiplicative Updating, Lee & Seung
        BPP,     ///< Block Principal Pivoting, Kim and Park
        HALS,    ///< Hierarchical Alternating Least Squares, Cichocki & Pan
        RANK2    ///< Rank2, Kuang and Park
    };

    ///
    /// The output file formats supported by the clustering code.
    ///
    enum OutputFormat
    {
        XML,  ///< Extensible Markup Language
        JSON  ///< JavaScript Object Notation
    };

    /// \name Initialization and Cleanup
    //@{
    ///
    /// \brief Perform initialization - call this function first.
    /// 
    /// Simply pass \c argc and \c argv from \c main to this function.
    ///
    /// \param[in] argc The number of command-line arguments.
    /// \param[in] argv The array of command-line argument strings.
    void Initialize(int& argc, char**& argv);

    ///
    /// \brief Determine whether the library has been initialized.
    /// \return A boolean indicating the initialization status.
    bool IsInitialized();

    ///
    /// \brief Perform cleanup - call this function last.
    void Finalize();
    //@}

    /// \name Versioning
    //@{
    ///
    /// \brief Returns the API major version.
    /// \return The major version number as an \c unsigned \c int.
    unsigned int GetMajorVersion();

    /// \brief Returns the API minor version.
    /// \return The minor version number as an \c unsigned \c int.
    unsigned int GetMinorVersion();

    /// \brief Returns the API patch level.
    /// \return The patch level as an \c unsigned \c int.
    unsigned int GetPatchLevel();

    /// \brief Returns the API version string.
    /// \return The version string formatted as \c major_version.minor_version.patch_level.
    std::string GetVersionString();
    //@}

    /// \name Common Functions
    //@{
    ///
    /// \brief Returns the precision with which output files will be written.
    /// \return The number of digits after the decimal point for elements of matrices \c W and \c H.
    unsigned int GetOutputPrecision();

    ///
    /// \brief Sets the precision with which output files will be written.
    /// \param[in] num_digits The desired number of digits after the decimal point.
    void SetOutputPrecision(const unsigned int num_digits = 6);

    ///
    /// \brief Returns the tolerance value used to terminate iterative NMF algorithms.
    /// \return The tolerance as a double-precision floating point number.
    /// 
    /// The NMF algorithms used by smallk are iterative.  At each iteration, a metric
    /// is computed that estimates the progress made towards a solution.  The metrics
    /// are designed to decrease as the computation progresses.  The 'tolerance' parameter
    /// is the cutoff value for this metric, meaning that whenever the metric decreases 
    /// to this value the iterations are terminated and the solution is said to have converged.
    /// There is a tradeoff between the maximum number of iterations and the tolerance.
    /// Smaller values of the tolerance generally require more iterations.
    double GetNmfTolerance();

    ///
    /// \brief Sets the cutoff tolerance for NMF iterations.
    /// \param[in] tol The tolerance value.
    void SetNmfTolerance(const double tol = 0.005);

    ///
    /// \brief Returns the maximum number of iterations allowed for either NMF or clustering.
    /// \return The iteration count as an unsigned integer.
    unsigned int GetMaxIter();

    ///
    /// \brief Sets the upper limit to the number of iterations for NMF and clustering.
    /// \param[in] max_iterations The maximum number of iterations to perform.
    void SetMaxIter(const unsigned int max_iterations = 5000);

    ///
    /// \brief Returns the minimum number of iterations performed by either the NMF or clustering code.
    /// \return The iteration count as an unsigned integer.
    unsigned int GetMinIter();

    ///
    /// \brief Sets the minimum number of iterations to be performed by either the NMF or clustering code.
    /// \param[in] min_iterations The minimum number of iterations to perform.
    void SetMinIter(const unsigned int min_iterations = 5);

    ///
    /// \brief Returns the maximum number of threads used for either the NMF or clustering computations.
    /// \return The thread count as an unsigned integer.
    unsigned int GetMaxThreads();

    ///
    /// \brief Sets an upper limit to the number of threads used by either the NMF or clustering code.
    /// \param[in] max_threads The maximum number of threads to use.
    ///
    /// The smallk code is multithreaded for improved performance.  The default upper limit 
    /// to the number of threads used by the code is determined by the hardware's capabilities.
    /// Use this function to override this limit and restrict the thread count to a smaller range.  
    /// Allowable values for max_threads are in the interval [1, hardware_max_threads].  
    void SetMaxThreads(const unsigned int max_threads);

    ///
    /// \brief Sets all params to default values and clears any loaded matrices.
    void Reset();

    ///
    /// \brief Seed the random number generator with the specified value.
    /// \param[in] seed An integer seed value.
    void SeedRNG(const int seed);

    ///
    /// \brief Load the matrix to be factored.
    /// \param[in] filepath The absolute or relative path to a matrix file.
    ///
    /// The smallk code can read matrix files in either CSV (comma separated value) 
    /// or MTX (MatrixMarket) format.  CSV files are assumed to contain dense
    /// matrices, i.e. the CSV files must contain an entry for each matrix element.
    /// MTX files are assumed to contain sparse matrices only.
    void LoadMatrix(const std::string& filepath);

    ///
    /// \brief Returns whether a CSV or MTX file has been successfully loaded.
    /// \return A boolean indicating the presence or absence of a loaded matrix.
    bool IsMatrixLoaded();

    ///
    /// \brief Returns the output directory where result files will be written.
    /// \return A string representing the output directory.
    std::string GetOutputDir();

    ///
    /// \brief Sets the directory into which result files will be written.
    /// \param[in] outdir The absolute or relative path to an output directory.
    ///
    /// The otuput directory must exist - the smallk code will not create it.
    void SetOutputDir(const std::string& outdir);

    //@}


    /// \name NMF
    //@{
    ///
    /// \brief Runs the selected NMF algorithm on the loaded matrix.
    /// \param[in] k The width of matrix W or the height of matrix H.
    /// \param[in] algorithm The desired NMF algorithm.
    /// \param[in] initfile_w A CSV file containing initial values for \c m x \c k matrix W.
    /// \param[in] initfile_h A CSV file containing initial values for \c k x \c n matrix H.
    /// This function performs the factorization \c A \c ~ \c W \c H, where
    /// matrix A is \c m x \c n, matrix W is \c m x \c k, and matrix H is
    /// \c k x \c n.  Matrix A must have already been loaded by a prior 
    /// call to \c LoadMatrix.
    void Nmf(const unsigned int k, 
             const Algorithm algorithm = Algorithm::BPP,
             const std::string& initfile_w = std::string(""),
             const std::string& initfile_h = std::string(""));

    ///
    /// \brief Returns a pointer to the buffer containing the NMF result matrix W.
    /// \param[out] ldim The leading dimension of the W-buffer.
    /// \param[out] height The height of the W matrix.
    /// \param[out] width The width of the W matrix.
    const double* LockedBufferW(unsigned int& ldim, 
                                unsigned int& height,
                                unsigned int& width);

    ///
    /// \brief Returns a pointer to the buffer containing the NMF result matrix H.
    /// \param[out] ldim The leading dimension of the H-buffer.
    /// \param[out] height The height of the H matrix.
    /// \param[out] width The width of the H matrix.
    const double* LockedBufferH(unsigned int& ldim,
                                unsigned int& height,
                                unsigned int& width);

    //@}

    /// \name Clustering
    //{@
    ///
    /// \brief Load the file containing the dictionary strings.
    /// \param[in] filepath The absolute or relative path to the dictionary file.
    /// Each row of the term-frequency matrix represents a 'term' in the dictionary.
    /// The dictionary file is the collection of all such terms, one per line.
    /// A 'term' can consist of multiple words.  The smallk code considers each individual
    /// line in the file a separate term.
    ///
    void LoadDictionary(const std::string& filepath);

    ///
    /// \brief Reuturns the maximum number of terms computed per node during clustering.
    /// \return The term count as an unsigned int.
    unsigned int GetMaxTerms();

    /// 
    /// \brief Set the maximum number of terms to compute per node during clustering.
    /// \param[in] max_terms The term count.
    void SetMaxTerms(const unsigned int max_terms = 5);

    ///
    /// \brief Returns the format used for writing clustering results.
    /// \returns \c OutputFormat::XML for XML format, or \c OutputFormat::JSON for JSON format.
    OutputFormat GetOutputFormat();

    ///
    /// \brief Sets the file format for clustering results.
    /// \param[in] format Either \c OutputFormat::XML for XML format, or 
    /// \c OutputFormat::JSON for JSON format.
    void SetOutputFormat(const OutputFormat format = OutputFormat::JSON);

    ///
    /// \brief Returns the tolerance for the NMF subproblems solved by HierNMF2.
    /// \return The tolerance as a double-precision value.
    double GetHierNmf2Tolerance();

    ///
    /// \brief Sets the tolerance for the NMF subproblems solved by HierNMF2.
    /// \param[in] tol The tolerance value, a double precision value on (0.0, 1.0).
    void SetHierNmf2Tolerance(const double tol = 0.0001);

    ///
    /// \brief Performs hierarchical clustering for the loaded matrix.
    /// \param[in] num_clusters The number of clusters (leaf nodes) to generate.
    /// 
    /// This function runs the HierNMF2 algorithm on the loaded matrix.  Two result files
    /// are written to the output directory, one containing the cluster labels for each
    /// document.  A value of -1 in this file represents an outlier, which is a document
    /// that was not assigned to a cluster.  The other output file contains information
    /// for each node in the factorization tree.
    void HierNmf2(const unsigned int num_clusters);

    ///
    /// \brief Performs hierarchical and flat clustering for the loaded matrix.
    /// \param[in] num_clusters The number of clusters to generate.
    ///
    /// This function performs the same actions as HierNmf2, but it also generates a flat
    /// clustering result, in addition the the hierarchical clustering result.  Two additional
    /// result files are written to the output directory, one for the flat clustering assignments,
    /// and anther for the flat clustering nodes.
    void HierNmf2WithFlat(const unsigned int num_clusters);

    //@}

} // namespace smallk
