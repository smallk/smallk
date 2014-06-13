#include <algorithm>
#include "random.hpp"
#include "index_set.hpp"
#include "fisher_yates_shuffle.hpp"

//-----------------------------------------------------------------------------
void RandomIndexSet(Random& rng, const unsigned int n,
                    const unsigned int index_set_size,
                    std::vector<unsigned int>& index_set)
{
    std::vector<unsigned int> all_indices(n);

    for (size_t q=0; q != n; ++q)
        all_indices[q] = q;

    FisherYatesShuffle(all_indices.begin(), all_indices.end(), rng);

    index_set.resize(index_set_size);
    for (size_t q=0; q != index_set_size; ++q)
        index_set[q] = all_indices[q];

    std::sort(index_set.begin(), index_set.end());
}

