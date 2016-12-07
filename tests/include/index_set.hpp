#pragma once

#include <vector>

class Random;

void RandomIndexSet(Random& rng,
                    const unsigned int n, // max number of indices
                    const unsigned int index_set_size, // number in the set
                    std::vector<unsigned int>& index_set);

