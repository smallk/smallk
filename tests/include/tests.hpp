#pragma once

#include <string>

bool TestRank2SystemSolve();
bool TestSparseGemm();
bool TestSparseGemmIndexed();
bool TestIndexedRank2Nmf();
bool TestDenseNmf();
bool TestOpenMP();
bool TestBpp(const std::string& data_dir);
bool TestParallelInit();
