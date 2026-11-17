
#ifndef BRILLHART_MORRISON_H
#define BRILLHART_MORRISON_H

#include <vector>
#include <cstddef>
#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;

int legendreSymbol(long long a, long long p);

long double calculateFactorBaseLimit(long long n);

std::vector<long long> buildFactorBase(long long n);

std::vector<long long> calculateContinuedFraction(long long n, int terms);

void calculateConvergents(
    const std::vector<long long>& coefficients,
    std::vector<cpp_int>& numerators,
    std::vector<cpp_int>& denominators
);

std::vector<cpp_int> calculateRelations(
    long long n,
    const std::vector<cpp_int>& numerators,
    const std::vector<cpp_int>& denominators
);

bool factorSmoothRelation(
    cpp_int value,
    const std::vector<long long>& factorBase,
    std::vector<int>& exponents
);

std::vector<std::vector<int>> buildParityMatrix(
    const std::vector<std::vector<int>>& exponents
);

std::vector<std::vector<int>> findDependencies(
    const std::vector<std::vector<int>>& matrix
);

long long findFactorFromDependencies(
    long long n,
    const std::vector<long long>& factorBase,
    const std::vector<cpp_int>& numerators,
    const std::vector<size_t>& smoothIndices,
    const std::vector<std::vector<int>>& smoothExponents,
    const std::vector<std::vector<int>>& dependencies
);

#endif
