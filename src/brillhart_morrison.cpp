
#include "brillhart_morrison.h"

#include <utility>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

static long long modLong(long long a, long long p)
{
    return (a % p + p) % p;
}

static cpp_int gcdBig(cpp_int a, cpp_int b)
{
    if (a < 0)
        a = -a;

    if (b < 0)
        b = -b;

    while (b != 0)
    {
        cpp_int r = a % b;
        a = b;
        b = r;
    }

    return a;
}

static cpp_int modBig(cpp_int a, const cpp_int& n)
{
    a %= n;

    if (a < 0)
        a += n;

    return a;
}

static cpp_int modPowBig(cpp_int base, int exponent, const cpp_int& modulus)
{
    cpp_int result = 1;
    base = modBig(base, modulus);

    while (exponent > 0)
    {
        if (exponent % 2 == 1)
            result = (result * base) % modulus;

        base = (base * base) % modulus;
        exponent /= 2;
    }

    return result;
}

int legendreSymbol(long long a, long long p)
{
    if (p <= 2 || p % 2 == 0)
        return 0;

    a = modLong(a, p);

    if (a == 0)
        return 0;

    long long result = 1;

    while (a != 0)
    {
        while (a % 2 == 0)
        {
            a /= 2;

            if (p % 8 == 3 || p % 8 == 5)
                result = -result;
        }

        std::swap(a, p);

        if (a % 4 == 3 && p % 4 == 3)
            result = -result;

        a %= p;
    }

    return p == 1 ? static_cast<int>(result) : 0;
}

long double calculateFactorBaseLimit(long long n)
{
    if (n <= 3)
        return 0;

    long double value = static_cast<long double>(n);
    long double a = 1.0L / sqrtl(2.0L);

    return expl(a * sqrtl(logl(value) * logl(logl(value))));
}

static bool isSmallPrime(long long n)
{
    if (n < 2)
        return false;

    if (n == 2)
        return true;

    if (n % 2 == 0)
        return false;

    for (long long i = 3; i <= n / i; i += 2)
    {
        if (n % i == 0)
            return false;
    }

    return true;
}

std::vector<long long> buildFactorBase(long long n)
{
    std::vector<long long> factorBase;

    if (n <= 1)
        return factorBase;

    factorBase.push_back(-1);
    factorBase.push_back(2);

    long double limit = calculateFactorBaseLimit(n);

    for (long long p = 3; p < limit; p += 2)
    {
        if (isSmallPrime(p) && legendreSymbol(n, p) == 1)
            factorBase.push_back(p);
    }

    return factorBase;
}

std::vector<long long> calculateContinuedFraction(long long n, int terms)
{
    std::vector<long long> coefficients;

    if (n <= 0 || terms <= 0)
        return coefficients;

    long long a0 = static_cast<long long>(
        sqrtl(static_cast<long double>(n))
    );

    if (a0 * a0 == n)
        return coefficients;

    long long m = 0;
    long long d = 1;
    long long a = a0;

    coefficients.push_back(a0);

    for (int i = 1; i < terms; i++)
    {
        __int128 nextM = static_cast<__int128>(d) * a - m;

        if (nextM < 0 ||
            nextM > std::numeric_limits<long long>::max())
            break;

        m = static_cast<long long>(nextM);

        __int128 numerator =
            static_cast<__int128>(n) -
            static_cast<__int128>(m) * m;

        if (d == 0 || numerator <= 0 || numerator % d != 0)
            break;

        __int128 nextD = numerator / d;

        if (nextD <= 0 ||
            nextD > std::numeric_limits<long long>::max())
            break;

        d = static_cast<long long>(nextD);

        __int128 nextA =
            (static_cast<__int128>(a0) + m) / d;

        if (nextA < 0 ||
            nextA > std::numeric_limits<long long>::max())
            break;

        a = static_cast<long long>(nextA);

        coefficients.push_back(a);
    }

    return coefficients;
}

void calculateConvergents(
    const std::vector<long long>& coefficients,
    std::vector<cpp_int>& numerators,
    std::vector<cpp_int>& denominators
)
{
    numerators.clear();
    denominators.clear();

    cpp_int pPrev2 = 0;
    cpp_int pPrev1 = 1;

    cpp_int qPrev2 = 1;
    cpp_int qPrev1 = 0;

    for (long long a : coefficients)
    {
        cpp_int p = cpp_int(a) * pPrev1 + pPrev2;
        cpp_int q = cpp_int(a) * qPrev1 + qPrev2;

        numerators.push_back(p);
        denominators.push_back(q);

        pPrev2 = pPrev1;
        pPrev1 = p;

        qPrev2 = qPrev1;
        qPrev1 = q;
    }
}

std::vector<cpp_int> calculateRelations(
    long long n,
    const std::vector<cpp_int>& numerators,
    const std::vector<cpp_int>& denominators
)
{
    std::vector<cpp_int> relations;

    if (n <= 0 || numerators.size() != denominators.size())
        return relations;

    for (size_t i = 0; i < numerators.size(); i++)
    {
        cpp_int p = numerators[i];
        cpp_int q = denominators[i];

        cpp_int value = p * p - cpp_int(n) * q * q;

        relations.push_back(value);
    }

    return relations;
}

bool factorSmoothRelation(
    cpp_int value,
    const std::vector<long long>& factorBase,
    std::vector<int>& exponents
)
{
    exponents.assign(factorBase.size(), 0);

    if (value == 0 || factorBase.size() < 2 || factorBase[0] != -1)
        return false;

    if (value < 0)
    {
        exponents[0] = 1;
        value = -value;
    }

    for (size_t i = 1; i < factorBase.size(); i++)
    {
        long long p = factorBase[i];

        if (p <= 1)
            continue;

        while (value % p == 0)
        {
            value /= p;
            exponents[i]++;
        }
    }

    return value == 1;
}

std::vector<std::vector<int>> buildParityMatrix(
    const std::vector<std::vector<int>>& exponents
)
{
    std::vector<std::vector<int>> matrix;

    if (exponents.empty())
        return matrix;

    size_t cols = exponents[0].size();

    for (const auto& relation : exponents)
    {
        if (relation.size() != cols)
            return {};

        std::vector<int> row;

        for (int exponent : relation)
            row.push_back((exponent % 2 + 2) % 2);

        matrix.push_back(row);
    }

    return matrix;
}

std::vector<std::vector<int>> findDependencies(
    const std::vector<std::vector<int>>& matrix
)
{
    std::vector<std::vector<int>> dependencies;

    if (matrix.empty())
        return dependencies;

    size_t rows = matrix.size();
    size_t cols = matrix[0].size();

    if (cols == 0)
        return dependencies;

    std::vector<std::vector<int>> a(
        cols, std::vector<int>(rows, 0)
    );

    for (size_t i = 0; i < rows; i++)
    {
        if (matrix[i].size() != cols)
            return {};

        for (size_t j = 0; j < cols; j++)
            a[j][i] = (matrix[i][j] % 2 + 2) % 2;
    }

    std::vector<size_t> pivotColumns;
    size_t pivotRow = 0;

    for (size_t col = 0; col < rows && pivotRow < cols; col++)
    {
        size_t selected = pivotRow;

        while (selected < cols && a[selected][col] == 0)
            selected++;

        if (selected == cols)
            continue;

        std::swap(a[pivotRow], a[selected]);

        for (size_t i = 0; i < cols; i++)
        {
            if (i != pivotRow && a[i][col] == 1)
            {
                for (size_t j = 0; j < rows; j++)
                    a[i][j] ^= a[pivotRow][j];
            }
        }

        pivotColumns.push_back(col);
        pivotRow++;
    }

    std::vector<bool> isPivot(rows, false);

    for (size_t col : pivotColumns)
        isPivot[col] = true;

    for (size_t freeCol = 0; freeCol < rows; freeCol++)
    {
        if (isPivot[freeCol])
            continue;

        std::vector<int> dependency(rows, 0);
        dependency[freeCol] = 1;

        for (size_t i = 0; i < pivotColumns.size(); i++)
            dependency[pivotColumns[i]] = a[i][freeCol];

        dependencies.push_back(dependency);
    }

    return dependencies;
}

long long findFactorFromDependencies(
    long long n,
    const std::vector<long long>& factorBase,
    const std::vector<cpp_int>& numerators,
    const std::vector<size_t>& smoothIndices,
    const std::vector<std::vector<int>>& smoothExponents,
    const std::vector<std::vector<int>>& dependencies
)
{
    if (n <= 1 ||
        factorBase.size() < 2 ||
        smoothIndices.size() != smoothExponents.size())
        return 0;

    cpp_int N = n;

    std::vector<std::vector<int>> candidates = dependencies;

    const size_t maxDependencies = 64;
    size_t limit = std::min(dependencies.size(), maxDependencies);

    for (size_t i = 0; i < limit; i++)
    {
        if (dependencies[i].size() != smoothExponents.size())
            continue;

        for (size_t j = i + 1; j < limit; j++)
        {
            if (dependencies[j].size() != smoothExponents.size())
                continue;

            std::vector<int> combined(dependencies[i].size(), 0);

            for (size_t k = 0; k < combined.size(); k++)
                combined[k] = dependencies[i][k] ^ dependencies[j][k];

            candidates.push_back(combined);
        }
    }

    for (const auto& dependency : candidates)
    {
        if (dependency.size() != smoothExponents.size())
            continue;

        std::vector<long long> exponentSums(factorBase.size(), 0);

        cpp_int X = 1;
        bool valid = true;
        bool selected = false;

        for (size_t i = 0; i < dependency.size(); i++)
        {
            if (dependency[i] == 0)
                continue;

            if (dependency[i] != 1 ||
                smoothExponents[i].size() != factorBase.size() ||
                smoothIndices[i] >= numerators.size())
            {
                valid = false;
                break;
            }

            selected = true;

            X = (X * numerators[smoothIndices[i]]) % N;

            for (size_t j = 0; j < factorBase.size(); j++)
                exponentSums[j] += smoothExponents[i][j];
        }

        if (!valid || !selected)
            continue;

        for (long long exponent : exponentSums)
        {
            if (exponent % 2 != 0)
            {
                valid = false;
                break;
            }
        }

        if (!valid)
            continue;

        cpp_int Y = 1;

        for (size_t j = 1; j < factorBase.size(); j++)
        {
            if (factorBase[j] <= 1)
                continue;

            long long exponent = exponentSums[j] / 2;

            if (exponent > 0)
            {
                cpp_int term = modPowBig(
                    cpp_int(factorBase[j]),
                    static_cast<int>(exponent),
                    N
                );

                Y = (Y * term) % N;
            }
        }

        X = modBig(X, N);
        Y = modBig(Y, N);

        cpp_int d1 = gcdBig(X - Y, N);

        if (d1 > 1 && d1 < N)
            return d1.convert_to<long long>();

        cpp_int d2 = gcdBig(X + Y, N);

        if (d2 > 1 && d2 < N)
            return d2.convert_to<long long>();
    }

    return 0;
}