
#include "brillhart_morrison.h"

#include <utility>
#include <cmath>
#include <vector>
#include <limits>

long long mod(long long a, long long p)
{
    return (a % p + p) % p;
}

int legendreSymbol(long long a, long long p)
{
    if (p <= 2 || p % 2 == 0)
        return 0;

    a = mod(a, p);

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

    return p == 1 ? result : 0;
}

long double calculateFactorBaseLimit(long long n)
{
    if (n <= 3)
        return 0;

    long double value = static_cast<long double>(n);
    long double a = 1.0L / sqrtl(2.0L);

    return expl(a * sqrtl(logl(value) * logl(logl(value))));
}

bool isSmallPrime(long long n)
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

    factorBase.push_back(-1);
    factorBase.push_back(2);

    long double limit = calculateFactorBaseLimit(n);

    for (long long p = 3; p < limit; p += 2)
    {
        if (isSmallPrime(p) && legendreSymbol(n, p) == 1)
        {
            factorBase.push_back(p);
        }
    }

    return factorBase;
}

std::vector<long long> calculateContinuedFraction(long long n, int terms)
{
    std::vector<long long> coefficients;

    if (n <= 0 || terms <= 0)
        return coefficients;

    long long a0 = static_cast<long long>(sqrtl(static_cast<long double>(n)));

    if (a0 * a0 == n)
        return coefficients;

    long long m = 0;
    long long d = 1;
    long long a = a0;

    coefficients.push_back(a0);

    for (int i = 1; i < terms; i++)
    {
        __int128 nextM = static_cast<__int128>(d) * a - m;

        if (nextM < 0 || nextM > std::numeric_limits<long long>::max())
            break;

        m = static_cast<long long>(nextM);

        __int128 numerator = static_cast<__int128>(n) -
                             static_cast<__int128>(m) * m;

        if (d == 0 || numerator <= 0 || numerator % d != 0)
            break;

        __int128 nextD = numerator / d;

        if (nextD <= 0 || nextD > std::numeric_limits<long long>::max())
            break;

        d = static_cast<long long>(nextD);

        __int128 nextA = (static_cast<__int128>(a0) + m) / d;

        if (nextA < 0 || nextA > std::numeric_limits<long long>::max())
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
