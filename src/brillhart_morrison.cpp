#include "brillhart_morrison.h"
#include <utility>
#include <cmath>

long long mod(long long a, long long p)
{
    return (a % p + p) % p;
}

int legendreSymbol(long long a, long long p)
{
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
    if (n <= 2)
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