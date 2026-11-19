
#include "pollard_rho.h"

#include <random>
#include <numeric>
#include <cstdlib>

using namespace std;

long long pollardRho(long long n)
{
    if (n < 4)
        return 0;

    if (n % 2 == 0)
        return 2;

    if (n % 3 == 0)
        return 3;

    mt19937_64 gen(random_device{}());
    uniform_int_distribution<long long> dist(2, n - 2);

    const int maxAttempts = 25;
    const int maxIterations = 100000;

    for (int attempt = 0; attempt < maxAttempts; ++attempt)
    {
        long long x = dist(gen);
        long long y = x;
        long long c = dist(gen);
        long long d = 1;

        auto f = [n, c](long long value)
        {
            return ((__int128)value * value + c) % n;
        };

        int iterations = 0;

        while (d == 1 && iterations < maxIterations)
        {
            x = f(x);
            y = f(f(y));

            d = gcd(abs(x - y), n);
            ++iterations;
        }

        if (d > 1 && d < n)
            return d;
    }

    return 0;
}
