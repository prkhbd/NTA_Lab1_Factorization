
#include "primality.h"

#include <iostream>
#include <random>
#include <numeric>
#include <cstdlib>

using namespace std;

long long modPow(long long a, long long d, long long n)
{
    long long result = 1;
    a %= n;

    while (d > 0)
    {
        if (d % 2 == 1)
        {
            result = (__int128)result * a % n;
        }

        a = (__int128)a * a % n;
        d /= 2;
    }

    return result;
}

bool isPrimeMillerRabin(long long n, int k)
{
    if (n < 2)
        return false;

    if (n == 2 || n == 3)
        return true;

    if (n % 2 == 0)
        return false;

    long long d = n - 1;
    int s = 0;

    while (d % 2 == 0)
    {
        d /= 2;
        s++;
    }

    random_device rd;
    mt19937_64 gen(rd());

    uniform_int_distribution<long long> dist(2, n - 2);

    for (int i = 0; i < k; i++)
    {
        long long a = dist(gen);

        if (gcd(a, n) != 1)
            return false;

        long long x = modPow(a, d, n);

        if (x == 1 || x == n - 1)
            continue;

        bool isStrongPseudoPrime = false;

        for (int r = 1; r < s; r++)
        {
            x = (__int128)x * x % n;

            if (x == n - 1)
            {
                isStrongPseudoPrime = true;
                break;
            }

            if (x == 1)
                return false;
        }

        if (!isStrongPseudoPrime)
            return false;
    }

    return true;
}

long long trialDivision(long long n)
{
    if (n < 2)
        return n;

    if (n % 2 == 0)
        return 2;

    for (long long i = 3; i <= n / i; i += 2)
    {
        if (n % i == 0)
            return i;
    }

    return n;
}

long long pollardRho(long long n)
{
    if (n < 2)
        return n;

    if (n % 2 == 0)
        return 2;

    if (n % 3 == 0)
        return 3;

    mt19937_64 gen(random_device{}());
    uniform_int_distribution<long long> dist(2, n - 2);

    while (true)
    {
        long long x = dist(gen);
        long long y = x;
        long long c = dist(gen);
        long long d = 1;

        auto f = [n, c](long long value)
        {
            return ((__int128)value * value + c) % n;
        };

        while (d == 1)
        {
            x = f(x);
            y = f(f(y));

            d = gcd(abs(x - y), n);
        }

        if (d != n)
            return d;
    }
}