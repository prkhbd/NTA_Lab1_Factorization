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
    if (k <= 0)
        return false;

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