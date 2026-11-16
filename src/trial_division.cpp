#include "trial_division.h"

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
