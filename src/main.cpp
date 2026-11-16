#include <iostream>
#include "primality.h"

using namespace std;

int main()
{
    long long n;

    cout << "Enter a number: ";
    cin >> n;

    long long divisor = trialDivision(n);

    if (divisor == n)
    {
        cout << n << " is prime." << endl;
    }
    else
    {
        cout << "Found divisor: " << divisor << endl;
    }

    return 0;
}