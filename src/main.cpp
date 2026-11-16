
#include <iostream>
#include "primality.h"

using namespace std;

int main()
{
    long long n;

    cout << "Enter a number: ";
    cin >> n;

    if (isPrimeMillerRabin(n, 20))
    {
        cout << n << " is probably prime." << endl;
    }
    else
    {
        cout << n << " is composite." << endl;
    }

    return 0;
}
