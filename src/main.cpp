#include <iostream>
#include <chrono>
#include "pollard_rho.h"

using namespace std;
using namespace chrono;

int main()
{
    long long n;

    cout << "Enter a number: ";
    cin >> n;

    auto start = high_resolution_clock::now();

    long long divisor = pollardRho(n);

    auto end = high_resolution_clock::now();

    duration<double, milli> elapsed = end - start;

    cout << "Found divisor: " << divisor << endl;
    cout << "Execution time: " << elapsed.count() << " ms" << endl;

    return 0;
}