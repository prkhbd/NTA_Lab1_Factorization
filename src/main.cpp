#include <iostream>
#include <vector>

#include "brillhart_morrison.h"

using namespace std;

int main()
{
    long long n = 691534156424661573;

    long double limit = calculateFactorBaseLimit(n);
    vector<long long> factorBase = buildFactorBase(n);

    cout << "Number: " << n << endl;
    cout << "Factor base limit: " << limit << endl;
    cout << "Factor base:" << endl;

    for (long long p : factorBase)
    {
        cout << p << " ";
    }

    cout << endl;
    cout << "Factor base size: " << factorBase.size() << endl;

    return 0;
}