
#include <iostream>
#include <vector>
#include <cstddef>

#include "brillhart_morrison.h"

using namespace std;

int main()
{
    long long n = 691534156424661573;

    int maxTerms = 5000;

    vector<long long> factorBase = buildFactorBase(n);

    vector<long long> coefficients =
        calculateContinuedFraction(n, maxTerms);

    vector<cpp_int> numerators;
    vector<cpp_int> denominators;

    calculateConvergents(coefficients, numerators, denominators);

    vector<cpp_int> relations =
        calculateRelations(n, numerators, denominators);

    vector<size_t> smoothIndices;
    vector<vector<int>> smoothExponents;

    size_t targetRelations = factorBase.size() + 20;

    for (size_t i = 0; i < relations.size(); i++)
    {
        vector<int> exponents;

        if (factorSmoothRelation(relations[i], factorBase, exponents))
        {
            smoothIndices.push_back(i);
            smoothExponents.push_back(exponents);

            if (smoothIndices.size() >= targetRelations)
                break;
        }
    }

    cout << "Number: " << n << endl;
    cout << "Factor base size: " << factorBase.size() << endl;
    cout << "Continued fraction terms generated: "
         << coefficients.size() << endl;
    cout << "Relations calculated: " << relations.size() << endl;
    cout << "Target B-smooth relations: " << targetRelations << endl;
    cout << "B-smooth relations found: " << smoothIndices.size() << endl;

    cout << endl;
    cout << "B-smooth relations:" << endl;

    for (size_t i = 0; i < smoothIndices.size(); i++)
    {
        size_t index = smoothIndices[i];

        cout << "R[" << index << "] = " << relations[index] << endl;
        cout << "Factorization: ";

        bool first = true;

        for (size_t j = 0; j < smoothExponents[i].size(); j++)
        {
            if (smoothExponents[i][j] != 0)
            {
                if (!first)
                    cout << " * ";

                cout << factorBase[j];

                if (smoothExponents[i][j] > 1)
                    cout << "^" << smoothExponents[i][j];

                first = false;
            }
        }

        cout << endl;
    }

    if (smoothIndices.size() < targetRelations)
    {
        cout << endl;
        cout << "The target number of B-smooth relations was not reached." << endl;
    }

    return 0;
}
