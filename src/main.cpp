
#include <iostream>
#include <vector>
#include <map>
#include <chrono>
#include <iomanip>
#include <string>
#include <functional>
#include <ctime>
#include <sstream>

#include "primality.h"
#include "trial_division.h"
#include "pollard_rho.h"
#include "brillhart_morrison.h"

using namespace std;
using namespace chrono;

const long long DEFAULT_NUMBER = 691534156424661573LL;

struct BMResult
{
    long long divisor = 0;
    size_t factorBaseSize = 0;
    size_t relationsCount = 0;
    size_t smoothCount = 0;
    size_t dependenciesCount = 0;
};

struct FactorEvent
{
    long long divisor;
    string method;
    string timestamp;
};

string getTimestamp()
{
    auto now = system_clock::now();
    auto time = system_clock::to_time_t(now);

    auto millis = duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()
    ).count() % 1000;

    tm localTime = *localtime(&time);

    stringstream ss;

    ss << put_time(&localTime, "%Y-%m-%d %H:%M:%S");
    ss << "." << setfill('0') << setw(3) << millis;

    return ss.str();
}

long long trialDivisionUpTo47(long long n)
{
    if (n % 2 == 0)
        return 2;

    for (long long i = 3; i <= 47 && i <= n / i; i += 2)
    {
        if (n % i == 0)
            return i;
    }

    return 0;
}

BMResult runBrillhartMorrison(long long n)
{
    BMResult result;

    if (n <= 3)
        return result;

    vector<long long> factorBase = buildFactorBase(n);

    result.factorBaseSize = factorBase.size();

    if (factorBase.size() < 2)
        return result;

    int maxTerms = 5000;

    vector<long long> coefficients =
        calculateContinuedFraction(n, maxTerms);

    vector<cpp_int> numerators;
    vector<cpp_int> denominators;

    calculateConvergents(
        coefficients,
        numerators,
        denominators
    );

    vector<cpp_int> relations =
        calculateRelations(n, numerators, denominators);

    result.relationsCount = relations.size();

    vector<size_t> smoothIndices;
    vector<vector<int>> smoothExponents;

    size_t targetRelations = factorBase.size() + 20;

    for (size_t i = 0;
         i < relations.size() && smoothIndices.size() < targetRelations;
         i++)
    {
        vector<int> exponents;

        if (factorSmoothRelation(
                relations[i],
                factorBase,
                exponents))
        {
            smoothIndices.push_back(i);
            smoothExponents.push_back(exponents);
        }
    }

    result.smoothCount = smoothIndices.size();

    if (smoothIndices.size() < targetRelations)
        return result;

    vector<vector<int>> parityMatrix =
        buildParityMatrix(smoothExponents);

    if (parityMatrix.empty())
        return result;

    vector<vector<int>> dependencies =
        findDependencies(parityMatrix);

    result.dependenciesCount = dependencies.size();

    result.divisor = findFactorFromDependencies(
        n,
        factorBase,
        numerators,
        smoothIndices,
        smoothExponents,
        dependencies
    );

    if (result.divisor <= 1 ||
        result.divisor >= n ||
        n % result.divisor != 0)
    {
        result.divisor = 0;
    }

    return result;
}

void printBMResult(long long n, const BMResult& result)
{
    cout << "Number: " << n << endl;
    cout << "Factor base size: "
         << result.factorBaseSize << endl;

    cout << "Relations calculated: "
         << result.relationsCount << endl;

    cout << "B-smooth relations found: "
         << result.smoothCount << endl;

    cout << "Linear dependencies found: "
         << result.dependenciesCount << endl;

    if (result.divisor != 0)
    {
        cout << "Brillhart-Morrison divisor: "
             << result.divisor << endl;

        cout << "Second factor: "
             << n / result.divisor << endl;
    }
    else
    {
        cout << "No non-trivial divisor found." << endl;
    }
}

bool factorRecursive(
    long long n,
    map<long long, int>& factors,
    vector<FactorEvent>& events
)
{
    if (n == 1)
        return true;

    if (isPrimeMillerRabin(n, 20))
    {
        factors[n]++;
        return true;
    }

    long long divisor = trialDivisionUpTo47(n);

    if (divisor > 1 && divisor < n)
    {
        events.push_back({
            divisor,
            "Trial division",
            getTimestamp()
        });

        return factorRecursive(divisor, factors, events) &&
               factorRecursive(n / divisor, factors, events);
    }

    divisor = pollardRho(n);

    if (divisor > 1 &&
        divisor < n &&
        n % divisor == 0)
    {
        events.push_back({
            divisor,
            "Pollard Rho",
            getTimestamp()
        });

        return factorRecursive(divisor, factors, events) &&
               factorRecursive(n / divisor, factors, events);
    }

    BMResult bmResult = runBrillhartMorrison(n);

    divisor = bmResult.divisor;

    if (divisor > 1 &&
        divisor < n &&
        n % divisor == 0)
    {
        events.push_back({
            divisor,
            "Brillhart-Morrison",
            getTimestamp()
        });

        return factorRecursive(divisor, factors, events) &&
               factorRecursive(n / divisor, factors, events);
    }

    return false;
}

void printCanonicalFactorization(
    const map<long long, int>& factors
)
{
    cout << "Canonical factorization: ";

    bool first = true;

    for (const auto& factor : factors)
    {
        if (!first)
            cout << " * ";

        cout << factor.first;

        if (factor.second > 1)
            cout << "^" << factor.second;

        first = false;
    }

    cout << endl;
}

void runGeneralFactorization(long long n)
{
    if (n <= 1)
    {
        cout << "Number must be greater than 1." << endl;
        return;
    }

    auto startTime = high_resolution_clock::now();
    string startTimestamp = getTimestamp();

    cout << "Canonical factorization" << endl;
    cout << "Number: " << n << endl;
    cout << "Start time: " << startTimestamp << endl;

    map<long long, int> factors;
    vector<FactorEvent> events;

    bool success = factorRecursive(n, factors, events);

    auto endTime = high_resolution_clock::now();
    string endTimestamp = getTimestamp();

    cout << endl;
    cout << "Found divisors:" << endl;

    if (events.empty())
    {
        cout << "No non-trivial divisors were required." << endl;
    }

    for (const auto& event : events)
    {
        cout << "Divisor: " << event.divisor << endl;
        cout << "Method: " << event.method << endl;
        cout << "Timestamp: " << event.timestamp << endl;
        cout << endl;
    }

    if (success)
    {
        printCanonicalFactorization(factors);
    }
    else
    {
        cout << "I cannot find the canonical factorization of this number :(" << endl;
    }

    cout << "Start time: " << startTimestamp << endl;
    cout << "End time: " << endTimestamp << endl;

    duration<double, milli> elapsed = endTime - startTime;

    cout << fixed << setprecision(6);
    cout << "Total execution time: "
         << elapsed.count() << " ms" << endl;
}

void runPrimalityTest(long long n)
{
    auto startTime = high_resolution_clock::now();

    bool prime = isPrimeMillerRabin(n, 20);

    auto endTime = high_resolution_clock::now();

    cout << "Miller-Rabin primality test" << endl;
    cout << "Number: " << n << endl;
    cout << "Result: "
         << (prime ? "Prime" : "Composite") << endl;

    duration<double, milli> elapsed = endTime - startTime;

    cout << fixed << setprecision(6);
    cout << "Execution time: "
         << elapsed.count() << " ms" << endl;
}

void runTrialDivisionTest(long long n)
{
    auto startTime = high_resolution_clock::now();

    long long divisor = trialDivision(n);

    auto endTime = high_resolution_clock::now();

    cout << "Trial division" << endl;
    cout << "Number: " << n << endl;

    if (divisor > 1 && divisor < n)
    {
        cout << "Divisor: " << divisor << endl;
        cout << "Second factor: " << n / divisor << endl;
    }
    else
    {
        cout << "No non-trivial divisor found." << endl;
    }

    duration<double, milli> elapsed = endTime - startTime;

    cout << fixed << setprecision(6);
    cout << "Execution time: "
         << elapsed.count() << " ms" << endl;
}

void runPollardTest(long long n)
{
    auto startTime = high_resolution_clock::now();

    long long divisor = pollardRho(n);

    auto endTime = high_resolution_clock::now();

    cout << "Pollard Rho factorization" << endl;
    cout << "Number: " << n << endl;

    if (divisor > 1 &&
        divisor < n &&
        n % divisor == 0)
    {
        cout << "Divisor: " << divisor << endl;
        cout << "Second factor: " << n / divisor << endl;
    }
    else
    {
        cout << "No non-trivial divisor found." << endl;
    }

    duration<double, milli> elapsed = endTime - startTime;

    cout << fixed << setprecision(6);
    cout << "Execution time: "
         << elapsed.count() << " ms" << endl;
}

void runBMTest(long long n)
{
    auto startTime = high_resolution_clock::now();

    BMResult result = runBrillhartMorrison(n);

    auto endTime = high_resolution_clock::now();

    cout << "Brillhart-Morrison factorization" << endl;

    printBMResult(n, result);

    duration<double, milli> elapsed = endTime - startTime;

    cout << fixed << setprecision(6);
    cout << "Execution time: "
         << elapsed.count() << " ms" << endl;
}

void runBenchmark()
{
    vector<long long> numbers = {
        3009182572376191LL,
        1021514194991569LL,
        4000852962116741LL,
        15196946347083LL,
        499664789704823LL,
        269322119833303LL,
        679321846483919LL,
        96267366284849LL,
        61333127792637LL,
        2485021628404193LL
    };

    cout << "Factorization benchmark" << endl;

    for (long long n : numbers)
    {
        cout << endl;
        cout << "================================" << endl;
        cout << "Number: " << n << endl;

        auto rhoStart = high_resolution_clock::now();

        long long rhoDivisor = pollardRho(n);

        auto rhoEnd = high_resolution_clock::now();

        duration<double, milli> rhoTime = rhoEnd - rhoStart;

        auto bmStart = high_resolution_clock::now();

        BMResult bmResult = runBrillhartMorrison(n);

        auto bmEnd = high_resolution_clock::now();

        duration<double, milli> bmTime = bmEnd - bmStart;

        cout << "Pollard Rho divisor: ";

        if (rhoDivisor > 1 &&
            rhoDivisor < n &&
            n % rhoDivisor == 0)
        {
            cout << rhoDivisor << endl;
        }
        else
        {
            cout << "Not found" << endl;
        }

        cout << "Pollard Rho time: "
             << fixed << setprecision(6)
             << rhoTime.count() << " ms" << endl;

        cout << "Brillhart-Morrison divisor: ";

        if (bmResult.divisor != 0)
            cout << bmResult.divisor << endl;
        else
            cout << "Not found" << endl;

        cout << "Brillhart-Morrison time: "
             << bmTime.count() << " ms" << endl;
    }
}

void runAllTests()
{
    cout << "TEST 1: MILLER-RABIN" << endl;
    runPrimalityTest(97);
    runPrimalityTest(91);

    cout << endl;
    cout << "TEST 2: TRIAL DIVISION" << endl;
    runTrialDivisionTest(8051);

    cout << endl;
    cout << "TEST 3: POLLARD RHO" << endl;
    runPollardTest(8051);

    cout << endl;
    cout << "TEST 4: BRILLHART-MORRISON" << endl;
    runBMTest(8051);

    cout << endl;
    cout << "TEST 5: CANONICAL FACTORIZATION" << endl;
    runGeneralFactorization(DEFAULT_NUMBER);

    cout << endl;
    cout << "TEST 6: BENCHMARK" << endl;
    runBenchmark();
}

int main(int argc, char* argv[])
{
    string mode = "--factor";
    long long n = DEFAULT_NUMBER;

    if (argc > 1)
    {
        mode = argv[1];
    }

    if (argc > 2)
    {
        try
        {
            size_t pos = 0;
            n = stoll(argv[2], &pos);

            if (pos != string(argv[2]).size())
            {
                cout << "Invalid input number." << endl;
                return 1;
            }
        }
        catch (...)
        {
            cout << "Invalid input number." << endl;
            return 1;
        }
    }

    if (mode == "--factor")
    {
        runGeneralFactorization(n);
    }
    else if (mode == "--test")
    {
        runPrimalityTest(n);
    }
    else if (mode == "--trial")
    {
        runTrialDivisionTest(n);
    }
    else if (mode == "--rho")
    {
        runPollardTest(n);
    }
    else if (mode == "--bm")
    {
        runBMTest(n);
    }
    else if (mode == "--benchmark")
    {
        runBenchmark();
    }
    else if (mode == "--all")
    {
        runAllTests();
    }
    else
    {
        cout << "Available modes:" << endl;
        cout << "--factor [number]" << endl;
        cout << "--test [number]" << endl;
        cout << "--trial [number]" << endl;
        cout << "--rho [number]" << endl;
        cout << "--bm [number]" << endl;
        cout << "--benchmark" << endl;
        cout << "--all" << endl;

        return 1;
    }

    return 0;
}
