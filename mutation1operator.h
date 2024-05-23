#include <iostream>
#include <vector>
#include <bits/stdc++.h>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <chrono>
using namespace std;

const double M1PI = acos(-1);

double randn() {
    double u1, u2, z, mean, stddev, randomNumber;

    // Generate two random numbers between 0 and 1
    u1 = static_cast<double>(rand()) / RAND_MAX;
    u2 = static_cast<double>(rand()) / RAND_MAX;

    // Use the Box-Muller transform to generate a normally distributed random number
    z = sqrt(-2.0 * log(u1)) * cos(2.0 * M1PI * u2);

    mean = 0.0;
    stddev = 1.0;
    randomNumber = mean + stddev * z;

    return randomNumber;
}

double levyFlight() {
    double sigmaNumerator, sigmaDenominator, sigma, levyFlightNumerator, levyFlightDenominator, levy;
    double betaL = 1.5;
    double randomU = randn();
    double randomV = randn();
    double absoluteV = abs(randomV);

    sigmaNumerator = (tgamma(1 + betaL)) * (sin((M1PI * betaL) / 2));
    sigmaDenominator = (tgamma((1 + betaL)/2)) * (betaL * (pow(2,((betaL - 1) / 2))));
    sigma = pow((sigmaNumerator / sigmaDenominator),(1 / betaL));

    levyFlightNumerator = 0.01 * (randomU * sigma);
    levyFlightDenominator = pow(absoluteV,(1 / betaL));

    levy = levyFlightNumerator / levyFlightDenominator;

    return levy;
}

void newEagle(vector<vector<vector<double>>> PhE, int X, int Xr1, int Xr2, int numFacilities, vector<vector<vector<double>>> &Xnew, double scalingFactor) {
    double eaglesX, eaglesY, finalEagleX, finalEagleY, levy, randomNum, x, y;

    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < numFacilities; j++) {
            levy = levyFlight();
            double S = static_cast<double>(rand()) / RAND_MAX; // random number j between 0 and 1

            eaglesX = PhE[Xr1][j][0] + PhE[X][j][0];
            eaglesY = PhE[Xr1][j][1] + PhE[X][j][1];

            finalEagleX = eaglesX - PhE[Xr2][j][0];
            finalEagleY = eaglesY - PhE[Xr2][j][1];

            x = scalingFactor * (finalEagleX) + (S * levy);
            y = scalingFactor * (finalEagleY) + (S * levy);

            Xnew[i][j][0] = y;
            Xnew[i][j][1] = x;
        }
    }
}

vector<vector<vector<double>>> mutation1Operator(int currentSolutionIndex, int bestSolutionIndex, vector<vector<vector<double>>> PhE, int PS1, int numFacilities, double scalingFactor) {
    vector<vector<vector<double>>> Xnew(1, vector<vector<double>>(numFacilities, vector<double>(2)));

    int Xr1, Xr2;

    auto currentTime = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    srand(static_cast<unsigned int>(currentTime));

    // random eagle Xr1
    do {
        Xr1 = rand() % PS1;
    }
    while (Xr1 == bestSolutionIndex || Xr1 == currentSolutionIndex); // make sure that the random eagle is not the same as best eagle

    do {
        Xr2 = rand() % PS1;
    }
    while (Xr2 == bestSolutionIndex || Xr2 == Xr1 || Xr2 == currentSolutionIndex); // make sure that the random eagle is not the same as best eagle and the first random eagle

    newEagle(PhE, bestSolutionIndex, Xr1, Xr2, numFacilities, Xnew, scalingFactor);

    return Xnew;
}