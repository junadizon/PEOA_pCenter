#include <iostream>
#include <vector>
#include <bits/stdc++.h>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <chrono>
using namespace std;

vector<vector<vector<double>>> movementOperator(int currentSolutionIndex, int bestSolutionIndex, vector<vector<vector<double>>> PhE, vector<vector<vector<double>>> ArchivedEagles, int PS1, int numFacilities, double scalingFactor, int generation) {
    vector<vector<vector<double>>> Xnew(1, vector<vector<double>>(numFacilities, vector<double>(2)));

    double NewX, d, NewY, dist, dx, x, dy, squaredx, squaredy, distance, expDistance, summation;
    double bestCurrentEaglesX, bestCurrentEaglesY, randomEaglesX, randomEaglesY, nearestEagleX, nearestEagleY, distanceX, distanceY;
    int Xnear, Xr1, Xarc, closestFacility;
    double minDistance;
    double minEuclideanDistance = numeric_limits<double>::max();
    vector<vector<double>> distances(PS1, vector<double>(numFacilities));

    vector<double> minDistances(PS1);

     // Compute minimum distances and nearest index
    for (int i = 0; i < PS1; ++i) {
        for (int j = 0; j < numFacilities; ++j) {
            dx = PhE[currentSolutionIndex][j][0] - PhE[i][j][0];
            dy = PhE[currentSolutionIndex][j][1] - PhE[i][j][1];
            double dist = sqrt(dx * dx + dy * dy);
            if (dist < minDistances[i]) {
                minDistances[i] = dist;
            }
        }
        if (minDistances[i] < minEuclideanDistance) {
            minEuclideanDistance = minDistances[i];
            Xnear = i;
        }
    }

    auto currentTime = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    srand(static_cast<unsigned int>(currentTime));

    do {
        Xr1 = rand() % PhE.size();
    }
    while (Xr1 == Xnear || Xr1 == bestSolutionIndex || Xr1 == currentSolutionIndex);

    if (generation == 1) {
        do {
            Xarc = rand() % ArchivedEagles.size();
        }
        while (Xarc == Xr1 || Xarc == bestSolutionIndex || Xarc == currentSolutionIndex || Xarc == Xnear);
    }
    else {
        Xarc = rand() % ArchivedEagles.size();
    }

    d = pow(minEuclideanDistance, 2);
    expDistance = exp(-d);

    for (int j = 0; j < numFacilities; j++) {
        bestCurrentEaglesX = PhE[bestSolutionIndex][j][0] - PhE[currentSolutionIndex][j][0];
        bestCurrentEaglesY = PhE[bestSolutionIndex][j][1] - PhE[currentSolutionIndex][j][1];

        randomEaglesX = PhE[Xr1][j][0] - ArchivedEagles[Xarc][j][0];
        randomEaglesY = PhE[Xr1][j][1] - ArchivedEagles[Xarc][j][1];

        nearestEagleX = PhE[Xnear][j][0] - PhE[currentSolutionIndex][j][0];
        nearestEagleY = PhE[Xnear][j][1] - PhE[currentSolutionIndex][j][1];

        distanceX = expDistance * nearestEagleX;
        distanceY = expDistance * nearestEagleY;

        NewX = PhE[currentSolutionIndex][j][0] + scalingFactor * (bestCurrentEaglesX + randomEaglesX + distanceX);
        NewY = PhE[currentSolutionIndex][j][1] + scalingFactor * (bestCurrentEaglesY + randomEaglesY + distanceY);

        Xnew[0][j][0] = NewX;
        Xnew[0][j][1] = NewY;

    }

    return Xnew;
}