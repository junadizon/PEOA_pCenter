#include <iostream>
#include <vector>
#include <bits/stdc++.h>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <chrono>
using namespace std;

void newEagle(vector<vector<vector<double>>> PhE, int Xhat, int bestSolutionIndex, double XmeanX, double XmeanY, int numFacilities, vector<vector<vector<double>>> &Xnew, double scalingFactor) {
    double eaglesX, eaglesY, finalEagleX, finalEagleY, levy, randomNum, x, y;

    for (int j = 0; j < numFacilities; j++) {
        double XCoordinate = PhE[Xhat][j][0] + PhE[bestSolutionIndex][j][0] - XmeanX;
        double YCoordinate = PhE[Xhat][j][1] + PhE[bestSolutionIndex][j][1] - XmeanY;

        Xnew[0][j][0] = scalingFactor * XCoordinate;
        Xnew[0][j][1] = scalingFactor * YCoordinate;
    }
}

vector<vector<vector<double>>> mutation2Operator(int currentSolutionIndex, int bestSolutionIndex, vector<vector<vector<double>>> PhE, int PS1, int numFacilities, double scalingFactor) {
    vector<vector<vector<double>>> Xnew(1, vector<vector<double>>(numFacilities, vector<double>(2)));

    int Xhat;
    double XmeanX, XmeanY, Xsum = 0.0, Ysum = 0.0, sumX = 0.0, sumY = 0.0;
    double SingleEagleMeanX, SingleEagleMeanY;

    auto currentTime = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    srand(static_cast<unsigned int>(currentTime));

    // random eagle inside the search space
    do {
        Xhat = rand() % PS1;
    }
    while (Xhat == bestSolutionIndex); // make sure that the random eagle is not the same as best eagle

    for (int i = 0; i < PS1; i++) {
        for (int j = 0; j < numFacilities; j++) {
            Xsum += PhE[i][j][0];
            Ysum += PhE[i][j][1];
        }
    }

    XmeanX = Xsum / (PS1 * numFacilities);
    XmeanY = Ysum / (PS1 * numFacilities);

    newEagle(PhE, Xhat, bestSolutionIndex, XmeanX, XmeanY, numFacilities, Xnew, scalingFactor);

    return Xnew;
}