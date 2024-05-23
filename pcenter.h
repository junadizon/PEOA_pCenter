#include <iostream>
#include <vector>
#include <string>
#include <bits/stdc++.h>
#include <cmath>
using namespace std;

void minj(vector<vector<double>> &minDistances, vector<vector<vector<double>>> distances, int sizeDemandPoints, int So, int numFacilities) {
    for (int i = 0; i < So; ++i) {
        for (int k = 0; k < sizeDemandPoints; ++k) {
            double minDistance = numeric_limits<double>::max();
            int closestFacility = -1;

            for (int j = 0; j < numFacilities; ++j) {
                double distance = distances[i][j][k];

                if (distance < minDistance) {
                    minDistance = distance;
                    closestFacility = j;
                }
            }
            minDistances[i][k] = minDistance;
        }
    }
}

vector<double> findMaxValue(vector<vector<double>> minDistances, int So, int sizeDemandPoints) {
    vector<double> maxDistances(So, numeric_limits<double>::min()); // Initialize with a small value

    for (int i = 0; i < So; ++i) { // So
        for (int k = 0; k < sizeDemandPoints; ++k) {
            double distance = minDistances[i][k];

            if (distance > maxDistances[i]) {
                maxDistances[i] = distance;
            }
        }
    }

    return maxDistances;
    maxDistances.clear();
}

void MinMaxDistances(vector<double> maxDistances, int &bestSolutionIndex) {
    double minDistance = 1000000;

    for (int i = 0; i < maxDistances.size(); i++) {
        if (maxDistances[i] < minDistance) {
            minDistance = maxDistances[i];
            bestSolutionIndex = i;
        }
    }
}

void pcenter(vector<vector<double>> demandPoints, vector<vector<vector<double>>> PhE, int numFacilities, int dimension, int So, double &bestEagle, int &bestSolutionIndex, vector<double> &maxDistance) {
    int n = demandPoints.size();
    double x, y, squaredx, squaredy, distance;

    vector<vector<vector<double>>> distances(So, vector<vector<double>>(numFacilities, vector<double>(n)));
    vector<vector<double>> minDistances(So, vector<double>(n));

    bestSolutionIndex = -1;

    for (int i = 0; i < So; i++) {
        for (int j = 0; j < numFacilities; j++) {
            for (int k = 0; k < n; k++) {
                x = demandPoints[k][0] - PhE[i][j][0];
                y = demandPoints[k][1] - PhE[i][j][1];

                squaredx = pow(x, 2);
                squaredy = pow(y, 2);

                distance = sqrt(squaredx + squaredy);
                distances[i][j][k]= distance;
            }
        }
    }

    minj(minDistances, distances, n, So, numFacilities);
    maxDistance = findMaxValue(minDistances, So, n);
    MinMaxDistances(maxDistance, bestSolutionIndex);

    bestEagle = maxDistance[bestSolutionIndex];

    distances.clear();
    minDistances.clear();
}