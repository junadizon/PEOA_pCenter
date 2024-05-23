#include <iostream>
#include <vector>
#include <bits/stdc++.h>
#include <algorithm>
#include <cmath>
using namespace std;

void localphase_pcenter(vector<vector<double>> demandPoints, vector<vector<vector<double>>> solution, int numFacilities, int dimension, int populationSize, double &bestSolution, int &bestSolutionIndex, vector<double> &maxDistance) {
    double x, y, squaredx, squaredy, distance, minDistance, largeValue = 1000000;
    int closestFacility;
    vector<double> maxDistances(populationSize, numeric_limits<double>::min());

    for (int i = 0; i < populationSize; i++) {
        for (int k = 0; k < demandPoints.size(); k++) {
            minDistance = numeric_limits<double>::max();
            closestFacility = -1;

            for (int j = 0; j < numFacilities; j++) {
                x = demandPoints[k][0] - solution[i][j][0];
                y = demandPoints[k][1] - solution[i][j][1];

                squaredx = pow(x, 2);
                squaredy = pow(y, 2);

                distance = sqrt(squaredx + squaredy);

                if (distance < minDistance) {
                    minDistance = distance;
                    closestFacility = j;
                }
            }

            if (minDistance > maxDistances[i]) {
                maxDistances[i] = distance;
            }
        }

        if (maxDistances[i] < largeValue) {
            largeValue = maxDistances[i];
            bestSolutionIndex = i;
        }
    }

    bestSolution = maxDistances[bestSolutionIndex];

    maxDistances.clear();
}

vector<vector<double>> XYminmax(vector<vector<double>> bestEagle, double radius, int numFacilities, int coordinate, string bound, int dim) {
    vector<vector<double>> boundXY(numFacilities, vector<double>(dim));

    if (bound == "lb") {
        for (int i = 0; i < numFacilities; i++) {
            boundXY[i][coordinate] = bestEagle[i][coordinate] - radius;
        }
    }
    else if (bound == "ub") {
        for (int i = 0; i < numFacilities; i++) {
            boundXY[i][coordinate] = bestEagle[i][coordinate] + radius;
        }
    }

    return boundXY;
}

void functionValueIPM(int localFoodSize, vector<vector<vector<double>>> randomValues, int numFacilities, int randValBestSolutionIndex, vector<vector<double>> lbX, vector<vector<double>> ubX, vector<vector<double>> lbY, vector<vector<double>> ubY, vector<vector<vector<double>>> &IPMSolutions) {
    double X, Y, Xlog1, Xlog2, Ylog1, Ylog2, Xsum, Ysum, t, quotient, functionValX, functionValY;

    for (int i = 0; i < localFoodSize; i++) {
        for (int j = 0; j < numFacilities; j++) {
            X = randomValues[randValBestSolutionIndex][j][0];
            Y = randomValues[randValBestSolutionIndex][j][1];

            Xlog1 = log(ubX[j][0] - randomValues[i][j][0]);
            Xlog2 = log(randomValues[i][j][0] + lbX[j][0]);

            Ylog1 = log(ubY[j][1] - randomValues[i][j][1]);
            Ylog2 = log(randomValues[i][j][1] + lbY[j][1]);

            Xsum = Xlog1 + Xlog2;
            Ysum = Ylog1 + Ylog2;

            t = i + 1;
            quotient = 1/t;

            functionValX = X - (quotient * Xsum);
            functionValY = Y - (quotient * Ysum);

            IPMSolutions[i][j][0] = functionValX;
            IPMSolutions[i][j][1] = functionValY;
        }
    }
}

void checkFuncVal(int localFoodSize, vector<vector<vector<double>>> IPMSolutions, vector<vector<double>> lbX, vector<vector<double>> ubX, vector<vector<double>> lbY, vector<vector<double>> ubY, int numFacilities, vector<vector<vector<double>>> &localFoods) {
    int Ybound, Xbound;

    for (int i = 0; i < localFoodSize; i++) {
        for (int j = 0; j < numFacilities; j++) {
            Xbound = 0, Ybound = 0;

            localFoods[i][j][0] = IPMSolutions[i][j][0];
            localFoods[i][j][1] = IPMSolutions[i][j][1];

            if (IPMSolutions[i][j][0] < lbX[j][0]) {
                localFoods[i][j][0] = lbX[j][0];

                Xbound = 1;
            }
            else if (IPMSolutions[i][j][0] > ubX[j][0]) {
                localFoods[i][j][0] = ubX[j][0];

                Xbound = 1;
            }
            // check the Y
            if (IPMSolutions[i][j][1] < lbY[j][1]) {
                localFoods[i][j][1] = lbY[j][1];

                Ybound = 1;
            }
            else if (IPMSolutions[i][j][1] > ubY[j][1]) {
                localFoods[i][j][1] = ubY[j][1];

                Ybound = 1;
            }
        }
    }
}

void localphase(double &bestFoodFitness, int &bestFoodSolutionIndex, vector<vector<vector<double>>> &localFoodSolutions, int numFacilities, int localFoodSize, vector<vector<double>> demandPoints, int dimension, double rho, vector<vector<double>> bestEagle, double lowerBoundX, double lowerBoundY, double upperBoundX, double upperBoundY, int dim) {
    double differenceX, differenceY, leastDifference;
    double radiusX, radiusY;
    double maxDistance;
    int bestFacilityIndex;
    vector<double> randomValuesXX, randomValuesYY, localFoodFitVals;
    vector<vector<vector<double>>> randomValues(localFoodSize, vector<vector<double>>(numFacilities, vector<double>(dim)));
    vector<vector<vector<double>>> IPMSolutions(localFoodSize, vector<vector<double>>(numFacilities, vector<double>(dim)));

    double bestRandValFitness = 0.0;
    int randValBestSolutionIndex;
    vector<double> randValFitnessValues;
    vector<double> onesVector(dimension);
    vector<vector<double>> lbXTerritoryXVector(numFacilities, vector<double>(dim));
    vector<vector<double>> ubXTerritoryXVector(numFacilities, vector<double>(dim));
    vector<vector<double>> lbYTerritoryYVector(numFacilities, vector<double>(dim));
    vector<vector<double>> ubYTerritoryYVector(numFacilities, vector<double>(dim));

    differenceX = upperBoundX - lowerBoundX;
    differenceY = upperBoundY - lowerBoundY;
    radiusX = max(rho * differenceX, 1.0);
    radiusY = max(rho * differenceY, 1.0);

    lbXTerritoryXVector = XYminmax(bestEagle, radiusX, numFacilities, 0, "lb", dim);
    ubXTerritoryXVector = XYminmax(bestEagle, radiusX, numFacilities, 0, "ub", dim);
    lbYTerritoryYVector = XYminmax(bestEagle, radiusY, numFacilities, 1, "lb", dim);
    ubYTerritoryYVector = XYminmax(bestEagle, radiusY, numFacilities, 1, "ub", dim);

    for (int i = 0; i < localFoodSize; i++) {
        for (int j = 0; j < numFacilities; j++) {
            randomValues[i][j][0] = lbXTerritoryXVector[j][0] + (static_cast<double>(std::rand()) / RAND_MAX) * (ubXTerritoryXVector[j][0] - lbXTerritoryXVector[j][0]);
            randomValues[i][j][1] = lbYTerritoryYVector[j][1] + (static_cast<double>(std::rand()) / RAND_MAX) * (ubYTerritoryYVector[j][1] - lbYTerritoryYVector[j][1]);
        }
    }

    localphase_pcenter(demandPoints, randomValues, numFacilities, dim, localFoodSize, bestRandValFitness, randValBestSolutionIndex, randValFitnessValues);

    functionValueIPM(localFoodSize, randomValues, numFacilities, randValBestSolutionIndex, lbXTerritoryXVector, ubXTerritoryXVector, lbYTerritoryYVector, ubYTerritoryYVector, IPMSolutions);
    checkFuncVal(localFoodSize, IPMSolutions, lbXTerritoryXVector, ubXTerritoryXVector, lbYTerritoryYVector, ubYTerritoryYVector, numFacilities, localFoodSolutions);
    localphase_pcenter(demandPoints, localFoodSolutions, numFacilities, dim, localFoodSize, bestFoodFitness, bestFoodSolutionIndex, localFoodFitVals);

    randomValuesXX.clear();
    randomValuesYY.clear();
    localFoodFitVals.clear();
    randomValues.clear();
    IPMSolutions.clear();
    randValFitnessValues.clear();
    onesVector.clear();
    lbXTerritoryXVector.clear();
    lbYTerritoryYVector.clear();
    ubXTerritoryXVector.clear();
    lbXTerritoryXVector.clear();
    bestEagle.clear();
}