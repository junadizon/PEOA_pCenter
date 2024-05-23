#include <iostream>
#include <vector>
#include <string>
#include <bits/stdc++.h>
using namespace std;

vector<vector<double>> lhsdesign(int rows, int cols) {
    vector<vector<double>> design(rows, vector<double>(cols));
    vector<int> indices(rows);

    for (int i = 0; i < rows; ++i) {
        indices[i] = i;
    }

    for (int j = 0; j < cols; ++j) {
        random_shuffle(indices.begin(), indices.end());

        for (int i = 0; i < rows; ++i) {
            double lower = static_cast<double>(indices[i]) / rows;
            double upper = static_cast<double>(indices[i] + 1) / rows;

            design[i][j] = lower + (upper - lower) * static_cast<double>(rand()) / RAND_MAX;
        }
    }

    return design;
}

void getBounds(vector<vector<double>> &demandPoints, double &lbx, double &ubx, double &lby, double &uby) {
    lbx = ubx = demandPoints[0][0];
    lby = uby = demandPoints[0][1];

    for (const auto& point : demandPoints) {
        if (point[0] < lbx) lbx = point[0];
        if (point[0] > ubx) ubx = point[0];
        if (point[1] < lby) lby = point[1];
        if (point[1] > uby) uby = point[1];
    }
}

void initialization(vector<vector<double>> demandPoints, vector<vector<vector<double>>> &PhE, int dset, int dimension, int So, int numFacilities, double &upperBoundX, double &upperBoundY, double &lowerBoundX, double &lowerBoundY) {
    string city;
    double yDifference, xDifference;
    vector<vector<double>> xlhsMatrix, ylhsMatrix;

    if (dset == 1) {
        city = "Davao City";
    }
    else if (dset == 2) {
        city = "Digos City";
    }
    else if (dset == 3) {
        city = "Tagum City";
    }

    getBounds(demandPoints, lowerBoundX, upperBoundX, lowerBoundY, upperBoundY);

    xDifference = upperBoundX - lowerBoundX;
    yDifference = upperBoundY - lowerBoundY;

    xlhsMatrix = lhsdesign(numFacilities, So);
    ylhsMatrix = lhsdesign(numFacilities, So);

    for (int i = 0; i < So; i++) {
        for (int j = 0; j < numFacilities; j++) {
            PhE[i][j][0] = lowerBoundX + xDifference * xlhsMatrix[j][i];
            PhE[i][j][1] = lowerBoundY + yDifference * ylhsMatrix[j][i];
        }
    }
}