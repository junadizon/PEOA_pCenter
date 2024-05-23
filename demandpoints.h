#include <iostream>
#include <vector>
#include <string>
#include <fstream>
using namespace std;

void dmndpnts(vector<vector<double>> &demandPoints, int &dset) {
    string Xline, Yline;
    ifstream readXFile, readYFile;

    if (dset == 1) {
        readXFile.open("Davao_Full_X.csv");
        readYFile.open("Davao_Full_Y.csv");
    }
    else if (dset == 2) {
        readXFile.open("Digos_X_full.csv");
        readYFile.open("Digos_Y_full.csv");
    }
    else if (dset == 3) {
        readXFile.open("TagumFullX.csv");
        readYFile.open("TagumFullY.csv");
    }


    if (readXFile.is_open() && readYFile.is_open()) {
        while (getline (readXFile, Xline) && getline (readYFile, Yline)) {
            if (Xline.size() >= 3 && Xline[0] == '\xEF' && Xline[1] == '\xBB' && Xline[2] == '\xBF') {
                Xline = Xline.substr(3);
            }

            if (Yline.size() >= 3 && Yline[0] == '\xEF' && Yline[1] == '\xBB' && Yline[2] == '\xBF') {
                Yline = Yline.substr(3);
            }

            vector<double> coord;
            coord.push_back(stod(Xline));
            coord.push_back(stod(Yline));
            demandPoints.push_back(coord);
        }
        readXFile.close();
        readYFile.close();
    }
    else {
        cout << "ERROR! FILE NOT FOUND!" << endl;
    }
}