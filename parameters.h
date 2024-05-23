#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
using namespace std;

void parameters(int &p, int &Smin, double &A, int &Nmax, double &V, double &M, int &H, int &So, int &Sloc, double &rho, int &dset, int &dim, double &ftrue, int &runs) {
    string line, str;
    ifstream params;

    params.open("parameters.txt");

    if (params.is_open()) {
        while (getline (params, line)) {
            stringstream ss(line);
            ss >> str;

            if (str == "p") {
                ss >> str;
                p = stoi(str);
            }
            else if (str == "Smin") {
                ss >> str;
                Smin = stoi(str);
            }
            else if (str == "A") {
                ss >> str;
                A = stod(str);
            }
            else if (str == "Nmax") {
                ss >> str;
                Nmax = stoi(str);
            }
            else if (str == "V") {
                ss >> str;
                V = stod(str);
            }
            else if (str == "M") {
                ss >> str;
                M = stod(str);
            }
            else if (str == "H") {
                ss >> str;
                H = stoi(str);
            }
            else if (str == "So") {
                ss >> str;
                So = stoi(str);
            }
            else if (str == "Sloc") {
                ss >> str;
                Sloc = stoi(str);
            }
            else if (str == "rho") {
                ss >> str;
                rho = stod(str);
            }
            else if (str == "dset") {
                ss >> str;
                dset = stoi(str);
            }
            else if (str == "dim") {
                ss >> str;
                dim = stoi(str);
            }
            else if (str == "ftrue") {
                ss >> str;
                ftrue = stoi(str);
            }
            else if (str == "runs") {
                ss >> str;
                runs = stoi(str);
            }
            else {
                cout << "Warning. Undefined parameter read.";
            }
        }
        params.close();
    }
    else {
        cout << "ERROR! FILE NOT FOUND!" << endl;
    }
}