#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <chrono>
#include <string>
#include <bits/stdc++.h>

#include "parameters.h" // reads parameters file
#include "demandpoints.h" // reads csv files
#include "initialization.h" // initialize eagle
#include "pcenter.h" // pcenter math model to get the best eagle
#include "localphase.h" // search the best food for the best eagle
#include "movementoperator.h" // movement operator for the global phase
#include "mutation1operator.h" // mutation 1 operator for the global phase
#include "mutation2operator.h" // mutation 2 operator for the global phase