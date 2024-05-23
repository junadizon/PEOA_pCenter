#include "headerfile.h"
using namespace std;
using namespace std::chrono;

const double PI = acos(-1);
void printEagle(vector<vector<vector<double>>> PhE, int solution, int numFacilities);
int operatorAssignment(double P1, double P2, double P3);
double ImprovementRate(vector<int> OperatorIndices, vector<double> OldFitVal, vector<double> NewFitVal);
void UpdateProbability(double R1, double R2, double R3, double &P1, double &P2, double &P3);
void removeRandomEagleArchive(vector<vector<vector<double>>> &ArchivedEagles);

/*
	Main function of PEOA
*/
vector<double> initBestFoodFit, initBestEagleFit;
vector<vector<double>> bestFoodFitTime, bestEagleFitTime, bestEagleFit, bestFoodFit;
vector<vector<double>> initBestEagleSolution, initBestFoodSolution, finalBestEagleSolution, finalBestFoodSolution, CurrentFood;
vector<vector<vector<double>>> initBestEagleSolutions, finalBestEagleSolutions, initBestFoodSolutions, finalBestFoodSolutions;
vector<vector<int>> countMovement, countMut1, countMut2;
int numFacilities;

void resizeEagle(int size, int numF, int dim, vector<vector<vector<double>>>& Eagle) {
    if (size < 0 || numF < 0 || dim < 0) {
        cerr << "Error: Negative size dimensions provided to resizeEagle." << endl;
        return; // Or handle more gracefully depending on your application needs
    }

    try {
        Eagle.resize(size, vector<vector<double>>(numF, vector<double>(dim)));
    } catch (const std::bad_alloc& e) {
        cerr << "Error: Memory allocation failed in resizeEagle: " << e.what() << endl;
    }
}

void mergeEagles(vector<vector<vector<double>>> Eagle, vector<vector<vector<double>>> ArchivedEagles) {
	ArchivedEagles.insert(ArchivedEagles.end(), Eagle.begin(), Eagle.end());
}

void removeMaxFit(vector<double>& fitnessValues, vector<vector<vector<double>>>& Eagle) {
	auto maxIndex = max_element(fitnessValues.begin(), fitnessValues.end());
	auto index = distance(fitnessValues.begin(), maxIndex);

	fitnessValues.erase(fitnessValues.begin() + index);
	Eagle.erase(Eagle.begin() + index);
}

int main() {
	srand(static_cast<unsigned>(time(0)));
	steady_clock::time_point start, stop1, stop2;
	duration<double> time1, time2;
	double timestamp1, timestamp2;

	vector<double> bestEagleFitRun,bestFoodFitRun, timeBE, timeBF;
	vector<int> counterM, counterM1, counterM2;

	double NewFitFood;
	int NewBestFoodSolIndex, counterMovement, counterMut1, counterMut2;

	// parameters
	int runs, currentRun = 0;
	int Dimension, minPopSize, maxFuncEval, memorySize, initPopSize, lclFoodSize, dset, dim;
	int bestSolutionIndex, bestFacilityIndex, best;
	double archiveRate, variance, initialMemory, rho, bestEaglefitness, ftrue;
	int bestSingleEagleSolIndex;
	int currentEvals; // number of function evaluation
	double upperBoundX, lowerBoundX, upperBoundY, lowerBoundY; // upper and lower bounds
	double bestEagleX, bestEagleY; // best eagle
	double bestFoodX, bestFoodY; // best eagle
	double memoryValue = 0.2;
	int generation = 0; // generation number
	int currentPopSize, histPos = 0;
	double S;	// entire eagle population size

	double P1 = 1.0 / 3.0, P2 = 1.0 / 3.0, P3 = 1.0 / 3.0; // Probability of Operators
	double R1 = 0.0, R2 = 0.0, R3 = 0.0; // Improve rate

	vector<double> goodF; // if there's improvement of the fitness values (in global phase)
	vector<double> fitnessValues, foodFitValues; // fitness val
	vector<vector<double>> demandPoints; // demand points

	vector<double> NewFitnessValues, NewFoodFitValues;
	double NewBestEaglefitness;
	int NewBestSolutionIndex;

	// read parameters file
	parameters(numFacilities, minPopSize, archiveRate, maxFuncEval, variance, initialMemory, memorySize, initPopSize, lclFoodSize, rho, dset, dim, ftrue, runs);

	Dimension = dim*numFacilities; // 2 coordinates * num of facilities

	double bestFoodFitness;
	int bestFoodSolutionIndex;

	// Archive F is the initial scaling factor
	// The initial setting of scaling factor
	vector<double> archiveF(memorySize);

	// External Archive of Eagles
	double ArchiveSize;
	vector<vector<vector<double>>> ArchivedEagles, CurrentEagle, Eagle, localFoodSolutions, newLocalFoodSolutions;

	dmndpnts(demandPoints, dset);

	for (currentRun = 1; currentRun <= runs; currentRun++) {
		cout << endl << endl << "Run: " << currentRun << endl;
		cout << "Current Function Evaluation: " << endl;

		currentPopSize = initPopSize;

		resizeEagle(initPopSize, numFacilities, dim, Eagle);
		resizeEagle(lclFoodSize, numFacilities, dim, localFoodSolutions);
		resizeEagle(lclFoodSize, numFacilities, dim, newLocalFoodSolutions);

		initBestFoodSolution.resize(numFacilities);
		initBestEagleSolution.resize(numFacilities);
		finalBestEagleSolution.resize(numFacilities);
		finalBestFoodSolution.resize(numFacilities);
		CurrentFood.resize(numFacilities);

		counterMovement = 0, counterMut1 = 0, counterMut2 = 0, generation = 0;

    	start = steady_clock::now();
		auto startIteration = chrono::steady_clock::now();
		currentEvals = 0;
		for (int i = 0; i < memorySize; i++) {
			archiveF[i] = initialMemory;
		}
		P1 = 1.0 / 3.0, P2 = 1.0 / 3.0, P3 = 1.0 / 3.0;
		R1 = 0.0, R2 = 0.0, R3 = 0.0;

		// -------------------------Initialization---------------------------
		initialization(demandPoints, Eagle, dset, dim, initPopSize, numFacilities, upperBoundX, upperBoundY, lowerBoundX, lowerBoundY);

		// ------------------------Objective Function------------------------
		pcenter(demandPoints, Eagle, numFacilities, dim, initPopSize, bestEaglefitness, bestSolutionIndex, fitnessValues);
		for (int j = 0; j < numFacilities; j++) {
			initBestEagleSolution[j].resize(2);
			initBestEagleSolution[j][0] = Eagle[bestSolutionIndex][j].at(0);
			initBestEagleSolution[j][1] = Eagle[bestSolutionIndex][j].at(1);
		}
		initBestEagleSolutions.push_back(initBestEagleSolution);
		currentEvals = currentEvals + initPopSize;

		// ------------------------Local Phase------------------------
		localphase(bestFoodFitness, bestFoodSolutionIndex, localFoodSolutions, numFacilities, lclFoodSize, demandPoints, Dimension, rho, initBestEagleSolution, lowerBoundX, lowerBoundY, upperBoundX, upperBoundY, dim);
		pcenter(demandPoints, localFoodSolutions, numFacilities, dim, lclFoodSize, bestFoodFitness, bestFoodSolutionIndex, foodFitValues);
		for (int j = 0; j < numFacilities; j++) {
			initBestFoodSolution[j].resize(2);
			initBestFoodSolution[j][0] = localFoodSolutions[bestFoodSolutionIndex][j].at(0);
			initBestFoodSolution[j][1] = localFoodSolutions[bestFoodSolutionIndex][j].at(1);
		}
		initBestFoodSolutions.push_back(initBestFoodSolution);
		bestFoodFitRun.push_back(bestFoodFitness);
		currentEvals = currentEvals + lclFoodSize;

		CurrentEagle = Eagle;
		CurrentFood = initBestFoodSolution;
		Eagle.clear();
		localFoodSolutions.clear();

		// ------------------------Start of the Loop------------------------
		while (currentEvals <= maxFuncEval)
		{
			generation++;
			cout << currentEvals << " ";

			// ------------------------Update Population Size------------------------
			int updatedPopSize = floor(initPopSize + (minPopSize - initPopSize) * static_cast<double>(currentEvals) / maxFuncEval);
			if (currentPopSize > updatedPopSize)
			{
				int reductionIndNum = currentPopSize - updatedPopSize;

				if ((currentPopSize - reductionIndNum) < minPopSize)
				{
					reductionIndNum = currentPopSize - minPopSize;
				}

				for (int r = 1; r < reductionIndNum; r++)
				{
					removeMaxFit(fitnessValues, CurrentEagle);
					currentPopSize = currentPopSize - 1;
					resizeEagle(currentPopSize, numFacilities, dim, CurrentEagle);
					fitnessValues.resize(currentPopSize);
				}

				ArchiveSize = round(archiveRate * currentPopSize);
				if (ArchivedEagles.size() > ArchiveSize)
				{
					random_device rd;
					mt19937 gen(rd());
					uniform_int_distribution<int> dist(0, ArchivedEagles.size() - 1);

					int randomIndex = dist(gen);
					ArchivedEagles.erase(ArchivedEagles.begin() + randomIndex);
				}
			}

			if (generation == 1)
			{
				for (int i = 0; i < currentPopSize; i++)
				{
					if (fitnessValues[i] == bestEaglefitness)
					{
						bestSolutionIndex = i;
					}
				}
			}
			else
			{
				for (int i = 0; i < currentPopSize; i++)
				{
					if (fitnessValues[i] == bestEaglefitness)
					{
						bestSolutionIndex = i;
					}
				}
			}

			vector<int> memoryRandomIndex(currentPopSize);
			vector<double> mu_sf(currentPopSize, initialMemory);

			for (int i = 0; i < currentPopSize; i++)
			{
				double randomValueF = static_cast<double>(rand()) / RAND_MAX;
				memoryRandomIndex[i] = ceil(memorySize * randomValueF);
			}

			vector<double> randomValuesF;
			randomValuesF.reserve(currentPopSize);

			vector<double> scalingFactor(currentPopSize);
			for (int i = 0; i < currentPopSize; i++)
			{
				double randomValue = static_cast<double>(rand()) / RAND_MAX;
				randomValuesF.push_back(randomValue);

				int index = memoryRandomIndex[i];
				double tempF;

				mu_sf[i] = archiveF[index];

				do
				{
					tempF = mu_sf[i] + variance * tan(PI * ((static_cast<double>(rand()) / RAND_MAX) - 0.5));

					if (tempF >= 1)
					{
						scalingFactor[i] = 1;
					}
					else if (tempF <= 0)
					{
						randomValuesF[i] = static_cast<double>(rand()) / RAND_MAX;
					}
					else
					{
						scalingFactor[i] = tempF;
					}
				} while (tempF <= 0);
			}

			// ------------------------Archive Eagles------------------------
			if (generation == 1)
			{
				ArchivedEagles = CurrentEagle;
			}
			else
			{
				mergeEagles(CurrentEagle, ArchivedEagles);
			}

			/*
			--------------------------Global Phase-----------------------
			(1) Movement Operator: pair-hunt strategy, more effective compared to hunting alone
			(2) Mutation I Operator: differing flight patterns of ph eagle, varying between straight glides and large circle
			(3) Mutation II Operator: perching behavior of the ph eagle -- watches its surroundings and looks out for prey at perch
			*/

			vector<vector<vector<double>>> movementOperatorXNew(1, vector<vector<double>>(numFacilities, vector<double>(dim)));
			vector<vector<vector<double>>> mutation1OperatorXNew(1, vector<vector<double>>(numFacilities, vector<double>(dim)));
			vector<vector<vector<double>>> mutation2OperatorXNew(1, vector<vector<double>>(numFacilities, vector<double>(dim)));
			vector<vector<vector<double>>> NewEagle(currentPopSize, vector<vector<double>>(numFacilities, vector<double>(dim)));

			vector<int> movementIndices, mutation1Indices, mutation2Indices;
			vector<int> assignOperator(currentPopSize);
			movementIndices.reserve(currentPopSize);
			mutation1Indices.reserve(currentPopSize);
			mutation2Indices.reserve(currentPopSize);

			for (int i = 0; i < currentPopSize; i++)
			{
				assignOperator[i] = operatorAssignment(P1, P2, P3);
				double currentScalingFactor = scalingFactor[i];

				if (assignOperator[i] == 1)
				{
					movementIndices.push_back(i);
					movementOperatorXNew = movementOperator(i, bestSolutionIndex, CurrentEagle, ArchivedEagles, currentPopSize, numFacilities, currentScalingFactor, generation);

					for (int j = 0; j < numFacilities; j++)
					{
						NewEagle[i][j][0] = movementOperatorXNew[0][j][0];
						NewEagle[i][j][1] = movementOperatorXNew[0][j][1];
					}
				}
				else if (assignOperator[i] == 2)
				{
					mutation1Indices.push_back(i);
					mutation1OperatorXNew = mutation1Operator(i, bestSolutionIndex, CurrentEagle, currentPopSize, numFacilities, currentScalingFactor);

					for (int j = 0; j < numFacilities; j++)
					{
						NewEagle[i][j][0] = mutation1OperatorXNew[0][j][0];
						NewEagle[i][j][1] = mutation1OperatorXNew[0][j][1];
					}
				}
				else if (assignOperator[i] == 3)
				{
					mutation2Indices.push_back(i);
					mutation2OperatorXNew = mutation2Operator(i, bestSolutionIndex, CurrentEagle, currentPopSize, numFacilities, currentScalingFactor);

					for (int j = 0; j < numFacilities; j++)
					{
						NewEagle[i][j][0] = mutation2OperatorXNew[0][j][0];
						NewEagle[i][j][1] = mutation2OperatorXNew[0][j][1];
					}
				}
			}

			// ------------------------p-Center for the New Eagles------------------------
			pcenter(demandPoints, NewEagle, numFacilities, dim, currentPopSize, NewBestEaglefitness, NewBestSolutionIndex, NewFitnessValues);
			// Check if there's improvement
			for (int i = 0; i < currentPopSize; i++)
			{
				if (NewFitnessValues[i] < fitnessValues[i])
				{
					fitnessValues[i] = NewFitnessValues[i];

					for (int j = 0; j < numFacilities; j++)
					{
						CurrentEagle[i][j][0] = NewEagle[i][j][0];
						CurrentEagle[i][j][1] = NewEagle[i][j][1];
					}

					if (assignOperator[i] == 1)
					{
						counterMovement++;
					}
					else if (assignOperator[i] == 2)
					{
						counterMut1++;
					}
					else if (assignOperator[i] == 3)
					{
						counterMut2++;
					}
				}
			}
			counterM.push_back(counterMovement);
			counterM1.push_back(counterMut1);
			counterM2.push_back(counterMut2);
			assignOperator.clear();

			currentEvals = currentEvals + currentPopSize;

			if (NewBestEaglefitness < bestEaglefitness)
			{
				bestEaglefitness = NewBestEaglefitness;
				bestSolutionIndex = NewBestSolutionIndex;

				for (int j = 0; j < numFacilities; j++)
				{
					finalBestEagleSolution[j].resize(2);
					finalBestEagleSolution[j][0] = NewEagle[NewBestSolutionIndex][j].at(0);
					finalBestEagleSolution[j][1] = NewEagle[NewBestSolutionIndex][j].at(1);
				}
			}
			else
			{
				for (int j = 0; j < numFacilities; j++)
				{
					finalBestEagleSolution[j].resize(2);
					finalBestEagleSolution[j][0] = CurrentEagle[bestSolutionIndex][j].at(0);
					finalBestEagleSolution[j][1] = CurrentEagle[bestSolutionIndex][j].at(1);
				}
			}
			/* ---------------------- STORE BEST EAGLE TIME ---------------------- */

			bestEagleFitRun.push_back(bestEaglefitness);

			stop1 = steady_clock::now();
			time1 = duration_cast<duration<double>>(stop1 - start);
			timestamp1 = time1.count();
			timeBE.push_back(timestamp1);

			// ------------------------Local Phase for the New Eagles------------------------
			localphase(NewFitFood, NewBestFoodSolIndex, newLocalFoodSolutions, numFacilities, lclFoodSize, demandPoints, Dimension, rho, finalBestEagleSolution, lowerBoundX, lowerBoundY, upperBoundX, upperBoundY, dim);
			pcenter(demandPoints, newLocalFoodSolutions, numFacilities, dim, lclFoodSize, NewFitFood, NewBestFoodSolIndex, NewFoodFitValues);

			// Check if there are improvements
			if (NewFitFood < bestFoodFitness)
			{
				bestFoodFitness = NewFitFood;
				bestFoodFitRun.push_back(bestFoodFitness);

				for (int j = 0; j < numFacilities; j++)
				{
					CurrentFood[j].at(0) = newLocalFoodSolutions[NewBestFoodSolIndex][j].at(0);
					CurrentFood[j].at(1) = newLocalFoodSolutions[NewBestFoodSolIndex][j].at(1);

					finalBestFoodSolution[j].resize(2);
					finalBestFoodSolution[j][0] = newLocalFoodSolutions[NewBestFoodSolIndex][j].at(0);
					finalBestFoodSolution[j][1] = newLocalFoodSolutions[NewBestFoodSolIndex][j].at(1);
				}
			}
			else {
				bestFoodFitRun.push_back(bestFoodFitness);
				for (int j = 0; j < numFacilities; j++)
				{
					finalBestFoodSolution[j].resize(2);
					finalBestFoodSolution[j][0] = CurrentFood[j].at(0);
					finalBestFoodSolution[j][1] = CurrentFood[j].at(1);
				}
			}
			currentEvals = currentEvals + lclFoodSize;

			/* ---------------------- STORE BEST FOOD TIME ---------------------- */

			stop2 = steady_clock::now();
			time2 = duration_cast<duration<double>>(stop2 - start);
			timestamp2 = time2.count();
			timeBF.push_back(timestamp2);

			/* ---------------------- UPDATE PROBABILITIES OF OPERATORS BASED ON IMPROVEMENT RATE ---------------------- */

			R1 = ImprovementRate(movementIndices, fitnessValues, NewFitnessValues);
			R2 = ImprovementRate(mutation1Indices, fitnessValues, NewFitnessValues);
			R3 = ImprovementRate(mutation2Indices, fitnessValues, NewFitnessValues);

			if (R1 == 0.0 || R2 == 0.0 || R3 == 0.0)
			{
				P1 = 1.0 / 3.0;
				P2 = 1.0 / 3.0;
				P3 = 1.0 / 3.0;
			}
			else
			{
				UpdateProbability(R1, R2, R3, P1, P2, P3);
			}

			/* ---------------------- UPDATE MEMORY ELEMENTS ---------------------- */

			double sumDiff = 0;
			vector<double> IndexImproved;
			vector<double> StoreDiff;
			for (int i = 0; i < currentPopSize; i++)
			{
				if (NewFitnessValues[i] < fitnessValues[i])
				{
					goodF.push_back(scalingFactor[i]);
					sumDiff += abs(fitnessValues[i] - NewFitnessValues[i]);
					StoreDiff.push_back(abs(fitnessValues[i] - NewFitnessValues[i]));

					IndexImproved.push_back(i);
				}
			}

			double sumNumerator, sumDenominator;
			if (IndexImproved.size() > 0)
			{
				for (int i = 0; i < IndexImproved.size(); i++)
				{
					double weightsDE = StoreDiff[i] / sumDiff;
					double Numerator = weightsDE * (pow(goodF[i], 2));
					double Denominator = weightsDE * goodF[i];

					sumNumerator += Numerator;
					sumDenominator += Denominator;
				}
				archiveF[histPos] = sumNumerator / sumDenominator;

				histPos++;

				if (histPos > memorySize)
				{
					histPos = 0;
				}
			}
			else
			{
				archiveF[histPos] = 0.5;
			}
		}

		for (int i = 0; i < timeBE.size(); i++) {
			if (i > 0 && bestEagleFitRun[i] == bestEagleFitRun[i - 1]) {
				timeBE[i] = timeBE[i - 1];
			}
		}

		for (int i = 0; i < bestFoodFitRun.size(); i++) {
			if (i > 0 && bestFoodFitRun[i] == bestFoodFitRun[i - 1]) {
				timeBF[i] = timeBF[i - 1];
			}
		}

		countMovement.push_back(counterM);
		countMut1.push_back(counterM1);
		countMut2.push_back(counterM2);
		bestEagleFitTime.push_back(timeBE);
		bestFoodFitTime.push_back(timeBF);
		bestEagleFit.push_back(bestEagleFitRun);
		bestFoodFit.push_back(bestFoodFitRun);
		finalBestEagleSolutions.push_back(finalBestEagleSolution);
		finalBestFoodSolutions.push_back(finalBestFoodSolution);

		timeBE.clear();
		timeBF.clear();
		bestEagleFitRun.clear();
		bestFoodFitRun.clear();
		counterM.clear();
		counterM1.clear();
		counterM2.clear();
		initBestEagleSolution.clear();
		initBestFoodSolution.clear();
		finalBestEagleSolution.clear();
		finalBestFoodSolution.clear();
		CurrentEagle.clear();
		CurrentFood.clear();
		ArchivedEagles.clear();
		fitnessValues.clear();
		foodFitValues.clear();
	}

	ofstream outfileInitBestEagleX("PEOA_InitBestEagleX.csv");
	ofstream outfileInitBestEagleY("PEOA_InitBestEagleY.csv");
	ofstream outfileFinalBestEagleX("PEOA_FinalBestEagleX.csv");
	ofstream outfileFinalBestEagleY("PEOA_FinalBestEagleY.csv");

	ofstream outfileInitBestFoodX("PEOA_InitBestFoodX.csv");
	ofstream outfileInitBestFoodY("PEOA_InitBestFoodY.csv");
	ofstream outfileFinalBestFoodX("PEOA_FinalBestFoodX.csv");
	ofstream outfileFinalBestFoodY("PEOA_FinalBestFoodY.csv");

	ofstream outfileBE("PEOA_BestEagleFit.csv");
	ofstream outfileBF("PEOA_BestFoodFit.csv");
	ofstream outfileTimeBE("PEOA_BETime.csv");
	ofstream outfileTimeBF("PEOA_BFTime.csv");

	ofstream outfileCounterM("PEOA_CounterM.csv");
	ofstream outfileCounterM1("PEOA_CounterM1.csv");
	ofstream outfileCounterM2("PEOA_CounterM2.csv");

	for (int i = 0; i < bestEagleFitTime.at(0).size(); i++)
	{
		for (int j = 0; j < bestEagleFitTime.size(); j++)
		{
			outfileBE << fixed << setprecision(4) << bestEagleFit[j][i];
			outfileTimeBE << fixed << setprecision(4) << bestEagleFitTime[j][i];
			outfileTimeBF << fixed << setprecision(4) << bestFoodFitTime[j][i];
			outfileCounterM << countMovement[j][i];
			outfileCounterM1 << countMut1[j][i];
			outfileCounterM2 << countMut2[j][i];

			if (j != bestEagleFitTime.size() - 1)
			{
				outfileBE << ",";
				outfileTimeBE << ",";
				outfileTimeBF << ",";
				outfileCounterM << ",";
				outfileCounterM1 << ",";
				outfileCounterM2 << ",";
			}
		}
		if (i != bestEagleFitTime.at(0).size() - 1)
		{
			outfileBE << "\n";
			outfileTimeBE << "\n";
			outfileTimeBF << "\n";
			outfileCounterM << "\n";
			outfileCounterM1 << "\n";
			outfileCounterM2 << "\n";
		}
	}

	for (int i = 0; i < bestEagleFitTime.at(0).size(); i++) {
		for (int j = 0; j < bestEagleFitTime.size(); j++) {
			outfileBF << fixed << setprecision(4) << bestFoodFit[j][i];

			if (j != bestEagleFitTime.size() - 1)
			{
				outfileBF << ",";
			}
		}
		if (i != bestEagleFitTime.at(0).size() - 1)
		{
			outfileBF << "\n";
		}
	}

	for (int i = 0; i < initBestEagleSolutions.size(); i++)
	{
		for (int j = 0; j < numFacilities; j++)
		{
			outfileInitBestEagleX << fixed << setprecision(4) << initBestEagleSolutions[i][j][0];
			outfileInitBestEagleY << fixed << setprecision(4) << initBestEagleSolutions[i][j][1];
			outfileInitBestFoodX << fixed << setprecision(4) << initBestFoodSolutions[i][j][0];
			outfileInitBestFoodY << fixed << setprecision(4) << initBestFoodSolutions[i][j][1];

			outfileFinalBestEagleX << fixed << setprecision(4) << finalBestEagleSolutions[i][j][0];
			outfileFinalBestEagleY << fixed << setprecision(4) << finalBestEagleSolutions[i][j][1];
			outfileFinalBestFoodX << fixed << setprecision(4) << finalBestFoodSolutions[i][j][0];
			outfileFinalBestFoodY << fixed << setprecision(4) << finalBestFoodSolutions[i][j][1];

			if (j != numFacilities - 1)
			{
				outfileInitBestEagleX << ",";
				outfileInitBestEagleY << ",";
				outfileInitBestFoodX << ",";
				outfileInitBestFoodY << ",";

				outfileFinalBestEagleX << ",";
				outfileFinalBestEagleY << ",";
				outfileFinalBestFoodX << ",";
				outfileFinalBestFoodY << ",";
			}
		}
		if (i != initBestEagleSolutions.size() - 1)
		{
			outfileInitBestEagleX << "\n";
			outfileInitBestEagleY << "\n";
			outfileInitBestFoodX << "\n";
			outfileInitBestFoodY << "\n";

			outfileFinalBestEagleX << "\n";
			outfileFinalBestEagleY << "\n";
			outfileFinalBestFoodX << "\n";
			outfileFinalBestFoodY << "\n";
		}
	}

	return 0;
}

int operatorAssignment(double P1, double P2, double P3) {
	double r = static_cast<double>(rand()) / RAND_MAX;
	int chosenOperator = 0;
	double mutate1Prob = P1 + P2;

	if (r <= P1) {
		chosenOperator = 1;
	}
	else if (r <= mutate1Prob && P1 < r) {
		chosenOperator = 2;
	}
	else if (r > mutate1Prob) {
		chosenOperator = 3;
	}

	return chosenOperator;
}

double ImprovementRate(vector<int> OperatorIndices, vector<double> OldFitVal, vector<double> NewFitVal) {
	double Numerator = 0.0, Denominator = 0.0;
	double ImprovedRate, FunctionValueDiff;
	int j;

	for (int i = 0; i < OperatorIndices.size(); i++) {
		j = OperatorIndices[i];
		FunctionValueDiff = OldFitVal[j] - NewFitVal[j];

		Denominator += OldFitVal[j];
		Numerator += max(FunctionValueDiff, 0.0);
	}

	if (Numerator == 0.0 && Denominator == 0.0) {
		ImprovedRate = 0.0;
	} else {
		ImprovedRate = Numerator / Denominator;
	}

	return ImprovedRate;
}

void UpdateProbability(double R1, double R2, double R3, double &P1, double &P2, double &P3) {
	double SumR;

	SumR = R1 + R2 + R3;

	P1 = max(0.1, min(0.9, R1 / SumR));
	P2 = max(0.1, min(0.9, R2 / SumR));
	P3 = max(0.1, min(0.9, R3 / SumR));
}

void removeRandomEagleArchive(vector<vector<vector<double>>> &ArchivedEagles) {
	if (!ArchivedEagles.empty()) {
		size_t randomIndex = rand() % ArchivedEagles.size();

		// Erase the eagle
		ArchivedEagles.erase(ArchivedEagles.begin() + randomIndex);
	}
}