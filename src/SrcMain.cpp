#include "SrcMain.h"
#include <iostream>
#include <random>
#include "TSP.h"
#include <fstream>
#include <algorithm>
#include <sstream>
#include <numeric>


/*
 * Helper function to read in all locations from inputfile
*/
std::vector<Location> ReadLocations(const std::string fileName)
{
	std::vector<Location> res;
	std::ifstream ifile(fileName);

	std::string line;
	if (ifile.is_open()) {
		while (static_cast<bool>(std::getline(ifile, line)))
		{
			//read in each line
			std::stringstream ss(line);
			std::string name;
			std::string latitude;
			std::string longitude;

			std::getline(ss, name, ',');
			std::getline(ss, latitude, ',');
			std::getline(ss, longitude, ',');

			//create the object
			Location toAdd;
			toAdd.mName = name;
			toAdd.mLatitude = std::stod(latitude);
			toAdd.mLongitude = std::stod(longitude);

			res.emplace_back(toAdd);
		}
	}

	ifile.close();

	return res;
}


/*
 * Function to initialize populations
*/ 
Population InitPolulations(const int popSize, const int numLocations, std::mt19937& RNG) 
{
	Population p; //contruct all populations

	std::vector<std::vector<int> > outerV(popSize);
	std::generate(outerV.begin(), outerV.end(), [&]() {
		std::vector<int> innerV(numLocations);
		std::iota(innerV.begin(), innerV.end(), 0);
		std::shuffle(innerV.begin() + 1, innerV.end(), RNG);
		p.mMembers.emplace_back(innerV);
		return innerV;
	});

	return p;
}

/*
 * Helper function to calculate the distance between two locations
*/
double CalcDistance(const Location& la, const Location& lb)
{
	const double convertDTR = 0.0174533;
	double lat1 = la.mLatitude*convertDTR;
	double lat2 = lb.mLatitude*convertDTR;
	double lon1 = la.mLongitude*convertDTR;
	double lon2 = lb.mLongitude*convertDTR;

	double dlon = lon2 - lon1;
	double dlat = lat2 - lat1;

	double partialA = std::pow((sin(dlat / 2)), 2) + cos(lat1) * cos(lat2) * pow((sin(dlon / 2)), 2);
	double partialC = 2 * atan2(sqrt(partialA), sqrt(1 - partialA));

	return 3961 * partialC;
}

/*
 * Helper function to calculate sum of distances between all pairs of location
*/
double CalcAllDistance(const std::vector<int>& individualPop, const std::vector<Location>& locations)
{	
	int numLocations = locations.size();
	std::vector<double> diff;
	std::adjacent_difference(individualPop.begin(), individualPop.end(), std::back_inserter(diff), [&locations](int a, int b) {
		Location la = locations[a];
		Location lb = locations[b];

		return CalcDistance(la,  lb);
	});

	//calculate the distance between last element and first one
	Location la = locations[individualPop[numLocations - 1]];
	Location lb = locations[individualPop[0]];
	double lastDist = CalcDistance(la, lb);
	diff[0] = lastDist;

	double fitness = std::accumulate(diff.begin(), diff.end(), 0.0, [](const double&a, const double& b) {
		return a + b;
	});

	return fitness;
}

/*
 * Function to compute the fitness numerics
 * Return vector of pairs, first = individual in the population, second = fitness of that individual
*/
std::vector<std::pair<int, double> > ComputeFitness(const Population& p, const int popSize, 
									const std::vector<Location>& locations, const int numLocations)
{
	std::vector<std::pair<int, double> > res;
	std::vector<int> helper(popSize);

	int index = 0;
	std::generate(helper.begin(), helper.end(), [&]() {
		std::vector<int> individualPop = p.mMembers[index]; //get each population
		double fitness = CalcAllDistance(individualPop, locations);
		res.emplace_back(std::make_pair(index, fitness));
		index++;
		return index;;
	});
	return res;
}

/*
 * Function to do Selection based on fitness value calculated
*/
std::vector<std::pair<int, int> > SelectPairs(const Population& p, const std::vector<std::pair<int, double> >& fitness, std::mt19937& RNG)
{
	const size_t popSize = p.mMembers.size();
	std::vector<std::pair<int, int> > res(popSize);

	std::vector<double> prob(popSize);
	std::generate(prob.begin(), prob.end(), [&popSize]() {
		return (1.0) / popSize; //equal probability
	});

	//multiply the first two by 6.0
	prob[fitness[0].first] = prob[fitness[0].first] * 6.0;
	prob[fitness[1].first] = prob[fitness[1].first] * 6.0;

	//multiply the reaminding top-half by 3.0
	for (unsigned i = 2; i <= (popSize / 2) - 1; ++i) {
		prob[fitness[i].first] = prob[fitness[i].first] * 3.0;
	}


	//find the total probability needed to renormalize
	double totalProb = std::accumulate(prob.begin(), prob.end(), 0.0, [](const double& a, const double& b) {
		return a + b;
	});

	//normalize the vector to have totalProb 1.0
	std::vector<double> normProb;
	std::transform(prob.begin(), prob.end(), std::back_inserter(normProb), [totalProb](const float& a) {
		return (a / totalProb);
	});

	std::generate(res.begin(), res.end(), [&](){
		unsigned parentOne;
		unsigned parentTwo;
		double sum = 0.0;
		std::uniform_real_distribution<double> d(0, 1);
		double randomNum = d(RNG);//generate an random number
		for (unsigned i = 0; i < popSize; ++i) 
		{
			sum += normProb[i];
			if (sum >= randomNum) {
				parentOne = i;
				break;
			}
		}

		sum = 0.0;
		double randomNum2 = d(RNG);
		for (unsigned i = 0; i < popSize; ++i)
		{
			sum += normProb[i];
			if (sum >= randomNum2) {
				parentTwo = i;
				break;
			}
		}
		return std::make_pair(parentOne, parentTwo);
	});

	return res;
}

/*
 * Function to do the crossover generation
*/
Population CrossoverGen(const Population& p, const std::vector<std::pair<int, int> >& selectedPairs, std::mt19937& RNG, const int& mutationChance)
{	
	Population newPopulation;
	unsigned int locationSize = p.mMembers[0].size();
	

	//populate new population with selected pairs
	std::transform(selectedPairs.begin(), selectedPairs.end(), std::back_inserter(newPopulation.mMembers), [&](const std::pair<int, int>& s){
		//two parent candiates
		std::vector<int> parentA = p.mMembers[s.first];
		std::vector<int> parentB = p.mMembers[s.second];

		std::uniform_int_distribution<int> d1(1, locationSize - 2);
		int crossIndex = d1(RNG);
		std::uniform_int_distribution<int> d2(0, 1);
		int toss = d2(RNG);

		std::vector<int> firstParent, secondParent;
		if (toss != 0)
		{
			firstParent = parentA;
			secondParent = parentB;
		}
		else
		{
			firstParent = parentB;
			secondParent = parentA;
		}

		std::vector<int> child;
		std::copy_n(firstParent.begin(), crossIndex + 1, std::back_inserter(child));

		//only copy elements are not in child now
		std::copy_if(secondParent.begin(), secondParent.end(), std::back_inserter(child), [&child](const int& a) {
			return (std::find(child.begin(), child.end(), a) == child.end());
		});

		std::uniform_real_distribution<double> d3(0, 1);

		double mutate = d3(RNG);
		if (mutate <= mutationChance / 100.0) {
			std::uniform_int_distribution<int> d4(1, locationSize - 1);
			int firstIndex = d4(RNG);
			int secondIndex = d4(RNG);
			std::swap(child[firstIndex], child[secondIndex]);
		}

		//insert the new crossed child to new population
		return child;
	});

	return newPopulation;
}

/*
 * All of the Output functions
*/


void OutputSelectedPairs(const std::vector<std::pair<int, int> > selected, std::ofstream& ofile)
{
	ofile << "SELECTED PAIRS:" << std::endl;
	for (const auto& i : selected)
	{
		ofile << "(" << i.first << "," << i.second << ")" << std::endl;
	}
}

void OutputFitness(const std::vector<std::pair<int, double> >& fitness, std::ofstream& ofile)
{
	ofile << "FITNESS:" << std::endl;
	for (const auto& i : fitness)
	{
		ofile << i.first << ":" << i.second << std::endl;
	}
}

void OutputRoute(const std::vector<int>& v, std::ofstream& ofile) 
{
	for (const auto& i : v) {
		if (i != v[v.size() - 1]) {
			ofile << i << ",";
		} 
		else //last element
		{
			ofile << i << std::endl;
		}
	}
}

void OutputPopulation(const Population& p, std::ofstream& ofile)
{
	for (unsigned i = 0; i < p.mMembers.size(); ++i)
	{
		OutputRoute(p.mMembers[i], ofile);
	}
}

void FinalOutput(std::ofstream& ofile, const Population& lastPopulation, const int& popSize, const std::vector<Location>& vLocations)
{
	//output the final route
	std::vector<std::pair<int, double> > fitnessLast;
	fitnessLast = ComputeFitness(lastPopulation, popSize, vLocations, vLocations.size());
	OutputFitness(fitnessLast, ofile);

	//sort the fitness vector, ascending
	std::sort(fitnessLast.begin(), fitnessLast.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
		return (a.second < b.second);
	});
	int bestRouteIndex = fitnessLast[0].first;
	std::vector<int> bestRoute = lastPopulation.mMembers[bestRouteIndex];

	ofile << "SOLUTION:" << std::endl;
	double totalDist = CalcAllDistance(bestRoute, vLocations);
	for (const auto& i : bestRoute)
	{
		ofile << vLocations[i].mName << std::endl;
	}
	ofile << vLocations[bestRoute[0]].mName << std::endl; //print out the starting point
	ofile << "DISTANCE: " << totalDist << " miles" << std::endl;
}

/*	
 * Main Training Function
*/
void TrainingIterations(std::ofstream& ofile, const Population& p, const std::vector<Location>& vLocations, const int& generations, 
						std::mt19937& randGen, const int& popSize, const int& mutationChance)
{
	Population lastPopulation = p;
	int counter = 0;
	while (counter < generations)
	{
		if (counter == 0)
		{
			ofile << "INITIAL POPULATION:" << std::endl;
		}
		else
		{
			ofile << "GENERATION: " << counter << std::endl;
		}
		OutputPopulation(lastPopulation, ofile);

		std::vector<std::pair<int, double> > fitness = ComputeFitness(lastPopulation, popSize, vLocations, vLocations.size());;
		OutputFitness(fitness, ofile);


		//sort the fitness to select pairs
		std::sort(fitness.begin(), fitness.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
			return (a.second < b.second);
		});
		std::vector<std::pair<int, int> > selected = SelectPairs(lastPopulation, fitness, randGen);;

		OutputSelectedPairs(selected, ofile);

		//Population newPopulation;
		Population newPopulation = CrossoverGen(lastPopulation, selected, randGen, mutationChance);
		lastPopulation = newPopulation;

		counter++;
	}

	//last iteration
	ofile << "GENERATION: " << counter << std::endl;
	OutputPopulation(lastPopulation, ofile);
	
	//output the final solution
	FinalOutput(ofile, lastPopulation, popSize, vLocations);
}


/*
 * Input Argument Processing and calling main training function
*/
void ProcessCommandArgs(int argc, const char* argv[])
{
	// check enough argument counts
	if (argc != 6) 
	{
		std::cerr << "Not Enought CMD arguments" << std::endl;
		return;
	}
	//read in CMD arguments
	std::string fileName = argv[1];
	int popSize = std::stoi(argv[2]);
	int generations = std::stoi(argv[3]);
	int mutationChance = std::stoi(argv[4]);
	int seed = std::stoi(argv[5]);
	std::ofstream ofile("log.txt");

	//create RNG
	std::mt19937 randGen(seed);

	//read in input locations
	std::vector<Location> vLocations = ReadLocations(fileName);; //all location data stored in this vector

	//create inital population
	Population p = InitPolulations(popSize, vLocations.size(), randGen);

	//main training function and output
	TrainingIterations(ofile, p, vLocations, generations, randGen, popSize, mutationChance);
}
