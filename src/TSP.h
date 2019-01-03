#pragma once
#include <string>
#include <vector>
#include <random>

struct Location
{
	std::string mName;
	double mLatitude = 0.0;
	double mLongitude = 0.0;
};

struct Population
{
	std::vector<std::vector<int>> mMembers;
};
