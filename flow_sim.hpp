#pragma once

#include <vector>
#include <functional>
#include <utility>

namespace flowsim
{

struct Flow
{
	int id;
	int src;
	int dst;
	std::vector<int> path;
	std::vector<int> links;
	double size;
	double tArr;
	double remaining;
	double tFinish;
};

using RouteFn = std::function<std::vector<int>(int src, int dst)>;

struct Config
{
	int N = 0;
	double linkCap = 1e9;
	RouteFn route;
};

struct Result
{
	double maxThp;
	double avgThp;
	double totalFCT;
	double avgHops;
	std::vector<std::pair<double, double>> thpCurve;
};

Result run(const Config& cfg, std::vector<Flow> flows);

}
