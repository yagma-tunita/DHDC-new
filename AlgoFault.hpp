#pragma once

#include "node.hpp"
#include <queue>
#include <vector>
#include <set>

std::queue<node> rans(int dim, int n);
float faultTorran(int dim, int n);
float checker(const std::vector<std::queue<node>>& paths, const std::queue<node>& randomNodes);
float diavg(const std::vector<std::queue<node>>& paths);

std::set<int> faultisis(int n, const std::vector<char>& faulty, unsigned seed);
std::set<int> faultisisConditional(int n, const std::vector<char>& faulty, unsigned seed);
std::set<int> faultisisG(int n, int g, const std::vector<char>& faulty, unsigned seed);

inline int diagOrdinary(int n) { return n; }
inline int diagConditional(int n) { return (2 * n - 1) < (n + 4) ? (2 * n - 1) : (n + 4); }
inline int diagG(int n, int g) { int a = 2 * n + g - 2, b = n + 3 * g + 1; return a < b ? a : b; }
