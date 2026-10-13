#pragma once

#include "node.hpp"
#include <queue>
#include <vector>
#include <cassert>

void print_path(std::queue<node>path);
void countAndPrintNodes(std::vector<std::vector<node>>& nodes);
bool isconnected(node u, node v);

std::queue<node> DHRouting(node u, node v);
std::vector<std::queue<node>> Disjointpath(node u, node v);
std::vector<std::queue<node>> StrictDisjointpath(node u, node v);
std::queue<node> BFS(node u, node v);
std::queue<node> BiDirBFS(node u, node v);
std::queue<node> AStar(node u, node v);
std::queue<node> Dijkstra(node u, node v);
