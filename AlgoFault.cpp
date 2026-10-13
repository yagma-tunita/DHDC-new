#include "AlgoFault.hpp"
#include "AlgoRouting.hpp"
#include <random>
#include <unordered_set>
#include <map>
#include <algorithm>

std::queue<node> rans(int dim, int n) 
{
	std::random_device rd;
	std::mt19937 gen(rd());

	int path_len = dim - 1;
	std::uniform_int_distribution<int> tag_dist(1, dim - 1);
	std::uniform_int_distribution<int> bool_dist(0, 1);
	std::uniform_int_distribution<int> bit_dist(0, 1);
	std::queue<node> q;

	for (int i = 0; i < n; ++i) 
	{
		int pathA = 0, pathB = 0;
		for (int j = 0; j < path_len; ++j) 
		{
			if (bit_dist(gen)) pathA |= (1 << j);
			if (bit_dist(gen)) pathB |= (1 << j);
		}
		int tag = tag_dist(gen);
		bool tail = bool_dist(gen) == 1;
		q.emplace(dim, pathA, pathB, tag, tail);
	}
	return q;
}

float faultTorran(int dim, int n) 
{
	const int TEST_TIMES = 20;
	float totalRate = 0.0f;

	for (int test = 0; test < TEST_TIMES; ++test) 
	{
		std::queue<node> nodes = rans(dim, 2);
		if (nodes.size() != 2) {
			continue;
		}
		node u = nodes.front();
		nodes.pop();
		node v = nodes.front();

		int retry = 0;
		while (u == v && retry < 5) {
			nodes = rans(dim, 2);
			u = nodes.front();
			nodes.pop();
			v = nodes.front();
			++retry;
		}
		if (u == v)continue;
		std::vector<std::queue<node>> paths = StrictDisjointpath(u, v);
		if (paths.empty())continue;
		std::queue<node> faultNodes = rans(dim, n);
		float rate = checker(paths, faultNodes);
		totalRate += rate;
	}

	return totalRate / static_cast<float>(TEST_TIMES);
}

float checker(const std::vector<std::queue<node>>& paths, const std::queue<node>& randomNodes) 
{
	if (paths.empty()) return 0.0f;

	std::unordered_set<node, NodeHash> randomSet;
	std::queue<node> randomCopy = randomNodes;
	while (!randomCopy.empty()) {
		randomSet.insert(randomCopy.front());
		randomCopy.pop();
	}

	int failCount = 0;
	int total = static_cast<int>(paths.size());

	for (const auto& pathQueue : paths) 
	{
		std::queue<node> pathCopy = pathQueue;
		bool failed = false;
		while (!pathCopy.empty()) 
		{
			const node& nd = pathCopy.front();
			if (randomSet.find(nd) != randomSet.end()) 
			{
				failed = true;
				break;
			}
			pathCopy.pop();
		}
		if (failed) 
		{
			++failCount;
		}
	}
	return static_cast<float>(failCount) / static_cast<float>(total);
}

float diavg(const std::vector<std::queue<node>>& paths) {
	if (paths.empty()) return 0.0f;

	size_t totalNodes = 0;
	for (const auto& pathQueue : paths) {
		totalNodes += pathQueue.size();
	}
	return static_cast<float>(totalNodes) / static_cast<float>(paths.size());
}

struct DHDCGraph {
	int n, m, half, N;
	std::vector<std::vector<int>> nbr;
	std::vector<std::vector<int>> ridx;
	DHDCGraph(int d) : n(d), m(d - 1), half(1 << (d - 1)), N(m * half * half * 2) {
		nbr.resize(N);
		ridx.resize(N);
		for (int id = 0; id < N; ++id) nbr[id] = neighbors(id);
		for (int u = 0; u < N; ++u) {
			ridx[u].resize(nbr[u].size());
			for (int k = 0; k < (int)nbr[u].size(); ++k) {
				int v = nbr[u][k];
				for (int k2 = 0; k2 < (int)nbr[v].size(); ++k2) if (nbr[v][k2] == u) { ridx[u][k] = k2; break; }
			}
		}
	}
	std::vector<int> neighbors(int id) const {
		int tag1 = id % m, rest = id / m;
		int tail = rest % 2; rest /= 2;
		int b = rest % half, a = rest / half;
		int tag = tag1 + 1;
		std::vector<int> nb;
		for (int t = 1; t <= m; ++t) if (t != tag) nb.push_back(((a * half + b) * 2 + tail) * m + (t - 1));
		nb.push_back(((a * half + b) * 2 + (1 - tail)) * m + tag1);
		if (tail == 0) { int a2 = a ^ (1 << (tag - 1)); nb.push_back(((a2 * half + b) * 2 + tail) * m + tag1); }
		else { int b2 = b ^ (1 << (tag - 1)); nb.push_back(((a * half + b2) * 2 + tail) * m + tag1); }
		return nb;
	}
};

namespace {
	using Res = std::vector<std::vector<char>>;
	Res makeSynd(const DHDCGraph& g, const std::vector<char>& faulty, std::mt19937& rng) {
		Res res(g.N);
		for (int u = 0; u < g.N; ++u) {
			res[u].resize(g.nbr[u].size());
			for (int k = 0; k < (int)g.nbr[u].size(); ++k) {
				int v = g.nbr[u][k];
				res[u][k] = faulty[u] ? (char)(rng() & 1) : (char)(faulty[v] ? 1 : 0);
			}
		}
		return res;
	}
	std::vector<int> mutual0Components(const DHDCGraph& g, const Res& res) {
		std::vector<int> par(g.N);
		for (int i = 0; i < g.N; ++i) par[i] = i;
		auto findRoot = [&](int x) { while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; } return x; };
		for (int u = 0; u < g.N; ++u) {
			for (int k = 0; k < (int)g.nbr[u].size(); ++k) {
				int v = g.nbr[u][k];
				if (u > v) continue;
				char ru = res[u][k];
				char rv = res[v][g.ridx[u][k]];
				if (ru == 0 && rv == 0) { int a = findRoot(u), b = findRoot(v); if (a != b) par[b] = a; }
			}
		}
		std::vector<int> comp(g.N);
		for (int i = 0; i < g.N; ++i) comp[i] = findRoot(i);
		return comp;
	}
}

static std::set<int> diagCore(int n, const std::vector<char>& faulty, unsigned seed, int gg, int defKind) {
	DHDCGraph g(n);
	std::mt19937 rng(seed);
	Res res = makeSynd(g, faulty, rng);
	std::vector<int> comp = mutual0Components(g, res);
	std::vector<int> size(g.N, 0);
	for (int i = 0; i < g.N; ++i) ++size[comp[i]];
	int best = -1, bestsz = -1;
	for (int i = 0; i < g.N; ++i) if (comp[i] == i && size[i] > bestsz) { bestsz = size[i]; best = i; }
	std::vector<std::vector<int>> cmembers(g.N);
	for (int i = 0; i < g.N; ++i) cmembers[comp[i]].push_back(i);
	std::vector<char> H(g.N, 0), F(g.N, 0);
	for (int i = 0; i < g.N; ++i) if (comp[i] == best) H[i] = 1; else F[i] = 1;

	bool changed = true;
	while (changed) {
		changed = false;
		for (int u = 0; u < g.N; ++u) if (H[u]) {
			for (int k = 0; k < (int)g.nbr[u].size(); ++k) {
				int v = g.nbr[u][k];
				char r = res[u][k];
				if (r == 0 && !H[v]) { H[v] = 1; F[v] = 0; changed = true; }
				else if (r == 1 && !F[v]) { H[v] = 0; F[v] = 1; changed = true; }
			}
		}
		for (int v = 0; v < g.N; ++v) if (F[v]) {
			for (int k = 0; k < (int)g.nbr[v].size(); ++k) {
				int u = g.nbr[v][k];
				char r = res[u][g.ridx[v][k]];
				if (r == 0 && !F[u]) { H[u] = 0; F[u] = 1; changed = true; }
			}
		}
		for (int r = 0; r < g.N; ++r) {
			const std::vector<int>& C = cmembers[r];
			if (C.empty()) continue;
			bool fullyF = true;
			for (int u : C) if (!F[u]) { fullyF = false; break; }
			if (!fullyF) continue;
			bool ok = true;
			for (int u : C) {
				for (int k = 0; k < (int)g.nbr[u].size(); ++k) {
					int v = g.nbr[u][k];
					if (comp[v] == comp[u]) continue;
					if (!F[v]) { ok = false; break; }
					if (res[u][k] != 1) { ok = false; break; }
				}
				if (!ok) break;
			}
			if (!ok) continue;
			if (defKind >= 1) {
				int need = (defKind == 1) ? 1 : gg;
				for (int u : C) {
					int c = 0;
					for (int v : g.nbr[u]) if (comp[v] == comp[u]) ++c;
					if (c < need) { ok = false; break; }
				}
			}
			if (ok) for (int u : C) { if (!H[u]) { H[u] = 1; F[u] = 0; changed = true; } }
		}
	}

	if (defKind == 1) {
		bool chg = true;
		while (chg) {
			chg = false;
			for (int v = 0; v < g.N; ++v) if (F[v]) {
				bool hasH = false;
				for (int u : g.nbr[v]) if (H[u]) { hasH = true; break; }
				if (!hasH) { H[v] = 1; F[v] = 0; chg = true; }
			}
		}
	}

	std::set<int> out;
	for (int i = 0; i < g.N; ++i) if (F[i]) out.insert(i);
	return out;
}

std::set<int> faultisis(int n, const std::vector<char>& faulty, unsigned seed) {
	return diagCore(n, faulty, seed, 0, 0);
}

std::set<int> faultisisConditional(int n, const std::vector<char>& faulty, unsigned seed) {
	return diagCore(n, faulty, seed, 1, 1);
}

std::set<int> faultisisG(int n, int gg, const std::vector<char>& faulty, unsigned seed) {
	return diagCore(n, faulty, seed, gg, 2);
}
