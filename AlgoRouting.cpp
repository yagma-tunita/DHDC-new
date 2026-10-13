#include "AlgoRouting.hpp"
#include <algorithm>
#include <string>
#include <optional>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <random>

void print_path(std::queue<node>path)
{
	while (!path.empty())
	{
		path.front().print_node();
		std::cout << "->";
		path.pop();
	}
	std::cout << "end";
}
void countAndPrintNodes(std::vector<std::vector<node>>& nodes)
{
	std::map<node, int> nodeCount;
	for (auto line : nodes)
		for (node n : line)
			nodeCount[n]++;
	for (auto& pair : nodeCount)
	{
		const_cast<node&>(pair.first).print_node();
		std::cout << "   " << pair.second << '\n';
	}
}


void DIF1(node u, node v, int& A, int& B)
{
	A = u.pathA ^ v.pathA;
	B = u.pathB ^ v.pathB;
}

void DIF2(node u, node v, int& A, int& B)
{
	A = u.pathA ^ v.pathA;
	B = u.pathB ^ v.pathB;
}

bool isconnected(node u, node v)
{
	int A, B;
	DIF2(u, v, A, B);
	if (A == 0 && B == 0)
	{
		if ((u.tail == v.tail && u.tag != v.tag) || (u.tail != v.tail && u.tag == v.tag)) return true;
	}
	
	if (u.tail == v.tail && u.tag == v.tag)
	{
		if (u.tail == 0 && A != 0 && (A & (A - 1)) == 0 && B == 0)
		{
			if (__builtin_ctz((unsigned)A) == u.tag - 1) return true;
		}

		if (u.tail == 1 && A == 0 && B != 0 && (B & (B - 1)) == 0)
		{
			if (__builtin_ctz((unsigned)B) == u.tag - 1) return true;
		}
	}

	return false;
}



std::queue<node> DHRouting(node u, node v)
{
	int A, B;
	int n = u.dimension, x = v.tag;
	DIF1(u, v, A, B);

	std::vector<int> TA, TB;
	for (int i = 0; i < n - 1; ++i)
	{
		if ((A >> i) & 1) TA.push_back(i + 1);
		if ((B >> i) & 1) TB.push_back(i + 1);
	}

	int t = u.tag;
	bool sameTail = (u.tail == v.tail);
	std::vector<int> T1 = (u.tail == false) ? TA : TB;
	std::vector<int> T2 = (u.tail == false) ? TB : TA;

	auto cliqueMoves = [](const std::vector<int>& S, int p, int q) {
		int k = (int)S.size();
		if (k == 0) return (p != q) ? 1 : 0;
		bool pIn = std::find(S.begin(), S.end(), p) != S.end();
		bool qIn = std::find(S.begin(), S.end(), q) != S.end();
		if (p == q) { if (pIn) return (k == 1) ? 0 : k; return k + 1; }
		if (pIn && qIn) return k - 1;
		if (pIn || qIn) return k;
		return k + 1;
	};
	auto cliqueOrder = [](const std::vector<int>& S, int p, int q) {
		std::vector<int> res;
		res.push_back(p);
		if (S.empty()) { if (p != q) res.push_back(q); return res; }
		bool pIn = std::find(S.begin(), S.end(), p) != S.end();
		bool qIn = std::find(S.begin(), S.end(), q) != S.end();
		if (p == q) {
			if (pIn) {
				bool any = false;
				for (int s : S) if (s != p) { res.push_back(s); any = true; }
				if (any) res.push_back(p);
			}
			else {
				for (int s : S) res.push_back(s);
				res.push_back(p);
			}
		}
		else if (pIn && qIn) {
			for (int s : S) if (s != p && s != q) res.push_back(s);
			res.push_back(q);
		}
		else if (pIn) {
			for (int s : S) if (s != p) res.push_back(s);
			res.push_back(q);
		}
		else if (qIn) {
			for (int s : S) if (s != q) res.push_back(s);
			res.push_back(q);
		}
		else {
			for (int s : S) res.push_back(s);
			res.push_back(q);
		}
		return res;
	};

	struct Phase { bool tail; std::vector<int> tags; };
	std::vector<Phase> phases;
	bool t1Empty = T1.empty(), t2Empty = T2.empty();
	const int INF = 1000000000;

	if (t1Empty && t2Empty)
		phases.push_back({ u.tail, cliqueOrder(std::vector<int>(), t, x) });
	else if (!t1Empty && t2Empty)
		phases.push_back({ u.tail, cliqueOrder(T1, t, x) });
	else if (t1Empty && !t2Empty)
	{
		phases.push_back({ u.tail, cliqueOrder(std::vector<int>(), t, t) });
		phases.push_back({ !u.tail, cliqueOrder(T2, t, x) });
		if (sameTail)
			phases.push_back({ u.tail, cliqueOrder(std::vector<int>(), x, x) });
	}
	else if (!sameTail)
	{
		int best = INF, bc = -1;
		for (int c = 1; c <= n - 1; ++c)
		{
			int cost = cliqueMoves(T1, t, c) + 1 + cliqueMoves(T2, c, x);
			if (cost < best) { best = cost; bc = c; }
		}
		phases.push_back({ u.tail, cliqueOrder(T1, t, bc) });
		phases.push_back({ !u.tail, cliqueOrder(T2, bc, x) });
	}
	else
	{
		auto cliqCost = [](int m, int p, int q, bool pIn, bool qIn) {
			if (m == 0) return (p != q) ? 1 : 0;
			if (p == q) return pIn ? ((m == 1) ? 0 : m) : (m + 1);
			if (pIn && qIn) return m - 1;
			if (pIn || qIn) return m;
			return m + 1;
		};

		int a = (int)T1.size();
		int best = INF, bc1 = -1, bc2 = -1;
		std::vector<int> bestA1, bestA2;

		for (int c1 = 1; c1 <= n - 1; ++c1)
			for (int c2 = 1; c2 <= n - 1; ++c2)
			{
				std::vector<int> special;
				for (int s : { t, c1, c2, x })
					if (std::find(T1.begin(), T1.end(), s) != T1.end()
						&& std::find(special.begin(), special.end(), s) == special.end())
						special.push_back(s);

				std::vector<int> freeTags;
				for (int s : T1)
					if (std::find(special.begin(), special.end(), s) == special.end())
						freeTags.push_back(s);

				int ks = (int)special.size();
				int r = (int)freeTags.size();

				for (int mask = 0; mask < (1 << ks); ++mask)
				{
					int nA = 0;
					bool tInA = false, c1InA = false, c2InC = false, xInC = false;
					for (int i = 0; i < ks; ++i)
					{
						int s = special[i];
						if ((mask >> i) & 1)
						{
							++nA;
							if (s == t) tInA = true;
							if (s == c1) c1InA = true;
						}
						else
						{
							if (s == c2) c2InC = true;
							if (s == x) xInC = true;
						}
					}

					int lo = nA, hi = nA + r;
					int cand[4] = { lo, hi, 1, a - 1 };
					for (int ci = 0; ci < 4; ++ci)
					{
						int m = cand[ci];
						if (m < lo || m > hi) continue;

						int cost = cliqCost(m, t, c1, tInA, c1InA)
							+ cliqCost(a - m, c2, x, c2InC, xInC)
							+ 1 + cliqueMoves(T2, c1, c2) + 1;

						if (cost < best)
						{
							best = cost;
							bc1 = c1;
							bc2 = c2;

							std::vector<int> A1, A2;
							for (int i = 0; i < ks; ++i)
								((mask >> i) & 1) ? A1.push_back(special[i]) : A2.push_back(special[i]);
							int need = m - nA;
							for (int s : freeTags)
							{
								if (need > 0) { A1.push_back(s); --need; }
								else A2.push_back(s);
							}
							bestA1 = A1;
							bestA2 = A2;
						}
					}
				}
			}

		phases.push_back({ u.tail, cliqueOrder(bestA1, t, bc1) });
		phases.push_back({ !u.tail, cliqueOrder(T2, bc1, bc2) });
		phases.push_back({ u.tail, cliqueOrder(bestA2, bc2, x) });
	}

	std::queue<node> path;
	path.push(u);
	auto setTag = [&](int k) { if (u.tag != k) { u.tag = k; path.push(u); } };
	for (size_t pi = 0; pi < phases.size(); ++pi)
	{
		if (u.tail != phases[pi].tail) { u.tail = phases[pi].tail; path.push(u); }
		for (int tg : phases[pi].tags)
		{
			setTag(tg);
			if (phases[pi].tail == false)
			{
				if ((A >> (tg - 1)) & 1) { u.pathA ^= (1 << (tg - 1)); path.push(u); A &= ~(1 << (tg - 1)); }
			}
			else
			{
				if ((B >> (tg - 1)) & 1) { u.pathB ^= (1 << (tg - 1)); path.push(u); B &= ~(1 << (tg - 1)); }
			}
		}
	}
	setTag(x);
	if (u.tail != v.tail) { u.tail = v.tail; path.push(u); }

	return path;
}


std::vector<std::queue<node>> StrictDisjointpath(node u, node v)
{
	int n = u.dimension;
	int m = n - 1, half = 1 << (n - 1);
	int N = m * half * half * 2;

	auto pack = [&](const node& nd) {
		int a = nd.pathA, b = nd.pathB;
		return ((a * half + b) * 2 + (nd.tail ? 1 : 0)) * m + (nd.tag - 1);
	};

	static std::map<int, std::vector<std::vector<int>>> cache;
	std::vector<std::vector<int>> adj;
	auto it = cache.find(n);
	if (it == cache.end())
	{
		adj.assign(N, {});
		for (int id = 0; id < N; ++id)
		{
			int tag1 = id % m, rest = id / m;
			int tail = rest % 2; rest /= 2;
			int b = rest % half, a = rest / half;
			int tag = tag1 + 1;
			for (int t = 1; t <= m; ++t)
				if (t != tag) adj[id].push_back(((a * half + b) * 2 + tail) * m + (t - 1));
			adj[id].push_back(((a * half + b) * 2 + (1 - tail)) * m + tag1);
			if (tail == 0) { int a2 = a ^ (1 << (tag - 1)); adj[id].push_back(((a2 * half + b) * 2 + tail) * m + tag1); }
			else { int b2 = b ^ (1 << (tag - 1)); adj[id].push_back(((a * half + b2) * 2 + tail) * m + tag1); }
		}
		cache[n] = adj;
	}
	else adj = it->second;

	int uid = pack(u), vid = pack(v);
	int V = 2 * N;
	std::vector<int> head(V, -1), to, cap, nxt;
	auto addEdge = [&](int x, int y) {
		to.push_back(y); cap.push_back(1); nxt.push_back(head[x]); head[x] = (int)to.size() - 1;
		to.push_back(x); cap.push_back(0); nxt.push_back(head[y]); head[y] = (int)to.size() - 1;
	};
	for (int x = 0; x < N; ++x) addEdge(2 * x, 2 * x + 1);
	for (int x = 0; x < N; ++x) for (int y : adj[x]) if (x < y) { addEdge(2 * x + 1, 2 * y); addEdge(2 * y + 1, 2 * x); }

	int S = 2 * uid + 1, T = 2 * vid;
	std::vector<int> parE(V, -1);
	while (true)
	{
		std::fill(parE.begin(), parE.end(), -1);
		std::queue<int> q; parE[S] = -2; q.push(S);
		while (!q.empty() && parE[T] < 0)
		{
			int x = q.front(); q.pop();
			for (int e = head[x]; e != -1; e = nxt[e])
				if (cap[e] > 0 && parE[to[e]] < 0) { parE[to[e]] = e; q.push(to[e]); }
		}
		if (parE[T] < 0) break;
		for (int x = T; x != S; ) { int e = parE[x]; cap[e] -= 1; cap[e ^ 1] += 1; x = to[e ^ 1]; }
	}

	auto unode = [&](int id) {
		int tag1 = id % m, rest = id / m;
		int tail = rest % 2; rest /= 2;
		int b = rest % half, a = rest / half;
		return node(n, a, b, tag1 + 1, tail != 0);
	};

	std::vector<std::queue<node>> result;
	for (int e = head[S]; e != -1; e = nxt[e])
	{
		if (cap[e] > 0) continue;
		std::vector<node> path;
		path.push_back(u);
		int cur = to[e];
		while (cur != T)
		{
			int orig = cur / 2;
			node nd = unode(orig);
			if (path.empty() || !(path.back() == nd)) path.push_back(nd);
			int nv = -1;
			for (int e2 = head[cur]; e2 != -1; e2 = nxt[e2])
				if (cap[e2] == 0 && (e2 % 2 == 0)) { nv = to[e2]; break; }
			if (nv < 0) break;
			cur = nv;
		}
		if (cur == T)
		{
			if (path.empty() || !(path.back() == v)) path.push_back(v);
			std::queue<node> q;
			for (auto& nd : path) q.push(nd);
			result.push_back(q);
		}
	}
	return result;
}


static void eq_path(std::queue<node>& c, node& p) { c.push(p); }
static void neTag(std::queue<node>& c, node& p, int x) { if (p.tag != x) { p.tag = x; c.push(p); } }
static void neTail(std::queue<node>& c, node& p, bool x) { if (p.tail != x) { p.tail = x; c.push(p); } }
static void flipTail(std::queue<node>& c, node& p) { p.tail = !p.tail; c.push(p); }
static void flipBit(std::queue<node>& c, node& p, int j) {
	if (p.tail == false) p.pathA ^= (1 << j);
	else p.pathB ^= (1 << j);
	c.push(p);
}

static void mod_path(node& p, std::set<int> S, std::queue<node>& c, int tip) {
	neTag(c, p, tip); flipBit(c, p, tip - 1);
	int pi = tip - 1; bool in = S.count(pi);
	for (int i : S) { if (i <= pi) continue; neTag(c, p, i + 1); flipBit(c, p, i); }
	for (int i : S) { if (i >= pi) break; neTag(c, p, i + 1); flipBit(c, p, i); }
	if (!in) { neTag(c, p, tip); flipBit(c, p, tip - 1); }
}

std::vector<std::queue<node>> Disjointpath(node u, node v)
{
	int n = u.dimension;
	std::set<int> A, B;
	int maskA = u.pathA ^ v.pathA, maskB = u.pathB ^ v.pathB;
	for (int i = 0; i < n - 1; ++i) { if ((maskA >> i) & 1) A.insert(i); }
	for (int i = 0; i < n - 1; ++i) { if ((maskB >> i) & 1) B.insert(i); }
	bool tu = u.tail, tv = v.tail;
	std::vector<std::queue<node>> result;

	if (A.empty() && B.empty()) {
		for (int i = 0; i <= n - 2; ++i) {
			node p = u; std::queue<node> c; eq_path(c, p);
			neTag(c, p, i + 1); neTail(c, p, tv); neTag(c, p, v.tag);
			if (tu != tv || u.tag != i + 1) result.push_back(c);
		}
		if (tu == tv) {
			{ node p = u; std::queue<node> c; eq_path(c, p); flipTail(c, p); neTag(c, p, v.tag); flipTail(c, p); result.push_back(c); }
			{ node p = u; std::queue<node> c; eq_path(c, p);
				flipBit(c, p, u.tag - 1); neTag(c, p, v.tag); flipBit(c, p, v.tag - 1); neTag(c, p, u.tag);
				flipBit(c, p, u.tag - 1); neTag(c, p, v.tag); flipBit(c, p, v.tag - 1); result.push_back(c); }
		}
		else {
			{ node p = u; std::queue<node> c; eq_path(c, p); int bu = u.tag - 1, bv = v.tag - 1;
				flipBit(c, p, bu); flipTail(c, p); neTag(c, p, v.tag); flipBit(c, p, bv); neTag(c, p, u.tag);
				flipTail(c, p); flipBit(c, p, bu); flipTail(c, p); neTag(c, p, v.tag); flipBit(c, p, bv); result.push_back(c); }
		}
	}
	else if (!A.empty() && !B.empty()) {
		for (int i = 0; i <= n - 2; ++i) {
			node p = u; std::queue<node> c; eq_path(c, p);
			std::set<int> S1 = (tu == false) ? A : B, S2 = (tu == false) ? B : A;
			mod_path(p, S1, c, i + 1); flipTail(c, p);
			mod_path(p, S2, c, p.tag); neTail(c, p, tv); neTag(c, p, v.tag);
			result.push_back(c);
		}
		{ node p = u; std::queue<node> c; eq_path(c, p);
			std::set<int> S1 = (tu == false) ? A : B, S2 = (tu == false) ? B : A;
			flipTail(c, p); mod_path(p, S2, c, p.tag); flipTail(c, p);
			if (S1.count(v.tag - 1)) { for (int i : S1) { if (i + 1 != v.tag) { neTag(c, p, i + 1); flipBit(c, p, i); } } neTag(c, p, v.tag); flipBit(c, p, v.tag - 1); }
			else { mod_path(p, S1, c, v.tag); }
			neTail(c, p, tv);
			result.push_back(c); }
	}
	else {
		bool lac = (tu == false && !A.empty()) || (tu == true && !B.empty());
		std::set<int> S1 = B.empty() ? A : B;
		for (int i = 0; i <= n - 2; ++i) {
			node p = u; std::queue<node> c; eq_path(c, p);
			if (lac == false) { neTag(c, p, i + 1); flipTail(c, p); }
			mod_path(p, S1, c, i + 1); neTail(c, p, tv); neTag(c, p, v.tag);
			result.push_back(c);
		}
		{ node p = u; std::queue<node> c; eq_path(c, p);
			if (lac == true) {
				flipTail(c, p); neTag(c, p, v.tag); flipBit(c, p, v.tag - 1);
				flipTail(c, p); mod_path(p, S1, c, p.tag);
				flipTail(c, p); neTag(c, p, v.tag); flipBit(c, p, v.tag - 1); neTail(c, p, tv);
			}
			else {
				flipBit(c, p, u.tag - 1); flipTail(c, p); mod_path(p, S1, c, u.tag); flipTail(c, p);
				std::set<int> R1 = { u.tag - 1 }; mod_path(p, R1, c, v.tag); neTail(c, p, tv);
			}
			result.push_back(c); }
	}
	return result;
}


std::queue<node> BFS(node u, node v) {
	if (u == v) {
		std::queue<node> path;
		path.push(u);
		return path;
	}

	std::queue<node> q;
	std::unordered_map<node, node, NodeHash> parent;
	std::unordered_set<node, NodeHash> visited;

	q.push(u);
	visited.insert(u);
	parent.emplace(u, u);
	bool found = false;
	int n = u.dimension;

	while (!q.empty() && !found) {
		node cur = q.front();
		q.pop();

		for (int new_tag = 1; new_tag <= n - 1; ++new_tag) {
			if (new_tag == cur.tag) continue;
			node nei(n, cur.pathA, cur.pathB, new_tag, cur.tail);
			if (visited.find(nei) == visited.end() && isconnected(cur, nei)) {
				parent.emplace(nei, cur);
				if (nei == v) { found = true; break; }
				visited.insert(nei);
				q.push(nei);
			}
		}
		if (found) break;

		{
			node nei(n, cur.pathA, cur.pathB, cur.tag, !cur.tail);
			if (visited.find(nei) == visited.end() && isconnected(cur, nei)) {
				parent.emplace(nei, cur);
				if (nei == v) { found = true; break; }
				visited.insert(nei);
				q.push(nei);
			}
		}
		if (found) break;

		if (cur.tail == 0) {
			node nei(n, cur.pathA ^ (1 << (cur.tag - 1)), cur.pathB, cur.tag, cur.tail);
			if (visited.find(nei) == visited.end() && isconnected(cur, nei)) {
				parent.emplace(nei, cur);
				if (nei == v) { found = true; break; }
				visited.insert(nei);
				q.push(nei);
			}
		}
		else {
			node nei(n, cur.pathA, cur.pathB ^ (1 << (cur.tag - 1)), cur.tag, cur.tail);
			if (visited.find(nei) == visited.end() && isconnected(cur, nei)) {
				parent.emplace(nei, cur);
				if (nei == v) { found = true; break; }
				visited.insert(nei);
				q.push(nei);
			}
		}
	}

	std::queue<node> path;
	if (found) {
		std::vector<node> rev;
		node cur = v;
		while (!(cur == u)) {
			rev.push_back(cur);
			cur = parent.at(cur);
		}
		rev.push_back(u);
		for (auto it = rev.rbegin(); it != rev.rend(); ++it)
			path.push(*it);
	}
	return path;
}

std::queue<node> BiDirBFS(node u, node v) {
	if (u == v) {
		std::queue<node> path;
		path.push(u);
		return path;
	}

	int n = u.dimension;
	std::queue<node> q1, q2;
	std::unordered_set<node, NodeHash> vis1, vis2;
	std::unordered_map<node, node, NodeHash> par1, par2;

	q1.push(u); vis1.insert(u); par1.emplace(u, u);
	q2.push(v); vis2.insert(v); par2.emplace(v, v);

	auto expand = [&](std::queue<node>& q,
		std::unordered_set<node, NodeHash>& vis,
		std::unordered_map<node, node, NodeHash>& par,
		const std::unordered_set<node, NodeHash>& other_vis) -> std::optional<node> {
			int sz = q.size();
			while (sz--) {
				node cur = q.front(); q.pop();

				for (int new_tag = 1; new_tag <= n - 1; ++new_tag) {
					if (new_tag == cur.tag) continue;
					node nei(n, cur.pathA, cur.pathB, new_tag, cur.tail);
					if (!isconnected(cur, nei)) continue;
					if (vis.find(nei) != vis.end()) continue;
					if (other_vis.find(nei) != other_vis.end()) return nei;
					vis.insert(nei);
					par.emplace(nei, cur);
					q.push(nei);
				}
				{
					node nei(n, cur.pathA, cur.pathB, cur.tag, !cur.tail);
					if (isconnected(cur, nei)) {
						if (vis.find(nei) != vis.end()) continue;
						if (other_vis.find(nei) != other_vis.end()) return nei;
						vis.insert(nei);
						par.emplace(nei, cur);
						q.push(nei);
					}
				}
				if (cur.tail == 0) {
					node nei(n, cur.pathA ^ (1 << (cur.tag - 1)), cur.pathB, cur.tag, cur.tail);
					if (isconnected(cur, nei)) {
						if (vis.find(nei) != vis.end()) continue;
						if (other_vis.find(nei) != other_vis.end()) return nei;
						vis.insert(nei);
						par.emplace(nei, cur);
						q.push(nei);
					}
				}
				else {
					node nei(n, cur.pathA, cur.pathB ^ (1 << (cur.tag - 1)), cur.tag, cur.tail);
					if (isconnected(cur, nei)) {
						if (vis.find(nei) != vis.end()) continue;
						if (other_vis.find(nei) != other_vis.end()) return nei;
						vis.insert(nei);
						par.emplace(nei, cur);
						q.push(nei);
					}
				}
			}
			return std::nullopt;
		};

	std::optional<node> meet_opt;
	while (!q1.empty() && !q2.empty()) {

		if (q1.size() <= q2.size()) {
			meet_opt = expand(q1, vis1, par1, vis2);
		}
		else {
			meet_opt = expand(q2, vis2, par2, vis1);
		}
		if (meet_opt.has_value()) break;
	}

	if (!meet_opt.has_value()) return std::queue<node>();

	node meet = meet_opt.value();

	std::vector<node> path;
	node cur = meet;
	while (!(cur == u)) {
		path.push_back(cur);
		cur = par1.at(cur);
	}
	path.push_back(u);
	std::reverse(path.begin(), path.end());

	cur = par2.at(meet);
	while (!(cur == v)) {
		path.push_back(cur);
		cur = par2.at(cur);
	}
	path.push_back(v);

	std::queue<node> result;
	for (const auto& nd : path) result.push(nd);
	return result;
}

std::queue<node> AStar(node u, node v) {
	if (u == v) {
		std::queue<node> path;
		path.push(u);
		return path;
	}

	int n = u.dimension;
	auto heuristic = [&](const node& a, const node& b) -> int {
		int diff = __builtin_popcount((unsigned)(a.pathA ^ b.pathA))
			+ __builtin_popcount((unsigned)(a.pathB ^ b.pathB));
		if (a.tail != b.tail) ++diff;
		if (a.tag != b.tag) ++diff;
		return diff;
		};

	struct QueueNode {
		node n;
		int g;
		int f;
		bool operator<(const QueueNode& other) const {
			return f > other.f;
		}
	};

	std::priority_queue<QueueNode> open;
	std::map<node, int, node::ncmp> g_score;
	std::map<node, node, node::ncmp> parent;
	std::set<node, node::ncmp> closed;

	g_score.emplace(u, 0);
	parent.emplace(u, u);
	open.push({ u, 0, heuristic(u, v) });

	while (!open.empty()) {
		QueueNode curNode = open.top();
		open.pop();
		node cur = curNode.n;

		if (cur == v) {
			std::queue<node> path;
			std::vector<node> rev;
			node c = v;
			while (!(c == u)) {
				rev.push_back(c);
				c = parent.at(c);
			}
			rev.push_back(u);
			for (auto it = rev.rbegin(); it != rev.rend(); ++it)
				path.push(*it);
			return path;
		}

		closed.insert(cur);
		int cur_g = g_score.at(cur);

		auto relax = [&](const node& nei) {
			if (!isconnected(cur, nei)) return;
			int new_g = cur_g + 1;
			auto it_g = g_score.find(nei);
			if (it_g != g_score.end() && new_g >= it_g->second) return;
			if (closed.find(nei) != closed.end() && new_g >= g_score.at(nei)) return; 
			if (it_g == g_score.end()) {
				g_score.emplace(nei, new_g);
			}
			else {
				it_g->second = new_g;
			}
			auto it_p = parent.find(nei);
			if (it_p == parent.end()) {
				parent.emplace(nei, cur);
			}
			else {
				it_p->second = cur;
			}
			int h = heuristic(nei, v);
			open.push({ nei, new_g, new_g + h });
			};

		for (int new_tag = 1; new_tag <= n - 1; ++new_tag) {
			if (new_tag == cur.tag) continue;
			node nei(n, cur.pathA, cur.pathB, new_tag, cur.tail);
			relax(nei);
		}
		{
			node nei(n, cur.pathA, cur.pathB, cur.tag, !cur.tail);
			relax(nei);
		}
		if (cur.tail == 0) {
			node nei(n, cur.pathA ^ (1 << (cur.tag - 1)), cur.pathB, cur.tag, cur.tail);
			relax(nei);
		}
		else {
			node nei(n, cur.pathA, cur.pathB ^ (1 << (cur.tag - 1)), cur.tag, cur.tail);
			relax(nei);
		}
	}

	return std::queue<node>();
}

std::queue<node> Dijkstra(node u, node v) {
	if (u == v) {
		std::queue<node> path;
		path.push(u);
		return path;
	}

	int n = u.dimension;

	struct QueueNode {
		node nd;
		int dist;
		bool operator<(const QueueNode& other) const {
			return dist > other.dist; 
		}
	};

	std::priority_queue<QueueNode> pq;
	std::map<node, int, node::ncmp> best;     
	std::map<node, node, node::ncmp> parent; 
	std::set<node, node::ncmp> settled;       

	best.emplace(u, 0);
	parent.emplace(u, u);
	pq.push({ u, 0 });

	while (!pq.empty()) {
		QueueNode curNode = pq.top();
		pq.pop();
		node cur = curNode.nd;
		int curDist = curNode.dist;

		if (settled.find(cur) != settled.end()) continue;
		if (curDist != best.at(cur)) continue;

		if (cur == v) break;

		settled.insert(cur);

		auto try_relax = [&](const node& nei) {
			if (!isconnected(cur, nei)) return;
			int newDist = curDist + 1;
			auto it = best.find(nei);
			if (it == best.end() || newDist < it->second) {
				if (it == best.end())
					best.emplace(nei, newDist);
				else
					it->second = newDist;
				auto pit = parent.find(nei);
				if (pit == parent.end())
					parent.emplace(nei, cur);
				else
					pit->second = cur;
				pq.push({ nei, newDist });
			}
			};

		for (int new_tag = 1; new_tag <= n - 1; ++new_tag) {
			if (new_tag == cur.tag) continue;
			node nei(n, cur.pathA, cur.pathB, new_tag, cur.tail);
			try_relax(nei);
		}
		{
			node nei(n, cur.pathA, cur.pathB, cur.tag, !cur.tail);
			try_relax(nei);
		}
		if (cur.tail == 0) {
			node nei(n, cur.pathA ^ (1 << (cur.tag - 1)), cur.pathB, cur.tag, cur.tail);
			try_relax(nei);
		}
		else {
			node nei(n, cur.pathA, cur.pathB ^ (1 << (cur.tag - 1)), cur.tag, cur.tail);
			try_relax(nei);
		}
	}

	std::queue<node> path;
	if (best.find(v) != best.end()) {
		std::vector<node> rev;
		node cur = v;
		while (!(cur == u)) {
			rev.push_back(cur);
			cur = parent.at(cur);
		}
		rev.push_back(u);
		for (auto it = rev.rbegin(); it != rev.rend(); ++it)
			path.push(*it);
	}
	return path;
}
