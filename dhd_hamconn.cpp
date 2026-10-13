#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <random>
using namespace std;

// ---------- hypercube Hamiltonian path (opposite parity) ----------
static vector<int> hamPathQ(int m, int X, int Y) {
	if (m == 1) return { X, Y };
	int msb = 1 << (m - 1);
	int xh = X >> (m - 1), yh = Y >> (m - 1);
	int xl = X & (msb - 1), yl = Y & (msb - 1);
	vector<int> res;
	if (xh != yh) {
		int Z = yl ^ 1;
		vector<int> P0 = hamPathQ(m - 1, xl, Z);
		vector<int> P1 = hamPathQ(m - 1, Z, yl);
		for (int x : P0) res.push_back(x | (xh * msb));
		for (int x : P1) res.push_back(x | (yh * msb));
	} else {
		vector<int> P = hamPathQ(m - 1, xl, yl);
		int U = P[0], V = P[1];
		vector<int> D = hamPathQ(m - 1, U, V);
		int om = (1 - xh) * msb;
		res.push_back(U | (xh * msb));
		for (int x : D) res.push_back(x | om);
		for (size_t i = 1; i < P.size(); ++i) res.push_back(P[i] | (xh * msb));
	}
	return res;
}

// ---------- DHDC ----------
struct DHDC {
	int n, m, half, N;
	DHDC(int d) : n(d), m(d - 1), half(1 << (d - 1)), N(m * half * half * 2) {}
	vector<int> nb(int id) const {
		int t = id % m, r = id / m, tail = r % 2; r /= 2;
		int b = r % half, a = r / half, tag = t + 1;
		vector<int> o;
		for (int x = 1; x <= m; ++x) if (x != tag) o.push_back(((a * half + b) * 2 + tail) * m + (x - 1));
		o.push_back(((a * half + b) * 2 + (1 - tail)) * m + t);
		if (tail == 0) { int a2 = a ^ (1 << (tag - 1)); o.push_back(((a2 * half + b) * 2 + tail) * m + t); }
		else { int b2 = b ^ (1 << (tag - 1)); o.push_back(((a * half + b2) * 2 + tail) * m + t); }
		return o;
	}
	void dec(int id, int& a, int& b, int& tail, int& tag) const {
		int t = id % m, r = id / m; tail = r % 2; r /= 2; b = r % half; a = r / half; tag = t + 1;
	}
	int enc(int a, int b, int tail, int tag) const { return ((a * half + b) * 2 + tail) * m + (tag - 1); }
	int col(int id) const { int t = id % m, r = id / m; r /= 2; int b = r % half, a = r / half; return a | (b << (n - 1)); }
};

// column K_{n-1} x K_2 Hamiltonian path from (tin,tagin) to (tout,tagout)
static bool colPath(int m, int tin, int tagin, int tout, int tagout, vector<pair<int,int>>& out) {
	int N = 2 * m;
	auto id = [&](int tail, int tag) { return tail * m + (tag - 1); };
	auto adj = [&](int tail, int tag) {
		vector<pair<int,int>> r;
		for (int t = 1; t <= m; ++t) if (t != tag) r.push_back({ tail, t });
		r.push_back({ 1 - tail, tag });
		return r;
	};
	vector<char> vis(N, 0);
	vector<pair<int,int>> path;
	bool ok = false;
	function<void(int,int,int)> dfs = [&](int tail, int tag, int cnt) {
		if (ok) return;
		if (cnt == N) { if (tail == tout && tag == tagout) ok = true; return; }
		vector<pair<int,pair<int,int>>> cand2;
		for (auto& p : adj(tail, tag)) if (!vis[id(p.first, p.second)]) {
			int d = 0; for (auto& q : adj(p.first, p.second)) if (!vis[id(q.first, q.second)]) ++d;
			cand2.push_back({ d, p });
		}
		sort(cand2.begin(), cand2.end());
		for (auto& c : cand2) { auto p = c.second; vis[id(p.first, p.second)] = 1; path.push_back(p); dfs(p.first, p.second, cnt + 1); if (ok) return; path.pop_back(); vis[id(p.first, p.second)] = 0; }
	};
	vis[id(tin, tagin)] = 1; path.push_back({ tin, tagin }); dfs(tin, tagin, 1);
	if (ok) out = path;
	return ok;
}

// DH_n Hamiltonian path u -> v (opposite column parity only)
static bool hamPathDH(DHDC& g, int u, int v, vector<int>& out) {
	int n = g.n, m = g.m, half = g.half;
	int cu = g.col(u), cv = g.col(v);
	if ((__builtin_popcount(cu) ^ __builtin_popcount(cv)) & 1) return false; // same parity not handled here
	vector<int> cols = hamPathQ(2 * (n - 1), cu, cv);
	int C = cols.size();
	auto bitOf = [&](int c1, int c2) { int d = c1 ^ c2, i = 0; while (!(d & 1)) { d >>= 1; ++i; } return i; };
	auto port = [&](int bit) { if (bit < n - 1) return make_pair(0, bit + 1); else return make_pair(1, bit - (n - 1) + 1); };
	int ut = 0, utag = 0, vt = 0, vtag = 0;
	g.dec(u, *(new int()), *(new int()), ut, utag);
	g.dec(v, *(new int()), *(new int()), vt, vtag);
	out.clear();
	for (int i = 0; i < C; ++i) {
		int a = cols[i] & (half - 1), b = (cols[i] >> (n - 1)) & (half - 1);
		int inBit = (i == 0) ? -1 : bitOf(cols[i - 1], cols[i]);
		int outBit = (i == C - 1) ? -1 : bitOf(cols[i], cols[i + 1]);
		auto pin = (i == 0) ? make_pair(ut, utag) : port(inBit);
		auto pout = (i == C - 1) ? make_pair(vt, vtag) : port(outBit);
		vector<pair<int,int>> cp;
		if (!colPath(m, pin.first, pin.second, pout.first, pout.second, cp)) return false;
		for (auto& p : cp) out.push_back(g.enc(a, b, p.first, p.second));
	}
	return true;
}

static bool verify(DHDC& g, int u, int v, const vector<int>& P) {
	int N = g.N;
	if ((int)P.size() != N) { cout << "  size " << P.size() << " != " << N << "\n"; return false; }
	if (P.front() != u || P.back() != v) { cout << "  endpoint mismatch\n"; return false; }
	vector<char> seen(N, 0);
	for (int x : P) { if (seen[x]) { cout << "  repeat\n"; return false; } seen[x] = 1; }
	for (size_t i = 0; i + 1 < P.size(); ++i) { bool adj = false; for (int w : g.nb(P[i])) if (w == P[i + 1]) adj = true; if (!adj) { cout << "  edge break at " << i << "\n"; return false; } }
	return true;
}

int main() {
	for (int n = 4; n <= 5; ++n) {
		DHDC g(n);
		mt19937 rng(5);
		int pass = 0, total = 0;
		for (int k = 0; k < 8; ++k) {
			int u = rng() % g.N, v = rng() % g.N; if (u == v) { --k; continue; }
			int cu = g.col(u), cv = g.col(v);
			if (((__builtin_popcount(cu) ^ __builtin_popcount(cv)) & 1)) { --k; continue; } // need opposite parity columns
			++total;
			vector<int> P;
			if (hamPathDH(g, u, v, P) && verify(g, u, v, P)) ++pass;
			else cout << "  n=" << n << " u=" << u << " v=" << v << " FAIL\n";
		}
		cout << "DH_" << n << " (N=" << g.N << "): opposite-parity Hamiltonian paths " << pass << "/" << total << "\n";
	}
	return 0;
}
