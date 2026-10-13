#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
using namespace std;

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
	int enc(int a, int b, int tail, int tag) const { return ((a * half + b) * 2 + tail) * m + (tag - 1); }
};

// column K_{n-1} x K_2: nodes (tail, tag); find Hamiltonian path from (tin,tagin) to (tout,tagout)
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
		vector<pair<int,int>> cand = adj(tail, tag);
		vector<pair<int,pair<int,int>>> cand2;
		for (auto& p : cand) if (!vis[id(p.first, p.second)]) {
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

int main() {
	for (int n = 4; n <= 6; ++n) {
		cout << "n=" << n << " start" << endl; cout.flush();
		DHDC g(n);
		int m = g.m, half = g.half, C = half * half;
		// binary reflected gray code over 2(n-1) bits: low n-1 = a, high n-1 = b
		vector<int> gray;
		for (int i = 0; i < C; ++i) gray.push_back(i ^ (i >> 1));
		auto bitOf = [&](int c1, int c2) { int d = c1 ^ c2; int i = 0; while (!(d & 1)) { d >>= 1; ++i; } return i; }; // bit index changed
		// port for a flip in dimension 'bit' (0..2n-3)
		auto port = [&](int bit) { if (bit < n - 1) return make_pair(0, bit + 1); else return make_pair(1, bit - (n - 1) + 1); };
		auto colA = [&](int c) { return c & (half - 1); };
		auto colB = [&](int c) { return (c >> (n - 1)) & (half - 1); };

		bool ok = true;
		vector<int> cycle;
		// for each column i, entry port = port(bit between gray[i-1], gray[i]); exit port = port(bit between gray[i], gray[i+1])
		for (int i = 0; i < C && ok; ++i) {
			int prev = gray[(i - 1 + C) % C], cur = gray[i], nxt = gray[(i + 1) % C];
			int inBit = bitOf(prev, cur), outBit = bitOf(cur, nxt);
			auto pin = port(inBit), pout = port(outBit);
			vector<pair<int,int>> cp;
			if (!colPath(m, pin.first, pin.second, pout.first, pout.second, cp)) { ok = false; cout << "  col " << i << " no path from (" << pin.first << "," << pin.second << ") to (" << pout.first << "," << pout.second << ")\n"; break; }
			int a = colA(cur), b = colB(cur);
			for (auto& p : cp) cycle.push_back(g.enc(a, b, p.first, p.second));
		}
		if (ok) {
			// verify cycle
			int N = g.N;
			vector<char> seen(N, 0);
			for (int x : cycle) { if (seen[x]) ok = false; seen[x] = 1; }
			if (ok) for (int i = 0; i + 1 < N; ++i) { bool adj = false; for (int v : g.nb(cycle[i])) if (v == cycle[i + 1]) adj = true; if (!adj) ok = false; }
			if (ok) { bool clos = false; for (int v : g.nb(cycle[N - 1])) if (v == cycle[0]) clos = true; ok = clos; }
		}
		cout << "DH_" << n << " (N=" << g.N << "): Hamiltonian cycle = " << (ok ? "YES" : "NO") << "\n";
	}
	return 0;
}
