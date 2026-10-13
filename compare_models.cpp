#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
#include <random>
#include <queue>
using namespace std;

struct DHDC {
	int n, m, half, N;
	DHDC(int d) : n(d), m(d - 1), half(1 << (d - 1)), N(m * half * half * 2) {}
	vector<int> nb(int id) const {
		int tag1 = id % m, rest = id / m, tail = rest % 2; rest /= 2;
		int b = rest % half, a = rest / half, tag = tag1 + 1;
		vector<int> o;
		for (int t = 1; t <= m; ++t) if (t != tag) o.push_back(((a * half + b) * 2 + tail) * m + (t - 1));
		o.push_back(((a * half + b) * 2 + (1 - tail)) * m + tag1);
		if (tail == 0) { int a2 = a ^ (1 << (tag - 1)); o.push_back(((a2 * half + b) * 2 + tail) * m + tag1); }
		else { int b2 = b ^ (1 << (tag - 1)); o.push_back(((a * half + b2) * 2 + tail) * m + tag1); }
		return o;
	}
};
struct HSDC {
	int n, N;
	HSDC(int d) : n(d), N(n * (1 << n)) {}
	vector<int> nb(int id) const {
		int y = id % n, u = id / n;
		vector<int> o;
		for (int yy = 0; yy < n; ++yy) if (yy != y) o.push_back(u * n + yy);
		o.push_back((u ^ (1 << y)) * n + y);
		return o;
	}
};
struct CSDC {
	int n, N;
	CSDC(int d) : n(d), N(n * (1 << n)) {}
	int cqNeighbor(int u, int x) const { int v = u ^ (1 << (x - 1)); for (int j = 0; j < (x - 1) / 2; ++j) if (u & (1 << (2 * j))) v ^= (1 << (2 * j + 1)); return v; }
	vector<int> nb(int id) const {
		int x = id % n, u = id / n;
		vector<int> o;
		for (int xx = 0; xx < n; ++xx) if (xx != x) o.push_back(u * n + xx);
		o.push_back(cqNeighbor(u, x + 1) * n + x);
		return o;
	}
};
struct CrossedCube {
	int n; vector<vector<int>> adj;
	CrossedCube(int d) : n(d) { build(); }
	static bool pairRel(int x0, int x1, int y0, int y1) { return x0 == y0 && ((x1 ^ y1) == x0); }
	static bool crossOk(int k, int u, int v) {
		if (k % 2 == 0) if (((u >> (k - 2)) & 1) != ((v >> (k - 2)) & 1)) return false;
		int lim = (k - 1) / 2;
		for (int i = 0; i < lim; ++i) { int x0 = (u >> (2 * i)) & 1, x1 = (u >> (2 * i + 1)) & 1, y0 = (v >> (2 * i)) & 1, y1 = (v >> (2 * i + 1)) & 1; if (!pairRel(x0, x1, y0, y1)) return false; }
		return true;
	}
	void build() {
		adj.assign(1 << n, {});
		function<void(int, vector<vector<int>>&)> rec = [&](int k, vector<vector<int>>& a) {
			if (k == 1) { a = { {1}, {0} }; return; }
			int h = 1 << (k - 1); vector<vector<int>> sub; rec(k - 1, sub); a.assign(2 * h, {});
			for (int u = 0; u < h; ++u) { for (int v : sub[u]) { a[u].push_back(v); a[u + h].push_back(v + h); } }
			for (int u = 0; u < h; ++u) for (int v = 0; v < h; ++v) if (crossOk(k, u, v)) { a[u].push_back(v + h); a[v + h].push_back(u); }
		};
		rec(n, adj);
	}
};
struct BCDC {
	int n, N; vector<vector<int>> adj; vector<pair<int,int>> edges;
	BCDC(int d) : n(d) { CrossedCube cq(d); for (int u = 0; u < (1 << n); ++u) for (int v : cq.adj[u]) if (u < v) edges.push_back({ u, v }); N = edges.size(); adj.assign(N, {}); for (int i = 0; i < N; ++i) for (int j = i + 1; j < N; ++j) { auto a = edges[i], b = edges[j]; if (a.first == b.first || a.first == b.second || a.second == b.first || a.second == b.second) { adj[i].push_back(j); adj[j].push_back(i); } } }
};

struct G { int N; vector<vector<int>> adj; G(int n_, vector<vector<int>> a) : N(n_), adj(move(a)) {} };

static bool bipartite(const G& g, int& cntOdd) {
	vector<int> c(g.N, -1); bool bip = true; cntOdd = 0;
	for (int s = 0; s < g.N; ++s) {
		if (c[s] >= 0) continue; c[s] = 0; queue<int> q; q.push(s);
		while (!q.empty()) { int u = q.front(); q.pop(); for (int v : g.adj[u]) { if (c[v] < 0) { c[v] = c[u] ^ 1; q.push(v); } else if (c[v] == c[u]) bip = false; } }
	}
	// count triangles
	for (int u = 0; u < g.N; ++u) for (int v : g.adj[u]) if (v > u) for (int w : g.adj[v]) if (w > v) for (int x : g.adj[u]) if (x == w) ++cntOdd;
	return bip;
}

static int connectivity(const G& g, int maxk) {
	for (int k = 1; k <= maxk; ++k) {
		int N = g.N; vector<int> idx(k); bool found = false;
		function<void(int,int)> rec = [&](int pos, int start) {
			if (found) return;
			if (pos == k) {
				vector<char> rem(N, 0); for (int x : idx) rem[x] = 1;
				int root = -1; for (int i = 0; i < N; ++i) if (!rem[i]) { root = i; break; }
				if (root < 0) { found = true; return; }
				vector<char> vis(N, 0); vector<int> st{ root }; vis[root] = 1;
				while (!st.empty()) { int c = st.back(); st.pop_back(); for (int v : g.adj[c]) if (!rem[v] && !vis[v]) { vis[v] = 1; st.push_back(v); } }
				for (int i = 0; i < N; ++i) if (!rem[i] && !vis[i]) { found = true; return; }
				return;
			}
			for (int x = start; x < N; ++x) { idx[pos] = x; rec(pos + 1, x + 1); }
		};
		rec(0, 0);
		if (found) return k;
	}
	return maxk;
}

static bool hamCycle(const G& g) {
	int N = g.N; vector<char> vis(N, 0); vector<int> path; bool ok = false;
	function<void(int,int)> dfs = [&](int u, int cnt) {
		if (ok) return;
		if (cnt == N) { for (int v : g.adj[u]) if (v == path[0]) { ok = true; return; } return; }
		vector<pair<int,int>> c; for (int v : g.adj[u]) if (!vis[v]) { int d = 0; for (int w : g.adj[v]) if (!vis[w]) ++d; c.push_back({ d, v }); }
		sort(c.begin(), c.end());
		for (auto& p : c) { int v = p.second; vis[v] = 1; path.push_back(v); dfs(v, cnt + 1); path.pop_back(); vis[v] = 0; if (ok) return; }
	};
	vis[0] = 1; path.push_back(0); dfs(0, 1); return ok;
}

static bool hamPath(const G& g, int s, int t) {
	int N = g.N; vector<char> vis(N, 0); vector<int> path; bool ok = false;
	function<void(int,int)> dfs = [&](int u, int cnt) {
		if (ok) return;
		if (cnt == N) { if (u == t) ok = true; return; }
		vector<pair<int,int>> c; for (int v : g.adj[u]) if (!vis[v]) { int d = 0; for (int w : g.adj[v]) if (!vis[w]) ++d; c.push_back({ d, v }); }
		sort(c.begin(), c.end());
		for (auto& p : c) { int v = p.second; vis[v] = 1; path.push_back(v); dfs(v, cnt + 1); path.pop_back(); vis[v] = 0; if (ok) return; }
	};
	vis[s] = 1; path.push_back(s); dfs(s, 1); return ok;
}

static int hamConnSample(const G& g, int samples) {
	mt19937 rng(7); int pass = 0, total = 0;
	while (total < samples) { int s = rng() % g.N, t = rng() % g.N; if (s == t) continue; ++total; if (hamPath(g, s, t)) ++pass; }
	return pass * 100 / total;
}

static G buildDHDC(int n) { DHDC d(n); vector<vector<int>> a(d.N); for (int i = 0; i < d.N; ++i) a[i] = d.nb(i); return G(d.N, move(a)); }
static G buildHSDC(int n) { HSDC h(n); vector<vector<int>> a(h.N); for (int i = 0; i < h.N; ++i) a[i] = h.nb(i); return G(h.N, move(a)); }
static G buildCSDC(int n) { CSDC c(n); vector<vector<int>> a(c.N); for (int i = 0; i < c.N; ++i) a[i] = c.nb(i); return G(c.N, move(a)); }
static G buildBCDC(int n) { BCDC b(n); return G(b.N, b.adj); }

int main() {
	for (int n = 3; n <= 4; ++n) {
		cout << "===== n=" << n << " =====\n";
		struct Row { const char* name; G g; int kappa; bool skipHam; };
		vector<Row> rows;
		rows.push_back({ "DHDC", buildDHDC(n), n, true });
		rows.push_back({ "HSDC", buildHSDC(n), n, false });
		rows.push_back({ "BCDC", buildBCDC(n), 2 * n - 2, false });
		rows.push_back({ "CSDC", buildCSDC(n), n, false });
		for (auto& r : rows) {
			const G& g = r.g;
			int deg = g.adj[0].size();
			int odd = 0; bool bip = bipartite(g, odd);
			string kap;
			if (g.N <= 40) kap = to_string(connectivity(g, r.kappa));
			else kap = to_string(r.kappa) + "(known)";
			string hc, hp;
			if (r.skipHam || g.N > 64) { hc = "known"; hp = "see note"; }
			else { hc = hamCycle(g) ? "Y" : "N"; hp = to_string(hamConnSample(g, 40)) + "%"; }
			cout << "  " << r.name << "  N=" << g.N << "  deg=" << deg << "  kappa=" << kap
				<< "  bipartite=" << (bip ? "Y" : "N") << "  hamCycle=" << hc << "  hamConn=" << hp << "\n";
		}
		cout.flush();
	}
	return 0;
}
