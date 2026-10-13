#include <iostream>
#include <vector>
#include <set>
#include <cstdint>
#include <functional>
using namespace std;

struct DH {
	int n, m, half, N;
	vector<vector<int>> adj;
	DH(int d) : n(d), m(d - 1), half(1 << (d - 1)), N(m * half * half * 2) { build(); }
	void build() {
		adj.assign(N, {});
		for (int id = 0; id < N; ++id) {
			int tag1 = id % m, rest = id / m;
			int tail = rest % 2; rest /= 2;
			int b = rest % half, a = rest / half;
			int tag = tag1 + 1;
			for (int t = 1; t <= m; ++t) if (t != tag) adj[id].push_back(((a * half + b) * 2 + tail) * m + (t - 1));
			adj[id].push_back(((a * half + b) * 2 + (1 - tail)) * m + tag1);
			if (tail == 0) { int a2 = a ^ (1 << (tag - 1)); adj[id].push_back(((a2 * half + b) * 2 + tail) * m + tag1); }
			else { int b2 = b ^ (1 << (tag - 1)); adj[id].push_back(((a * half + b2) * 2 + tail) * m + tag1); }
		}
	}
	int enc(int a, int b, int t, int tag) const { return ((a * half + b) * 2 + t) * m + (tag - 1); }
};

static int deg(const DH& g, int u, const set<int>& S) { int c = 0; for (int v : g.adj[u]) if (S.count(v)) ++c; return c; }
static bool indist(const DH& g, const set<int>& F1, const set<int>& F2) {
	set<int> U = F1; U.insert(F2.begin(), F2.end());
	for (int x = 0; x < g.N; ++x) if (F1.count(x) != F2.count(x)) for (int y : g.adj[x]) if (!U.count(y)) return false;
	return true;
}
static bool goodA(const DH& g, const set<int>& F) {
	for (int v = 0; v < g.N; ++v) { int c = 0; for (int u : g.adj[v]) if (!F.count(u)) ++c; if (c < 1) return false; }
	return true;
}
static bool goodB(const DH& g, const set<int>& F, int gg) {
	for (int v = 0; v < g.N; ++v) if (!F.count(v)) { int c = 0; for (int u : g.adj[v]) if (!F.count(u)) ++c; if (c < gg) return false; }
	return true;
}

// n=3: full exhaustive over S of size d
static int fullExhaustive(const DH& g, int d, int gg, int defKind) {
	int best = 1 << 30;
	vector<int> idx(d);
	function<void(int, int)> rec = [&](int pos, int start) {
		if (pos == d) {
			set<int> S(idx.begin(), idx.end());
			set<int> Ns; for (int x : S) for (int y : g.adj[x]) Ns.insert(y);
			set<int> C; for (int y : Ns) if (!S.count(y)) C.insert(y);
			int csz = (int)C.size();
			for (uint64_t mA = 0; mA < (1ull << d); ++mA) {
				int na = __builtin_popcountll(mA), nb = d - na;
				int mx = csz + (na > nb ? na : nb);
				if (mx >= best) continue;
				set<int> A, B;
				for (int i = 0; i < d; ++i) { if (mA >> i & 1) A.insert(idx[i]); else B.insert(idx[i]); }
				if (defKind != 0 && (A.empty() || B.empty())) continue;
				set<int> F1 = C; F1.insert(A.begin(), A.end());
				set<int> F2 = C; F2.insert(B.begin(), B.end());
				if (F1 == F2) continue;
				if (!indist(g, F1, F2)) continue;
				if (defKind == 1 && (!goodA(g, F1) || !goodA(g, F2))) continue;
				if (defKind == 2 && (!goodB(g, F1, gg) || !goodB(g, F2, gg))) continue;
				best = mx;
			}
			return;
		}
		for (int x = start; x < g.N; ++x) { idx[pos] = x; rec(pos + 1, x + 1); }
	};
	rec(0, 0);
	return best;
}

// n>=4: exhaustive over A,B within column(0,0) (extremal is within a column)
static int columnExhaustive(const DH& g, int gg, int defKind) {
	vector<int> col;
	for (int t = 0; t < 2; ++t) for (int tag = 1; tag <= g.n - 1; ++tag) col.push_back(g.enc(0, 0, t, tag));
	int L = (int)col.size();
	int best = 1 << 30;
	for (int mA = 0; mA < (1 << L); ++mA) {
		set<int> A; for (int i = 0; i < L; ++i) if (mA >> i & 1) A.insert(col[i]);
		if (A.empty()) continue;
		int rest = ((1 << L) - 1) ^ mA;
		for (int mB = rest; mB; mB = (mB - 1) & rest) {
			set<int> B; for (int i = 0; i < L; ++i) if (mB >> i & 1) B.insert(col[i]);
			if (B.empty()) continue;
			bool ok = true;
			if (defKind == 2) { for (int x : A) if (deg(g, x, A) < gg) { ok = false; break; } if (ok) for (int y : B) if (deg(g, y, B) < gg) { ok = false; break; } }
			else { for (int x : A) if (deg(g, x, A) < 1 || deg(g, x, B) < 1) { ok = false; break; } if (ok) for (int y : B) if (deg(g, y, B) < 1 || deg(g, y, A) < 1) { ok = false; break; } }
			if (!ok) continue;
			set<int> S = A; S.insert(B.begin(), B.end());
			set<int> Ns; for (int x : S) for (int y : g.adj[x]) Ns.insert(y);
			set<int> C; for (int y : Ns) if (!S.count(y)) C.insert(y);
			set<int> F1 = C; F1.insert(A.begin(), A.end());
			set<int> F2 = C; F2.insert(B.begin(), B.end());
			int mx = max((int)F1.size(), (int)F2.size());
			if (mx >= best) continue;
			if (!indist(g, F1, F2)) continue;
			if (defKind == 1 && (!goodA(g, F1) || !goodA(g, F2))) continue;
			if (defKind == 2 && (!goodB(g, F1, gg) || !goodB(g, F2, gg))) continue;
			best = mx;
		}
	}
	return best;
}

static int constructRung(const DH& g, int gg) {
	int gs = gg + 1;
	set<int> A, B;
	for (int tag = 1; tag <= gs; ++tag) { A.insert(g.enc(0, 0, 0, tag)); B.insert(g.enc(0, 0, 1, tag)); }
	set<int> S = A; S.insert(B.begin(), B.end());
	set<int> Ns; for (int x : S) for (int y : g.adj[x]) Ns.insert(y);
	set<int> C; for (int y : Ns) if (!S.count(y)) C.insert(y);
	set<int> F1 = C; F1.insert(A.begin(), A.end());
	set<int> F2 = C; F2.insert(B.begin(), B.end());
	if (!indist(g, F1, F2) || !goodB(g, F1, gg) || !goodB(g, F2, gg)) return -1;
	return (int)F1.size();
}
static int constructSameHalf(const DH& g, int gg) {
	int gs = gg + 1;
	set<int> A, B;
	for (int tag = 1; tag <= gs; ++tag) A.insert(g.enc(0, 0, 0, tag));
	for (int tag = gs + 1; tag <= 2 * gs; ++tag) B.insert(g.enc(0, 0, 0, tag));
	set<int> S = A; S.insert(B.begin(), B.end());
	set<int> Ns; for (int x : S) for (int y : g.adj[x]) Ns.insert(y);
	set<int> C; for (int y : Ns) if (!S.count(y)) C.insert(y);
	set<int> F1 = C; F1.insert(A.begin(), A.end());
	set<int> F2 = C; F2.insert(B.begin(), B.end());
	if (!indist(g, F1, F2) || !goodB(g, F1, gg) || !goodB(g, F2, gg)) return -1;
	return (int)F1.size();
}

int main() {
	DH g3(3);
	int t0 = 1 << 30, tc = 1 << 30, t1 = 1 << 30;
	for (int d = 1; d <= 4; ++d) {
		t0 = min(t0, fullExhaustive(g3, d, 0, 0));
		tc = min(tc, fullExhaustive(g3, d, 0, 1));
		t1 = min(t1, fullExhaustive(g3, d, 1, 2));
	}
	cout << "n=3 full-exhaustive:  t0(ordinary)=" << (t0 - 1) << " (formula n=" << 3 << ")  "
		<< "tc(defA)=" << (tc - 1) << " (min(2n-1,n+4)=" << min(2 * 3 - 1, 3 + 4) << ")  "
		<< "t1(defB)=" << (t1 - 1) << " (min(2n+g-2,n+3g+1)=" << min(2 * 3 + 1 - 2, 3 + 3 * 1 + 1) << ")\n";

	cout << "--- n>=4: column-exhaustive lower bound + construction upper bound ---\n";
	for (int n = 4; n <= 7; ++n) {
		DH g(n);
		int tcE = columnExhaustive(g, 0, 1);
		int tcF = min(2 * n - 1, n + 4);
		cout << "n=" << n << " defA: exhaustive tc=" << (tcE - 1) << "  formula min(2n-1,n+4)=" << tcF << "  " << (tcE - 1 == tcF ? "OK" : "MISMATCH") << "\n";
		for (int gg = 1; gg <= n - 2; ++gg) {
			int e = columnExhaustive(g, gg, 2);
			int f = min(2 * n + gg - 2, n + 3 * gg + 1);
			int r = constructRung(g, gg);
			int s = (2 * gg + 2 <= n - 1) ? constructSameHalf(g, gg) : -1;
			bool match = (e - 1 == f);
			cout << "n=" << n << " g=" << gg << "  exhaustive t_g=" << (e - 1) << "  formula=" << f
				<< "  rung=" << r << " sameHalf=" << s << "  " << (match ? "OK" : "MISMATCH") << "\n";
		}
	}
	return 0;
}
