#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>
using namespace std;

// Hamiltonian path in Q_m from X to Y (X != Y, parity(X) != parity(Y)).
static vector<int> hamPathQ(int m, int X, int Y) {
	if (m == 1) return { X, Y };
	int msb = 1 << (m - 1);
	int xh = X >> (m - 1), yh = Y >> (m - 1);
	int xl = X & (msb - 1), yl = Y & (msb - 1);
	vector<int> res;
	if (xh != yh) {
		// X in Q^{xh}, Y in Q^{yh}; parity(xl)==parity(yl).
		int Z = yl ^ 1;             // opposite parity to xl (and yl)
		vector<int> P0 = hamPathQ(m - 1, xl, Z);
		vector<int> P1 = hamPathQ(m - 1, Z, yl);
		for (int x : P0) res.push_back(x | (xh * msb));
		for (int x : P1) res.push_back(x | (yh * msb));
	} else {
		// same sub-cube Q^{xh}; detour through Q^{1-xh} at the first edge.
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

// verify path is a valid Hamiltonian path in Q_m
static bool verify(int m, const vector<int>& P, int X, int Y) {
	int N = 1 << m;
	if ((int)P.size() != N) return false;
	if (P.front() != X || P.back() != Y) return false;
	vector<char> seen(N, 0);
	for (int x : P) { if (seen[x]) return false; seen[x] = 1; }
	for (size_t i = 0; i + 1 < P.size(); ++i) {
		int d = P[i] ^ P[i + 1];
		if (d & (d - 1)) return false; // not single-bit
	}
	return true;
}

int main() {
	for (int m = 2; m <= 6; ++m) {
		int N = 1 << m;
		int pass = 0, total = 0, fail = 0;
		for (int X = 0; X < N; ++X) for (int Y = 0; Y < N; ++Y) {
			if (X == Y) continue;
			if ((__builtin_popcount(X) ^ __builtin_popcount(Y)) & 1) { // opposite parity
				++total;
				vector<int> P = hamPathQ(m, X, Y);
				if (verify(m, P, X, Y)) ++pass; else { ++fail; if (fail <= 2) cout << "  m=" << m << " X=" << X << " Y=" << Y << " FAIL\n"; }
			}
		}
		cout << "Q_" << m << ": opposite-parity pairs " << pass << "/" << total << " (fail=" << fail << ")\n";
	}
	return 0;
}
