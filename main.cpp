#include "node.hpp"
#include "AlgoRouting.hpp"
#include "flow_sim.hpp"

#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <queue>
#include <set>
#include <cstdlib>

using namespace flowsim;

// DHDC ������������������ id <-> node ˫��ӳ�� + ·��
struct DhdicTopology
{
	int n;
	int N;
	int m;
	long long half;

	DhdicTopology(int dim) : n(dim)
	{
		m = n - 1;
		half = 1LL << (n - 1);
		N = (int)((long long)m * half * half * 2);
	}

	node idToNode(int id) const
	{
		long long b = id % half;
		long long t = id / half;
		long long a = t % half;
		long long r = t / half;
		int tail = (int)(r / m);
		int tag = (int)(r % m) + 1;
		int pathA = 0, pathB = 0;
		for (int i = 0; i < n - 1; ++i)
		{
			pathA |= (int)((a >> i) & 1) << (n - 2 - i);
			pathB |= (int)((b >> i) & 1) << (n - 2 - i);
		}
		return node(n, pathA, pathB, tag, tail != 0);
	}

	int nodeToId(const node& nd) const
	{
		long long a = 0, b = 0;
		for (int i = 0; i < n - 1; ++i)
		{
			a = (a << 1) | ((nd.pathA >> i) & 1);
			b = (b << 1) | ((nd.pathB >> i) & 1);
		}
		long long tail = nd.tail ? 1 : 0;
		long long tag = nd.tag - 1;
		return (int)(b + half * (a + half * (tag + m * tail)));
	}

	std::vector<int> route(int src, int dst) const
	{
		node u = idToNode(src);
		node v = idToNode(dst);
		std::queue<node> q = DHRouting(u, v);
		std::vector<int> path;
		int prev = -1;
		while (!q.empty())
		{
			int id = nodeToId(q.front());
			q.pop();
			if (id != prev) path.push_back(id);
			prev = id;
		}
		return path;
	}
};

// ���������log-uniform ����С + ���ɵ���
std::vector<Flow> genRandom(int N, int numFlows, double minBytes, double maxBytes,
	double duration, std::mt19937& rng)
{
	std::vector<Flow> flows;
	flows.reserve(numFlows);
	std::uniform_int_distribution<int> uid(0, N - 1);
	std::uniform_real_distribution<double> u01(0.0, 1.0);
	std::exponential_distribution<double> interArr((double)numFlows / duration);
	double t = 0.0;
	double totalBytes = 0.0;
	for (int i = 0; i < numFlows; ++i)
	{
		Flow f;
		f.id = i;
		f.src = uid(rng);
		do { f.dst = uid(rng); } while (f.dst == f.src);
		double s = std::exp(std::log(minBytes) + u01(rng) * (std::log(maxBytes) - std::log(minBytes)));
		f.size = s * 8.0;
		f.tArr = t;
		t += interArr(rng);
		totalBytes += f.size;
		flows.push_back(f);
	}
	std::cout << "[random] flows=" << numFlows << " total=" << totalBytes / 8.0 / 1e12 << " TB\n";
	return flows;
}

// all-to-all�����ܶȲ�����ͬ��������
std::vector<Flow> genAllToAll(int N, int numFlows, double minBytes, double maxBytes,
	double duration, std::mt19937& rng)
{
	std::vector<Flow> flows;
	flows.reserve(numFlows);
	std::uniform_int_distribution<int> uid(0, N - 1);
	std::uniform_real_distribution<double> u01(0.0, 1.0);
	std::exponential_distribution<double> interArr((double)numFlows / duration);
	std::set<long long> seen;
	double t = 0.0;
	double totalBytes = 0.0;
	for (int i = 0; i < numFlows; ++i)
	{
		Flow f;
		f.id = i;
		long long key;
		do
		{
			f.src = uid(rng);
			f.dst = uid(rng);
			int lo = f.src < f.dst ? f.src : f.dst;
			int hi = f.src < f.dst ? f.dst : f.src;
			key = ((long long)lo << 32) | (unsigned int)hi;
		} while (f.src == f.dst || seen.count(key));
		seen.insert(key);
		double s = std::exp(std::log(minBytes) + u01(rng) * (std::log(maxBytes) - std::log(minBytes)));
		f.size = s * 8.0;
		f.tArr = t;
		t += interArr(rng);
		totalBytes += f.size;
		flows.push_back(f);
	}
	std::cout << "[alltoall] flows=" << numFlows << " total=" << totalBytes / 8.0 / 1e12 << " TB\n";
	return flows;
}

void report(const char* name, const Result& r)
{
	std::cout << "=== " << name << " ===\n";
	std::cout << "  maxThp   = " << r.maxThp / 1e9 << " Gbps\n";
	std::cout << "  avgThp   = " << r.avgThp / 1e9 << " Gbps\n";
	std::cout << "  totalFCT = " << r.totalFCT << " s\n";
	std::cout << "  avgHops  = " << r.avgHops << "\n";
}

void dumpCurve(const char* path, const Result& r, int buckets = 200)
{
	if (r.thpCurve.empty()) return;
	double t0 = r.thpCurve.front().first;
	double t1 = r.thpCurve.back().first;
	double span = t1 - t0;
	if (span <= 0.0) return;
	std::vector<double> sum(buckets, 0.0);
	std::vector<long long> cnt(buckets, 0);
	for (auto& p : r.thpCurve)
	{
		int b = (int)((p.first - t0) / span * buckets);
		if (b >= buckets) b = buckets - 1;
		if (b < 0) b = 0;
		sum[b] += p.second;
		cnt[b]++;
	}
	std::ofstream of(path);
	of << "time,thpGbps\n";
	for (int b = 0; b < buckets; ++b)
	{
		double avg = cnt[b] ? sum[b] / cnt[b] : 0.0;
		double tt = t0 + span * (b + 0.5) / buckets;
		of << tt << "," << (avg / 1e9) << "\n";
	}
	std::cout << "curve -> " << path << "\n";
}

int main(int argc, char** argv)
{
	int dim = 6;
	int numFlows = 80000;
	if (argc > 1) dim = std::atoi(argv[1]);
	if (argc > 2) numFlows = std::atoi(argv[2]);

	DhdicTopology topo(dim);
	std::cout << "DHDC n=" << dim << ", servers=" << topo.N << "\n";

	Config cfg;
	cfg.N = topo.N;
	cfg.linkCap = 1e9;
	cfg.route = [&](int s, int t) { return topo.route(s, t); };

	std::mt19937 rng(12345u);

	{
		auto flows = genRandom(topo.N, numFlows, 1e3, 1e9, 100.0, rng);
		Result r = run(cfg, std::move(flows));
		report("random", r);
		dumpCurve("thp_random.csv", r);
	}

	{
		int n2 = numFlows / 8 > 1000 ? numFlows / 8 : 1000;
		auto flows = genAllToAll(topo.N, n2, 1e3, 1e9, 50.0, rng);
		Result r = run(cfg, std::move(flows));
		report("all-to-all", r);
		dumpCurve("thp_alltoall.csv", r);
	}

	return 0;
}
