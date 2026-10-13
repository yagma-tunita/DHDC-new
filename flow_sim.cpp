#include "flow_sim.hpp"

#include <algorithm>
#include <unordered_map>
#include <limits>
#include <cmath>
#include <queue>
#include <functional>

namespace flowsim
{

namespace
{

struct Link
{
	double cap;
	std::vector<int> flows;
};

// 精确 max-min 公平带宽分配（水填充）
// rate[f] 为输出，只对 active[f]==1 的流赋值
void maxMinFair(const std::vector<Link>& links, std::vector<double>& rate,
	const std::vector<char>& active, const std::vector<std::vector<int>>& flowLinks)
{
	int F = (int)rate.size();
	int L = (int)links.size();

	std::vector<char> satFlow(F, 0);
	std::vector<int> cnt(L, 0);
	int remaining = 0;
	for (int f = 0; f < F; ++f)
	{
		rate[f] = 0.0;
		if (!active[f]) continue;
		++remaining;
		for (int lid : flowLinks[f]) ++cnt[lid];
	}
	if (remaining == 0) return;

	using P = std::pair<double, int>;
	std::priority_queue<P, std::vector<P>, std::greater<P>> heap;
	for (int l = 0; l < L; ++l)
		if (cnt[l] > 0) heap.push({ links[l].cap / cnt[l], l });

	while (remaining > 0 && !heap.empty())
	{
		P top = heap.top();
		heap.pop();
		int l = top.second;
		if (cnt[l] <= 0) continue;
		double cur = links[l].cap / cnt[l];
		if (top.first != cur) continue;   // 过期条目（cnt 只减不增，fair 单调上升）

		for (int f : links[l].flows)
		{
			if (satFlow[f]) continue;
			rate[f] = cur;
			satFlow[f] = 1;
			--remaining;
			for (int lid2 : flowLinks[f])
			{
				if (lid2 == l) continue;
				if (cnt[lid2] > 0)
				{
					--cnt[lid2];
					if (cnt[lid2] > 0) heap.push({ links[lid2].cap / cnt[lid2], lid2 });
				}
			}
		}
		cnt[l] = 0;
	}
}

} // namespace

Result run(const Config& cfg, std::vector<Flow> flows)
{
	int F = (int)flows.size();
	Result res;
	res.maxThp = 0.0;
	res.avgThp = 0.0;
	res.totalFCT = 0.0;
	res.avgHops = 0.0;

	if (F == 0) return res;

	// 1. 路由所有流，建立边 -> 链路映射
	std::unordered_map<long long, int> edge2link;
	std::vector<Link> links;
	std::vector<std::vector<int>> flowLinks(F);
	double totalHops = 0.0;

	for (int f = 0; f < F; ++f)
	{
		flows[f].id = f;
		flows[f].path = cfg.route(flows[f].src, flows[f].dst);
		int hops = (int)flows[f].path.size() - 1;
		if (hops < 0) hops = 0;
		totalHops += (double)hops;
		for (size_t i = 0; i + 1 < flows[f].path.size(); ++i)
		{
			int a = flows[f].path[i];
			int b = flows[f].path[i + 1];
			if (a == b) continue;
			int lo = a < b ? a : b;
			int hi = a < b ? b : a;
			long long key = ((long long)lo << 32) | (unsigned int)hi;
			int lid;
			auto it = edge2link.find(key);
			if (it == edge2link.end())
			{
				lid = (int)links.size();
				Link lk;
				lk.cap = cfg.linkCap;
				links.push_back(lk);
				edge2link[key] = lid;
			}
			else
			{
				lid = it->second;
			}
			flowLinks[f].push_back(lid);
		}
	}
	res.avgHops = totalHops / (double)F;

	// 2. 按到达时刻排序
	std::vector<Flow*> order;
	order.reserve(F);
	for (int f = 0; f < F; ++f)
		if (flows[f].path.size() >= 2) order.push_back(&flows[f]);
	if (order.empty()) return res;

	std::sort(order.begin(), order.end(),
		[](const Flow* a, const Flow* b) { return a->tArr < b->tArr; });

	// 3. 事件循环
	std::vector<double> rate(F, 0.0);
	std::vector<char> active(F, 0);

	double now = 0.0;
	size_t idx = 0;
	int activeCount = 0;
	double totalBytes = 0.0;
	double minArr = order.front()->tArr;
	double maxFinish = minArr;
	const double epsAbs = 1e-6;    // 忽略 < 1e-6 比特的传输残留
	const double epsTime = 1e-12;  // 忽略 < 1ps 的剩余传输时间
	const double epsArr = 1e-9;    // 到达时刻比较容差

	while (idx < order.size() || activeCount > 0)
	{
		double activeRate = 0.0;
		double tComp = std::numeric_limits<double>::infinity();
		for (int f = 0; f < F; ++f)
		{
			if (!active[f]) continue;
			activeRate += rate[f];
			if (rate[f] > 0.0)
			{
				double rt = flows[f].remaining / rate[f];
				if (rt < tComp) tComp = rt;
			}
		}

		double tArr = (idx < order.size()) ? order[idx]->tArr : std::numeric_limits<double>::infinity();
		double tArrGap = (tArr == std::numeric_limits<double>::infinity()) ? std::numeric_limits<double>::infinity() : (tArr - now);
		double dt = (tComp < tArrGap) ? tComp : tArrGap;
		if (!(dt > 0.0)) dt = 0.0;
		if (dt == std::numeric_limits<double>::infinity()) break;

		for (int f = 0; f < F; ++f)
			if (active[f]) flows[f].remaining -= rate[f] * dt;

		totalBytes += activeRate * dt;
		if (activeRate > res.maxThp) res.maxThp = activeRate;
		res.thpCurve.push_back({ now, activeRate });
		now += dt;

		for (int f = 0; f < F; ++f)
		{
			if (!active[f]) continue;
			bool done = flows[f].remaining <= epsAbs;
			if (!done && rate[f] > 0.0 && flows[f].remaining / rate[f] <= epsTime) done = true;
			if (!done) continue;
			flows[f].tFinish = now;
			if (now > maxFinish) maxFinish = now;
			active[f] = 0;
			--activeCount;
			for (int lid : flowLinks[f])
			{
				auto& vec = links[lid].flows;
				vec.erase(std::remove(vec.begin(), vec.end(), f), vec.end());
			}
		}

		while (idx < order.size() && order[idx]->tArr <= now + epsArr)
		{
			int f = order[idx]->id;
			if (!active[f])
			{
				active[f] = 1;
				++activeCount;
				flows[f].remaining = flows[f].size;
				for (int lid : flowLinks[f]) links[lid].flows.push_back(f);
			}
			++idx;
		}

		maxMinFair(links, rate, active, flowLinks);
	}

	double totalTime = maxFinish - minArr;
	if (totalTime > 0.0) res.avgThp = totalBytes / totalTime;
	res.totalFCT = totalTime;

	return res;
}

} // namespace flowsim
