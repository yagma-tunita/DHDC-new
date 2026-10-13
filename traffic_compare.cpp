#include "flow_sim.hpp"
#include <iostream>
#include <vector>
#include <functional>
#include <queue>
#include <unordered_map>
#include <random>
#include <algorithm>
using namespace std;
using namespace flowsim;

// ---------- four models as integer adjacency ----------
struct DHDC { int n,m,half,N; DHDC(int d):n(d),m(d-1),half(1<<(d-1)),N(m*half*half*2){}
	vector<int> nb(int id)const{int t=id%m,r=id/m,tail=r%2;r/=2;int b=r%half,a=r/half,tag=t+1;vector<int>o;
	for(int x=1;x<=m;++x)if(x!=tag)o.push_back(((a*half+b)*2+tail)*m+(x-1));
	o.push_back(((a*half+b)*2+(1-tail))*m+t); if(tail==0){int a2=a^(1<<(tag-1));o.push_back(((a2*half+b)*2+tail)*m+t);}else{int b2=b^(1<<(tag-1));o.push_back(((a*half+b2)*2+tail)*m+t);}return o;} };
struct HSDC { int n,N; HSDC(int d):n(d),N(n*(1<<n)){}
	vector<int> nb(int id)const{int y=id%n,u=id/n;vector<int>o;for(int yy=0;yy<n;++yy)if(yy!=y)o.push_back(u*n+yy);o.push_back((u^(1<<y))*n+y);return o;} };
struct CSDC { int n,N; CSDC(int d):n(d),N(n*(1<<n)){}
	int cq(int u,int x)const{int v=u^(1<<(x-1));for(int j=0;j<(x-1)/2;++j)if(u&(1<<(2*j)))v^=(1<<(2*j+1));return v;}
	vector<int> nb(int id)const{int x=id%n,u=id/n;vector<int>o;for(int xx=0;xx<n;++xx)if(xx!=x)o.push_back(u*n+xx);o.push_back(cq(u,x+1)*n+x);return o;} };
struct CrossedCube { int n; vector<vector<int>> adj; CrossedCube(int d):n(d){build();}
	static bool pr(int x0,int x1,int y0,int y1){return x0==y0&&((x1^y1)==x0);}
	static bool cok(int k,int u,int v){if(k%2==0)if(((u>>(k-2))&1)!=((v>>(k-2))&1))return false;int lim=(k-1)/2;for(int i=0;i<lim;++i){int x0=(u>>(2*i))&1,x1=(u>>(2*i+1))&1,y0=(v>>(2*i))&1,y1=(v>>(2*i+1))&1;if(!pr(x0,x1,y0,y1))return false;}return true;}
	void build(){adj.assign(1<<n,{});function<void(int,vector<vector<int>>&)>rec=[&](int k,vector<vector<int>>&a){if(k==1){a={{1},{0}};return;}int h=1<<(k-1);vector<vector<int>>s;rec(k-1,s);a.assign(2*h,{});for(int u=0;u<h;++u){for(int v:s[u]){a[u].push_back(v);a[u+h].push_back(v+h);}}for(int u=0;u<h;++u)for(int v=0;v<h;++v)if(cok(k,u,v)){a[u].push_back(v+h);a[v+h].push_back(u);}};rec(n,adj);} };
struct BCDC { int n,N; vector<vector<int>> adj; BCDC(int d):n(d){CrossedCube cq(d);vector<pair<int,int>>e;for(int u=0;u<(1<<n);++u)for(int v:cq.adj[u])if(u<v)e.push_back({u,v});N=e.size();adj.assign(N,{});for(int i=0;i<N;++i)for(int j=i+1;j<N;++j){auto a=e[i],b=e[j];if(a.first==b.first||a.first==b.second||a.second==b.first||a.second==b.second){adj[i].push_back(j);adj[j].push_back(i);}}} };

struct G { int N; vector<vector<int>> adj; G(int n_,vector<vector<int>>a):N(n_),adj(move(a)){} };

static G buildDHDC(int n){DHDC d(n);vector<vector<int>>a(d.N);for(int i=0;i<d.N;++i)a[i]=d.nb(i);return G(d.N,move(a));}
static G buildHSDC(int n){HSDC h(n);vector<vector<int>>a(h.N);for(int i=0;i<h.N;++i)a[i]=h.nb(i);return G(h.N,move(a));}
static G buildCSDC(int n){CSDC c(n);vector<vector<int>>a(c.N);for(int i=0;i<c.N;++i)a[i]=c.nb(i);return G(c.N,move(a));}
static G buildBCDC(int n){BCDC b(n);return G(b.N,b.adj);}

static vector<int> bfsRoute(const G& g, int s, int t) {
	vector<int> par(g.N, -1); queue<int> q; q.push(s); par[s] = s;
	while (!q.empty()) { int u = q.front(); q.pop(); if (u == t) break; for (int v : g.adj[u]) if (par[v] < 0) { par[v] = u; q.push(v); } }
	vector<int> path; if (par[t] < 0) return path;
	for (int c = t; ; c = par[c]) { path.push_back(c); if (c == s) break; }
	reverse(path.begin(), path.end()); return path;
}

// max link load under all-to-all shortest-path routing (congestion metric)
static long long maxLinkLoadAllToAll(const G& g) {
	unordered_map<long long, int> load;
	for (int s = 0; s < g.N; ++s) {
		for (int t = 0; t < g.N; ++t) {
			if (s == t) continue;
			vector<int> p = bfsRoute(g, s, t);
			for (size_t i = 0; i + 1 < p.size(); ++i) {
				int a = p[i], b = p[i + 1]; if (a > b) swap(a, b);
				load[((long long)a << 32) | (unsigned)b]++;
			}
		}
	}
	long long mx = 0; for (auto& kv : load) mx = max(mx, (long long)kv.second);
	return mx;
}

static double avgHops(const G& g, int samples) {
	mt19937 rng(11); double tot = 0;
	for (int i = 0; i < samples; ++i) {
		int s = rng() % g.N, t = rng() % g.N; if (s == t) { --i; continue; }
		tot += bfsRoute(g, s, t).size() - 1;
	}
	return tot / samples;
}

int main() {
	for (int n = 3; n <= 4; ++n) {
		cout << "===== n=" << n << " =====\n";
		struct Row { const char* name; G g; };
		vector<Row> rows = { { "DHDC", buildDHDC(n) }, { "HSDC", buildHSDC(n) }, { "BCDC", buildBCDC(n) }, { "CSDC", buildCSDC(n) } };
		for (auto& r : rows) {
			const G& g = r.g;
			double ah = avgHops(g, 2000);
			long long mll = maxLinkLoadAllToAll(g);
			// throughput via flow_sim (random flows)
			Config cfg; cfg.N = g.N; cfg.linkCap = 1e9;
			cfg.route = [&](int s, int t) { return bfsRoute(g, s, t); };
			mt19937 rng(123);
			vector<Flow> flows; int num = (g.N <= 64 ? 2000 : 2000);
			uniform_int_distribution<int> uid(0, g.N - 1);
			exponential_distribution<double> inter((double)num / 100.0);
			double tt = 0;
			for (int i = 0; i < num; ++i) { Flow f; f.id = i; f.src = uid(rng); do { f.dst = uid(rng); } while (f.dst == f.src); f.size = 1e6; f.tArr = tt; tt += inter(rng); flows.push_back(f); }
			Result res = run(cfg, std::move(flows));
			cout << "  " << r.name << "  N=" << g.N << "  avgHops=" << ah
				<< "  maxLinkLoad(all-to-all)=" << mll
				<< "  maxThp=" << res.maxThp / 1e9 << "Gbps  avgThp=" << res.avgThp / 1e9 << "Gbps  FCT=" << res.totalFCT << "s\n";
		}
		cout.flush();
	}
	return 0;
}
