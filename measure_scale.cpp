#include "flow_sim.hpp"
#include <iostream>
#include <vector>
#include <functional>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <cstdio>
#include <cmath>
using namespace std;
using namespace flowsim;

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
static vector<int> bfsRoute(const G& g,int s,int t){vector<int>par(g.N,-1);queue<int>q;q.push(s);par[s]=s;while(!q.empty()){int u=q.front();q.pop();if(u==t)break;for(int v:g.adj[u])if(par[v]<0){par[v]=u;q.push(v);}}vector<int>p;if(par[t]<0)return p;for(int c=t;;c=par[c]){p.push_back(c);if(c==s)break;}reverse(p.begin(),p.end());return p;}

static void measure(const G& g, double& avgHops, double& loadFactor) {
	mt19937 rng(12345u);
	int samples = 10000;
	double tot = 0;
	unordered_map<long long,int> load;
	for (int k = 0; k < samples; ++k) {
		int s = rng() % g.N, t = rng() % g.N; if (s == t) { --k; continue; }
		vector<int> p = bfsRoute(g, s, t);
		tot += p.size() - 1;
		for (size_t i = 0; i + 1 < p.size(); ++i) { int a = p[i], b = p[i+1]; if (a > b) swap(a,b); load[((long long)a<<32)|(unsigned)b]++; }
	}
	avgHops = tot / samples;
	long long mll = 0; for (auto& kv : load) mll = max(mll, (long long)kv.second);
	double avgLinkLoad = tot / (double)load.size();
	loadFactor = mll / avgLinkLoad;
}

static void throughput(const G& g, double& maxThp, double& avgThp) {
	Config cfg; cfg.N = g.N; cfg.linkCap = 1e9;
	cfg.route = [&](int s, int t) { return bfsRoute(g, s, t); };
	mt19937 rng(6789u);
	vector<Flow> flows;
	int F = 300;
	for (int i = 0; i < F; ++i) {
		Flow f; f.id = i; f.src = rng() % g.N; do { f.dst = rng() % g.N; } while (f.dst == f.src);
		f.size = 1e8; f.tArr = 0; flows.push_back(f);
	}
	Result r = run(cfg, std::move(flows));
	maxThp = r.maxThp / 1e9;
	avgThp = r.avgThp / 1e9;
}

int main() {
	// target server scales and the n for each model (chosen so N is closest)
	struct S { int target; int dn, hn, cn, bn; };
	vector<S> scales = { {64,3,4,4,5}, {384,4,6,6,7}, {2048,5,8,8,9}, {10240,6,10,10,11} };
	printf("ScaleN,Model,n,N,Degree,Diameter,avgHops,loadFactor,maxThpGbps,avgThpGbps\n");
	for (auto& sc : scales) {
		struct Row { const char* name; int n; int deg; int diam; function<G(int)> build; };
		vector<Row> rows = {
			{"DHDC", sc.dn, sc.dn, 4*sc.dn-4, buildDHDC},
			{"HSDC", sc.hn, sc.hn, 2*sc.hn, buildHSDC},
			{"CSDC", sc.cn, sc.cn, 2*((sc.cn+1)/2)+1, buildCSDC},
			{"BCDC", sc.bn, 2*sc.bn-2, (sc.bn+1)/2+1, buildBCDC},
		};
		for (auto& r : rows) {
			G g = r.build(r.n);
			double ah, lf; measure(g, ah, lf);
			double mx, av; throughput(g, mx, av);
			printf("%d,%s,%d,%d,%d,%d,%.2f,%.2f,%.2f,%.2f\n", sc.target, r.name, r.n, g.N, r.deg, r.diam, ah, lf, mx, av);
			fflush(stdout);
		}
	}
	return 0;
}
