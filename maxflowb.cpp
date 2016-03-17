#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <cstdio>
#include <vector>

using namespace std;

const int N = 100 +  5;
const int oo = 0x3f3f3f3f;

struct Edge {
  int from, to, cap, flow;
};

struct Dinic {
  int n, m, s, t;
  int dis[N], cur[N], que[N << 1];
  bool vis[N];
  vector <Edge> edges;
  vector <int> G[N];

  void add(int from, int to, int cap) {
	edges.push_back((Edge) {from, to, cap, 0});
	edges.push_back((Edge) {to, from, 0, 0});
	m = edges.size();
	G[from].push_back(m - 2);
	G[to].push_back(m - 1);
  }

  bool bfs() {
	int head = 1, tail = 1;

	memset(vis, false, sizeof vis);
	dis[s] = 0; vis[s] = true; que[head] = s;
	while(head <= tail) {
	  int x = que[head];

	  for(int i = 0; i < (signed) G[x].size(); ++ i) {
		Edge &e = edges[G[x][i]];

		if(!vis[e.to] && e.cap > e.flow) {
		  vis[e.to] = true;
		  dis[e.to] = dis[x] + 1;
		  que[++ tail] = e.to;
		}
	  }
	  ++ head;
	}
	return vis[t];
  }

  int dfs(int x, int a) {
	if(x == t || a == 0) return a;

	int flw = 0, f;

	for(int &i = cur[x]; i < (signed) G[x].size(); ++ i) {
	  Edge &e = edges[G[x][i]];
	  
	  if(dis[e.to] == dis[x] + 1 && (f = dfs(e.to, min(a, e.cap - e.flow))) > 0) {
		e.flow += f; edges[G[x][i] ^ 1].flow -= f; flw += f; a -= f;
		if(a == 0) break;
	  }
	}
	return flw;
  }

  int MaxFlow(int s, int t) {
	this->s = s; this->t = t;

	int flw = 0;

	while(bfs()) {
	  memset(cur, 0, sizeof cur);
	  flw += dfs(s, oo);
	}
	return flw;
  }
}net;

int n, M[N];

int main() {
#ifndef ONLINE_JUDGE
  freopen("maxflowb.in", "r", stdin);
  freopen("maxflowb.out", "w", stdout);
#endif

  int Up, Down;
  
  scanf("%d", &n); net.n = n + 1;
  for(int i = 1; i <= n; ++ i) {
	for(int j = 1; j <= n; ++ j) {
	  scanf("%d%d", &Down, &Up);
	  M[i] -= Down; M[j] += Down;
	  net.add(i, j, Up - Down);
	}
  }
  net.add(n, 1, oo);
  for(int i = 1; i <= n; ++ i) {
	if(M[i] > 0) {
	  net.add(0, i, M[i]);
	}
	else if(M[i] < 0){
	  net.add(i, n + 1, -M[i]);
	}
  }
  net.MaxFlow(0, n + 1);
  printf("%d\n", net.MaxFlow(1, n));
  
#ifndef ONLINE_JUDGE
  fclose(stdin); fclose(stdout);
#endif
  return 0;
}

