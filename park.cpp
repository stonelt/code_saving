//WC2013 糖果公园
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <stack>

using namespace std;

const int N = 100000 + 5;
typedef long long ll;

inline int read() {
  int x = 0;
  char c = getchar();

  while(!isdigit(c)) c = getchar();
  while(isdigit(c)) {
	x = x * 10 + c - '0';
	c = getchar();
  }
  return x;
}

int buf[30];

inline void output(ll x) {
  int p = 0;

  buf[0] = 0;
  if(x == 0) p ++;
  else {
	while(x) {
	  buf[p ++] = x % 10;
	  x /= 10;
	}
  }
  for(int j = p - 1; j >= 0; -- j)
	putchar('0' + buf[j]);
}

ll ans[N], tot = 0;
int n, m, q, size, cnt, tk, tim;
int head[N], fa[N], f[N][18], seq[N], v[N], w[N];
int ckn[N], depth[N], son[N], pre[N], c[N];
bool flag[N]; int vis[N];
stack <int> s;

struct Edge {
  int from, to, next;
}edges[N << 1];
struct Query {
  int x, v, pre;
}C[N];
struct query {
  int id, u, v, t;
  bool operator < (const query &k) const {
	if(ckn[u] == ckn[k.u] && ckn[v] == ckn[k.v]) return t < k.t;
	else if(ckn[u] == ckn[k.u]) return ckn[v] < ckn[k.v];
	else return ckn[u] < ckn[k.u];
  }
}Q[N];

void insert(int from, int to) {
  ++ cnt;
  edges[cnt].from = from; edges[cnt].to = to;
  edges[cnt].next = head[from]; head[from] = cnt;
}

void dfs(int u, int f) {
  fa[u] = f; seq[u] = ++ tim;
  for(int i = head[u]; i; i = edges[i].next) {
	int v = edges[i].to;

	if(v != f) {
	  depth[v] = depth[u] + 1;
	  dfs(v, u);
	  son[u] += son[v];
	  if(son[u] >= size) {
		++ tk;
		son[u] = 0;
		while(!s.empty()) {
		  int cur = s.top(); s.pop();
		  ckn[cur] = tk;
		}
	  }
	}
  }
  s.push(u); son[u] ++;
}

void prepare() {
  memset(f, -1, sizeof f);
  for(int i = 1; i <= n; ++ i) f[i][0] = fa[i];
  for(int j = 1; (1 << j) <= n; ++ j) {
	for(int i = 1; i <= n; ++ i) {
	  if(f[i][j - 1] != -1) {
		f[i][j] = f[f[i][j - 1]][j - 1];
	  }
	}
  }
}

int lca(int a, int b) {
  int i;

  if(depth[a] < depth[b]) swap(a, b);
  for(i = 0; (1 << i) <= depth[a]; ++ i);
  -- i;
  for(int j = i; j >= 0; -- j) {
	if(depth[a] - depth[b] >= (1 << j))
	  a = f[a][j];
  }
  if(a == b) return a;
  for(int j = i; j >= 0; -- j) {
	if(f[a][j] != -1 && f[a][j] != f[b][j]) {
	  a = f[a][j];
	  b = f[b][j];
	}
  }
  return f[a][0];
}

void rever(int u) {
  if(flag[u]) {
	flag[u] = false;
	tot -= 1LL * w[vis[c[u]]] * v[c[u]];
	vis[c[u]] --;
  }
  else {
	flag[u] = true;
	vis[c[u]] ++;
	tot += 1LL * w[vis[c[u]]] * v[c[u]];
  }
}

void Getpath(int u, int fi) {
  while(u != fi) {
	rever(u);
	u = fa[u];
  }
}

void change(int x, int v) {
  if(flag[x]) {
	rever(x); c[x] = v; rever(x);
  }
  else c[x] = v;
}

void solve(int nl, int nr, int &zl, int &zr, int ts) {
  int Lca;

  Lca = lca(nl, zl);
  Getpath(nl, Lca); Getpath(zl, Lca);
  Lca = lca(nr, zr);
  Getpath(nr, Lca); Getpath(zr, Lca);
  Lca = lca(nl, nr);
  rever(Lca); ans[Q[ts].id] = tot; rever(Lca);
  zl = nl; zr = nr;
}

#define ONLINE_JUDGE

int main() {
  //int __size__ =  64 <<20;  //这里谁<<20位，就是多少M的栈
  //char*__p__ =(char*)malloc(__size__)+ __size__;
  //__asm__("movl %0, %%esp\n"::"r"(__p__));

#ifndef ONLINE_JUDGE
  freopen("park.in", "r", stdin);
  freopen("park.out", "w", stdout);
#endif

  int x, y;

  n = read(); m = read(); q = read();
  for(int i = 1; i <= m; ++ i) v[i] = read();
  for(int i = 1; i <= n; ++ i) w[i] = read();
  for(int i = 1; i < n; ++ i) {
	x = read(); y = read();
	insert(x, y); insert(y, x);
  }
  for(int i = 1; i <= n; ++ i) pre[i] = c[i] = read();
  depth[1] = 1; //size = (int) pow((double) n, (double) 2 / 3);
  size = 1500;
  dfs(1, -1);
  if(!s.empty()) {
	++ tk;
	while(!s.empty()) {
	  int cur = s.top(); s.pop();
	  ckn[cur] = tk;
	}
  }
  prepare();

  int type, t1 = 0, t2 = 0;
  
  for(int i = 1; i <= q; ++ i) {
	type = read();
	if(!type) {
	  ++ t1;
	  C[t1].x = read(); C[t1].v = read(); C[t1].pre = pre[C[t1].x];
	  pre[C[t1].x] = C[t1].v;
	}
	else {
	  ++ t2;
	  Q[t2].u = read(); Q[t2].v = read(); Q[t2].id = t2; Q[t2].t = t1;
	  if(seq[Q[t2].u] > seq[Q[t2].v]) swap(Q[t2].u, Q[t2].v);
	}
  }

  int L = 1, R = 1;

  sort(Q + 1, Q + t2 + 1);
  for(int i = 1; i <= Q[1].t; ++ i) {
	change(C[i].x, C[i].v);
  }
  solve(Q[1].u, Q[1].v, L, R, 1);
  for(int i = 2; i <= t2; ++ i) {
	for(int t = Q[i - 1].t + 1; t <= Q[i].t; ++ t)
	  change(C[t].x, C[t].v);
	for(int t = Q[i - 1].t; t > Q[i].t; -- t)
	  change(C[t].x, C[t].pre);
	solve(Q[i].u, Q[i].v, L, R, i);
  }
  for(int i = 1; i <= t2; ++ i) {
	output(ans[i]);
	puts("");
  }
  
#ifndef ONLINE_JUDGE
  fclose(stdin); fclose(stdout);
#endif
  return 0;
}
