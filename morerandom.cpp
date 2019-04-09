#include <bits/stdc++.h>
using namespace std;

int n0 = 200, n1 = 6; //input nodes, output nodes
double sig = 4.9; //sigmoid multiplier
queue<int> q; //heap queue for network query - node identifiers
map<int, pair<int, double> > m; //heap map for network query - node identifiers to remaining edge dependencies, activation
map<int, int> rm; //heap map for phenotype node index based on identifier
vector<double> v, rs; //query, result vectors

inline double sgm(double x) {return 1 / (1 + exp(-sig * x));}

struct edg {
    int n, a, b; //identifier, source, destination
    double w; //weight
    bool e; //enabled
};

struct net {
    bool E = 0; //fitness evaluated
    double F = 0; //fitness value
    vector<int> node; //nodes
    vector<edg> edge; //edges - implicitly always sorted by identifier
    vector<vector<edg> > adj, ret; //forward adjacency list, reverse adjacency list

    void reg() { //regenerate network from edge list
        int i;
        adj.clear();
        ret.clear();
        for (i = 0; i < node.size(); ++i) rm[node[i]] = i;
        for (i = 0; i < edge.size(); ++i) adj[rm[edge[i].a]].push_back(edge[i]), ret[rm[edge[i].b]].push_back(edge[i]);
    }

    void query() { //run network
        assert(v.size() == n0);
        int i, j;
        for (i = 0; i < node.size(); ++i) m[node[i]] = make_pair(ret[i].size(), i < n0 + 1 ? (i < n0 ? v[i] : 1) : 0);
        for (i = 0; i < n0 + 1; ++i) q.push(i);
        while (!q.empty()) {
            i = q.front();
            for (j = 0; j < adj[i].size(); ++j) {
                if (adj[i][j].e) m[adj[i][j].b].second += sgm(m[i].second) * adj[i][j].w;
                if (!--m[adj[i][j].b].first) q.push(adj[i][j].b);
            }
            q.pop();
        }
        for (i = n0 + 1; i < n0 + n1 + 1; ++i) rs.push_back(sgm(m[node[i]].second));
        m.clear();
    }
};

net N;
edg eg;

int main() {
    int M, L, x, i;
    double F;
    FILE* src = fopen("net.txt", "r"), *in, *out, *flag;
    fscanf(src, "%d %d %d %d", &n0, &n1, &M, &L);
    for (i = 0; i < M; ++i) fscanf(src, "%d", &x), N.node.push_back(x);
    for (i = 0; i < L; ++i) fscanf(src, "%d %d %d %lf %d", &eg.n, &eg.a, &eg.b, &eg.w, &eg.e), N.edge.push_back(eg);
    fclose(src);
    N.reg();
    while (1) {
        if (flag = fopen("end.txt", "r")) {
            fscanf(flag, "%ld", &F);
            fclose(flag);
            remove("end.txt");
            break;
        }
        if (in = fopen("input.txt", "r")) {
            for (i = 0; i < n0; ++i) scanf("%d", &x), v.push_back(x);
            v.push_back(1000000000);
            fclose(in);
            remove("input.txt");
            N.query();
            out = fopen("output.txt", "w");
            for (i = 0; i < rs.size(); ++i) fprintf(out, "%lf%c", rs[i], i < rs.size() - 1 ? ' ' : '\n');
            fclose(out);
            rs.clear();
        }
        this_thread::sleep_for(chrono::microseconds(2000));
    }
    return 0;
}
