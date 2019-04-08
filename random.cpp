#include <bits/stdc++.h>
using namespace std;

int T = 300, N = 150, inn = 0, inw = 0, ins = 0; //epochs, total population size, innovation numbers - node, edge, identifier of species
int n0 = 200, n1 = 6; //input nodes, output nodes
int lar = 5; //threshold for large species - unchanged copy of champion
int kil = 15; //max generations w/o improvement for species
double sig = 4.9; //sigmoid multiplier
double cw1 = 1, cw2 = 1, cw3 = 0.4, delt = 3; //excess, disjoint, matching weight, threshold - compatibility constants
double mut = 0.8, rel = 0.9; //weight mutation rate, relative mutation vs absolute reassignment
double mug = 0.25, msg = 0.5, muh = 0.7, msh = 0.6; //{intensity, standard deviation} of relative mutation, absolute mutation
double dsd = 0.75, dtd = 0.05, dte = 0.05; //gene deactivation inheritance, gene deactivation mutation rate, gene activation mutation rate
double adn = 0.1, adw = 0.1; //edge-split node addition mutation, edge addition mutation
double sur = 0.7, isp = 0.005; //proportion of average fitness threshold for survival, interspecies mating rate
double pur = 0.25; //probability of descent with mutation only without crossover
queue<int> q; //heap queue for network query - node identifiers
map<int, pair<int, double> > m; //heap map for network query - node identifiers to remaining edge dependencies, activation
map<int, int> rm; //heap map for phenotype node index based on identifier
map<int, int> idn; //duplication removal map for node addition - identifier of split edge
map<pair<int, int>, int> idl; //duplication removal map for edge addition - identifier of two nodes

inline double sgm(double x) {return 1 / (1 + exp(-sig * x));}
inline long long rng(long long n) {long long i = 0, j; for (j = 0; j < 6; ++j) i *= RAND_MAX, i += rand(), i %= n; return n;}
inline double rdn() {return ((double) rng(1000000001)) / 1000000000;}
default_random_engine nrmgn;
normal_distribution<double> dsg(0, msg), dsh(0, msh);

struct edg {
    int n, a, b; //identifier, source, destination
    double w; //weight
    bool e; //enabled

    void mwt() { //mutate weight
        /*if (rdn() < rel) w += dsg(nrmgn) * mug;
        else w = dsh(nrmgn) * muh;*/
        if (rdn() < rel) w += (rdn() * 2 - 1) * mug;
        else w = (rdn() * 2 - 1) * muh;
        if (e) {if (rdn() < dtd) e = 0;}
        else {if (rdn() < dte) e = 1;}
    }
};

struct net {
    int L; //layers
    double F; //fitness value
    vector<int> node; //nodes
    vector<edg> edge; //edges - implicitly always sorted by identifier
    vector<vector<edg> > adj, ret; //forward adjacency list, reverse adjacency list

    void init() { //initialize one network
        int i;
        L = 2;
        for (i = 0; i < n0 + 1 + n1; ++i) node.push_back(inn++); //extra 1 for bias
        adj.resize(n0 + 1 + n1);
        ret.resize(n0 + 1 + n1);
    }

    void muw() {for (int i = 0; i < edge.size(); ++i) if (rdn() < mut) edge[i].mwt();} //mutate all edge weights

    void reg() { //regenerate network from edge list
        int i;
        adj.clear();
        ret.clear();
        for (i = 0; i < node.size(); ++i) rm[node[i]] = i;
        for (i = 0; i < edge.size(); ++i) adj[rm[edge[i].a]].push_back(edge[i]), ret[rm[edge[i].b]].push_back(edge[i]);
    }

    bool trv(int a, int b) { //check if b is accessible from a - prevent cycles
        int i, j;
        queue<int> q;
        vector<bool> vt;
        vt.assign(node.size(), 0);
        q.push(a);
        while (!q.empty() && q.front() != b) {
            i = q.front();
            q.pop();
            if (vt[i]) continue;
            vt[i] = 1;
            for (j = 0; j < adj[i].size(); ++j) if (!vt[adj[i][j].b]) q.push(adj[i][j].b);
        }
        return !q.empty();
    }

    bool dhe(int a, int b) { //check if there already exists a direct edge from a to b
        for (int i = 0; i < adj[a].size(); ++i) if (adj[a][i].b == b) return 1;
        return 0;
    }

    vector<double> query(vector<double> v) { //run network
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
        vector<double> rs;
        for (i = n0 + 1; i < n0 + n1 + 1; ++i) rs.push_back(sgm(m[node[i]].second));
        m.clear();
        return rs;
    }
};

struct invnetcmp {inline bool operator() (const net& a, const net& b) {return (a.F > b.F);}};

net mft(net _n) { //mutate
    int a, b;
    net n = _n;
    n.muw();
    if (rdn() < adn) { //new node + 2 edges
        edg* eg = &n.edge[rng(n.edge.size())];
        a = eg->a;
        b = eg->b;
        int c = (idn.find(eg->n) != idn.end() ? idn[eg->n] : idn[eg->n] = inn++);
        n.node.push_back(c);
        edg ea, eb;
        ea.n = (idl.find(make_pair(a, c)) != idl.end() ? idl[make_pair(a, c)] : idl[make_pair(a, c)] = inw++);
        ea.a = a;
        ea.b = c;
        ea.e = eg->e;
        ea.w = 1;
        eb.n = (idl.find(make_pair(c, b)) != idl.end() ? idl[make_pair(c, b)] : idl[make_pair(c, b)] = inw++);
        eb.a = c;
        eb.b = b;
        eb.e = eg->e;
        eb.w = eg->w;
        if (ea.n > eb.n) swap(ea, eb);
        n.edge.push_back(ea);
        n.edge.push_back(eb);
        eg->e = 0;
    }
    if (rdn() < adw) { //new edge
        a = b = -1;
        while (n.dhe(a, b)) {while (a == b) a = n.node[rng(n.node.size())], b = n.node[rng(n.node.size())]; if (n.trv(b, a)) swap(a, b);} //new, acyclic edge
        edg eg;
        eg.n = (idl.find(make_pair(a, b)) != idl.end() ? idl[make_pair(a, b)] : idl[make_pair(a, b)] = inw++);
        eg.a = a;
        eg.b = b;
        eg.e = 1;
        eg.w = (rdn() * 2 - 1) * muh;
    }
}

net ncr(net a, net b) { //crossbreed

}

struct spec {
    int n, ag, ls; //identifier, species age, time since last improvement
    double mf; //max fitness so far
    vector<net> pop; //species subpopulation

    net rep() {return pop[0];} //species representative for containment test

    spec sbp(int m) { //best m of species subpopulation
        spec spr = {.n = n, .ag = ag, .ls = ls, .mf = mf, .pop = pop};
        sort(spr.pop.begin(), spr.pop.end(), invnetcmp());
        while (spr.pop.size() > m) spr.pop.pop_back();
        return spr;
    }
};

vector<net> ppn[2]; //population
vector<spec> spc[2]; //species
net nnet; //empty net
spec nspec; //empty species

inline void init() {
    //TODO READ INPUT - JAVA INTERFACE
    //TODO GENERATE NETS FROM DEFAULT OR INPUT
}

inline void epoch(int tI) {
    //TODO EVALUATE NETS
    //TODO SURVIVAL & ELIMINATION
    //TODO REGEN POPULATION
    //TODO REFILL SPECIES


    idn.clear();
    idl.clear();
}

int main() {
    int tI, i;
    srand(0);
    init();
    for (tI = 0; tI < T; ++tI) epoch(tI);
    //TODO OUTPUT RESULT
}
