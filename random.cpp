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
double dsi = 0.75, dtd = 0.05, dte = 0.05; //gene deactivation inheritance, gene deactivation mutation rate, gene activation mutation rate
double adn = 0.06, adw = 0.16; //rate of edge-split node addition mutation, edge addition mutation
double sur = 0.7, isp = 0.005; //proportion of immediate survival, interspecies mating rate
double pur = 0.25, pum = 0.04, pud = 0.1; //probability of descent with mutation absent crossover, neither given former, descent with crossover absent mutation
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

struct edgcmp {inline bool operator() (const edg& a, const edg& b) {return (a.n < b.n);}};

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

    int dme(edg i) { //check if an edge is {disjoint, matched, excess} -> {-1, 0, 1}, sort nodes first
        int j = lower_bound(edge.begin(), edge.end(), i, edgcmp()) - edge.begin();
        if (j < edge.size() && i.n == edge[j].n) return 0;
        j = lower_bound(node.begin(), node.end(), i.a) - node.begin();
        if (j + 1 > node.size() || i.a != node[j]) return -1;
        j = lower_bound(node.begin(), node.end(), i.b) - node.begin();
        if (j + 1 > node.size() || i.b != node[j]) return -1;
        return 1;
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
        ea.n = (idl.find(make_pair(c, b)) != idl.end() ? idl[make_pair(c, b)] : idl[make_pair(c, b)] = inw++);
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
    return n;
}

net ncr(net a, net b) { //crossbreed
    int i, j;
    if (a.F < b.F) swap(a, b);
    for (i = j = 0; i < a.edge.size(); ++i) {
        while (j < b.edge.size() && b.edge[j].n < a.edge[i].n) ++j;
        if (a.edge[i].n == b.edge[i].n) {
            a.edge[i].w = rng(2) ? b.edge[j].w : a.edge[i].w;
            a.edge[i].e = (a.edge[i].e && b.edge[i].e) ? 1 : (rdn() > dsi);
        }
    }
    if (rdn() > pud) a = mft(a);
    return a;
}

struct spec {
    int n, ag, ls; //identifier, species age, time since last improvement
    double F, mf; //current fitness, max fitness so far
    vector<net> pop; //species subpopulation
    net rep;

    spec sbp(int m, int r) { //bottleneck to best m and regenerate to r
        int a, b;
        sort(pop.begin(), pop.end(), invnetcmp());
        spec spr = {.n = n, .ag = ag, .ls = ls, .F = 0, .mf = mf, .pop = pop};
        if (r >= lar) spr.pop.push_back(pop[0]);
        while (spr.pop.size() < r) {
            if (rdn() < pur) spr.pop.push_back(rdn() < pum ? pop[rng(m)] : mft(pop[rng(m)]));
            else {
                a = rng(m);
                b = rng(m);
                net pnl = ncr(pop[a], pop[b]);
                spr.pop.push_back(rdn() < pud ? pnl : mft(pnl));
            }
        }
        return spr;
    }

    double cmf() {
        for (int i = 0; i < pop.size(); ++i) F = max(F, pop[i].F);
        F /= pop.size();
    }

    static bool eq(net a, net b) { //check for same species
        int c, mt = 0, i;
        double dst = 0, mdw = 0;
        if (a.edge.size() < b.edge.size()) swap(a, b);
        sort(a.node.begin(), a.node.end());
        sort(b.node.begin(), b.node.end());
        for (i = 0; i < a.edge.size(); ++i) {
            c = b.dme(a.edge[i]);
            if (c == -1) dst += cw2;
            if (!c) {++mt; mdw += fabs(a.edge[i].w - lower_bound(b.edge.begin(), b.edge.end(), a.edge[i], edgcmp())->w);}
            if (c == 1) dst += cw1;
        }
        for (i = 0; i < b.edge.size(); ++i) {
            c = a.dme(b.edge[i]);
            if (c == -1) dst += cw2;
            if (c == 1) dst += cw1;
        }
        dst /= a.edge.size();
        dst += mdw * cw3 / mt;
        return dst < delt;
    }
};

vector<net> ppn[2], npv; //population
vector<spec> spc[2]; //species
net nnet; //empty net
spec nspec = {.n = -1, .ag = 0, .ls = 0, .F = 0, .mf = 0, .pop = npv, .rep = nnet}; //empty species

inline void init() {
    //TODO READ INPUT - JAVA INTERFACE
    //TODO GENERATE NETS FROM DEFAULT OR INPUT
}

inline void epoch(int tI) {
    int np = 0, sf = 0, a, b, i, j;
    double sum = 0;
    //TODO EVALUATE NETS
    //TODO SPECIATION
    vector<pair<int, int> > cbp; //bottleneck and regeneration
    for (i = 0; i < N; ++i) {
        for (j = 0; j < spc[0].size(); ++j) if (spec::eq(ppn[0][i], spc[0][j].rep)) {spc[0][j].pop.push_back(ppn[0][i]); break;}
        if (j + 1 > spc[0].size()) {
            spc[0].push_back(nspec);
            spc[0].back().n = ins++;
            spc[0].back().pop.push_back(ppn[0][i]);
            spc[0].back().rep = ppn[0][i];
        }
    }
    for (i = 0; i < N; ++i) np += rdn() > isp;
    for (i = 0; i < spc[0].size(); ++i) spc[0][i].cmf(), sum += spc[0][i].F;
    for (i = 0; i < spc[0].size(); ++i) cbp.push_back(make_pair(max((int) (spc[0][i].pop.size() * sur), 1), (N - np) * spc[0][i].F / sum)), sf += (N - np) * spc[0][i].F / sum;
    for (i = sf; i < N - np; ++i) ++cbp[rng(spc[0].size())].second;
    for (i = 0; i < spc[0].size(); ++i) if (cbp[i].second) {

    }
    for (i = 0; i < np; ++i) {

    }



    //TODO FINISH REGEN
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
