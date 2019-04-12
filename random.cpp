#include <bits/stdc++.h>
#include <windows.h>
#include <io.h>
using namespace std;

//remember to use -U__STRICT_ANSI__ compiler flag
HANDLE H;

int T = 300, N = 150, inn = 0, inw = 0, ins = 0; //epochs, total population size, innovation numbers - node, edge, identifier of species
int n0 = 407, n1 = 6; //input nodes, output nodes
int lar = 12, mup = 400; //threshold for large species - unchanged copy of champion, passes of attempted mutation
double stg = 4; //stagnation time constant
double sig = 4.9; //sigmoid multiplier
double cw1 = 1, cw2 = 1, cw3 = 0.4, delt = 3; //excess, disjoint, matching weight, threshold - compatibility constants
double mut = 0.8, rel = 0.9; //weight mutation rate, relative mutation vs absolute reassignment
double mug = 0.25, msg = 0.5, muh = 0.7, msh = 0.6; //{intensity, standard deviation} of relative mutation, absolute mutation
double dsi = 0.75, dtd = 0.05, dte = 0.05; //gene deactivation inheritance, gene deactivation mutation rate, gene activation mutation rate
double adn = 0.06, adw = 0.16; //rate of edge-split node addition mutation, edge addition mutation
double sur = 0.7, isp = 0.05; //proportion of immediate survival, interspecies mating rate
double pur = 0.25, pum = 0.04, pud = 0.1; //probability of descent with mutation absent crossover, neither given former, descent with crossover absent mutation
queue<int> q; //heap queue for network query - node identifiers
vector<pair<int, double> > m; //heap vector for network query - node identifiers to remaining edge dependencies, activation
map<int, int> rm; //heap map for phenotype node index based on identifier
map<int, int> idn; //duplication removal map for node addition - identifier of split edge
map<pair<int, int>, int> idl; //duplication removal map for edge addition - identifier of two nodes

inline bool eq(double a, double b) {return fabs(a - b) < 0.0001;}
inline double sgm(double x) {return 1 / (1 + exp(-sig * x));}
inline long long rng(long long n) {long long i = 0, j; for (j = 0; j < 6; ++j) i *= RAND_MAX, i += rand(), i %= n; return i;}
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

vector<double> qv, esr;

struct net {
    bool E = 0; //fitness evaluated
    double F = 0; //fitness value
    vector<int> node; //nodes
    vector<edg> edge; //edges - implicitly always sorted by identifier
    vector<vector<edg> > adj, ret; //forward adjacency list, reverse adjacency list

    void init() { //initialize one network
        for (int i = 0; i < n0 + 1 + n1; ++i) node.push_back(inn++); //extra 1 for bias
        adj.resize(n0 + 1 + n1);
        ret.resize(n0 + 1 + n1);
    }

    void muw() {for (int i = 0; i < edge.size(); ++i) if (rdn() < mut) edge[i].mwt();} //mutate all edge weights

    void reg() { //regenerate network from edge list
        int i;
        adj.clear();
        ret.clear();
        adj.resize(node.size());
        ret.resize(node.size());
        for (i = 0; i < node.size(); ++i) rm[node[i]] = i;
        for (i = 0; i < edge.size(); ++i) adj[rm[edge[i].a]].push_back(edge[i]), ret[rm[edge[i].b]].push_back(edge[i]);
    }

    bool trv(int a, int b) { //check if b is accessible from a - prevent cycles
        //printf("traversal check: %d -> %d\n", a, b);
        if (a < n0 + 1) return 1;
        if (n0 < b && b < n0 + 1 + n1) return 1;
        int i, j;
        queue<int> q;
        vector<bool> vt;
        vt.assign(node.size(), 0);
        for (i = 0; i < node.size(); ++i) rm[node[i]] = i;
        q.push(a);
        while (!q.empty() && q.front() != b) {
            i = q.front();
            q.pop();
            if (vt[i]) continue;
            vt[i] = 1;
            for (j = 0; j < adj[i].size(); ++j) if (!vt[rm[adj[i][j].b]]) q.push(rm[adj[i][j].b]);
        }
        rm.clear();
        return !q.empty();
    }

    bool dhe(int a, int b) { //check if there already exists a direct edge from a to b
        //printf("direct edge check: %d %d\n", a, b);
        for (int i = 0; i < adj[a].size(); ++i) if (adj[a][i].b == node[b]) return 1;
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

    void query() { //run network
        assert(qv.size() == n0);
        int i, j;
        for (i = 0; i < node.size(); ++i) rm[node[i]] = i, m.push_back(make_pair(ret[i].size(), i < n0 + 1 ? (i < n0 ? qv[i] : 1000000000) : 0));
        for (i = 0; i < n0 + 1; ++i) if (!m[i].first) q.push(i);
        while (!q.empty()) {
            i = q.front();
            for (j = 0; j < adj[i].size(); ++j) {
                if (adj[i][j].e) m[rm[adj[i][j].b]].second += sgm(m[i].second) * adj[i][j].w;
                if (!--m[rm[adj[i][j].b]].first) q.push(rm[adj[i][j].b]);
            }
            q.pop();
        }
        for (i = n0 + 1; i < n0 + n1 + 1; ++i) esr.push_back(sgm(m[i].second));
        m.clear();
        qv.clear();
        rm.clear();
    }
};

struct invnetcmp {inline bool operator() (const net& a, const net& b) {return (a.F > b.F);}};

net mft(net _n) { //mutate
    int a, b, zf, i;
    net n = _n;
    n.muw();
    for (i = 0; i < mup; ++i) {
        if (rdn() < adn && !n.edge.empty()) { //new node + 2 edges
            //printf("mutating node\n");
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
            //printf("!? %d %d: %d -> %d, %d -> %d\n", ea.n, eb.n, ea.a, ea.b, eb.a, eb.b);
            eg->e = 0;
        } else if (rdn() < adw) { //new edge
            //printf("mutating edge\n");
            //printf("nodes: %d\n", n.node.size());
            a = b = -1;
            for (zf = 0; zf < 10 && ((a < 0 && b < 0) || n.trv(b, a) || n.dhe(a, b)); ++zf) {do {a = rng(n.node.size()); b = rng(n.node.size());} while (a == b); if (n.trv(b, a)) swap(a, b);} //new, acyclic edge
            if (zf > 9) continue;
            a = n.node[a];
            b = n.node[b];
            edg eg;
            eg.n = (idl.find(make_pair(a, b)) != idl.end() ? idl[make_pair(a, b)] : (idl[make_pair(a, b)] = inw++));
            eg.a = a;
            eg.b = b;
            eg.e = 1;
            eg.w = (rdn() * 2 - 1) * muh;
            //printf("!!! %d | inw %d: %d -> %d\n", eg.n, inw, a, b);
            n.edge.push_back(eg);
        }
        n.reg();
    }
    return n;
}

net ncr(net a, net b) { //crossbreed
    //printf("crossbreeding\n");
    int i, j;
    if (a.F < b.F) swap(a, b);
    for (i = j = 0; i < a.edge.size(); ++i) {
        while (j < b.edge.size() && b.edge[j].n < a.edge[i].n) ++j;
        if (j < b.edge.size() && a.edge[i].n == b.edge[j].n) {
            a.edge[i].w = rng(2) ? b.edge[j].w : a.edge[i].w;
            a.edge[i].e = (a.edge[i].e && b.edge[j].e) ? 1 : (rdn() > dsi);
        }
    }
    if (rdn() > pud) a = mft(a);
    a.reg();
    return a;
}

struct spec {
    int n = -1, ag = 0, ls = 0; //identifier, species age, time since last improvement
    double F = 0, mf = -1000000000; //current fitness, max fitness so far
    vector<net> pop; //species subpopulation
    net rep;

    spec sbp(int m, int r) { //bottleneck to best m and regenerate to r
        int a, b;
        sort(pop.begin(), pop.end(), invnetcmp());
        spec spr = *this;
        spr.pop.clear();
        if (r >= lar) {printf("sp %d: descending best unmutated\n"); spr.pop.push_back(pop[0]);}
        while (spr.pop.size() < r) {
            printf("sp %d: descending %d of %d\n", n, spr.pop.size(), r);

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
        if (ag > stg * 3) F /= ag / stg - 2;
    }

    static bool eq(net a, net b) { //check for same species
        int c, mt = 0, i;
        double dst = 0, mdw = 0;
        if (a.edge.size() < b.edge.size()) swap(a, b);
        if (a.edge.empty()) return 1;
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
net iet, nnet; //initial net, empty net
spec nspec;

void augm() {
    int np = 0, gl = 0, gm = 0, h, sf = 0, a, b, i, j, k;
    double sum = 0;
    spec tsp;
    net na, nb;
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
    for (i = 0; i < spc[0].size(); ++i) if (!spc[0][i].pop.empty()) spc[1].push_back(spc[0][i]);
    swap(spc[0], spc[1]);
    spc[1].clear();
    printf("%d species\n", spc[0].size());
    for (i = 0; i < spc[0].size(); ++i) printf("sp %d: %d\n", i, spc[0][i].pop.size());

    if (spc[0].size() > 1) for (i = 0; i < N; ++i) np += rdn() < isp;

    printf("%d normal %d interspecies:\n", N - np, np);
    for (i = 0; i < N; ++i) printf("net %d: score %f\n", i, ppn[0][i].F);

    for (i = 0; i < spc[0].size(); ++i) spc[0][i].cmf(), sum += spc[0][i].F;
    for (i = 0; i < spc[0].size(); ++i) cbp.push_back(make_pair(max((int) (spc[0][i].pop.size() * sur), 1), (N - np) * spc[0][i].F / sum)), gl += max((int) (spc[0][i].pop.size() * sur), 1), sf += (N - np) * spc[0][i].F / sum;
    for (i = sf; i < N - np; ++i) ++cbp[rng(spc[0].size())].second;

    for (i = 0; i < spc[0].size(); ++i) printf("species %d: %d bottleneck %d final\n", i, cbp[i].first, cbp[i].second);

    for (i = 0; i < spc[0].size(); ++i) if (cbp[i].second) {
        tsp = spc[0][i].sbp(cbp[i].first, cbp[i].second);

        printf("species %d descended successfully\n", i);

        for (j = 0; j < tsp.pop.size(); ++j) ppn[1].push_back(tsp.pop[j]), ppn[1].back().reg();
        spc[1].push_back(spc[0][i]);

        printf("species %d: %d genomes\n", i, spc[1].back().pop.size());

        spc[1].back().rep = spc[1].back().pop[(int) rng(spc[1].back().pop.size())];
        spc[1].back().pop.clear();
        ++spc[1].back().ag;
        if (spc[1].back().mf < spc[1].back().F && !eq(spc[1].back().mf, spc[1].back().F)) {spc[1].back().mf = spc[1].back().F; spc[1].back().ls = 0;}
        else ++spc[1].back().ls;
        spc[1].back().F = 0;
    }

    for (i = 0; i < np; ++i) {
        printf("interspecies %d / %d\n", i, np);
        a = rng(gl);
        for (j = h = 0; h < a + 1; ++j) h += cbp[j].first;
        --j;

        //printf("?! %d / %d -> %d: %d / %d\n", a, gl, j, a - h + cbp[j].first, cbp[j].first);

        na = spc[0][j].pop[a - h + cbp[j].first];
        gm = gl - cbp[j].first;
        b = rng(gm);
        for (k = h = 0; h < b + 1; ++k) if (j != k) h += cbp[k].first;
        --k;

        //printf("!? %d / %d -> %d: %d / %d\n", b, gm, k, b - h + cbp[k].first, cbp[k].first);

        nb = spc[0][k].pop[b - h + cbp[k].first];
        ppn[1].push_back(rdn() < pud ? ncr(na, nb) : mft(ncr(na, nb)));
        ppn[1].back().reg();
    }

    printf("%d new species\n", spc[1].size());
    for (i = 0; i < spc[1].size(); ++i) printf("sp %d: %d\n", i, spc[1][i].pop.size());

    swap(spc[0], spc[1]);
    spc[1].clear();
    swap(ppn[0], ppn[1]);
    for (i = 0; i < N; ++i) ppn[0][i].E = ppn[0][i].F = 0;
    ppn[1].clear();
    idn.clear();
    idl.clear();
}

inline void init() {
    int M, L, n, a, b, e, al0, al1, i, j;
    double w;
    bool cut; //custom nets

    printf("getting src\n");

    FILE* srf = fopen("src.txt", "r");
    fscanf(srf, "%d %d %d %d", &N, &T, &n0, &n1);
    fscanf(srf, "%d %d %d %d", &inn, &inw, &ins, &cut);
    if (cut) {
        for (i = 0; i < N; ++i) {
            ppn[0].push_back(nnet);
            fscanf(srf, "%d %d %d %lf", &M, &L, &ppn[0].back().E, &ppn[0].back().F);
            for (j = 0; j < M; ++j) fscanf(srf, "%d", &n), ppn[0].back().node.push_back(n);
            for (j = 0; j < L; ++j) {
                fscanf(srf, "%d %d %d %lf %d", &n, &a, &b, &w, &e);
                edg eg = {.n = n, .a = a, .b = b, .w = w, .e = e};
                ppn[0].back().edge.push_back(eg);
            }
            ppn[0].back().reg();
        }
    } else {
        iet.init();
        for (i = 0; i < N; ++i) ppn[0].push_back(iet), ppn[0][i].F = 1;
        augm();
    }
    fclose(srf);

    printf("src received\n");
}

inline void store(string stg) {
    cout << "storing to " << stg << endl;

    int i, j;
    FILE* res = fopen(stg.c_str(), "w");
    fprintf(res, "%d %d %d %d\n", N, T, n0, n1);
    fprintf(res, "%d %d %d %d\n", inn, inw, ins, 1);
    for (i = 0; i < N; ++i) {
        ppn[0].push_back(nnet);
        fprintf(res, "\n%d %d %d %f\n", ppn[0][i].node.size(), ppn[0][i].edge.size(), ppn[0][i].E, ppn[0][i].F);
        for (j = 0; j < ppn[0][i].node.size(); ++j) fprintf(res, "%d%c", ppn[0][i].node[j], j < ppn[0][i].node.size() - 1 ? ' ' : '\n');
        for (j = 0; j < ppn[0][i].edge.size(); ++j) fprintf(res, "%d %d %d %f %d\n", ppn[0][i].edge[j].n, ppn[0][i].edge[j].a, ppn[0][i].edge[j].b, ppn[0][i].edge[j].w, (int) ppn[0][i].edge[j].e);
    }
    fclose(res);
}

inline void eval(int n) {
    int x, tL = 0, i;
    FILE* flag, *in, *out;
    flag = fopen("run.txt", "w");
    fclose(flag);
    while (1) {
        if ((H = CreateFileA("end.txt", GENERIC_READ | GENERIC_WRITE, 0, NULL, OPEN_EXISTING, FILE_FLAG_DELETE_ON_CLOSE, NULL)) != INVALID_HANDLE_VALUE) {
            in = _fdopen(_open_osfhandle((int) H, 0), "r");
            fscanf(in, "%lf", &ppn[0][n].F);
            printf("score %f\n", ppn[0][n].F);
            ppn[0][n].E = 1;
            fclose(in);
            CloseHandle(H);
            break;
        }
        if ((H = CreateFileA("input.txt", GENERIC_READ | GENERIC_WRITE, 0, NULL, OPEN_EXISTING, FILE_FLAG_DELETE_ON_CLOSE, NULL)) != INVALID_HANDLE_VALUE) {
            in = _fdopen(_open_osfhandle((int) H, 0), "r");
            printf("input %d\n", tL++);

            for (i = 0; i < n0; ++i) fscanf(in, "%d", &x), qv.push_back(x);
            fclose(in);
            CloseHandle(H);
            ppn[0][n].query();
            while (!(out = fopen("output.txt", "w"))) this_thread::sleep_for(chrono::microseconds(1000));
            for (i = 0; i < esr.size(); ++i) fprintf(out, "%f%c", esr[i], i < esr.size() - 1 ? ' ' : '\n');
            fclose(out);
            esr.clear();
        }
        this_thread::sleep_for(chrono::microseconds(1000));
    }
}

char buf[1000]; //epoch process buffer

inline void epoch(int tI) {
    printf("\nepoch %d\n", tI);
    for (int i = 0; i < N; ++i) if (!ppn[0][i].E) {
        printf("net %d\n", i);
        eval(i);
    }
    store("intermediates/epoch " + to_string(tI) + " end.txt");
    if (tI + 1 == T) store("res.txt");
    printf("epoch %d evaluated\n", tI);
    augm();
    store("intermediates/epoch " + to_string(tI + 1) + " begin.txt");
}

int main() {
    SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);
    int tI = 0, i;
    FILE* flag;
    //srand(0);
    srand(time(0));
    while (1) {
        if (flag = fopen("init.txt", "r")) {
            fclose(flag);
            remove("init.txt");
            printf("init flag received\n");
            break;
        }
        this_thread::sleep_for(chrono::milliseconds(1000));
    }
    init();
    while (tI < T) epoch(tI++);
    flag = fopen("kys.txt", "w");
    fclose(flag);
    return 0;
}
