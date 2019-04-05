#include <bits/stdc++.h>
using namespace std;

int N = 150, inv; //total population size, innovation number
int lar = 15; //threshold for large population - unchanged copy of champion, new mutation constants
int kil = 15; //max generations w/o improvement for species
double sig = 4.9; //sigmoid multiplier
double cw1 = 1, cw2 = 1, cw3 = 0.4, delt = 3; //excess, disjoint, matching weight, threshold - compatibility constants
double mut = 0.8, rel = 0.9; //point mutation rate, relative mutation vs random reassignment
double dsd = 0.75; //gene deactivation inheritance
double pur = 0.25; //probability of descent with mutation only
vector<vector<

vector<int> spec; //species



inline double sgm(double x) {return 1 / (1 + exp(-sig * x));}

inline void init(int _N) {
    N = _N;
}

int main() {
    int i;
    init(150);
}
