#include <bits/stdc++.h>
using namespace std;

int A[8];
double B[8], r;

int main() {
    int N, i, j;
    for (i = 0; i < 8; ++i) A[i] = __builtin_popcount(i);
    FILE* fl = fopen("init.txt", "w");
    fprintf(fl, "%d %d\n", 3, 1);
    fclose(fl);

    printf("init.txt created\n");

    this_thread::sleep_for(chrono::milliseconds(2000));
    while (1) {
        if (fl = fopen("kys.txt", "r")) {
            fclose(fl);
            remove("kys.txt");
            break;
        }
        if (fl = fopen("run.txt", "r")) {
            printf("tester beginning ops\n");

            fclose(fl);
            remove("run.txt");
            for (i = 0; i < 8; ++i) {
                printf("attempting %d\n", i);

                fl = fopen("input.txt", "w");
                for (j = 0; j < 3; ++j) fprintf(fl, "%d%c", i & (1 << j), i < 7 ? ' ' : '\n');
                fclose(fl);
                this_thread::sleep_for(chrono::milliseconds(100));
                if (!(fl = fopen("output.txt", "r"))) {printf("FAILURE\n"); while(1);}
                fscanf(fl, "%lf", &B[i]);
                fclose(fl);
                remove("output.txt");

                printf("finished %d\n", i);
            }
            for (i = r = 0; i < 8; ++i) r += (A[i] - B[i]) * (A[i] - B[i]);
            fl = fopen("end.txt", "w");
            fprintf(fl, "%lf", r);
            fclose(fl);

            printf("tester ending ops\n");
        }
        this_thread::sleep_for(chrono::milliseconds(100));
    }
    return 0;
}
