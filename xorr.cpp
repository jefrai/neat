#include <bits/stdc++.h>
#include <windows.h>
#include <io.h>
using namespace std;

//remember to use -U__STRICT_ANSI__ compiler flag

int A[8];
double B[8], r;

int main() {
    int N, i, j;
    for (i = 0; i < 8; ++i) A[i] = __builtin_popcount(i) % 2;
    HANDLE H;
    FILE* fl = _fdopen(_open_osfhandle((int) (H = CreateFileA("init.txt", GENERIC_READ | GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, 0, NULL)), 0), "w");
    fclose(fl);
    CloseHandle(H);

    printf("init.txt created\n");

    this_thread::sleep_for(chrono::milliseconds(2000));
    while (1) {
        if (fl = fopen("kys.txt", "r")) {
            fclose(fl);
            remove("kys.txt");
            break;
        }
        if (fl = fopen("run.txt", "r")) {
            printf("running net\n");

            fclose(fl);
            while (!remove("run.txt"));
            for (i = 0; i < 8; ++i) {
                printf("input %d\n", i);

                while (!(fl = fopen("input.txt", "w"))) this_thread::sleep_for(chrono::microseconds(1000));
                for (j = 0; j < 3; ++j) fprintf(fl, "%d%c", !!(i & (1 << j)), i < 7 ? ' ' : '\n');
                fclose(fl);
                this_thread::sleep_for(chrono::milliseconds(100));
                while (!(fl = fopen("output.txt", "r"))) printf("FAILURE\n"), this_thread::sleep_for(chrono::milliseconds(30));
                fscanf(fl, "%lf", &B[i]);
                fclose(fl);
                remove("output.txt");
            }
            for (i = r = 0; i < 8; ++i) r += 1 - sqrt(fabs(A[i] - B[i]));
            while (!(fl = fopen("end.txt", "w"))) this_thread::sleep_for(chrono::microseconds(1000));
            printf("score %f\n", r);
            fprintf(fl, "%f", r);
            fclose(fl);

            printf("finished net\n");
        }
        this_thread::sleep_for(chrono::milliseconds(10));
    }
    return 0;
}
