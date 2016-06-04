#include "Matrix.h"
#include <cstdio>
#include <cstring>
#include <algorithm>
using namespace std;

void Matrix::init() {
    scanf("%d", &N);
    for(int i = 0; i< N; i++) {
        for(int j = 0; j < N; j++)
            scanf("%lf", &a[i][j]);
        scanf("%lf", &b[i]);
    }
}

void Matrix::showN() {
    printf("N = %d\n", N);
}

void Matrix::showA() {
    printf("A:\n");
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            printf((j == N-1)?"%.3f\n":"%.3f ", a[i][j]);
}

void Matrix::showb() {
    printf("b:\n");
    for(int i = 0; i < N; i++) printf("%.3f\n", b[i]);
}

void Matrix::showx() {
    printf("x:\n");
    for(int i = 0; i < N; i++) printf("x%d = %.3f\n", i, x[i]);
}

void Matrix::show(){
    printf("--\nN = %d:\n", N);
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++) printf("%.3f ", a[i][j]);
        printf("%.3f\n", b[i]);
    }

}

int Matrix::GuassInOrder() {
    for(int i = 0; i < N - 1; i++) {
        for(int j = i+1; j < N; j++) {
            double aji = a[j][i];
            for(int k = i; k < N; k++) {
                if(a[i][i] == 0) return 0 * printf("div by 0 in elimination.\n");
                a[j][k] -= a[i][k] * aji / a[i][i];
            }
            b[j] -= b[i] * aji / a[i][i];
        }
    }

    for(int i = N - 1; i > -1; i--) {
        if(a[i][i] == 0) return 0 * printf("div by 0 when get solution.\n");
        x[i] = b[i] / a[i][i];
        for(int j = i-1; j > -1; j--) b[j] -= a[j][i] * x[i];
    }
    showx();
}

void Matrix::swapRow(int _x, int _y){
    if(_x == _y) return ;
    for(int i = 0; i < N; i++) swap(a[_x][i], a[_y][i]);
    swap(b[_x], b[_y]);
}

int Matrix::GuassColumnPrincipleComponent(){
   for(int i = 0; i < N - 1; i++) {

        /* Change rows */
        double maxx = -1.0;
        int rowId;
        for(int j = i; j < N; j++)
            if(abs(a[j][i]) > maxx){
                maxx = abs(a[j][i]);
                rowId = j;
            }
        swapRow(i, rowId);

        /* calculation */
        for(int j = i+1; j < N; j++) {
            double aji = a[j][i];
            for(int k = i; k < N; k++) {
                if(a[i][i] == 0) return 0 * printf("div by 0 in elimination.\n");
                a[j][k] -= a[i][k] * aji / a[i][i];
            }
            b[j] -= b[i] * aji / a[i][i];
        }
        //show();
    }

    for(int i = N - 1 ; i > -1; i--) {
        if(a[i][i] == 0) return 0 * printf("div by 0 when get solution.\n");
        x[i] = b[i] / a[i][i];
        for(int j = i-1; j > -1; j--) b[j] -= a[j][i] * x[i];
    }
    showx();
}

int Matrix::SquareRoot(){
    double g[100][100];
    memset(g, 0, sizeof(g));
    for(int i = 0; i < N; i++){

        /* g[i][0]~g[i][i-1] */
        for(int j = 0; j < i; j++){
            double sum = 0;
            for(int k = 0; k < j; k++) sum += g[i][k]*g[j][k];
            if(g[j][j] == 0) return 0 * printf("div by 0: G[%d][%d] = 0.\n", j ,j);
            g[i][j] = (a[i][j] - sum) / g[j][j];
        }

        /* g[i][i] */
        double sum = 0;
        for(int j = 0; j < i; j++) sum += g[i][j] * g[i][j];
        if(a[i][i] - sum < 0) return 0 * printf("sqrt()'s parameter is negative.\n");
        g[i][i] = sqrt(a[i][i] - sum);
    }

    printf("G:\n");
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            printf(j == N-1 ? "%.3f\n" : "%.3f ", g[i][j]);

    /* Solve G y = b */
    double y[100];
    for(int i = 0; i < N; i++){
        if(g[i][i] == 0) return 0 * printf("div by 0 when solve Gy = b .\n");
        y[i] = b[i] / g[i][i];
        for(int j = i + 1; j < N; j++) b[j] -= g[j][i] * y[i];
    }
    for(int i = 0; i < N; i++) printf("y%d = %.3f\n", i, y[i]);

    /* Reverse G */
    for(int i = 0; i < N; i++)
    for(int j = i+1; j < N; j++){
        g[i][j] = g[j][i];
        g[j][i] = 0;
    }

    /* Solve G' x = y */
    for(int i = N - 1; i > -1; i--){
        if(g[i][i] == 0) return 0 * printf("div by 0 when solve G'x = y .\n");
        x[i] = y[i] / g[i][i];
        for(int j = i - 1; j > -1; j--) y[j] -= g[j][i] * x[i];
    }

    showx();
}

int Matrix::Chasing(){
    double alpha[100], beta[100], gamma[100];
    for(int i = 1; i < N; i++) gamma[i] = a[i][i-1];
    for(int i = 0; i < N; i++){
        alpha[i] = (i == 0 ? a[i][i] : a[i][i] - gamma[i] * beta[i-1]);
        if(i < N-1){
            if(alpha[i] == 0) return 0 * printf("div by 0: alpha[%d] = 0\n", i);
            beta[i] = a[i][i+1] / alpha[i];
        }
    }
    //for(int i = 0; i < N; i++) printf("%.3f\n",beta[i]);

    /* Solve Ty = b */
    double y[100];
    for(int i = 0; i < N; i++){
        if(alpha[i] == 0) return 0 * printf("div by 0 when solve Ty = b : alpha[%d] = 0\n", i);
        y[i] = b[i] / alpha[i];
        if(i < N-1) b[i+1] -= gamma[i+1] * y[i];
    }
    for(int i = 0; i < N; i++) printf("y%d = %.3f\n", i, y[i]);

    /* Solve Mx = y */
    for(int i = N-1; i > -1; i--){
        x[i] = y[i];
        if(i > 0) y[i-1] -= beta[i-1] * x[i];
    }
    showx();
}
