//Elliptic Partial Differential Equations(타원형 편미분방정식)
//SOR 방법 사용


/*

    U
*/

#include <stdio.h>
#include <math.h>

#define EPS 1e-5
#define MAX_ITER 100

int main() {
    int i, j ,count;
    int sect_x, sect_y;
    double dU;
    double U[20][20];
    double rel_coeff, fxy, norm;
    double x_0, x_n, y_0, y_n;
    double bound_a, bound_b, bound_c, bound_d;
    double dx, dy, sqr_dx, sqr_dy, sqr_sum, dx_dy;

    printf("Write the section No. for x & y direction : ");
    scanf("%d%d", &sect_x, &sect_y);
    printf("Input positions(x, y) for boundary condition : ");
    scanf("%lf%lf%lf%lf", &x_0, &y_0, &x_n, &y_n);          //격자 범위 입력
    printf("Input function value at every position : ");
    scanf("%lf%lf%lf%lf", &bound_a, &bound_b, &bound_c, &bound_d);  //경계조건 입력

    fxy = 0.0;
    rel_coeff = 1.3333333;

    dx = (x_n-x_0) / sect_x;
    dy = (y_n-y_0) / sect_y;
    sqr_dx = dx*dx;
    sqr_dy = dy*dy;
    sqr_sum = 2.0 * (sqr_dx + sqr_dy);
    dx_dy = 1.0 / sqr_sum;

    for(i = 1; i < sect_y; i++){
        U[0][i] = bound_a;
        U[sect_x][i] = bound_b;
    }
    for(i = 1; i < sect_x; i++){
        U[i][0] = bound_c;
        U[i][sect_y] = bound_d;
        for(j = 1; j < sect_y; j++){   
            U[i][j] = 0.0;
        }
    }

    count = 0;
    do{
        norm = 0.0;
        for(i = 1; i < sect_x; i++){
            for(j = 1; j < sect_y; j++){
                dU = dx_dy * ( (U[i-1][j] + U[i+1][j]) * sqr_dy + (U[i][j-1] + U[i][j+1]) * sqr_dx - sqr_dx * sqr_dy * fxy / sqr_sum) - U[i][j];
                U[i][j] += (rel_coeff * dU);
                norm += fabs(dU*rel_coeff);
            }
        }

        printf("    %10d    %18.7le %18.7le \n", ++count, norm, dU);
    }while((count<MAX_ITER)&&(norm>EPS));

    printf("Number of iteration = %5d\n", count);
    printf("At each mesh of x & y, value as U(x,y)\n");

    printf("         2         3         4        5        6\n");
    for(i = 1; i < sect_x; i++){
        printf("%3d", i + 1);
        for(j = 1; j < sect_y; j++){
            printf("%10.4lf", U[i][j]);
        }
        printf("\n");
    }
    
    return 0;
}