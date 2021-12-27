#include <stdio.h>
#include <math.h>

/*
exam
y1' = 4x + 3y1                <1식>

h =  0.05
y1(0) = 3
a ~ b = 0 ~ 1
i = (b - a) / h = (1.0 - 0.0) / 0.5 = 4

*/

double f(int i, double x, double y[])
{
    switch(i){
        case 1 : return(4 * x + 3 * y[i]);
        //case 2 : return(4 - 0.3 * y[2] - 0.1 * y[1]);
    }
    return -1;
}

void Run_kut(int n, double x, double y[], double h);

int main()
{
    int i, n;
    double x, y[3], h;

    printf("Input Data : n, h, x, y10, y20\n");
    scanf("%d%lf%lf%lf%lf",&n, &h, &x, &y[1], &y[2]);

    printf("No  x       y[1]    y[2]\n");
    
    for(i = 1; i <= 20; i++) // 4 는 (0 ~ 1 길이) / (H = 0.05 증분) 의 결과임
    {
        Run_kut(n, x, y, h);
        x += h;
        printf("%2d %6.3lf  %6.3lf  %6.3lf\n",i, x, y[1], y[2]);
    }

    return 0;
}

void Run_kut(int n, double x, double y[], double h)
{
    int i;
    double k1[3], k2[3], k3[3], k4[3];
    double y_temp1[3], y_temp2[3], y_temp3[3];

    for(int i = 1; i <= n; i++)
    {
        k1[i] = h * f(i, x, y);
        y_temp1[i] = y[i] + k1[i] / 2.0;
    }
    for(int i = 1; i <= n; i++)
    {
        k2[i] = h * f(i, x + h / 2.0, y_temp1);
        y_temp2[i] = y[i] + k2[i] / 2.0;
    }
    for(int i = 1; i <= n; i++)
    {
        k3[i] = h * f(i, x + h / 2.0, y_temp2);
        y_temp3[i] = y[i] + k3[i];
    }
    for(int i = 1; i <= n; i++)
    {
        k4[i] = h * f(i , x + h, y_temp3);
        y[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
    }
}

/* 
4계 룽게쿠타법 공식 : 

y(x + h) = y(x) + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4)

이때 x y<- 처음엔 초기값 다음부턴 이전 루프에서 계산된 값이 들어가게된다 /
k1 = h * f(x , y)
k2 = h * f(x + 1/2 * h, y + 1/2 * k1)
k3 = h * f(x + 1/2 * h, y + 1/2 * k2)
k4 = h * f(x + h, x + k3)
*/