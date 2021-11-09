#include <stdio.h>
#include <math.h>

double f(int i, double x, double y[])
{
    switch(i){
        case 1 : return(5.0*y[2] - 10.0*y[1]);
        case 2 : return(10.0*y[1] - 5.0*y[2]);
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

    printf("No  x       y[0]    y[1]\n");
    for(i = 1;i <= 10;i ++)
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

    for(int i=1;i<=n;i++)   //k1
    {
        k1[i] = h * f(i, x, y);
        y_temp1[i] = y[i] + k1[i]/2.0;
    }
    for(int i=1;i<=n;i++)   //k2
    {
        k2[i] = h * f(i, x + h/2.0, y_temp1);
        y_temp2[i] = y[i] + k2[i]/2.0;
    }
    for(int i=1;i<=n;i++)   //k3
    {
        k3[i] = h * f(i, x + h/2.0, y_temp2);
        y_temp3[i] = y[i] + k1[i]/2.0;
    }
    for(int i=1;i<=n;i++)   //k4
    {
        k4[i] = h * f(i, x + h, y_temp3);
        y[i] += (k1[i] + k2[i] * 2.0 + k3[i] * 2.0 + k4[i])/6.0;
    }

}